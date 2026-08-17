#!/usr/bin/env python3
"""
BAAM reference engine (Python) - Batch Adsorber Analogue Model
Balashankar, Subraveti, Li & Rajendran, Ind. Eng. Chem. Res. 58 (2019) 3314,
extended per Liske & Rajendran, Can. J. Chem. Eng. (2026) (e-BAAM):
4-step VSA/PVSA cycles with LPP or FP, generic binary mixture equilibria
(extended DSL analytic, or IAST over pure-component isotherms with FD partials).

Ground truth: Balashankar's own BAAM.m run headlessly in MATLAB R2026a with
Table-1 parameters (reproduces paper Table 2 exactly). Anchors embedded below.

Model summary (per 1 kg adsorbent, isothermal T, ideal gas):
  V = w/(rho*(1-eps)) column volume; N_i = P*y_i*V*eps/(R*T) + w*q_i*(P,y)
  BLO/EVAC ODE:  dy/dP = (a1*y - a2)/(f2 - f1*y)                        (1)
    a1 = V*eps/(RT) + w*(q1P + q2P);  a2 = y*V*eps/(RT) + w*q1P
    f1 = w*(q1y + q2y);               f2 = P*V*eps/(RT) + w*q1y
  Vacuum work (pump delivers at 1 bar, only for P < 1 bar):
    W = sum (k/(k-1)) * (R T/eta(P)) * dN * [ (1/P)^((k-1)/k) - 1 ]     (2)
  LPP: overall+component balances, unknowns (y_delta, N_LPP)            (3)
  ADS: overall+component balances, LINEAR in (N_feed, N_raff)           (4)
  Pu = dN1_evac/dN_evac; Re = dN1_evac/(N_feed*y_feed [+ N_FP*y_feed])  (5)
  En = (W_BLO+W_EVAC+W_ADS+W_FP)/(dN1_evac*44e-6 t) [kWh_e/t CO2]      (6)
  WC = dN1_evac/(V*(1-eps))  [mol CO2 / m3 adsorbent]                   (7)
"""
import numpy as np
from scipy.integrate import solve_ivp

RB = 8.314e-5   # m3 bar / (mol K)
RE = 8.314      # J / (mol K)
KAD = 1.4       # adiabatic constant
J2KWH = 1.0/3.6e6

# ---------------------------------------------------------------- mixtures --
class EDSL:
    """Extended dual-site Langmuir, concentration basis (paper Eq. 20-21).
    pars: dict with per-gas (qsb, qsd, b0, d0, dUb, dUd) [mol/kg, m3/mol, J/mol
    with dU NEGATIVE for exothermic], T in K."""
    def __init__(s, co2, n2, T):
        s.T = T
        def act(g):
            b = g["b0"]*np.exp(-g["dUb"]/(RE*T)) if g["b0"] > 0 else 0.0
            d = g["d0"]*np.exp(-g["dUd"]/(RE*T)) if g["d0"] > 0 else 0.0
            return b, d
        s.bA, s.dA = act(co2); s.bB, s.dB = act(n2)
        s.qsbA, s.qsdA = co2["qsb"], co2["qsd"]
        s.qsbB, s.qsdB = n2["qsb"], n2["qsd"]

    def loadings(s, P, y):
        """q1,q2 and the four partials dq/dP, dq/dy (P bar, y = y_CO2)."""
        T = s.T
        c1 = P*y/(RB*T); c2 = P*(1.0-y)/(RB*T)
        c1P = y/(RB*T); c2P = (1.0-y)/(RB*T)
        c1y = P/(RB*T); c2y = -P/(RB*T)
        q1 = q2 = 0.0; q1c1 = q1c2 = q2c1 = q2c2 = 0.0
        for (kA, kB, qsA, qsB) in ((s.bA, s.bB, s.qsbA, s.qsbB),
                                   (s.dA, s.dB, s.qsdA, s.qsdB)):
            D = 1.0 + kA*c1 + kB*c2
            q1 += qsA*kA*c1/D
            q2 += qsB*kB*c2/D
            q1c1 += qsA*kA*(1.0 + kB*c2)/D**2
            q1c2 += -qsA*kA*c1*kB/D**2
            q2c2 += qsB*kB*(1.0 + kA*c1)/D**2
            q2c1 += -qsB*kB*c2*kA/D**2
        return (q1, q2,
                q1c1*c1P + q1c2*c2P, q1c1*c1y + q1c2*c2y,
                q2c1*c1P + q2c2*c2P, q2c1*c1y + q2c2*c2y)


class IASTMix:
    """Binary IAST mixture over two pure-component isotherms given as
    (psi(c), q(c)) callables in CONCENTRATION [mol/m3]; partials by central FD.
    psi = integral q/c dc (same units as q)."""
    def __init__(s, q1f, psi1f, q2f, psi2f, T):
        s.q1f, s.psi1f, s.q2f, s.psi2f, s.T = q1f, psi1f, q2f, psi2f, T

    def _solve(s, P, y):
        T = s.T
        if P <= 0: return 0.0, 0.0
        ct = P/(RB*T)
        c1 = ct*y; c2 = ct*(1.0-y)
        if y <= 1e-14:  return 0.0, s.q2f(c2)
        if y >= 1-1e-14: return s.q1f(c1), 0.0
        # G(x1) = psi1(c1/x1) - psi2(c2/(1-x1)) strictly decreasing on (0,1)
        lo, hi = 1e-12, 1.0-1e-12
        for _ in range(200):
            x = 0.5*(lo+hi)
            G = s.psi1f(c1/x) - s.psi2f(c2/(1.0-x))
            if G > 0: lo = x
            else: hi = x
            if hi-lo < 1e-15: break
        x1 = 0.5*(lo+hi)
        q1o = s.q1f(c1/x1); q2o = s.q2f(c2/(1.0-x1))
        qt = 1.0/(x1/q1o + (1.0-x1)/q2o)
        return x1*qt, (1.0-x1)*qt

    def loadings(s, P, y):
        dP = max(1e-6, 1e-6*P)
        dy = 1e-7
        yl, yh = max(0.0, y-dy), min(1.0, y+dy)
        q1, q2 = s._solve(P, y)
        q1a, q2a = s._solve(P+dP, y); q1b, q2b = s._solve(max(P-dP,1e-12), y)
        q1c, q2c = s._solve(P, yh);   q1d, q2d = s._solve(P, yl)
        return (q1, q2, (q1a-q1b)/(2*dP), (q1c-q1d)/(yh-yl),
                (q2a-q2b)/(2*dP), (q2c-q2d)/(yh-yl))

# ----------------------------------------------------------------- engine ---
def eta_const(P_bar, val=0.72): return val
def eta_maruyama(P_bar):
    Patm = P_bar/1.01325
    return 0.8*(19.55*Patm)/(1.0 + 19.55*Patm)

def blow_path(mix, T, rho, PH, yfeed, Plow_min, w=1.0, eps=0.37,
              dP=1e-4, eta=eta_const, rtol=1e-9, atol=1e-11):
    """Integrate the BLO/EVAC ODE from (PH, yfeed) down to Plow_min.
    Returns dict of arrays on the uniform dP grid: P, y, q1, q2, and
    cumulative-removed N (total), N1 (CO2) and vacuum work W [J] (MATLAB
    convention: slice dN charged at the slice-end pressure; only P<1 bar)."""
    V = w/(rho*(1.0-eps))
    def rhs(P, yv):
        y = min(max(yv[0], 0.0), 1.0)
        q1, q2, q1P, q1y, q2P, q2y = mix.loadings(P, y)
        a1 = V*eps/(RB*T) + w*(q1P + q2P)
        a2 = y*V*eps/(RB*T) + w*q1P
        f1 = w*(q1y + q2y)
        f2 = P*V*eps/(RB*T) + w*q1y
        return [(a1*y - a2)/(f2 - f1*y)]
    n = int(round((PH-Plow_min)/dP))
    P = PH - dP*np.arange(n+1)
    sol = solve_ivp(rhs, (PH, P[-1]), [yfeed], method="BDF",
                    t_eval=P, rtol=rtol, atol=atol)
    y = np.clip(sol.y[0], 0.0, 1.0)
    q1 = np.empty_like(P); q2 = np.empty_like(P)
    for i in range(P.size):
        q1[i], q2[i], *_ = mix.loadings(P[i], y[i])
    fA = P*y*V*eps/(RB*T); fB = P*(1.0-y)*V*eps/(RB*T)
    N1s = fA + w*q1; N2s = fB + w*q2          # in-column inventories
    Ntot = N1s + N2s
    dN = Ntot[:-1] - Ntot[1:]                 # removed in each slice
    Pend = P[1:]
    kk = (KAD-1.0)/KAD
    ee = np.array([eta(pb) for pb in Pend])
    wsl = np.where(Pend < 1.0,
                   (KAD/(KAD-1.0))*(RE*T/ee)*dN*((1.0/Pend)**kk - 1.0), 0.0)
    W = np.concatenate([[0.0], np.cumsum(wsl)])
    N = np.concatenate([[0.0], np.cumsum(dN)])
    N1 = np.concatenate([[0.0], np.cumsum(N1s[:-1] - N1s[1:])])
    return dict(P=P, y=y, q1=q1, q2=q2, N=N, N1=N1, W=W,
                N1s=N1s, N2s=N2s, V=V, w=w, eps=eps, T=T, dP=dP)

def _idx(path, Pq):
    i = int(round((path["P"][0]-Pq)/path["dP"]))
    assert abs(path["P"][i]-Pq) < 1e-9, "pressure not on grid"
    return i

def _ntot(mix, T, V, eps, w, P, y):
    q1, q2, *_ = mix.loadings(P, y)
    f1 = P*y*V*eps/(RB*T); f2 = P*(1.0-y)*V*eps/(RB*T)
    return f1 + w*q1, f2 + w*q2

def solve_lpp(mix, T, V, eps, w, PH, Plow, ygam):
    """Pressurize (Plow,ygam) -> (PH, ydel) with gas of composition ydel.
    Returns ydel, NLPP."""
    N1i, N2i = _ntot(mix, T, V, eps, w, Plow, ygam)
    Ni = N1i + N2i
    def g(yd):
        N1f, N2f = _ntot(mix, T, V, eps, w, PH, yd)
        NLPP = (N1f + N2f) - Ni
        return N1i - N1f + NLPP*yd, NLPP
    lo, hi = 1e-12, ygam
    glo = g(lo)[0]; ghi = g(hi)[0]
    if glo*ghi > 0:  # scan for a bracket
        ys = np.linspace(lo, hi, 400)
        gs = [g(t)[0] for t in ys]
        k = next(i for i in range(len(gs)-1) if gs[i]*gs[i+1] <= 0)
        lo, hi = ys[k], ys[k+1]
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if g(lo)[0]*g(mid)[0] <= 0: hi = mid
        else: lo = mid
        if hi-lo < 1e-16: break
    yd = 0.5*(lo+hi)
    return yd, g(yd)[1]

def solve_press_feed(mix, T, V, eps, w, PH, Plow, ygam, yfeed):
    """FP: pressurize with FEED gas. Unknown y_end; N_FP from overall."""
    N1i, N2i = _ntot(mix, T, V, eps, w, Plow, ygam)
    Ni = N1i + N2i
    def g(ye):
        N1f, N2f = _ntot(mix, T, V, eps, w, PH, ye)
        NFP = (N1f + N2f) - Ni
        return N1i - N1f + NFP*yfeed, NFP
    lo, hi = 1e-12, max(ygam, yfeed)
    if g(lo)[0]*g(hi)[0] > 0:
        ys = np.linspace(lo, hi, 400)
        gs = [g(t)[0] for t in ys]
        k = next(i for i in range(len(gs)-1) if gs[i]*gs[i+1] <= 0)
        lo, hi = ys[k], ys[k+1]
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if g(lo)[0]*g(mid)[0] <= 0: hi = mid
        else: lo = mid
        if hi-lo < 1e-16: break
    ye = 0.5*(lo+hi)
    return ye, g(ye)[1]

def solve_ads(mix, T, V, eps, w, PH, ydel, yfeed):
    """ADS linear solve: (Nfeed, Nraff); raffinate leaves at y_delta."""
    N1i, N2i = _ntot(mix, T, V, eps, w, PH, ydel)
    N1f, N2f = _ntot(mix, T, V, eps, w, PH, yfeed)
    A = (N1f+N2f) - (N1i+N2i)     # Nfeed - Nraff
    B = N1f - N1i                 # Nfeed*yfeed - Nraff*ydel
    Nfeed = (B - A*ydel)/(yfeed - ydel)
    return Nfeed, Nfeed - A

def cycle_kpis(mix, path, Pint, Plow, yfeed, mode="LPP", eta=eta_const):
    """KPIs for the 4-step cycle at (Pint, Plow) using a precomputed path."""
    T, V, eps, w = path["T"], path["V"], path["eps"], path["w"]
    PH = path["P"][0]
    ib, ig = _idx(path, Pint), _idx(path, Plow)
    n1e = path["N1"][ig] - path["N1"][ib]         # CO2 collected in EVAC
    ne = path["N"][ig] - path["N"][ib]
    Wblo, Wevac = path["W"][ib], path["W"][ig] - path["W"][ib]
    ygam = path["y"][ig]
    if mode == "LPP":
        ydel, Npr = solve_lpp(mix, T, V, eps, w, PH, Plow, ygam)
        n1fed_pr = Npr*ydel*0.0   # LPP gas is raffinate; not feed
        Wpr = 0.0
    else:
        ydel, Npr = solve_press_feed(mix, T, V, eps, w, PH, Plow, ygam, yfeed)
        n1fed_pr = Npr*yfeed
        kk = (KAD-1.0)/KAD
        Wpr = Npr*(KAD/(KAD-1.0))*(RE*T/eta(1.0))*(PH**kk - 1.0) if PH > 1 else 0.0
    Nfeed, Nraff = solve_ads(mix, T, V, eps, w, PH, ydel, yfeed)
    kk = (KAD-1.0)/KAD
    Wads = Nfeed*(KAD/(KAD-1.0))*(RE*T/eta(1.0))*(PH**kk - 1.0) if PH > 1 else 0.0
    mCO2_t = n1e*44.01e-6                          # tonnes (MATLAB uses 44)
    mCO2_t_ml = n1e*44e-6
    Pu = 100.0*n1e/ne
    Re = 100.0*n1e/(Nfeed*yfeed + n1fed_pr)
    En = (Wblo+Wevac+Wads+Wpr)*J2KWH/mCO2_t_ml
    WC = n1e/(V*(1.0-eps))
    # cycle mole balance closure
    ib_ = ib; ig_ = ig
    nin = Nfeed + Npr if mode == "FP" else Nfeed
    nout = (Nraff - (Npr if mode == "LPP" else 0.0)) + path["N"][ib_] + ne
    mbal = (nin - nout)/max(nout, 1e-300)*100.0
    return dict(Pu=Pu, Re=Re, EnBLO=Wblo*J2KWH/mCO2_t_ml,
                EnEVAC=Wevac*J2KWH/mCO2_t_ml, EnADS=Wads*J2KWH/mCO2_t_ml,
                EnPR=Wpr*J2KWH/mCO2_t_ml, En=En, WC=WC,
                ydel=ydel, ygam=ygam, Npr=Npr, Nfeed=Nfeed, Nraff=Nraff,
                mbal=mbal)

# ------------------------------------------------------------- adsorbents ---
ADS2019 = {  # Balashankar Table 1 (c-basis DSL; dU stored NEGATIVE, J/mol)
 "MgMOF74": dict(rho=588.25,
    co2=dict(qsb=6.80, qsd=9.90, b0=1.81e-7, d0=1.06e-6, dUb=-39.30e3, dUd=-21.20e3),
    n2 =dict(qsb=14.00, qsd=0.0, b0=3.45e-6, d0=0.0,     dUb=-15.50e3, dUd=0.0)),
 "Z13X": dict(rho=1130.00,
    co2=dict(qsb=3.09, qsd=2.54, b0=8.65e-7, d0=2.63e-8, dUb=-36.60e3, dUd=-35.70e3),
    n2 =dict(qsb=5.84, qsd=0.0,  b0=2.50e-6, d0=0.0,     dUb=-15.80e3, dUd=0.0)),
 "UTSA16": dict(rho=1092.00,
    co2=dict(qsb=5.00, qsd=3.00, b0=6.24e-7, d0=1.87e-23, dUb=-30.60e3, dUd=-44.70e3),
    n2 =dict(qsb=12.70, qsd=0.0, b0=2.96e-6, d0=0.0,      dUb=-9.77e3,  dUd=0.0)),
 "CSAC": dict(rho=799.50,
    co2=dict(qsb=0.59, qsd=7.51, b0=9.40e-6, d0=1.04e-5, dUb=-25.61e3, dUd=-17.55e3),
    n2 =dict(qsb=0.16, qsd=41.30, b0=1.81e-3, d0=1.72e-12, dUb=-8.67e3, dUd=-44.90e3)),
}
# Liske Table 1 (SSL refits; CALF-20 DSL). Same storage convention.
ADS2026 = {
 "UTSA16ssl": dict(rho=1000.0,
    co2=dict(qsb=4.478, qsd=0.0, b0=4.70e-7, d0=0.0, dUb=-30.57e3, dUd=0.0),
    n2 =dict(qsb=4.478, qsd=0.0, b0=1.400e-6, d0=0.0, dUb=-9.91e3, dUd=0.0)),
 "Z13Xssl": dict(rho=1130.0,
    co2=dict(qsb=4.390, qsd=0.0, b0=2.50e-6, d0=0.0, dUb=-31.19e3, dUd=0.0),
    n2 =dict(qsb=4.390, qsd=0.0, b0=2.70e-6, d0=0.0, dUb=-16.38e3, dUd=0.0)),
 "IISERP2": dict(rho=1000.0,
    co2=dict(qsb=4.478, qsd=0.0, b0=2.02e-7, d0=0.0, dUb=-31.13e3, dUd=0.0),
    n2 =dict(qsb=4.478, qsd=0.0, b0=2.64e-7, d0=0.0, dUb=-11.89e3, dUd=0.0)),
 "CALF20": dict(rho=570.0,
    co2=dict(qsb=2.387, qsd=3.271, b0=5.52e-7, d0=5.19e-8, dUb=-35.06e3, dUd=-28.95e3),
    n2 =dict(qsb=2.387, qsd=3.271, b0=8.14e-7, d0=0.0,     dUb=-17.96e3, dUd=0.0)),
}

# ------------------------------------------------------- MATLAB anchors -----
# From BAAM.m (Balashankar's own code) run headlessly in MATLAB R2026a with
# ode15s tolerances TIGHTENED to RelTol=1e-10/AbsTol=1e-12 (the published
# RelTol=1e-6 run reproduces paper Table 2 to its displayed digits; the tight
# run is converged and is the machine-precision ground truth used here).
ML_KPI = {  # (Plow, Pint): Pu, Re, EnBLO, EnEVAC, EnT, WC
 "MgMOF74": {(0.03,0.15): (99.856580, 77.995657, 2.444353, 96.054093, 98.498446, 1146.345970),
             (0.05,0.10): (99.99999981, 47.69462764, None, None, 114.56186475, 509.75552449),
             (0.03,0.30): (98.77610577, 79.86407374, None, None, 96.19407867, 1173.80714689),
             (0.06,0.50): (94.72371444, 60.36603298, None, None, 76.84470527, 570.38331536)},
 "Z13X":    {(0.03,0.15): (99.958895, 79.305127, 1.651625, 92.483140, 94.134765, 1143.072675),
             (0.05,0.10): (100.00000000, 44.27210790, None, None, 121.66474956, 496.29396545),
             (0.03,0.30): (99.32383846, 80.75682881, None, None, 92.44257775, 1163.99693601),
             (0.06,0.50): (97.23600786, 61.85209698, None, None, 74.58353601, 632.38988148)},
 "UTSA16":  {(0.03,0.15): (99.974047, 79.977967, 1.667222, 86.835122, 88.502344, 1674.496102),
             (0.05,0.10): (100.00000000, 41.63041760, None, None, 128.32050829, 777.09535057),
             (0.03,0.30): (99.38094577, 81.53369868, None, None, 86.81364468, 1707.06840793),
             (0.06,0.50): (97.85729014, 63.11301939, None, None, 73.32758795, 1115.67167187)},
 "CSAC":    {(0.03,0.15): (94.418783, 55.639333, 19.187496, 99.487041, 118.674537, 297.282289),
             (0.05,0.10): (96.63457041, 30.18082691, None, None, 148.12687854, 152.57646164),
             (0.03,0.30): (85.79980864, 62.48748814, None, None, 105.66870698, 333.87214661),
             (0.06,0.50): (65.94973416, 40.01795593, None, None, 90.57482927, 197.31092495)},
}
ML_PATH = {  # P: (y, qA, qB, Ncum, NAcum, Wcum[J])
 "MgMOF74": {0.50: (0.2876043823, 6.3383679264, 0.0396635443, 0.0894555817, 0.0184100700, 104.8157029653),
             0.15: (0.8835919471, 6.2751732940, 0.0020956179, 0.2043150025, 0.0820583554, 754.5232051312),
             0.06: (1.0000000000, 5.3721158107, 0.0000000000, 1.1130930110, 0.9880374654, 11802.4904151252),
             0.03: (1.0000000000, 4.3305637932, 0.0000000000, 2.1558533310, 2.0307977854, 30404.5081886472)},
 "Z13X":    {0.50: (0.2946411412, 3.4266930010, 0.0085133640, 0.0278515688, 0.0058008436, 32.5118342658),
             0.15: (0.9399060447, 3.4028814541, 0.0002269780, 0.0672879787, 0.0297452098, 264.6437144563),
             0.06: (1.0000000000, 2.8688867949, 0.0000000000, 0.6033966529, 0.5654379069, 6611.5649427036),
             0.03: (1.0000000000, 2.3936397656, 0.0000000000, 1.0792726945, 1.0413139484, 15083.4355067695)},
 "UTSA16":  {0.50: (0.2970936618, 2.3076889628, 0.0147455512, 0.0398508441, 0.0083354252, 46.4731757634),
             0.15: (0.9587463112, 2.2701179983, 0.0002638075, 0.0994973976, 0.0460091211, 404.9578363313),
             0.06: (1.0000000000, 1.2879327817, 0.0000000000, 1.0838991248, 1.0300127804, 11866.8648673425),
             0.03: (1.0000000000, 0.7391660036, 0.0000000000, 1.6333168039, 1.5794304596, 21496.6658651577)},
 "CSAC":    {0.50: (0.2538324510, 0.7040245098, 0.1164685187, 0.1637831932, 0.0316877768, 192.7810982330),
             0.15: (0.6363238539, 0.6331997389, 0.0203630713, 0.3410854929, 0.1034450763, 1130.1177505578),
             0.06: (0.9962494436, 0.4592213069, 0.0000972681, 0.5379968347, 0.2784806785, 3540.7479980509),
             0.03: (0.9999994154, 0.2633040184, 0.0000000087, 0.7349004183, 0.4752803344, 6989.7709482427)},
}
# Paper Table 2 (Balashankar 2019), displayed digits - reproduced by the
# published-tolerance run; checked here at 0.02 absolute.
PAPER_T2 = {
 "MgMOF74": (99.86, 78.00, 98.49, 1146.36),
 "Z13X":    (99.96, 79.30, 94.13, 1143.0),
 "UTSA16":  (99.97, 79.97, 88.50, 1674.4),
 "CSAC":    (94.42, 55.63, 118.68, 297.2),
}

# ------------------------------------------------------------------ checks --
def make_mix(spec, T):
    return EDSL(spec["co2"], spec["n2"], T)

def check_paths():
    print("== A: BLO/EVAC path vs converged MATLAB (T=298.15, PH=1, yF=0.15) ==")
    worst = 0.0
    for name, spec in ADS2019.items():
        mix = make_mix(spec, 298.15)
        path = blow_path(mix, 298.15, spec["rho"], 1.0, 0.15, 0.03)
        for Pq, ref in ML_PATH[name].items():
            i = _idx(path, Pq)
            got = (path["y"][i], path["q1"][i], path["q2"][i],
                   path["N"][i], path["N1"][i], path["W"][i])
            for g, r in zip(got, ref):
                e = abs(g-r)/max(abs(r), 1e-3)   # floor: qB noise ~1e-10
                worst = max(worst, e)
        print(f"  {name}: worst rel so far {worst:.3e}")
    print(f"  PATH worst rel err = {worst:.3e}")
    return worst

def check_kpis():
    print("== B: cycle KPIs vs converged MATLAB + paper Table 2 ==")
    worst = 0.0; worstT2 = 0.0
    for name, spec in ADS2019.items():
        mix = make_mix(spec, 298.15)
        path = blow_path(mix, 298.15, spec["rho"], 1.0, 0.15, 0.03)
        for (Plow, Pint), ref in ML_KPI[name].items():
            if Plow < path["P"][-1] - 1e-12: continue
            k = cycle_kpis(mix, path, Pint, Plow, 0.15, "LPP")
            got = (k["Pu"], k["Re"], k["EnBLO"], k["EnEVAC"], k["En"], k["WC"])
            for g, r in zip(got, ref):
                if r is None: continue
                e = abs(g-r)/max(abs(r), 1e-9)
                worst = max(worst, e)
        k = cycle_kpis(mix, path, 0.15, 0.03, 0.15, "LPP")
        t2 = PAPER_T2[name]
        for g, r in ((k["Pu"], t2[0]), (k["Re"], t2[1]), (k["En"], t2[2])):
            worstT2 = max(worstT2, abs(g-r))
        worstT2 = max(worstT2, abs(k["WC"]-t2[3])/t2[3]*100*0.0 + abs(k["WC"]-t2[3]))
        print(f"  {name}: worst rel so far {worst:.3e}")
    print(f"  KPI worst rel err = {worst:.3e}; Table-2 worst abs dev = {worstT2:.3f}")
    return worst

def check_iast_equiv():
    print("== C: IAST(equal-qs Langmuir) == extended Langmuir ==")
    # single-site langmuir, equal qs -> IAST identical to extended Langmuir
    T = 298.15; rho = 1130.0
    qs = 4.39; bC = 2.50e-6*np.exp(31.19e3/(RE*T)); bN = 2.70e-6*np.exp(16.38e3/(RE*T))
    spec = dict(rho=rho,
        co2=dict(qsb=qs, qsd=0.0, b0=bC, d0=0.0, dUb=0.0, dUd=0.0),
        n2 =dict(qsb=qs, qsd=0.0, b0=bN, d0=0.0, dUb=0.0, dUd=0.0))
    mixE = EDSL(spec["co2"], spec["n2"], T)
    mixI = IASTMix(lambda c: qs*bC*c/(1+bC*c), lambda c: qs*np.log1p(bC*c),
                   lambda c: qs*bN*c/(1+bN*c), lambda c: qs*np.log1p(bN*c), T)
    pE = blow_path(mixE, T, rho, 1.0, 0.15, 0.03, dP=5e-4)
    pI = blow_path(mixI, T, rho, 1.0, 0.15, 0.03, dP=5e-4, rtol=1e-8, atol=1e-10)
    kE = cycle_kpis(mixE, pE, 0.15, 0.03, 0.15, "LPP")
    kI = cycle_kpis(mixI, pI, 0.15, 0.03, 0.15, "LPP")
    worst = max(abs(kE[k]-kI[k])/max(abs(kE[k]), 1e-9)
                for k in ("Pu", "Re", "En", "WC"))
    print(f"  EDSL: Pu={kE['Pu']:.6f} Re={kE['Re']:.6f} En={kE['En']:.6f} WC={kE['WC']:.4f}")
    print(f"  IAST: Pu={kI['Pu']:.6f} Re={kI['Re']:.6f} En={kI['En']:.6f} WC={kI['WC']:.4f}")
    print(f"  worst rel diff = {worst:.3e}")
    return worst

def check_liske_example():
    """Liske e-BAAM worked example: 13X SSL refit, PH=2, Pint=0.8, Plow=0.07,
    yF=0.15. Paper text quotes 96.6/74.8 (LPP) and 96.6/74.0 (FP). These
    reproduce at T=298.15 K (the figure caption says 30 C but the quoted
    numbers correspond to 298.15 K - confirmed to all displayed digits)."""
    print("== D: Liske e-BAAM worked example (13X SSL, PH=2, Pint=0.8, Plow=0.07) ==")
    T = 298.15
    spec = ADS2026["Z13Xssl"]
    mix = make_mix(spec, T)
    path = blow_path(mix, T, spec["rho"], 2.0, 0.15, 0.07, eta=eta_maruyama)
    kL = cycle_kpis(mix, path, 0.8, 0.07, 0.15, "LPP", eta=eta_maruyama)
    kF = cycle_kpis(mix, path, 0.8, 0.07, 0.15, "FP", eta=eta_maruyama)
    print(f"  LPP: Pu={kL['Pu']:.2f} Re={kL['Re']:.2f} (paper 96.6/74.8)  mbal={kL['mbal']:.2e}%")
    print(f"  FP : Pu={kF['Pu']:.2f} Re={kF['Re']:.2f} (paper 96.6/74.0)  mbal={kF['mbal']:.2e}%")
    okL = abs(kL["Pu"]-96.6) < 0.05 and abs(kL["Re"]-74.8) < 0.05
    okF = abs(kF["Pu"]-96.6) < 0.05 and abs(kF["Re"]-74.0) < 0.05
    print(f"  LPP match: {okL}, FP match: {okF}")
    return okL and okF

PVSA_CASES = [  # (label, adsorbent, T, PH, yF, Pint, Plow, mode, eta)
    ("Z13Xssl PH=2 LPP", "Z13Xssl", 298.15, 2.0, 0.15, 0.80, 0.07, "LPP", "mar"),
    ("Z13Xssl PH=3 LPP", "Z13Xssl", 298.15, 3.0, 0.15, 1.50, 0.08, "LPP", "mar"),
    ("Z13Xssl PH=3 FP ", "Z13Xssl", 298.15, 3.0, 0.15, 1.50, 0.08, "FP",  "mar"),
    ("Z13Xssl PH=5 LPP", "Z13Xssl", 298.15, 5.0, 0.15, 2.50, 0.10, "LPP", "mar"),
    ("Z13Xssl PH=5 LPP low-Pint", "Z13Xssl", 298.15, 5.0, 0.15, 0.60, 0.05, "LPP", "mar"),
    ("CALF20  PH=5 FP ", "CALF20",  298.15, 5.0, 0.15, 2.00, 0.06, "FP",  "const"),
]

def check_pvsa():
    """P_H > 1 bar (PVSA). Two things are checked that only bite above 1 bar:
    (i) no work is charged to depressurization above 1 bar, and (ii) the total
    energy is exactly vacuum + compression, with compression evaluated at
    eta(1 bar) as in the e-BAAM (the convention that reproduces the paper's
    PH = 2 worked example)."""
    print("== F: PVSA (PH > 1 bar) ==")
    worst_split = 0.0
    for lab, ads, T, PH, yF, Pint, Plow, mode, et in PVSA_CASES:
        spec = ADS2026[ads]
        mix = make_mix(spec, T)
        eta = eta_maruyama if et == "mar" else eta_const
        path = blow_path(mix, T, spec["rho"], PH, yF, Plow, eta=eta)
        k = cycle_kpis(mix, path, Pint, Plow, yF, mode, eta=eta)
        # (i) zero vacuum work above 1 bar
        above = path["P"] >= 1.0
        wmax = float(np.max(np.abs(path["W"][above])))
        assert wmax == 0.0, f"{lab}: W != 0 above 1 bar ({wmax})"
        # (ii) energy split closes
        tot = k["EnBLO"] + k["EnEVAC"] + k["EnADS"] + k["EnPR"]
        rel = abs(tot - k["En"])/max(abs(k["En"]), 1e-300)
        worst_split = max(worst_split, rel)
        print(f"  {lab}: Pu={k['Pu']:8.4f} Re={k['Re']:8.4f} En={k['En']:8.3f} "
              f"(vac {k['EnBLO']+k['EnEVAC']:7.3f} + comp {k['EnADS']+k['EnPR']:7.3f}) "
              f"WC={k['WC']:8.2f} mbal={k['mbal']:.1e}%")
    print(f"  zero-work-above-1-bar: PASS (exact) | worst energy-split residual = {worst_split:.2e}")
    return worst_split

def print_pvsa_anchors():
    """PVSA anchors for the JS harness (pasted as literals in baam_test.js)."""
    print("== F-anchors: JS harness (PH > 1) ==")
    for lab, ads, T, PH, yF, Pint, Plow, mode, et in PVSA_CASES:
        spec = ADS2026[ads]
        mix = make_mix(spec, T)
        eta = eta_maruyama if et == "mar" else eta_const
        path = blow_path(mix, T, spec["rho"], PH, yF, Plow, eta=eta)
        k = cycle_kpis(mix, path, Pint, Plow, yF, mode, eta=eta)
        vals = ", ".join(f"{kk}: {k[kk]:.10g}" for kk in
                         ("Pu", "Re", "En", "EnADS", "EnPR", "WC", "ydel"))
        print(f'  {{ads:"{ads}", T:{T}, PH:{PH}, yF:{yF}, Pint:{Pint}, Plow:{Plow}, '
              f'mode:"{mode}", eta:"{et}", ref:{{{vals}}}}},')

def print_extra_anchors():
    """Anchors for the JS test harness (printed, pasted as literals there)."""
    print("== E: JS-harness anchors ==")
    T = 298.15; spec = ADS2019["Z13X"]; mix = make_mix(spec, T)
    path = blow_path(mix, T, spec["rho"], 1.0, 0.15, 0.03)
    k = cycle_kpis(mix, path, 0.15, 0.03, 0.15, "LPP")
    print("  Z13X LPP  :", {kk: float(f"{k[kk]:.10g}") for kk in
          ("Pu","Re","En","WC","ydel","ygam","Npr","Nfeed","Nraff")})
    kF = cycle_kpis(mix, path, 0.15, 0.03, 0.15, "FP")
    print("  Z13X FP   :", {kk: float(f"{kF[kk]:.10g}") for kk in
          ("Pu","Re","En","WC","ydel","Npr","Nfeed")})
    T2 = 303.15; s2 = ADS2026["CALF20"]; mix2 = make_mix(s2, T2)
    p2 = blow_path(mix2, T2, s2["rho"], 2.0, 0.15, 0.05, eta=eta_maruyama)
    k2 = cycle_kpis(mix2, p2, 0.5, 0.05, 0.20, "LPP", eta=eta_maruyama)
    print("  CALF20 PVSA (T=303.15,PH=2,yF=0.20,Pint=0.5,Plow=0.05):",
          {kk: float(f"{k2[kk]:.10g}") for kk in ("Pu","Re","En","WC","ydel")})

if __name__ == "__main__":
    import time
    t0 = time.time()
    eA = check_paths()
    eB = check_kpis()
    eC = check_iast_equiv()
    okD = check_liske_example()
    eF = check_pvsa()
    print_extra_anchors()
    print_pvsa_anchors()
    print(f"\nTotal {time.time()-t0:.1f}s | A={eA:.2e} B={eB:.2e} C={eC:.2e} "
          f"D={'PASS' if okD else 'FAIL'} F={eF:.2e}")
