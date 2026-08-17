#!/usr/bin/env python3
"""Independent Python/scipy reference for the DAC 0-D cycle engine.

Same model as dac_engine.js (Glaser-style 0-D TVSA with modified-Toth /
GAB equilibria and Simplified-DAC energy forms), but an independent
implementation: scipy BDF at tight tolerances, with the flow/power
quadratures carried as augmented ODE states (exact, no trapezoid error).

State y = [N, yc, yw, qc, qw, T, Tw] + 9 quadrature states
          [feedC, feedW, feedN, outC, outW, outN, Wfan, Wvac, Qheat]
('out' means the pump/vent outflow of the current step).
"""
import json
import numpy as np
from scipy.integrate import solve_ivp

R = 8.314
GAM = 1.4
MW_C = 0.04401
CPG = {"c": 37.1, "w": 33.6, "n": 29.1}
CPA_C, CPA_W = 50.0, 75.0
Y_CO2_AIR = 400e-6

PRESETS = {
    "TMCM-41": dict(toth=dict(qs0=2.7, chi=0, T0=298, b0=0.38776, dH0=90000, w0=0.3123, alpha=0.3472),
                    gab=None, wet=None, dHc=90000, dHw=0,
                    rho_p=930, eb=0.40, ep=0.40, cps=1000, kC=2e-4, kW=0.2),
    "TPMS": dict(toth=dict(qs0=1.0, chi=0, T0=298, b0=2.0123, dH0=90000, w0=0.2696, alpha=0.08069),
                 gab=None, wet=None, dHc=90000, dHw=0,
                 rho_p=1250, eb=0.40, ep=0.40, cps=900, kC=2e-4, kW=0.2),
    "NbOFFIVE": dict(toth=dict(qs0=2.22, chi=0, T0=273, b0=0.17567, dH0=50000, w0=1.166, alpha=-0.4937),
                     gab=dict(Cm=6.159, CG=36.81, K=0.4063),
                     wet=dict(psi=0, beta=-0.11707, xmax=0.75), dHc=50000, dHw=45036,
                     rho_p=1173.6, eb=0.40, ep=0.40, cps=1000, kC=2e-4, kW=0.2),
    "SIFSIX": dict(toth=dict(qs0=3.0, chi=0, T0=273, b0=0.04572, dH0=52000, w0=0.5455, alpha=1.592),
                   gab=dict(Cm=0.4067, CG=93.13, K=0.7859),
                   wet=dict(psi=0, beta=-0.6772, xmax=0.90), dHc=52000, dHw=45036,
                   rho_p=786, eb=0.40, ep=0.40, cps=1000, kC=2e-4, kW=0.2),
    "APDES": dict(toth=dict(qs0=2.2, chi=0, T0=296, b0=0.373, dH0=60000, w0=0.4247, alpha=-0.4921),
                  gab=dict(Cm=36.48, CG=0.1489, K=0.5751),
                  wet=dict(psi=0.00958, beta=3.448, xmax=1.0), dHc=60000, dHw=43800,
                  rho_p=61, eb=0.092, ep=0.9614, cps=2070, kC=2e-4, kW=0.2),
    "Lewatit": dict(toth=dict(qs0=3.93873034, chi=0, T0=303.15, b0=20112.4573/1e5, dH0=70999,
                              w0=0.26906407, alpha=0.04501535),
                    gab=dict(Cm=3.8640, CG=1.7950, K=0.76830),
                    wet=dict(psi=4.70e-3, beta=0.569, xmax=1.0), dHc=70999, dHw=43800,
                    rho_p=760, eb=0.37, ep=0.35, cps=1580, kC=4e-3, kW=0.2),
}


def psat(T):
    return 10.0 ** (4.6543 - 1435.264 / (T - 64.848)) * 1e5  # Pa


def qstarW(x, S):
    g = S["gab"]
    if g is None:
        return 0.0
    a = min(max(x, 0.0), 0.999)
    Ka = g["K"] * a
    return g["Cm"] * g["CG"] * Ka / ((1 - Ka) * (1 + (g["CG"] - 1) * Ka))


def qstarC(p, T, qw, S, pw=0.0):
    """Modified-Toth CO2 loading. Water coupling is EITHER S['wet']
    (Stampi-Bombelli modulation by the water loading qw) OR S['cmp']
    (extended/competitive Toth: the water partial pressure pw [Pa] joins
    the shared denominator, b_w Arrhenius on the Toth T0)."""
    t = S["toth"]
    qs = t["qs0"] * np.exp(t["chi"] * (1 - T / t["T0"]))
    b = t["b0"] * np.exp((t["dH0"] / (R * t["T0"])) * (t["T0"] / T - 1))
    w = t["w0"] + t["alpha"] * (1 - t["T0"] / T)
    wet = S["wet"]
    if wet is not None and qw > 0:
        qwE = min(qw, wet["qwClamp"]) if wet.get("qwClamp") is not None else qw
        den = 1 - wet["psi"] * qwE
        qs = qs / den if den > 0.05 else qs / 0.05
        b = max(b * (1 + wet["beta"] * qwE), 0.0)
    f = b * max(p, 0.0)
    if f <= 0:
        return 0.0
    g = f
    cmp_ = S.get("cmp")
    if cmp_ is not None and cmp_["bw0"] > 0:
        bw = cmp_["bw0"] * np.exp((cmp_["dHw"] / (R * t["T0"])) * (t["T0"] / T - 1))
        g = f + bw * max(pw, 0.0)
    return qs * f / (1 + g ** w) ** (1 / w)


class Model:
    def __init__(self, u):
        # "Custom" carries its full spec in u["iso"], same shape as a preset
        src = u["iso"] if (u["sorbent"] == "Custom" and u.get("iso")) \
            else PRESETS[u["sorbent"]]
        S = {k: (dict(v) if isinstance(v, dict) else v) for k, v in src.items()}
        self.u = u

        def pick(key, default):
            v = u.get(key)
            return v if v is not None else default

        # sorbent properties may be overridden by the user (defaults = preset)
        eb, ep = pick("eb", S["eb"]), pick("ep", S["ep"])
        self.eb = eb
        self.epsT = eb + (1 - eb) * ep
        self.ms = pick("rhop", S["rho_p"]) * (1 - eb)
        self.VC = self.epsT
        self.Vsolid = 1 - eb
        self.cps = pick("cps", S["cps"])
        self.dHc = pick("dHc", S["dHc"])
        self.dHw = pick("dHw", S["dHw"]) if S["gab"] else 0.0
        # the Toth temperature coefficient IS the heat of adsorption unless
        # the user unlinks them (see engine comment)
        toth = dict(S["toth"])
        if not u.get("unlinkDH") and self.dHc != S["dHc"]:
            toth["dH0"] = self.dHc
        wet = S["wet"]
        if wet is not None:
            wet = dict(wet)
            wet["qwClamp"] = qstarW(wet["xmax"], S) if (wet["xmax"] < 1 and S["gab"]) else None
        cmp_ = S.get("cmp")
        if cmp_ is not None and cmp_["bw0"] > 0:
            cmp_ = dict(bw0=cmp_["bw0"], dHw=pick("dHw", S.get("dHw") or 0.0))
        else:
            cmp_ = None
        self.S = {"toth": toth, "gab": S["gab"], "wet": wet, "cmp": cmp_}
        self.inertW = S["gab"] is None
        pw = (u["RH"] / 100) * psat(u["Tamb"])
        ywA = pw / u["Pamb"]
        ycA = Y_CO2_AIR * (1 - ywA)
        self.yAir = dict(c=ycA, w=ywA, n=1 - ycA - ywA)
        self.ndotFan = u["Pamb"] * (u["Vfan"] / 3600) / (R * u["Tamb"])
        self.cpaC = u.get("cpaC", CPA_C)
        self.cpaW = u.get("cpaW", CPA_W)
        # structure carried with the sorbent: (m_str/m_sorb)*cp_str [J/kg/K]
        self.cpStr = u.get("mstr", 0.0) * u.get("cpstr", 0.0)


TREF_LDF = 293.15  # [K] reference T for the LDF coefficients (20 C)


def k_ldf(kref, Ea, T):
    """Arrhenius LDF: k(T) = k_ref exp[-Ea/R (1/T - 1/Tref)], Ea in J/mol."""
    if Ea:
        return kref * np.exp(-(Ea / R) * (1.0 / T - 1.0 / TREF_LDF))
    return kref


def rhs(M, mode, y):
    u = M.u
    N, yc, yw, qc, qw, T, Tw = y[:7]
    yn = 1 - yc - yw
    P = N * R * T / M.VC

    xRH = max(yw, 0.0) * P / psat(T)
    qwS = 0.0 if M.inertW else qstarW(xRH, M.S)
    qcS = qstarC(max(yc, 0.0) * P, T, max(qw, 0.0), M.S, max(yw, 0.0) * P)
    rc = k_ldf(u["kC"], u.get("EaC", 0.0), T) * (qcS - qc) * M.ms
    rw = 0.0 if M.inertW else k_ldf(u["kW"], u.get("EaW", 0.0), T) * (qwS - qw) * M.ms
    rSum = rc + rw

    Cgas = N * (yc * CPG["c"] + yw * CPG["w"] + yn * CPG["n"])
    Ctot = M.ms * (M.cps + M.cpStr) + Cgas + M.ms * (qc * M.cpaC + qw * M.cpaW)
    Qext = u["Ua"] * (Tw - T)
    Hads = M.dHc * rc + M.dHw * rw

    dTw = (u["Tdes"] - Tw) / u["tauHeat"] if mode == "heat" else (u["Tamb"] - Tw) / u["tauCool"]
    ndotI = ndotO = 0.0

    if mode == "ads":
        ndotI = M.ndotFan
        cpgI = CPG["c"] * M.yAir["c"] + CPG["w"] * M.yAir["w"] + CPG["n"] * M.yAir["n"]
        dT = (ndotI * cpgI * (u["Tamb"] - T) + Hads + Qext) / Ctot
        dN = -(N / T) * dT
        ndotO = ndotI - rSum - dN
        dyc = (ndotI * (M.yAir["c"] - yc) - rc + yc * rSum) / N
        dyw = (ndotI * (M.yAir["w"] - yw) - rw + yw * rSum) / N
    elif mode == "cool":
        cpgI = CPG["c"] * M.yAir["c"] + CPG["w"] * M.yAir["w"] + CPG["n"] * M.yAir["n"]
        A = Ctot + (N / T) * cpgI * (u["Tamb"] - T)
        dT = (rSum * cpgI * (u["Tamb"] - T) + Hads + Qext) / A
        dN = -(N / T) * dT
        ndotI = rSum + dN
        if ndotI < 0:
            ndotI = 0.0
            dT = (Hads + Qext) / Ctot
            dN = -(N / T) * dT
            ndotO = max(-(dN + rSum), 0.0)
        dyc = (ndotI * (M.yAir["c"] - yc) - rc + yc * rSum) / N
        dyw = (ndotI * (M.yAir["w"] - yw) - rw + yw * rSum) / N
    elif mode == "evac":
        ndotO = (u["Spump"] / 3600) * N / M.VC
        dN = -ndotO - rSum
        dT = (Hads + Qext + R * T * dN) / (Ctot - R * N)
        dyc = (-rc + yc * rSum) / N
        dyw = (-rw + yw * rSum) / N
    elif mode == "heat":
        dT = (Hads + Qext) / Ctot
        dN = -(N / T) * dT
        ndotO = max(-rSum - dN, 0.0)
        dyc = (-rc + yc * rSum) / N
        dyw = (-rw + yw * rSum) / N
    else:
        raise ValueError(mode)

    # quadrature states
    # isentropic blower over dPfan (same form as dac.energyModel)  [W]
    Wfan = ((GAM / (GAM - 1)) * (u["Pamb"] / u["etaFan"]) *
            (((u["Pamb"] + u["dPfan"]) / u["Pamb"]) ** ((GAM - 1) / GAM) - 1) *
            (u["Vfan"] / 3600)) if mode == "ads" else 0.0
    Wvac = 0.0
    if mode in ("evac", "heat") and ndotO > 0 and P < u["Pamb"]:
        if u["pumpModel"] == "adia":
            Wvac = (1 / u["etaVac"]) * (GAM / (GAM - 1)) * ndotO * R * T * \
                   ((u["Pamb"] / P) ** ((GAM - 1) / GAM) - 1)
        else:
            Wvac = ndotO * R * u["Tamb"] * np.log(u["Pamb"] / P) / u["etaVac"]
    dQheat = max(Qext, 0.0) if mode == "heat" else 0.0

    return [dN, dyc, dyw, rc / M.ms, (0.0 if M.inertW else rw / M.ms), dT, dTw,
            ndotI * M.yAir["c"], ndotI * M.yAir["w"], ndotI * M.yAir["n"],
            ndotO * yc, ndotO * yw, ndotO * yn,
            Wfan, Wvac, dQheat]


def integrate(M, mode, x0, tspan, event_pregen=False):
    y0 = list(x0) + [0.0] * 9
    ev = None
    if event_pregen:
        def ev(t, y):
            return y[0] * R * y[5] / M.VC - M.u["Pregen"]
        ev.terminal = True
        ev.direction = -1
    sol = solve_ivp(lambda t, y: rhs(M, mode, y), tspan, y0, method="BDF",
                    rtol=1e-10, atol=[1e-12, 1e-14, 1e-13, 1e-12, 1e-12, 1e-9, 1e-9]
                    + [1e-12] * 9, events=[ev] if ev else None, dense_output=False)
    if not sol.success:
        raise RuntimeError(f"{mode}: {sol.message}")
    yend = sol.y[:, -1]
    tend = sol.t[-1]
    stopped = bool(event_pregen and sol.t_events[0].size > 0)
    return yend[:7], yend[7:], tend, stopped


def run_cycle(M, x0):
    u = M.u
    Q = dict(feedC=0, feedW=0, feedN=0, prodC=0, prodW=0, prodN=0,
             ventC=0, ventW=0, ventN=0, Wfan=0, Wvac=0, Qheat=0, Qwall=0)

    def addQ(q, outkey):
        Q["feedC"] += q[0]; Q["feedW"] += q[1]; Q["feedN"] += q[2]
        Q[outkey + "C"] += q[3]; Q[outkey + "W"] += q[4]; Q[outkey + "N"] += q[5]
        Q["Wfan"] += q[6]; Q["Wvac"] += q[7]; Q["Qheat"] += q[8]

    x, q1, t1, _ = integrate(M, "ads", x0, (0, u["t_ads"]))
    addQ(q1, "vent")
    x_ads_end = list(x)
    tguard = max(4 * 3600, u["t_ads"])
    x, q2, t2, stopped = integrate(M, "evac", x, (0, tguard), event_pregen=True)
    addQ(q2, "vent")
    x = list(x)
    Tw_before_heat = x[6]
    x[0] = u["Pregen"] * M.VC / (R * x[5])
    x3_0 = x
    x, q3, t3, _ = integrate(M, "heat", x3_0, (0, u["t_heat"]))
    addQ(q3, "prod")
    marks = dict(qcAds=x_ads_end[3], qwAds=x_ads_end[4], Tads=x_ads_end[5],
                 qcDes=x[3], qwDes=x[4], Tdes=x[5])
    Q["Qwall"] += u["wallCp"] * max(x[6] - Tw_before_heat, 0.0)
    # instantaneous repressurisation
    Nadd = max((u["Pamb"] - u["Pregen"]) * M.VC / (R * x[5]), 0.0)
    N3 = x[0]
    Nn = N3 + Nadd
    yc4 = (N3 * x[1] + Nadd * M.yAir["c"]) / Nn
    yw4 = (N3 * x[2] + Nadd * M.yAir["w"]) / Nn
    Q["feedC"] += Nadd * M.yAir["c"]; Q["feedW"] += Nadd * M.yAir["w"]; Q["feedN"] += Nadd * M.yAir["n"]
    x4 = [Nn, yc4, yw4, x[3], x[4], x[5], x[6]]
    x, q4, t4, _ = integrate(M, "cool", x4, (0, u["t_cool"]))
    addQ(q4, "vent")
    t_cycle = t1 + t2 + t3 + t4
    return (list(x), Q, t_cycle, dict(ads=t1, evac=t2, heat=t3, cool=t4),
            stopped, marks)


def initial_state(M):
    u = M.u
    N0 = u["Pamb"] * M.VC / (R * u["Tamb"])
    qw0 = 0.0 if M.inertW else qstarW(u["RH"] / 100, M.S)
    qc0 = qstarC(M.yAir["c"] * u["Pamb"], u["Tamb"], qw0, M.S,
                 M.yAir["w"] * u["Pamb"])
    return [N0, M.yAir["c"], M.yAir["w"], qc0, qw0, u["Tamb"], u["Tamb"]]


def run_to_css(u, tolCSS=1e-5, maxCyc=60):
    M = Model(u)
    x = initial_state(M)
    scale = [1, 1e-4, 1e-2, 1, 1, 100, 100]
    ncyc, delta = 0, np.inf
    for _ in range(maxCyc):
        x1, Q, tc, times, stopped, _mk = run_cycle(M, x)
        ncyc += 1
        delta = max(abs(a - b) / (abs(b) + s * 1e-3) for a, b, s in zip(x1, x, scale))
        x = x1
        if delta < tolCSS:
            break
    x1, Q, tc, times, stopped, marks = run_cycle(M, x)   # KPI cycle at CSS
    kpi = compute_kpis(M, Q, tc, times, stopped)
    kpi["simp"] = simplified_energy(M, Q, marks)
    mb = mole_balance(M, x, x1, Q)
    return dict(M=M, cycles=ncyc, delta=delta, state=x1, Q=Q, kpi=kpi, mb=mb)


def compute_kpis(M, Q, tc, times, stopped):
    u = M.u
    mCO2 = Q["prodC"] * MW_C
    puDry = Q["prodC"] / max(Q["prodC"] + Q["prodN"], 1e-30)
    re = Q["prodC"] / max(Q["feedC"], 1e-30)
    Eth = (Q["Qheat"] + Q["Qwall"]) / max(mCO2, 1e-30) / 1e6
    Efan = Q["Wfan"] / max(mCO2, 1e-30) / 1e6
    Evac = Q["Wvac"] / max(mCO2, 1e-30) / 1e6
    prC = mCO2 / tc * 3600
    return dict(purityDry=puDry * 100, recovery=re * 100, Eth=Eth, Efan=Efan,
                Evac=Evac, prContact=prC, prSorb=prC / (1 - M.eb),
                waterPerCO2=Q["prodW"] / max(Q["prodC"], 1e-30),
                mCO2cycle=mCO2, t_cycle=tc, t_evac=times["evac"], evacReached=stopped)


def simplified_energy(M, Q, mk):
    """dac.energyModel forms evaluated with this cycle's own working
    capacity, temperature swing and capture fraction (see engine comment)."""
    u = M.u
    dqc = mk["qcAds"] - mk["qcDes"]
    dqw = mk["qwAds"] - mk["qwDes"]
    dT = mk["Tdes"] - mk["Tads"]
    if not (dqc > 1e-9 and dT > 1e-9):
        return None
    thermalMass = M.cps + M.cpStr + dqc * M.cpaC + dqw * M.cpaW
    Esens = thermalMass * dT / dqc
    if u["wallCp"] > 0:
        Esens += u["wallCp"] * dT / (M.ms * dqc)
    Ec = M.dHc
    Ew = (dqw / dqc) * M.dHw
    wAir = (GAM / (GAM - 1)) * (u["Pamb"] / u["etaFan"]) * \
        (((u["Pamb"] + u["dPfan"]) / u["Pamb"]) ** ((GAM - 1) / GAM) - 1)
    cCO2air = M.yAir["c"] * u["Pamb"] / (R * u["Tamb"])
    nBlown = M.ndotFan * M.yAir["c"] * u["t_ads"]
    etaCap = Q["prodC"] / max(nBlown, 1e-30)
    Eblow = wAir / (cCO2air * max(etaCap, 1e-12))
    Evac = (R * u["Tamb"] / u["etaVac"]) * float(np.log(u["Pamb"] / u["Pregen"])) * (1 + dqw / dqc)
    conv = 1 / (1e6 * MW_C)
    return dict(Esens=Esens * conv, EdesC=Ec * conv, EdesW=Ew * conv,
                Eth=(Esens + Ec + Ew) * conv, Efan=Eblow * conv, Evac=Evac * conv,
                Eel=(Eblow + Evac) * conv,
                Etot=(Esens + Ec + Ew + Eblow + Evac) * conv,
                dqC=dqc, dqW=dqw, dT=dT, etaCap=etaCap * 100)


def mole_balance(M, x0, x1, Q):
    invC0 = x0[0] * x0[1] + x0[3] * M.ms; invC1 = x1[0] * x1[1] + x1[3] * M.ms
    invW0 = x0[0] * x0[2] + x0[4] * M.ms; invW1 = x1[0] * x1[2] + x1[4] * M.ms
    nb0 = x0[0] * (1 - x0[1] - x0[2]); nb1 = x1[0] * (1 - x1[1] - x1[2])
    return dict(
        C=(Q["feedC"] - Q["ventC"] - Q["prodC"] - (invC1 - invC0)) / max(Q["feedC"], 1e-30),
        W=(Q["feedW"] - Q["ventW"] - Q["prodW"] - (invW1 - invW0)) / max(Q["feedW"], 1e-30),
        N=(Q["feedN"] - Q["ventN"] - Q["prodN"] - (nb1 - nb0)) / max(Q["feedN"], 1e-30))


BASE = dict(sorbent="Lewatit", Tamb=293.15, RH=50, Tdes=393.15,
            Pregen=0.20e5, Pamb=101325.0,
            t_ads=3 * 3600, t_heat=1.5 * 3600, t_cool=0.75 * 3600,
            Vfan=20000.0, Spump=50.0, dPfan=1013.0, etaFan=0.70, etaVac=0.50,
            pumpModel="iso", tauHeat=300.0, tauCool=300.0, Ua=20000.0,
            wallCp=0.0, kC=4e-3, kW=0.2)

ANCHOR_CASES = [
    ("Lewatit-base", dict(BASE)),
    ("Lewatit-adia", dict(BASE, pumpModel="adia")),
    ("APDES-base", dict(BASE, sorbent="APDES", kC=2e-4, Tdes=373.15, Pregen=0.10e5)),
    ("TMCM-dry", dict(BASE, sorbent="TMCM-41", kC=2e-4, RH=30.0)),
    ("SIFSIX-hum", dict(BASE, sorbent="SIFSIX", kC=2e-4, RH=70.0, Tdes=373.15)),
    ("NbOFFIVE-clamp", dict(BASE, sorbent="NbOFFIVE", kC=2e-4, RH=85.0)),
    ("Lewatit-hotwet", dict(BASE, Tamb=308.15, RH=80.0, Tdes=403.15, Pregen=0.05e5,
                            Spump=200.0, wallCp=2e5)),
    ("Lewatit-struct", dict(BASE, mstr=2.0, cpstr=710.0)),
    # user-edited sorbent properties: dHc drives the Toth dH0 as well
    ("Lewatit-props", dict(BASE, eb=0.30, ep=0.45, rhop=900.0, cps=1200.0,
                           dHc=80000.0, dHw=48000.0)),
    # same dHc, but unlinked from the isotherm (energy balance only)
    ("Lewatit-unlinkDH", dict(BASE, dHc=80000.0, unlinkDH=True)),
    # Arrhenius LDF temperature dependence (k sliders = value at 20 C):
    # wet sorbent, both coefficients T-dependent
    ("Lewatit-Ea", dict(BASE, EaC=30000.0, EaW=20000.0)),
    # dry sorbent (H2O inert), CO2 coefficient only
    ("APDES-Ea", dict(BASE, sorbent="APDES", kC=2e-4, Tdes=373.15,
                      Pregen=0.10e5, EaC=45000.0)),
]


def custom_iso(src, **over):
    """Deep-copy a preset as a Custom iso spec, with overrides."""
    S = {k: (dict(v) if isinstance(v, dict) else v)
         for k, v in PRESETS[src].items()}
    S.update(over)
    return S


ANCHOR_CASES += [
    # competitive (shared-denominator) extended Toth, GAB water:
    # Lewatit numbers with the psi/beta modulation replaced by b_w
    ("Custom-cmp-wet", dict(BASE, sorbent="Custom",
                            iso=custom_iso("Lewatit", wet=None,
                                           cmp=dict(bw0=2e-3)))),
    # competition with INERT water: pw still acts through the gas phase;
    # dHw passed as user input so b_w stays T-dependent
    ("Custom-cmp-dry", dict(BASE, sorbent="Custom", kC=2e-4, Tdes=373.15,
                            Pregen=0.10e5, dHw=43800.0,
                            iso=custom_iso("APDES", cmp=dict(bw0=2e-3)))),
]


def main():
    out = {}
    for name, u in ANCHOR_CASES:
        r = run_to_css(u)
        k = r["kpi"]
        out[name] = dict(u=u, cycles=r["cycles"],
                         purityDry=k["purityDry"], recovery=k["recovery"],
                         Eth=k["Eth"], Efan=k["Efan"], Evac=k["Evac"],
                         prContact=k["prContact"], waterPerCO2=k["waterPerCO2"],
                         mCO2cycle=k["mCO2cycle"], t_evac=k["t_evac"],
                         t_cycle=k["t_cycle"],
                         simp={kk: float(vv) for kk, vv in k["simp"].items()}
                         if k["simp"] else None,
                         state=[float(v) for v in r["state"]])
        mb = r["mb"]
        print(f"{name:16s} cyc={r['cycles']:2d} Pu={k['purityDry']:8.4f} Re={k['recovery']:8.4f} "
              f"Eth={k['Eth']:9.4f} Efan={k['Efan']:7.4f} Evac={k['Evac']:8.5f} "
              f"Pr={k['prContact']:8.5f} w/c={k['waterPerCO2']:7.4f} "
              f"mb={max(abs(mb['C']), abs(mb['W']) if not np.isnan(mb['W']) else 0, abs(mb['N'])):.2e}")
    with open("dac_anchors.json", "w") as f:
        json.dump(out, f, indent=1)
    print("anchors written to dac_anchors.json")


if __name__ == "__main__":
    main()
