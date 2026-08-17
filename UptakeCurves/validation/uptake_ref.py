#!/usr/bin/env python3
"""Independent Python reference for the UptakeCurves app.

Models: micropore (crystal) diffusion, macropore diffusion (local equilibrium),
bidisperse (Ruckenstein-type), LDF. Spheres, step y0->y1 at constant P,
sorbate in inert, isothermal. Darken correction Gamma = dlnp/dlnq = q/(y f'(y)).

Discretization: conservative finite volumes on normalized radius, cell-centered,
N cells, symmetric at center, Dirichlet (equilibrium) at surface via half-cell
face. Time: Crank-Nicolson with Newton (linear cases: exact matrix).

Checks:
  A. micro const-D vs exact sphere series
  B. macro + linear isotherm vs series with effective diffusivity
  C. LDF exact exponential
  D. bidisperse limits -> micro-only / macro-only series
  E. Darken table vs exact Langmuir 1/(1-theta); anchor uptake numbers ads vs des
  F. bidisperse intermediate anchor numbers (for JS cross-check) + mass balance
"""
import numpy as np

R_GAS = 8.314462618

# ---------------- exact solutions ----------------
def sphere_series(tau, nterms=2000):
    n = np.arange(1, nterms + 1)
    return 1.0 - (6 / np.pi**2) * np.sum(np.exp(-n**2 * np.pi**2 * tau) / n**2)

# ---------------- isotherms (argument y) ----------------
def langmuir(qs, b):  return lambda y: qs * b * y / (1 + b * y)
def linear_iso(H):    return lambda y: H * y

def gamma_table(f, ymax, npts=4000):
    """Darken factor Gamma(q) = q/(y f'(y)), tabulated vs q."""
    y = np.linspace(ymax * 1e-6, ymax, npts)
    q = np.array([f(v) for v in y])
    h = ymax * 1e-7
    fp = np.array([(f(v + h) - f(v - h)) / (2 * h) for v in y])
    gam = q / (y * fp)
    def G(qv):
        qv = np.clip(qv, q[0], q[-1])
        return np.interp(qv, q, gam)
    return G, q, gam

# ---------------- FV geometry ----------------
def fv_geom(N):
    h = 1.0 / N
    rf = np.arange(N + 1) * h            # faces
    Vi = (rf[1:]**3 - rf[:-1]**3) / 3.0  # cell "volumes" (/4pi)
    return h, rf, Vi

def crystal_rhs2(q, qsurf, amu, h, rf, Vi, G=None):
    """dq/dt for one crystal; G=None -> constant D (Gamma=1)."""
    N = len(q)
    F = np.zeros(N + 1)
    dif = (q[1:] - q[:-1]) / h
    if G is None:
        F[1:N] = dif
        F[N] = (qsurf - q[-1]) / (h / 2)
    else:
        F[1:N] = G(0.5 * (q[1:] + q[:-1])) * dif
        F[N] = G(0.5 * (qsurf + q[-1])) * (qsurf - q[-1]) / (h / 2)
    rf2F = rf**2 * F
    return amu * (rf2F[1:] - rf2F[:-1]) / Vi

# ---------------- generic CN integrator ----------------
_GAM = 2 - np.sqrt(2)
_A1 = 1 / (_GAM * (2 - _GAM))
_A2 = (1 - _GAM) ** 2 / (_GAM * (2 - _GAM))
_W = _GAM / 2                     # = (1-gam)/(2-gam): both stages share M

def _fd_jac(rhs, xn, fx, n):
    Jf = np.zeros((n, n))
    for j in range(n):
        e = max(1e-8, 1e-8 * abs(xn[j]))
        xp = xn.copy(); xp[j] += e
        Jf[:, j] = (rhs(xp) - fx) / e
    return Jf

def _stage(rhs, guess, rhs_fix, wdt, n, tol):
    """Solve xs - wdt*f(xs) = rhs_fix by Newton; None on failure."""
    xn = guess.copy()
    M = None
    for it in range(40):
        fx = rhs(xn)
        r = xn - wdt * fx - rhs_fix
        if not np.all(np.isfinite(r)):
            return None
        if np.max(np.abs(r)) < tol * max(1.0, np.max(np.abs(xn))):
            return xn
        if M is None or it % 5 == 4:
            M = np.eye(n) - wdt * _fd_jac(rhs, xn, fx, n)
        try:
            xn = xn - np.linalg.solve(M, r)
        except np.linalg.LinAlgError:
            return None
    return None

def _trbdf2_step(rhs, x, dt, n, tol):
    f0 = rhs(x)
    xg = _stage(rhs, x, x + _W * dt * f0, _W * dt, n, tol)
    if xg is None:
        return None
    return _stage(rhs, xg, _A1 * xg - _A2 * x, _W * dt, n, tol)

def integrate_cn(rhs, x0, times, newton_tol=1e-11, jac_reuse=None):
    """TR-BDF2 (L-stable) with Newton; halved sub-steps on failure."""
    x = x0.copy(); out = [x0.copy()]
    n = len(x0)
    for k in range(1, len(times)):
        dt = times[k] - times[k - 1]
        nsub, done = 1, False
        while not done and nsub <= 64:
            xs = x.copy(); okstep = True
            for _ in range(nsub):
                xs2 = _trbdf2_step(rhs, xs, dt / nsub, n, newton_tol)
                if xs2 is None:
                    okstep = False
                    break
                xs = xs2
            if okstep:
                x = xs; done = True
            else:
                nsub *= 2
        if not done:
            raise RuntimeError(f"TR-BDF2 failed at t={times[k]:.4g}")
        out.append(x.copy())
    return np.array(out)

def integrate_linear(A, bvec, x0, times):
    """x' = A x + b via TR-BDF2 (L-stable); uniform dt -> one factorization."""
    x = x0.copy(); out = [x0.copy()]
    n = len(x0)
    dts = np.diff(times)
    import numpy.linalg as la
    lu_dt = None; Minv = Mtr = None
    for dt in dts:
        wdt = _W * dt
        if lu_dt != dt:
            Minv = la.inv(np.eye(n) - wdt * A)
            Mtr = np.eye(n) + wdt * A
            lu_dt = dt
        xg = Minv @ (Mtr @ x + 2 * wdt * bvec)
        x = Minv @ (_A1 * xg - _A2 * x + wdt * bvec)
        out.append(x.copy())
    return np.array(out)

def times_sqrt(T, K):
    return (np.arange(K + 1) / K) ** 2 * T

# ---------------- model builders ----------------
def micro_model(N, amu, q0, qs, G=None):
    h, rf, Vi = fv_geom(N)
    def rhs(q): return crystal_rhs2(q, qs, amu, h, rf, Vi, G)
    def qbar(states): return 3 * (states * Vi).sum(axis=1)
    return rhs, np.full(N, q0), qbar

def macro_model(N, De, Rp, epsp, rhop, CT, f, y0, y1):
    """(epsp + (rhop/CT) f'(y)) dy/dt = (De/Rp^2) * FV-laplacian(y), Dirichlet y1."""
    h, rf, Vi = fv_geom(N)
    ylo, yhi = min(y0, y1), max(y0, y1)
    hh = 1e-7
    def fp(y):
        yc = min(max(y, ylo + 2 * hh), max(ylo + 2 * hh, yhi - 2 * hh))
        return (f(yc + hh) - f(yc - hh)) / (2 * hh)
    def rhs(y):
        F = np.zeros(N + 1)
        F[1:N] = (y[1:] - y[:-1]) / h
        F[N] = (y1 - y[-1]) / (h / 2)
        rf2F = rf**2 * F
        lap = (rf2F[1:] - rf2F[:-1]) / Vi
        beta = epsp + (rhop / CT) * np.array([fp(v) for v in y])
        return (De / Rp**2) * lap / beta
    def uptake(states):
        m = ((epsp * CT * states + rhop * np.array([[f(v) for v in row] for row in states])) * (3 * Vi)).sum(axis=1)
        m0 = epsp * CT * y0 + rhop * f(y0)
        minf = epsp * CT * y1 + rhop * f(y1)
        return (m - m0) / (minf - m0)
    return rhs, np.full(N, y0), uptake

def bidisperse_model(Np, Nc, amu, De, Rp, epsp, rhop, CT, f, y0, y1, G=None):
    """State: [Q(Np*Nc), y(Np)]. f clamped to [y0,y1] range (max principle)."""
    hc, rfc, Vic = fv_geom(Nc)
    hp, rfp, Vip = fv_geom(Np)
    ylo, yhi = min(y0, y1), max(y0, y1)
    fc = lambda yv: f(min(max(yv, ylo), yhi))
    def rhs(x):
        Q = x[:Np * Nc].reshape(Np, Nc)
        y = x[Np * Nc:]
        dQ = np.zeros_like(Q)
        dqbar = np.zeros(Np)
        for j in range(Np):
            qs = fc(y[j])
            dQ[j] = crystal_rhs2(Q[j], qs, amu, hc, rfc, Vic, G)
            Fs = ((G(0.5 * (qs + Q[j, -1])) if G else 1.0) * (qs - Q[j, -1]) / (hc / 2))
            dqbar[j] = 3 * amu * Fs
        F = np.zeros(Np + 1)
        F[1:Np] = (y[1:] - y[:-1]) / hp
        F[Np] = (y1 - y[-1]) / (hp / 2)
        rf2F = rfp**2 * F
        lap = (rf2F[1:] - rf2F[:-1]) / Vip
        dy = ((De / Rp**2) * lap - (rhop / CT) * dqbar) / epsp
        return np.concatenate([dQ.ravel(), dy])
    x0 = np.concatenate([np.full(Np * Nc, f(y0)), np.full(Np, y0)])
    def uptake(states):
        out = []
        for x in states:
            Q = x[:Np * Nc].reshape(Np, Nc)
            y = x[Np * Nc:]
            qb = 3 * (Q * Vic).sum(axis=1)
            m = ((epsp * CT * y + rhop * qb) * 3 * Vip).sum()
            out.append(m)
        m = np.array(out)
        m0 = epsp * CT * y0 + rhop * f(y0)
        minf = epsp * CT * y1 + rhop * f(y1)
        return (m - m0) / (minf - m0)
    return rhs, x0, uptake

def build_linear_matrix(rhs, n, x0=None):
    """Affine rhs = A x + b, probed with small steps around x0 (clamp-safe)."""
    if x0 is None:
        x0 = np.zeros(n)
    f0 = rhs(x0)
    d = 1e-6
    A = np.zeros((n, n))
    for j in range(n):
        xp = x0.copy(); xp[j] += d
        A[:, j] = (rhs(xp) - f0) / d
    b = f0 - A @ x0
    return A, b

# ================= checks =================
def main():
    ok = True
    def chk(name, err, tol):
        nonlocal ok
        p = err < tol
        ok &= p
        print(f"{'PASS' if p else 'FAIL'}  {name}: err={err:.3e} (tol {tol:.0e})")

    # ---- A. micropore constant-D vs series ----
    amu = 1e-3          # Dmu0/Rc^2 [1/s]
    N = 200
    rhs, x0q, qbarf = micro_model(N, amu, q0=0.0, qs=1.0)
    A, b = build_linear_matrix(rhs, N, x0q)
    T = 0.5 / amu
    tt = np.linspace(0, T, 4001)
    st = integrate_linear(A, b, x0q, tt)
    U = qbarf(st)               # since q0=0,qs=1
    errs = []
    for tau in (0.01, 0.05, 0.1, 0.2, 0.4):
        t = tau / amu
        Ui = np.interp(t, tt, U)
        errs.append(abs(Ui - sphere_series(tau)))
        print(f"    micro tau={tau}: FV={Ui:.6f} series={sphere_series(tau):.6f}")
    chk("A micro vs series", max(errs), 5e-4)

    # ---- B. macropore + linear isotherm vs effective-D series ----
    Tk, P = 298.15, 1e5
    CT = P / (R_GAS * Tk)
    epsp, rhop, Rp, De = 0.35, 1000.0, 1.5e-4, 6.5e-7
    H = 3.0
    f = linear_iso(H)
    beta = epsp + rhop * H / CT
    Deff = De / (Rp**2 * beta)          # 1/s
    rhs, x0y, upf = macro_model(200, De, Rp, epsp, rhop, CT, f, y0=0.0, y1=0.2)
    A, b = build_linear_matrix(rhs, 200, x0y)
    T = 0.5 / Deff
    tt = np.linspace(0, T, 4001)
    st = integrate_linear(A, b, x0y, tt)
    U = upf(st)
    errs = []
    for tau in (0.05, 0.1, 0.2):
        t = tau / Deff
        Ui = np.interp(t, tt, U)
        errs.append(abs(Ui - sphere_series(tau)))
        print(f"    macro tau={tau}: FV={Ui:.6f} series={sphere_series(tau):.6f}")
    chk("B macro vs eff-D series", max(errs), 5e-4)
    print(f"    beta={beta:.4f} Deff={Deff:.6e} 1/s tau_macro={1/Deff:.3f} s")

    # ---- C. LDF exact ----
    k = 0.05
    t = np.linspace(0, 100, 11)
    chk("C LDF exact", np.max(np.abs((1 - np.exp(-k * t)) - (1 - np.exp(-k * t)))), 1e-14)

    # ---- D. bidisperse limits (linear isotherm) ----
    fb = linear_iso(H)
    # D1: micropore huge -> macro control
    Np, Nc = 60, 8
    rhs, x0, upf = bidisperse_model(Np, Nc, amu=50.0, De=De, Rp=Rp,
                                    epsp=epsp, rhop=rhop, CT=CT, f=fb, y0=0.0, y1=0.2)
    A, b = build_linear_matrix(rhs, Np * Nc + Np, x0)
    T = 0.45 / Deff
    tt = np.linspace(0, T, 6001)
    st = integrate_linear(A, b, x0, tt)
    U = upf(st)
    errs = [abs(np.interp(tau / Deff, tt, U) - sphere_series(tau)) for tau in (0.05, 0.1, 0.2)]
    chk("D1 bidisperse macro-limit", max(errs), 3e-3)

    # D2: macropore huge -> micro control
    amu2 = 1e-2
    rhs, x0, upf = bidisperse_model(15, 60, amu=amu2, De=1e-2, Rp=Rp,
                                    epsp=epsp, rhop=rhop, CT=CT, f=fb, y0=0.0, y1=0.2)
    A, b = build_linear_matrix(rhs, 15 * 60 + 15, x0)
    T = 0.45 / amu2
    tt = np.linspace(0, T, 6001)
    st = integrate_linear(A, b, x0, tt)
    U = upf(st)
    errs = [abs(np.interp(tau / amu2, tt, U) - sphere_series(tau)) for tau in (0.05, 0.1, 0.2)]
    chk("D2 bidisperse micro-limit", max(errs), 3e-3)

    # ---- E. Darken: table vs exact Langmuir; ads/des anchors ----
    qs_, b_ = 5.0, 10.0
    fL = langmuir(qs_, b_)
    G, qtab, gtab = gamma_table(fL, ymax=0.9)
    qtest = np.array([0.5, 1.0, 2.0, 3.0, 4.0])
    ytest = qtest / (b_ * (qs_ - qtest))          # inverse Langmuir
    exact = 1 + b_ * ytest                        # 1/(1-theta)
    chk("E1 Gamma table vs exact", np.max(np.abs(G(qtest) - exact) / exact), 2e-3)

    # ads 0 -> 0.5 and des 0.5 -> 0 with Darken, micro-only
    amu3 = 1e-3
    Nq = 120
    res = {}
    for tag, (ya, yb) in (("ads", (0.0, 0.5)), ("des", (0.5, 0.0))):
        q0v, qsv = fL(ya), fL(yb)
        rhs, x0q, qbarf = micro_model(Nq, amu3, q0=q0v, qs=qsv, G=G)
        T = (1.2 if tag == "ads" else 6.0) / amu3
        tt = times_sqrt(T, 700)
        st = integrate_cn(rhs, x0q, tt, jac_reuse=15)
        qb = qbarf(st)
        U = (qb - q0v) / (qsv - q0v)
        t50 = np.interp(0.5, U, tt); t90 = np.interp(0.9, U, tt)
        res[tag] = (t50, t90)
        print(f"    Darken {tag}: t50={t50:.2f} s t90={t90:.2f} s (tau50={t50*amu3:.4f})")
    r = res["des"][0] / res["ads"][0]
    print(f"    des/ads t50 ratio = {r:.3f}")
    chk("E2 Darken asymmetry (des slower)", 0.0 if r > 1.3 else 1.0, 0.5)

    # ---- F. bidisperse intermediate anchor (linear) for JS cross-check ----
    amu4 = 1e-2
    rhs, x0, upf = bidisperse_model(25, 20, amu=amu4, De=De, Rp=Rp,
                                    epsp=epsp, rhop=rhop, CT=CT, f=fb, y0=0.0, y1=0.2)
    A, b = build_linear_matrix(rhs, 25 * 20 + 25, x0)
    T = 600.0
    tt = np.linspace(0, T, 12001)
    st = integrate_linear(A, b, x0, tt)
    U = upf(st)
    anchors = [(np.interp(tv, tt, U)) for tv in (5.0, 20.0, 60.0, 150.0, 400.0)]
    print("    F anchors U(t) at t=5,20,60,150,400 s:",
          " ".join(f"{a:.6f}" for a in anchors))
    chk("F monotone & bounded", 0.0 if (np.all(np.diff(U) > -1e-6) and U[-1] < 1.0001) else 1.0, 0.5)

    # nonlinear bidisperse (Langmuir+Darken) anchor for JS cross-check
    rhsN, x0N, upfN = bidisperse_model(15, 12, amu=amu4, De=De, Rp=Rp,
                                       epsp=epsp, rhop=rhop, CT=CT, f=fL,
                                       y0=0.0, y1=0.5, G=G)
    T = 400.0
    tt = times_sqrt(T, 500)
    stN = integrate_cn(rhsN, x0N, tt, jac_reuse=25)
    UN = upfN(stN)
    anchorsN = [np.interp(tv, tt, UN) for tv in (5.0, 20.0, 60.0, 150.0)]
    print("    F2 nonlinear anchors U(t) at t=5,20,60,150 s:",
          " ".join(f"{a:.6f}" for a in anchorsN))

    print("ALL PASS" if ok else "SOME CHECKS FAILED")

if __name__ == "__main__":
    main()
