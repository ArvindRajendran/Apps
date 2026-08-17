#!/usr/bin/env python3
"""Reference for the Isotherm Fitter app, using scipy.optimize.least_squares.

Generates seeded synthetic datasets (printed as JS literals for the jsc
harness), fits with scipy (trf, native bounds), prints recovered parameters
and weighted-SSE/AICc anchors. The JS engine (own LM, logit-bounded
transform, multi-start) must find the same minima.

Checks:
  A. zero-noise recovery: Langmuir, Sips, Toth exact to ~1e-8
  B. noisy Langmuir + DSL fits (seed 7): anchor params + AICc ranking
     (DSL must NOT beat Langmuir on Langmuir data after AICc penalty)
  C. weighting matters: relative vs absolute give different optima (anchor both)
  D. global van 't Hoff fit, 3 temperatures, DSL: recover qs, b0, U
  E. AICc formula spot value
"""
import numpy as np
from scipy.optimize import least_squares

R = 8.314462618

# ---------------- models (must match JS exactly) ----------------
def langmuir(th, p, T=None):
    qs, b = th
    return qs * b * p / (1 + b * p)

def dsl(th, p, T=None):
    qs1, b1, qs2, b2 = th
    return qs1 * b1 * p / (1 + b1 * p) + qs2 * b2 * p / (1 + b2 * p)

def sips(th, p, T=None):
    qs, b, n = th
    t = np.power(np.maximum(b * p, 0), n)
    return qs * t / (1 + t)

def toth(th, p, T=None):
    qs, b, t = th
    return qs * b * p / np.power(1 + np.power(np.maximum(b * p, 0), t), 1 / t)

def dsl_vh(th, p, T):
    # global van 't Hoff: qs1, b01, U1[J/mol], qs2, b02, U2
    qs1, b01, U1, qs2, b02, U2 = th
    b1 = b01 * np.exp(U1 / (R * T))
    b2 = b02 * np.exp(U2 / (R * T))
    return qs1 * b1 * p / (1 + b1 * p) + qs2 * b2 * p / (1 + b2 * p)

def weights(q, mode):
    qf = 1e-3 * np.max(q)
    if mode == "abs": return np.ones_like(q)
    if mode == "rel": return 1 / np.maximum(q, qf)
    return 1 / np.sqrt(np.maximum(q, qf))     # sqrt

def fit(model, th0, lo, hi, p, q, T=None, mode="rel"):
    w = weights(q, mode)
    res = least_squares(lambda th: w * (model(th, p, T) - q), th0,
                        bounds=(lo, hi), method="trf", xtol=1e-14, ftol=1e-14)
    ssew = float(np.sum(res.fun ** 2))
    return res.x, ssew

def aicc(ssew, N, k):
    if N - k - 1 <= 0: return np.inf
    return N * np.log(ssew / N) + 2 * k + 2 * k * (k + 1) / (N - k - 1)

def js_array(name, arr, prec=10):
    vals = ", ".join(f"{v:.{prec}g}" for v in arr)
    return f"const {name} = [{vals}];"

def main():
    ok = True
    def chk(name, err, tol):
        nonlocal ok
        p = err < tol
        ok &= p
        print(f"{'PASS' if p else 'FAIL'}  {name}: err={err:.3e} (tol {tol:.0e})")

    P = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.6, 1.0, 2.0, 4.0,
                  7.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 80.0, 100.0])

    # ---- A. zero-noise recovery ----
    for name, model, thT, th0, lo, hi in [
        ("langmuir", langmuir, [4.0, 8.0], [2.0, 1.0], [1e-6]*2, [1e3]*2),
        ("sips", sips, [5.0, 2.0, 0.7], [3.0, 1.0, 1.0], [1e-6, 1e-6, 0.1], [1e3, 1e3, 3.0]),
        ("toth", toth, [4.0, 8.0, 0.5], [2.0, 1.0, 1.0], [1e-6, 1e-6, 0.05], [1e3, 1e3, 5.0]),
    ]:
        q = model(thT, P)
        th, _ = fit(model, th0, lo, hi, P, q)
        err = np.max(np.abs(np.array(th) / np.array(thT) - 1))
        chk(f"A zero-noise {name}", err, 1e-6)

    # ---- B. noisy Langmuir data (seed 7, 3% relative noise) ----
    rng = np.random.RandomState(7)
    qL = langmuir([4.0, 8.0], P) * (1 + 0.03 * rng.randn(len(P)))
    print(js_array("PB_p", P))
    print(js_array("PB_q", qL))
    thL, sseL = fit(langmuir, [2.0, 1.0], [1e-6]*2, [1e3]*2, P, qL)
    thD, sseD = fit(dsl, [2.0, 10.0, 2.0, 0.5], [1e-6]*4, [1e3]*4, P, qL)
    aL = aicc(sseL, len(P), 2)
    aD = aicc(sseD, len(P), 4)
    print(f"    B langmuir fit: qs={thL[0]:.8f} b={thL[1]:.8f} ssew={sseL:.10e} AICc={aL:.6f}")
    print(f"    B dsl fit:      ssew={sseD:.10e} AICc={aD:.6f}")
    chk("B AICc prefers true (simpler) model", 0.0 if aL < aD else 1.0, 0.5)

    # ---- C. weighting changes the optimum (Toth data, 5% noise, seed 11) ----
    rng = np.random.RandomState(11)
    qT = toth([4.0, 8.0, 0.45], P) * (1 + 0.05 * rng.randn(len(P)))
    print(js_array("PC_q", qT))
    thR, _ = fit(toth, [2.0, 1.0, 1.0], [1e-6, 1e-6, 0.05], [1e3, 1e3, 5.0], P, qT, mode="rel")
    thA, _ = fit(toth, [2.0, 1.0, 1.0], [1e-6, 1e-6, 0.05], [1e3, 1e3, 5.0], P, qT, mode="abs")
    print(f"    C rel: qs={thR[0]:.8f} b={thR[1]:.8f} t={thR[2]:.8f}")
    print(f"    C abs: qs={thA[0]:.8f} b={thA[1]:.8f} t={thA[2]:.8f}")
    chk("C weighting matters (optima differ)",
        0.0 if np.max(np.abs(thR / thA - 1)) > 1e-3 else 1.0, 0.5)

    # ---- D. global van 't Hoff DSL, 3 temperatures, 2% noise seed 3 ----
    thTrue = [3.0, 2.466e-4, 3.0e4, 4.0, 9.058e-3, 1.2e4]
    # b1(298)=44.6, b2(298)=1.15 — CO2-zeolite-ish
    Ts = [273.15, 298.15, 323.15]
    rng = np.random.RandomState(3)
    Pm, Qm, Tm = [], [], []
    for T in Ts:
        q = dsl_vh(thTrue, P, T) * (1 + 0.02 * rng.randn(len(P)))
        Pm.extend(P); Qm.extend(q); Tm.extend([T] * len(P))
    Pm, Qm, Tm = map(np.array, (Pm, Qm, Tm))
    print(js_array("PD_q", Qm))
    th0 = [2.0, 1e-3, 2e4, 2.0, 1e-3, 2e4]
    lo = [1e-6, 1e-12, 1e3, 1e-6, 1e-12, 1e3]
    hi = [1e3, 1e3, 1.5e5, 1e3, 1e3, 1.5e5]
    thG, sseG = fit(dsl_vh, th0, lo, hi, Pm, Qm, T=Tm)
    # site order can swap: sort by U descending
    def site_sort(th):
        s = [(th[2], th[0], th[1]), (th[5], th[3], th[4])]
        s.sort(reverse=True)
        return [s[0][1], s[0][2], s[0][0], s[1][1], s[1][2], s[1][0]]
    gs = site_sort(thG)
    ts = site_sort(thTrue)
    print(f"    D global: qs1={gs[0]:.6f} b01={gs[1]:.6e} U1={gs[2]:.2f}"
          f" qs2={gs[3]:.6f} b02={gs[4]:.6e} U2={gs[5]:.2f} ssew={sseG:.8e}")
    err = max(abs(gs[0]/ts[0]-1), abs(gs[3]/ts[3]-1), abs(gs[2]/ts[2]-1), abs(gs[5]/ts[5]-1))
    chk("D global vh recovery (qs, U within 10%)", err, 0.1)

    # ---- E. AICc spot value ----
    chk("E AICc formula", abs(aicc(1.0, 20, 3) - (20*np.log(1/20) + 6 + 24/16)), 1e-12)

    print("ALL PASS" if ok else "SOME CHECKS FAILED")

if __name__ == "__main__":
    main()
