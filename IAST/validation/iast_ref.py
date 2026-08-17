#!/usr/bin/env python3
"""Independent Python reference for the IAST app.

Binary IAST (Myers & Prausnitz): psi_1(p1o) = psi_2(p2o), P*y_i = x_i*p_i_o,
1/qT = sum x_i/q_i(p_i_o).

ALGORITHM (deliberately different from the JS engine): bisection on p1o with
x1 = P*y1/p1o, p2o = P*y2/(1-x1), root of G(p1o) = psi1(p1o) - psi2(p2o).
G is strictly increasing (psi1 up, p2o down) -> unique root. Only forward psi
evaluations; no inversion tables (the JS engine bisects on psi with tabulated
inverses -> algorithmic diversity for the cross-check).

Checks:
  A. equal-qs binary Langmuir == extended (competitive) Langmuir, machine prec
  B. pure-component limits y1 -> 0, 1
  C. Henry limit: selectivity(P->0) -> H1/H2
  D. closure residual sum(x)-1 at solution
  E. Toth(t=1) numerical-psi == Langmuir analytic-psi IAST
  F. BET closed-form psi == numerical quadrature
  G. cross-check anchors (DSL vs Sips pair) for the JS engine
  H. pyIAST cross-validation (skipped gracefully if not installed)
"""
import math

# ---------------- adaptive Simpson in ln p for numerical psi ----------------
def _simp(g, a, b, fa, fm, fb, tol, depth):
    m = 0.5 * (a + b)
    lm, rm = 0.5 * (a + m), 0.5 * (m + b)
    flm, frm = g(lm), g(rm)
    left = (m - a) / 6 * (fa + 4 * flm + fm)
    right = (b - m) / 6 * (fm + 4 * frm + fb)
    whole = (b - a) / 6 * (fa + 4 * fm + fb)
    if depth <= 0 or abs(left + right - whole) < 15 * tol:
        return left + right + (left + right - whole) / 15
    return (_simp(g, a, m, fa, flm, fm, tol / 2, depth - 1) +
            _simp(g, m, b, fm, frm, fb, tol / 2, depth - 1))

def psi_numeric(f, p, p0_frac=1e-10):
    """psi(p) = int_0^p q/p' dp'. Tail below p0 via local power-law exponent."""
    if p <= 0:
        return 0.0
    p0 = p * p0_frac
    q0, qh = f(p0), f(p0 / 2)
    alpha = math.log(q0 / qh) / math.log(2) if qh > 0 else 1.0
    tail = q0 / alpha if alpha > 0 else 0.0
    g = lambda u: f(math.exp(u))
    a, b = math.log(p0), math.log(p)
    fa, fb = g(a), g(b)
    fm = g(0.5 * (a + b))
    scale = max(abs(fb), 1e-300)
    return tail + _simp(g, a, b, fa, fm, fb, 1e-12 * scale * (b - a), 48)

# ---------------- isotherm models: (f, psi_analytic_or_None, henry) ---------
def langmuir(qs, b):
    return (lambda p: qs * b * p / (1 + b * p),
            lambda p: qs * math.log(1 + b * p), qs * b)

def dsl(qs1, b1, qs2, b2):
    return (lambda p: qs1 * b1 * p / (1 + b1 * p) + qs2 * b2 * p / (1 + b2 * p),
            lambda p: qs1 * math.log(1 + b1 * p) + qs2 * math.log(1 + b2 * p),
            qs1 * b1 + qs2 * b2)

def sips(qs, b, n):
    return (lambda p: qs * (b * p) ** n / (1 + (b * p) ** n),
            lambda p: qs / n * math.log(1 + (b * p) ** n),
            qs * b if n == 1 else (math.inf if n < 1 else 0.0))

def toth(qs, b, t):
    return (lambda p: qs * b * p / (1 + (b * p) ** t) ** (1 / t), None, qs * b)

def bet(qs, bs, bl):
    return (lambda p: qs * bs * p / ((1 - bl * p) * (1 - bl * p + bs * p)),
            lambda p: qs * math.log((1 + (bs - bl) * p) / (1 - bl * p)), qs * bs)

def make_psi(iso):
    f, psiA, _ = iso
    return psiA if psiA is not None else (lambda p: psi_numeric(f, p))

# ---------------- IAST solve: bisection on p1o ----------------
def iast(P, y1, iso1, iso2, tol=1e-13, pcap=1e10):
    f1, _, _ = iso1
    f2, _, _ = iso2
    if y1 <= 0:
        return dict(x1=0.0, q1=0.0, q2=f2(P), qT=f2(P), p1o=0.0, p2o=P)
    if y1 >= 1:
        return dict(x1=1.0, q1=f1(P), q2=0.0, qT=f1(P), p1o=P, p2o=0.0)
    psi1, psi2 = make_psi(iso1), make_psi(iso2)
    y2 = 1 - y1
    def G(p1o):
        x1 = P * y1 / p1o
        p2o = P * y2 / (1 - x1)
        return psi1(p1o) - psi2(p2o)
    lo = P * y1 * (1 + 1e-12)        # x1 -> 1-: p2o -> inf: G < 0
    hi = P * y1 * 2
    while G(hi) < 0 and hi < pcap:
        hi *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if G(mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol * hi:
            break
    p1o = 0.5 * (lo + hi)
    x1 = P * y1 / p1o
    p2o = P * y2 / (1 - x1)
    q1o, q2o = f1(p1o), f2(p2o)
    qT = 1.0 / (x1 / q1o + (1 - x1) / q2o)
    return dict(x1=x1, q1=x1 * qT, q2=(1 - x1) * qT, qT=qT, p1o=p1o, p2o=p2o)

def selectivity(r, y1):
    return (r["x1"] / y1) / ((1 - r["x1"]) / (1 - y1))

# ================= checks =================
def main():
    ok = True
    def chk(name, err, tol):
        nonlocal ok
        p = err < tol
        ok &= p
        print(f"{'PASS' if p else 'FAIL'}  {name}: err={err:.3e} (tol {tol:.0e})")

    # ---- A. equal-qs Langmuir pair == extended Langmuir (exact IAST result)
    qs, b1, b2 = 4.0, 8.0, 1.5
    isoA1, isoA2 = langmuir(qs, b1), langmuir(qs, b2)
    errA = 0.0
    for P in (0.01, 0.1, 1.0, 10.0, 100.0):
        for y1 in (0.05, 0.3, 0.5, 0.8, 0.95):
            r = iast(P, y1, isoA1, isoA2)
            d = 1 + b1 * P * y1 + b2 * P * (1 - y1)
            q1el = qs * b1 * P * y1 / d
            q2el = qs * b2 * P * (1 - y1) / d
            errA = max(errA, abs(r["q1"] - q1el), abs(r["q2"] - q2el))
    chk("A IAST == extended Langmuir (equal qs)", errA, 1e-9)

    # ---- B. pure limits
    iso1 = dsl(3.0, 50.0, 4.0, 1.0)
    iso2 = sips(5.0, 2.0, 1.0)          # n=1: Langmuir-equivalent, Henry ok
    r = iast(1.0, 1 - 1e-12, iso1, iso2)
    errB1 = abs(r["q1"] - iso1[0](1.0))
    r = iast(1.0, 1e-12, iso1, iso2)
    errB2 = abs(r["q2"] - iso2[0](1.0))
    chk("B pure-component limits", max(errB1, errB2), 1e-6)

    # ---- C. Henry limit selectivity -> H1/H2
    H1, H2 = iso1[2], iso2[2]
    r = iast(1e-9, 0.5, iso1, iso2)
    errC = abs(selectivity(r, 0.5) / (H1 / H2) - 1)
    chk("C Henry-limit selectivity", errC, 1e-5)
    print(f"    H1/H2 = {H1/H2:.6f}  S(P=1e-9) = {selectivity(r,0.5):.6f}")

    # ---- D. closure residual
    r = iast(5.0, 0.37, iso1, iso2)
    x2 = 5.0 * 0.63 / r["p2o"]
    chk("D closure sum(x)-1", abs(r["x1"] + x2 - 1), 1e-10)

    # ---- E. Toth(t=1) numerical psi == Langmuir analytic
    isoT = toth(4.0, 8.0, 1.0)
    isoL = langmuir(4.0, 8.0)
    errE = 0.0
    for P in (0.1, 1.0, 10.0):
        for y1 in (0.2, 0.5, 0.8):
            rT = iast(P, y1, isoT, isoA2)
            rL = iast(P, y1, isoL, isoA2)
            errE = max(errE, abs(rT["q1"] - rL["q1"]), abs(rT["q2"] - rL["q2"]))
    chk("E Toth(t=1) quadrature == analytic", errE, 1e-7)

    # ---- F. BET psi closed form vs quadrature
    fB, psiB, _ = bet(2.0, 30.0, 0.8)
    errF = max(abs(psi_numeric(fB, p) / psiB(p) - 1) for p in (0.01, 0.1, 0.5, 1.0, 1.2))
    chk("F BET psi closed form vs quadrature", errF, 1e-8)

    # ---- G. cross-check anchors for JS (DSL comp1, Sips n=0.85 comp2)
    isoG1 = dsl(3.0, 50.0, 4.0, 1.0)
    isoG2 = sips(5.0, 2.0, 0.85)
    print("    G anchors: P y1 -> q1 q2 x1 p1o p2o")
    for (P, y1) in ((0.1, 0.3), (1.0, 0.5), (10.0, 0.15), (100.0, 0.85)):
        r = iast(P, y1, isoG1, isoG2)
        print(f"    G {P:g} {y1:g} -> {r['q1']:.8f} {r['q2']:.8f} "
              f"{r['x1']:.8f} {r['p1o']:.6e} {r['p2o']:.6e}")

    # ---- H. pyIAST cross-validation (graceful skip)
    try:
        import numpy as np
        import pandas as pd
        import pyiast
        def fit_lang(qs_, b_):
            P = np.logspace(-3, 2, 60)
            df = pd.DataFrame({"P": P, "L": qs_ * b_ * P / (1 + b_ * P)})
            return pyiast.ModelIsotherm(df, loading_key="L", pressure_key="P",
                                        model="Langmuir",
                                        param_guess={"M": qs_, "K": b_})
        i1, i2 = fit_lang(3.0, 20.0), fit_lang(5.0, 0.7)
        print(f"    pyIAST fitted params: {i1.params} {i2.params}")
        errH = 0.0
        ours1, ours2 = langmuir(3.0, 20.0), langmuir(5.0, 0.7)
        for (P, y1) in ((0.5, 0.2), (2.0, 0.5), (20.0, 0.8)):
            qpy = pyiast.iast(np.array([P * y1, P * (1 - y1)]), [i1, i2],
                              verboseflag=False)
            r = iast(P, y1, ours1, ours2)
            errH = max(errH, abs(qpy[0] - r["q1"]) / r["q1"],
                       abs(qpy[1] - r["q2"]) / r["q2"])
        chk("H pyIAST cross-validation (rel)", errH, 1e-4)
    except Exception as e:
        print(f"SKIP  H pyIAST cross-validation ({type(e).__name__}: {e})")

    print("ALL PASS" if ok else "SOME CHECKS FAILED")

if __name__ == "__main__":
    main()
