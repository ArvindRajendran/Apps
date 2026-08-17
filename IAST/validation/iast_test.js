/* jsc test harness for iast_engine.js — run: jsc iast_engine.js iast_test.js
   Anchors from iast_ref.py (2026-08-11 run, all PASS incl. pyIAST at 1.7e-13). */
"use strict";
const E = globalThis.IAST;
let allPass = true;
function chk(name, err, tol) {
  const p = err < tol;
  allPass = allPass && p;
  print((p ? "PASS" : "FAIL") + "  " + name + ": err=" + err.toExponential(3) + " (tol " + tol + ")");
}
const mk = (id, pvals) => E.makeComponent({ id, pvals }, 1e10);

/* ---- A. equal-qs Langmuir pair == extended Langmuir ---- */
{
  const c1 = mk("langmuir", { qs: 4, b: 8 }), c2 = mk("langmuir", { qs: 4, b: 1.5 });
  let err = 0;
  for (const P of [0.01, 0.1, 1, 10, 100])
    for (const y1 of [0.05, 0.3, 0.5, 0.8, 0.95]) {
      const r = E.iastSolve(P, y1, c1, c2);
      const el = E.extendedLangmuir(P, y1, { qs: 4, b: 8 }, { qs: 4, b: 1.5 });
      err = Math.max(err, Math.abs(r.q1 - el.q1), Math.abs(r.q2 - el.q2));
    }
  chk("A IAST == extended Langmuir (equal qs)", err, 1e-9);
}

/* ---- B. pure limits ---- */
{
  const c1 = mk("dsl", { qs1: 3, b1: 50, qs2: 4, b2: 1 });
  const c2 = mk("sips", { qs: 5, b: 2, n: 1 });
  const r1 = E.iastSolve(1, 1 - 1e-12, c1, c2);
  const r0 = E.iastSolve(1, 1e-12, c1, c2);
  chk("B pure-component limits",
      Math.max(Math.abs(r1.q1 - c1.f(1)), Math.abs(r0.q2 - c2.f(1))), 1e-6);
}

/* ---- C. Henry-limit selectivity -> H1/H2 = 15.4 ---- */
{
  const c1 = mk("dsl", { qs1: 3, b1: 50, qs2: 4, b2: 1 });
  const c2 = mk("sips", { qs: 5, b: 2, n: 1 });
  const r = E.iastSolve(1e-9, 0.5, c1, c2);
  const S = E.iastSelectivity(r, 0.5);
  print("    S(P=1e-9) = " + S.toFixed(6) + " (want 15.4)");
  chk("C Henry-limit selectivity", Math.abs(S / 15.4 - 1), 1e-5);
}

/* ---- D. closure residual ---- */
{
  const c1 = mk("dsl", { qs1: 3, b1: 50, qs2: 4, b2: 1 });
  const c2 = mk("sips", { qs: 5, b: 2, n: 1 });
  const r = E.iastSolve(5, 0.37, c1, c2);
  const x2 = 5 * 0.63 / r.p2o;
  chk("D closure sum(x)-1", Math.abs(r.x1 + x2 - 1), 1e-10);
}

/* ---- E. Toth(t=1) table-psi == Langmuir analytic ---- */
{
  const cT = mk("toth", { qs: 4, b: 8, t: 1 });
  const cL = mk("langmuir", { qs: 4, b: 8 });
  const c2 = mk("langmuir", { qs: 4, b: 1.5 });
  let err = 0;
  for (const P of [0.1, 1, 10])
    for (const y1 of [0.2, 0.5, 0.8]) {
      const rT = E.iastSolve(P, y1, cT, c2);
      const rL = E.iastSolve(P, y1, cL, c2);
      err = Math.max(err, Math.abs(rT.q1 - rL.q1) / rL.q1, Math.abs(rT.q2 - rL.q2) / rL.q2);
    }
  chk("E Toth(t=1) table == analytic (rel)", err, 1e-6);
}

/* ---- F. BET psi closed form vs table ----
   Table path is only load-bearing for pole-free models (Toth/custom); BET in
   the app uses analytic psi. Check away from the pole (p <= 80% of 1/bl). */
{
  const cB = mk("bet", { qs: 2, bs: 30, bl: 0.8 });
  const tab = E.buildPsiTable(cB.f, 1 / 0.8 * (1 - 1e-9));
  let err = 0;
  for (const p of [0.01, 0.1, 0.5, 1.0])
    err = Math.max(err, Math.abs(tab.eval(p) / cB.psi(p) - 1));
  chk("F BET psi closed form vs table (rel)", err, 2e-6);
}

/* ---- G. cross anchors vs Python (DSL / Sips n=0.85) ---- */
{
  const c1 = mk("dsl", { qs1: 3, b1: 50, qs2: 4, b2: 1 });
  const c2 = mk("sips", { qs: 5, b: 2, n: 0.85 });
  const PY = [   // P, y1, q1, q2, x1, p1o, p2o  (iast_ref.py)
    [0.1, 0.3, 1.69484730, 0.41146444, 0.80465169, 3.728321e-2, 3.583343e-1],
    [1.0, 0.5, 3.86052802, 0.47496954, 0.89044636, 5.615161e-1, 4.563974e0],
    [10.0, 0.15, 3.36197662, 2.02236558, 0.62439877, 2.402311e0, 2.263038e1],
    [100.0, 0.85, 6.89883293, 0.03895783, 0.99438469, 8.548000e1, 2.671270e3]
  ];
  let err = 0;
  for (const [P, y1, q1, q2, x1, p1o, p2o] of PY) {
    const r = E.iastSolve(P, y1, c1, c2);
    err = Math.max(err,
      Math.abs(r.q1 / q1 - 1), Math.abs(r.q2 / q2 - 1), Math.abs(r.x1 / x1 - 1),
      Math.abs(r.p1o / p1o - 1), Math.abs(r.p2o / p2o - 1));
    print("    G P=" + P + " y1=" + y1 + ": q1=" + r.q1.toFixed(8) +
          " q2=" + r.q2.toFixed(8) + " x1=" + r.x1.toFixed(8));
  }
  chk("G cross anchors vs python (rel)", err, 2e-6);
}

/* ---- H. sweep + xy sanity, timing ---- */
{
  const c1 = mk("dsl", { qs1: 3, b1: 50, qs2: 4, b2: 1 });
  const c2 = mk("sips", { qs: 5, b: 2, n: 0.85 });
  const t0 = Date.now();
  const sw = E.iastSweep(c1, c2, 0.01, 100, 61, 61);
  const ms = Date.now() - t0;
  let mono = true;                // q1 increasing in y1 at fixed P
  for (let i = 0; i < sw.NP; i++)
    for (let j = 1; j < sw.NY; j++)
      if (sw.q1[i * sw.NY + j] < sw.q1[i * sw.NY + j - 1] - 1e-9) mono = false;
  const xy = E.xyCurve(c1, c2, 1.0);
  print("    sweep 61x61 in " + ms + " ms · extMax=" + sw.extMax.toExponential(2) +
        " · x(y=0.5)=" + xy.x[60].toFixed(6) + " · azeotropes: " + xy.az.length);
  chk("H sweep sane (q1 monotone in y1, finite)", mono && isFinite(sw.extMax) ? 0 : 1, 0.5);
}

/* ---- I. selectivity reversal (IAST cannot make composition azeotropes) ----
   S = p2o/p1o depends on psi only, so S=1 forces p1o=p2o=P: the x-y curve can
   touch the diagonal only at P = p* where psi1(p*)=psi2(p*). Checks:
   (a) p* found for a crossing pair; (b) S(P=p*, any y) = 1; (c) no interior
   x=y crossing at P != p*; (d) S same sign of (S-1) across compositions. */
{
  const c1 = mk("langmuir", { qs: 2.0, b: 20 });   // strong site, small capacity
  const c2 = mk("langmuir", { qs: 6.0, b: 0.8 });  // weak site, big capacity
  const rev = E.findSReversal(c1, c2, 1e-6, 1e6).ps;
  print("    I S-reversal pressures: [" + rev.map(v => v.toFixed(5)).join(", ") + "]");
  let err = rev.length === 1 ? 0 : 1;
  if (rev.length === 1) {
    const pstar = rev[0];
    for (const y1 of [0.2, 0.5, 0.8]) {
      const S = E.iastSelectivity(E.iastSolve(pstar, y1, c1, c2), y1);
      err = Math.max(err, Math.abs(S - 1));
    }
    const xy50 = E.xyCurve(c1, c2, 50.0);
    const xy1 = E.xyCurve(c1, c2, 1.0);
    if (xy50.az.length !== 0 || xy1.az.length !== 0) err = Math.max(err, 1);
    const Slow = E.iastSelectivity(E.iastSolve(1.0, 0.5, c1, c2), 0.5);
    const Shigh = E.iastSelectivity(E.iastSolve(50.0, 0.5, c1, c2), 0.5);
    print("    I p*=" + pstar.toFixed(4) + " bar · S(1 bar)=" + Slow.toFixed(3) +
          " · S(50 bar)=" + Shigh.toFixed(3));
    if (!(Slow > 1 && Shigh < 1)) err = Math.max(err, 1);
  }
  chk("I selectivity reversal at psi-crossing", err, 1e-6);
}

/* ---- J. degenerate pair (identical isotherms) ---- */
{
  const c1 = mk("quad", { qs: 3, b1: 5, b2: 10 });
  const c2 = mk("quad", { qs: 3, b1: 5, b2: 10 });
  const rr = E.findSReversal(c1, c2, 1e-6, 1e6);
  const r = E.iastSolve(2.0, 0.3, c1, c2);
  print("    J degenerate=" + rr.degenerate + " ps=" + rr.ps.length +
        " x1(y=0.3)=" + r.x1.toFixed(6));
  chk("J identical pair -> degenerate, x=y",
      (rr.degenerate && rr.ps.length === 0) ? Math.abs(r.x1 - 0.3) : 1, 1e-9);
}

print(allPass ? "ALL PASS" : "SOME CHECKS FAILED");
