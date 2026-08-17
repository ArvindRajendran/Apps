/* jsc test harness for uptake_engine.js — run: jsc uptake_engine.js uptake_test.js
   Mirrors uptake_ref.py checks; PYREF anchors filled from the Python run. */
"use strict";
const U = globalThis.UPTAKE;
let allPass = true;
function chk(name, err, tol) {
  const p = err < tol;
  allPass = allPass && p;
  print((p ? "PASS" : "FAIL") + "  " + name + ": err=" + err.toExponential(3) + " (tol " + tol + ")");
}
function sphereSeries(tau) {
  let s = 0;
  for (let n = 1; n <= 2000; n++) s += Math.exp(-n * n * Math.PI * Math.PI * tau) / (n * n);
  return 1 - (6 / (Math.PI * Math.PI)) * s;
}
function linspace(a, b, n) {
  const o = new Float64Array(n);
  for (let i = 0; i < n; i++) o[i] = a + (b - a) * i / (n - 1);
  return o;
}
function interp1(xq, xs, ys) {
  let lo = 0, hi = xs.length - 1;
  if (xq <= xs[0]) return ys[0];
  if (xq >= xs[hi]) return ys[hi];
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= xq) lo = m; else hi = m; }
  const w = (xq - xs[lo]) / (xs[hi] - xs[lo]);
  return ys[lo] + w * (ys[hi] - ys[lo]);
}

const Tk = 298.15, P = 1e5, CT = P / (8.314462618 * Tk);
const epsp = 0.35, rhop = 1000.0, Rp = 1.5e-4;
// De = 6.5e-7 built from Dm, DK, taup, epsp: choose Dm=DK=2*Dpore, Dpore=De*taup/epsp
const taup = 3.0;
const Dpore = 6.5e-7 * taup / epsp;
const Dm = 2 * Dpore, DK = 2 * Dpore; // Bosanquet gives Dpore back
const H = 3.0;

const base = { T: Tk, P, epsp, rhop, Rp, Dm, DK, taup, y0: 0.0, y1: 0.2,
               darken: false, kldf: 0.05 };

/* ---- A. micropore constant-D vs series ---- */
{
  const amu = 1e-3, Nc = 200;
  const cfg = { ...base, model: "micro", iso: { id: "linear", pvals: { H } },
                Dmu0: amu, Rc: 1.0, Nc, Np: 10 };
  const model = U.buildModel(cfg);
  const tt = linspace(0, 0.5 / amu, 4001);
  const res = U.integrateCN(model, tt, {});
  let maxe = 0;
  for (const tau of [0.01, 0.05, 0.1, 0.2, 0.4]) {
    const Ui = interp1(tau / amu, tt, res.uptake);
    maxe = Math.max(maxe, Math.abs(Ui - sphereSeries(tau)));
    print("    micro tau=" + tau + ": JS=" + Ui.toFixed(6) + " series=" + sphereSeries(tau).toFixed(6));
  }
  chk("A micro vs series", maxe, 5e-4);
}

/* ---- B. macropore linear vs effective-D series ---- */
{
  const beta = epsp + rhop * H / CT;
  const Deff = 6.5e-7 / (Rp * Rp * beta);
  const cfg = { ...base, model: "macro", iso: { id: "linear", pvals: { H } },
                Dmu0: 1e-12, Rc: 1e-6, Nc: 10, Np: 200 };
  const model = U.buildModel(cfg);
  print("    De check: " + model.De.toExponential(4) + " (want 6.5e-7)");
  const tt = linspace(0, 0.5 / Deff, 4001);
  const res = U.integrateCN(model, tt, {});
  let maxe = 0;
  for (const tau of [0.05, 0.1, 0.2]) {
    const Ui = interp1(tau / Deff, tt, res.uptake);
    maxe = Math.max(maxe, Math.abs(Ui - sphereSeries(tau)));
  }
  chk("B macro vs eff-D series", maxe, 5e-4);
}

/* ---- D1/D2 bidisperse limits ---- */
{
  const beta = epsp + rhop * H / CT;
  const Deff = 6.5e-7 / (Rp * Rp * beta);
  let cfg = { ...base, model: "bidisperse", iso: { id: "linear", pvals: { H } },
              Dmu0: 50.0, Rc: 1.0, Nc: 8, Np: 60 };
  let model = U.buildModel(cfg);
  let tt = linspace(0, 0.45 / Deff, 6001);
  let res = U.integrateCN(model, tt, {});
  let maxe = 0;
  for (const tau of [0.05, 0.1, 0.2])
    maxe = Math.max(maxe, Math.abs(interp1(tau / Deff, tt, res.uptake) - sphereSeries(tau)));
  chk("D1 bidisperse macro-limit", maxe, 3e-3);

  const amu2 = 1e-2;
  cfg = { ...base, model: "bidisperse", iso: { id: "linear", pvals: { H } },
          Dmu0: amu2, Rc: 1.0, Nc: 60, Np: 15,
          Dm: 2e-2 * taup / epsp * 2, DK: 2e-2 * taup / epsp * 2 }; // De=1e-2... set below
  // easier: pick Dm,DK so De=1e-2: Dpore=De*taup/epsp
  const Dp2 = 1e-2 * taup / epsp;
  cfg.Dm = 2 * Dp2; cfg.DK = 2 * Dp2;
  model = U.buildModel(cfg);
  tt = linspace(0, 0.45 / amu2, 6001);
  res = U.integrateCN(model, tt, {});
  maxe = 0;
  for (const tau of [0.05, 0.1, 0.2])
    maxe = Math.max(maxe, Math.abs(interp1(tau / amu2, tt, res.uptake) - sphereSeries(tau)));
  chk("D2 bidisperse micro-limit", maxe, 5e-3);
}

/* ---- E1. Gamma table vs exact Langmuir ---- */
{
  const qs = 5.0, b = 10.0;
  const f = U.makeIsotherm({ id: "langmuir", pvals: { qs, b } });
  const gt = U.gammaTable(f, 0, 0.9);
  let maxe = 0;
  for (const q of [0.5, 1.0, 2.0, 3.0, 4.0]) {
    const y = q / (b * (qs - q));
    const exact = 1 + b * y;
    maxe = Math.max(maxe, Math.abs(gt.G(q) - exact) / exact);
  }
  chk("E1 Gamma table vs exact", maxe, 2e-3);
}

/* ---- E2. Darken Langmuir ads/des t50 vs Python ---- */
{
  const amu3 = 1e-3, Nq = 120;
  const PYREF = { ads: [11.24, 59.69], des: [19.97, 140.85] };  // from uptake_ref.py
  for (const [tag, y0v, y1v, Tmul] of [["ads", 0.0, 0.5, 1.2], ["des", 0.5, 0.0, 6.0]]) {
    const cfg = { ...base, model: "micro", iso: { id: "langmuir", pvals: { qs: 5, b: 10 } },
                  darken: true, y0: y0v, y1: y1v, Dmu0: amu3, Rc: 1.0, Nc: Nq, Np: 10 };
    const model = U.buildModel(cfg);
    const tt = U.timesSqrt(Tmul / amu3, 700);
    const res = U.integrateCN(model, tt, {});
    const t50 = interp1(0.5, res.uptake, tt), t90 = interp1(0.9, res.uptake, tt);
    print("    Darken " + tag + ": t50=" + t50.toFixed(3) + " t90=" + t90.toFixed(3) +
          (isNaN(PYREF[tag][0]) ? "" :
           "  (py " + PYREF[tag][0].toFixed(3) + "/" + PYREF[tag][1].toFixed(3) + ")"));
    if (!isNaN(PYREF[tag][0]))
      chk("E2 " + tag + " t50/t90 vs python",
          Math.max(Math.abs(t50 - PYREF[tag][0]) / PYREF[tag][0],
                   Math.abs(t90 - PYREF[tag][1]) / PYREF[tag][1]), 5e-3);
  }
}

/* ---- F. bidisperse anchors vs Python ---- */
{
  const amu4 = 1e-2;
  const PYF = [0.599232, 0.912529, 0.998233, 1.000000, 1.000000];   // uptake_ref.py
  const cfg = { ...base, model: "bidisperse", iso: { id: "linear", pvals: { H } },
                Dmu0: amu4, Rc: 1.0, Nc: 20, Np: 25 };
  const model = U.buildModel(cfg);
  const tt = linspace(0, 600, 12001);
  const res = U.integrateCN(model, tt, {});
  const tsF = [5, 20, 60, 150, 400];
  let maxe = 0, msg = "    F anchors JS:";
  tsF.forEach((tv, i) => {
    const Ui = interp1(tv, tt, res.uptake);
    msg += " " + Ui.toFixed(6);
    if (!isNaN(PYF[i])) maxe = Math.max(maxe, Math.abs(Ui - PYF[i]));
  });
  print(msg);
  if (!isNaN(PYF[0])) chk("F bidisperse linear vs python", maxe, 1.5e-3);

  const PYF2 = [0.841206, 0.999919, 1.000000, 1.000000];       // uptake_ref.py
  const cfg2 = { ...base, model: "bidisperse", iso: { id: "langmuir", pvals: { qs: 5, b: 10 } },
                 darken: true, y1: 0.5, Dmu0: amu4, Rc: 1.0, Nc: 12, Np: 15 };
  const model2 = U.buildModel(cfg2);
  const tt2 = U.timesSqrt(400, 500);
  const res2 = U.integrateCN(model2, tt2, {});
  let maxe2 = 0, msg2 = "    F2 anchors JS:";
  [5, 20, 60, 150].forEach((tv, i) => {
    const Ui = interp1(tv, tt2, res2.uptake);
    msg2 += " " + Ui.toFixed(6);
    if (!isNaN(PYF2[i])) maxe2 = Math.max(maxe2, Math.abs(Ui - PYF2[i]));
  });
  print(msg2);
  if (!isNaN(PYF2[0])) chk("F2 bidisperse nonlinear vs python", maxe2, 2e-3);
}

/* ---- G. solveUptake smoke test (app defaults) ---- */
{
  const cfg = { model: "bidisperse", iso: { id: "langmuir", pvals: { qs: 5, b: 10 } },
                T: Tk, P, y0: 0, y1: 0.5, darken: true,
                Dmu0: 5e-14, Rc: 1e-6, Dm: 1e-5, DK: 1.26e-5, taup: 3, epsp: 0.35,
                rhop: 1000, Rp: 1.5e-4, kldf: 0.05, Nc: 15, Np: 20, K: 240 };
  const t0 = Date.now();
  const sol = U.solveUptake(cfg, null);
  const ms = Date.now() - t0;
  const ok = sol.U[0] === 0 && Math.abs(sol.U[sol.U.length - 1] - 1) < 0.05 &&
             sol.snaps && sol.snaps.length > 30;
  print("    smoke: T=" + sol.T.toFixed(2) + "s U_end=" + sol.U[sol.U.length-1].toFixed(4) +
        " snaps=" + sol.snaps.length + " (" + ms + " ms)");
  chk("G solveUptake smoke", ok ? 0 : 1, 0.5);
}

print(allPass ? "ALL PASS" : "SOME CHECKS FAILED");
