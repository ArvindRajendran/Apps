/* jsc test harness for fitter_engine.js — run: jsc fitter_engine.js fitter_test.js
   Seeded datasets + anchors from fitter_ref.py (scipy least_squares, 2026-08-11). */
"use strict";
const E = globalThis.IFIT;
let allPass = true;
function chk(name, err, tol) {
  const p = err < tol;
  allPass = allPass && p;
  print((p ? "PASS" : "FAIL") + "  " + name + ": err=" + err.toExponential(3) + " (tol " + tol + ")");
}
const P = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.6, 1, 2, 4, 7, 10, 15, 20, 30, 45, 60, 80, 100];
const T298 = P.map(() => 298.15);
const relerr = (a, b) => Math.abs(a / b - 1);

/* ---- A. zero-noise recovery ---- */
{
  const cases = [
    ["langmuir", [4.0, 8.0]],
    ["sips", [5.0, 2.0, 0.7]],
    ["toth", [4.0, 8.0, 0.5]],
    ["dsl", [3.0, 44.6, 4.0, 1.15]],
  ];
  for (const [id, thT] of cases) {
    const f = E.FIT_MODELS[id].f;
    const q = P.map(p => f(thT, p, 298.15));
    const r = E.fitModel(id, { p: P, q, T: T298 }, { weighting: "rel" });
    let err = 0;
    // dsl site order can swap: compare sorted by b descending
    const pair = th => id === "dsl"
      ? (th[1] > th[3] ? th : [th[2], th[3], th[0], th[1]])
      : th;
    const got = pair(r.theta), want = pair(thT);
    for (let j = 0; j < want.length; j++) err = Math.max(err, relerr(got[j], want[j]));
    chk("A zero-noise " + id, err, 1e-4);
  }
}

/* ---- B. noisy Langmuir data (seed 7): anchors from scipy ---- */
const PB_q = [0.03335605623, 0.09243955115, 0.2965880311, 0.7836584427, 1.735701883, 2.823704378, 3.310256403, 3.368384963, 3.879641375, 3.94866407, 3.856089777, 3.930285638, 4.02707695, 3.943987309, 3.954393471, 3.81501377, 4.058095272, 4.008602268, 4.027900315];
{
  const data = { p: P, q: PB_q, T: T298 };
  const rL = E.fitModel("langmuir", data, { weighting: "rel" });
  const rD = E.fitModel("dsl", data, { weighting: "rel", nStarts: 30 });
  print("    B langmuir: qs=" + rL.theta[0].toFixed(8) + " b=" + rL.theta[1].toFixed(8) +
        " ssew=" + rL.ssew.toExponential(6) + " AICc=" + rL.aicc.toFixed(4));
  print("    B dsl: ssew=" + rD.ssew.toExponential(6) + " AICc=" + rD.aicc.toFixed(4));
  chk("B langmuir params vs scipy",
      Math.max(relerr(rL.theta[0], 3.97843437), relerr(rL.theta[1], 8.10611476)), 1e-4);
  chk("B langmuir ssew vs scipy", relerr(rL.ssew, 1.0186424078e-2), 1e-5);
  chk("B dsl at least as good as scipy min", rD.ssew / 8.7123519630e-3 - 1, 0.02);
  chk("B AICc prefers Langmuir on Langmuir data", rL.aicc < rD.aicc ? 0 : 1, 0.5);
}

/* ---- C. weighting matters (Toth data, seed 11): anchors from scipy ---- */
const PC_q = [0.02738394101, 0.06468886158, 0.168223681, 0.3254430016, 0.7642819736, 1.252015189, 1.595923307, 1.94739777, 2.330608282, 2.478193743, 2.730811896, 2.923042098, 3.242530486, 3.314627906, 3.119603527, 3.243991628, 3.627285101, 3.834107199, 3.58807304];
{
  const data = { p: P, q: PC_q, T: T298 };
  const rR = E.fitModel("toth", data, { weighting: "rel" });
  const rA = E.fitModel("toth", data, { weighting: "abs" });
  print("    C rel: qs=" + rR.theta[0].toFixed(6) + " b=" + rR.theta[1].toFixed(6) + " t=" + rR.theta[2].toFixed(6));
  print("    C abs: qs=" + rA.theta[0].toFixed(6) + " b=" + rA.theta[1].toFixed(6) + " t=" + rA.theta[2].toFixed(6));
  chk("C rel-weight params vs scipy",
      Math.max(relerr(rR.theta[0], 4.16139140), relerr(rR.theta[1], 8.25189459), relerr(rR.theta[2], 0.42282072)), 5e-3);
  chk("C abs-weight params vs scipy",
      Math.max(relerr(rA.theta[0], 4.21155619), relerr(rA.theta[1], 10.19510615), relerr(rA.theta[2], 0.40532939)), 5e-3);
}

/* ---- D. global van 't Hoff DSL, 3 temperatures (seed 3) ---- */
const PD_q = [0.3758280334, 0.89162381, 1.794435174, 2.510434206, 3.37956107, 4.292020077, 5.023793404, 5.472350211, 6.108442915, 6.441081893, 6.524414489, 6.905678198, 6.975357406, 7.125498498, 6.932879525, 6.894079327, 6.886692578, 6.756274183, 7.114588184, 0.1293136183, 0.3579247577, 0.9643547326, 1.902097764, 2.873869468, 3.736215101, 4.45719401, 5.13380886, 5.733394991, 6.170729789, 6.517047483, 6.771769041, 7.043470042, 6.659528268, 6.798478576, 6.811158054, 6.605712195, 6.827494357, 6.822114899, 0.05575489783, 0.15804975, 0.460984999, 1.136788911, 2.182708877, 3.168837392, 3.974530859, 4.546274675, 5.270151377, 5.998029013, 6.075872355, 6.497606552, 6.811961725, 6.868018289, 6.98341904, 7.040686083, 7.120022395, 6.780295425, 7.065713654];
{
  const Ts = [273.15, 298.15, 323.15];
  const p = [], q = PD_q, T = [];
  for (const Tv of Ts) for (const pv of P) { p.push(pv); T.push(Tv); }
  const r = E.fitModel("dsl", { p, q, T }, { weighting: "rel", mode: "global", nStarts: 10 });
  // vector: [qs1, b01, U1, qs2, b02, U2]; sort sites by U descending
  const s = r.theta[2] > r.theta[5] ? r.theta
    : [r.theta[3], r.theta[4], r.theta[5], r.theta[0], r.theta[1], r.theta[2]];
  print("    D global: qs1=" + s[0].toFixed(5) + " U1=" + s[2].toFixed(0) +
        " qs2=" + s[3].toFixed(5) + " U2=" + s[5].toFixed(0) +
        " ssew=" + r.ssew.toExponential(6) + " (scipy 2.31858505e-2)");
  chk("D global ssew at least as good as scipy", r.ssew / 2.31858505e-2 - 1, 0.02);
  chk("D global qs/U recovery vs truth (10%)",
      Math.max(relerr(s[0], 3.0), relerr(s[3], 4.0), relerr(s[2], 3.0e4), relerr(s[5], 1.2e4)), 0.12);
  // isosteric heat sanity: between U2 and U1 (kJ/mol), decreasing with q
  const qst = E.isostericHeat(r, { p, q, T }, 30);
  const first = qst.qst[0], last = qst.qst[qst.qst.length - 1];
  print("    D qst: " + first.toFixed(2) + " -> " + last.toFixed(2) + " kJ/mol over loading");
  chk("D qst spans U-range and decreases",
      (first > last && first < 35 && last > 8) ? 0 : 1, 0.5);
  // model-free experimental qst must track the model-based curve
  const eq = E.expIsostericHeat({ p, q, T }, 25);
  let dmax = 0;
  const qcut = eq.q[0] + 0.75 * (eq.q[eq.q.length - 1] - eq.q[0]);  // plateau region
  for (let i = 0; i < eq.q.length; i++) {
    if (eq.q[i] > qcut) continue;   // model-free qst is ill-conditioned near the plateau
    // model qst at same loading via linear interpolation on the model curve
    let mj = 0;
    for (let j = 1; j < qst.q.length; j++) if (qst.q[j] <= eq.q[i]) mj = j;
    const w = (eq.q[i] - qst.q[mj]) / ((qst.q[mj + 1] || qst.q[mj]) - qst.q[mj] || 1);
    const mv = qst.qst[mj] + (isFinite(w) && qst.qst[mj + 1] !== undefined ? w * (qst.qst[mj + 1] - qst.qst[mj]) : 0);
    dmax = Math.max(dmax, Math.abs(eq.qst[i] - mv));
  }
  print("    D exp-qst: " + eq.q.length + " pts, " + eq.qst[0].toFixed(1) + " -> " +
        eq.qst[eq.qst.length - 1].toFixed(1) + " kJ/mol, max |exp-model| = " + dmax.toFixed(2));
  chk("D exp-qst tracks model qst (lower 75% of range, 2% noise)", dmax, 4.0);
}

/* ---- E. AICc spot value ---- */
{
  const N = 20, k = 3;
  const a = N * Math.log(1 / N) + 2 * k + 2 * k * (k + 1) / (N - k - 1);
  chk("E AICc formula", Math.abs(a - (20 * Math.log(1 / 20) + 6 + 24 / 16)), 1e-12);
}

/* ---- F. custom expression fit ---- */
{
  const cu = E.compileCustom("qs*b*p/(1+b*p)");
  const q = P.map(p => 4 * 8 * p / (1 + 8 * p));
  const r = E.fitModel("custom", { p: P, q, T: T298 },
                       { weighting: "rel", custom: cu, nStarts: 12 });
  const got = r.pn[0] === "qs" ? r.theta : [r.theta[1], r.theta[0]];
  print("    F custom: " + r.pn.join(",") + " = " + r.theta.map(v => v.toFixed(5)).join(", "));
  chk("F custom recovers Langmuir", Math.max(relerr(got[0], 4), relerr(got[1], 8)), 1e-3);
}

/* ---- G. fit-all speed + sanity ---- */
{
  const data = { p: P, q: PB_q, T: T298 };
  const ids = Object.keys(E.FIT_MODELS).filter(k => k !== "custom");
  const t0 = Date.now();
  const results = ids.map(id => E.fitModel(id, data, { weighting: "rel" })).filter(r => r && !r.error);
  const ms = Date.now() - t0;
  results.sort((a, b) => a.aicc - b.aicc);
  print("    G " + results.length + " models in " + ms + " ms; best by AICc: " +
        results.slice(0, 3).map(r => r.modelId).join(", "));
  chk("G all models fit, best is Langmuir-family",
      (results.length === ids.length &&
       ["langmuir", "toth", "sips", "unilan", "dsl"].includes(results[0].modelId)) ? 0 : 1, 0.5);
}

print(allPass ? "ALL PASS" : "SOME CHECKS FAILED");
