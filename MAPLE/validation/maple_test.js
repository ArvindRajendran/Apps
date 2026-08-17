/* maple_test.js — jsc harness for the MAPLE JS engine.
 * Run: jsc maple_weights.js maple_engine.js maple_anchors.js
 *          maple_labels_LPP.js maple_labels_FP.js maple_paper_pts.js maple_test.js
 * Section A: engine vs MATLAB network outputs (200 LHS rows + 12 named
 *            preset x condition anchors, per cycle, 5 KPIs) — expect ~1e-12.
 * Section B: R2 of engine predictions vs the detailed-model labels over
 *            ALL labelled rows; must match MATLAB's R2 to 1e-6 and stay
 *            above the published 0.99 claim. NOTE the labels file has
 *            swapped Pu/Re column headers (documented in README) — the
 *            mapping below uses the corrected semantics.
 * Section C: guard sanity + isotherm identities + NSGA-II smoke test.
 * Section D: engine vs the 216 PUBLISHED optimum points of the Limits
 *            paper SI (detailed-model KPIs at MAPLE-Opt conditions,
 *            curves E5/P5, E6/P6, E7/P7) — agreement at the surrogate's
 *            own accuracy (per-cent level), the external anchor.
 */
"use strict";
var NPASS = 0, NFAIL = 0, WORST = { v: 0, id: "" };
function check(id, got, want, tol) {
  var den = Math.max(Math.abs(want), 1e-10);
  var rel = Math.abs(got - want) / den;
  if (rel > WORST.v && isFinite(rel)) { WORST.v = rel; WORST.id = id; }
  if (rel <= tol && isFinite(got)) { NPASS++; }
  else { NFAIL++; print("FAIL " + id + ": got " + got + " want " + want + " rel " + rel); }
}
function ok(id, cond) {
  if (cond) NPASS++; else { NFAIL++; print("FAIL " + id); }
}

var CYCLES = ["LPP", "FP"];

/* ---------- Section A: vs MATLAB net outputs ---------- */
CYCLES.forEach(function (cyc) {
  var A = MAPLE_ANCHORS[cyc], S = A.samples;
  var n = S.inputs.length;
  for (var i = 0; i < n; i++) {
    var r = S.inputs[i];
    var u = { qsat: r[0], b0c: r[1], b0n: r[2], dUc: r[3], dUn: r[4], rho: r[5],
              y: r[6], tads: r[7], PH: r[8], PI: r[9], PL: r[10], vF: r[11] };
    var k = MAPLE.evalKPI(cyc, u);
    check(cyc + ".s" + i + ".Pu", k.Pu, S.outputs.purity[i], 1e-10);
    check(cyc + ".s" + i + ".Re", k.Re, S.outputs.recovery[i], 1e-10);
    check(cyc + ".s" + i + ".En", k.En, S.outputs.energy[i], 1e-10);
    check(cyc + ".s" + i + ".Pr", k.Prod, S.outputs.productivity[i], 1e-10);
    check(cyc + ".s" + i + ".E100", k.En100, S.outputs.energy100[i], 1e-10);
  }
  A.presets.forEach(function (P) {
    var r = P.inputs;
    var u = { qsat: r[0], b0c: r[1], b0n: r[2], dUc: r[3], dUn: r[4], rho: r[5],
              y: r[6], tads: r[7], PH: r[8], PI: r[9], PL: r[10], vF: r[11] };
    var k = MAPLE.evalKPI(cyc, u);
    check(cyc + "." + P.name + ".Pu", k.Pu, P.purity, 1e-10);
    check(cyc + "." + P.name + ".Re", k.Re, P.recovery, 1e-10);
    check(cyc + "." + P.name + ".En", k.En, P.energy, 1e-10);
    check(cyc + "." + P.name + ".Pr", k.Prod, P.productivity, 1e-10);
    check(cyc + "." + P.name + ".E100", k.En100, P.energy100, 1e-10);
  });
});

/* ---------- Section B: R2 vs detailed-model labels ---------- */
var LBL = { LPP: MAPLE_LABELS_LPP, FP: MAPLE_LABELS_FP };
CYCLES.forEach(function (cyc) {
  var D = LBL[cyc], rows = D.rows, n = rows.length;
  var names = ["Pu", "Re", "En", "Prod", "En100"];
  // label columns: index 13 carries the swapped "Re" header but IS
  // purity; index 12 carries "Pu" but IS recovery
  var li = [13, 12, 14, 15, 16];
  var sum = [0,0,0,0,0], se = [0,0,0,0,0], cnt = [0,0,0,0,0], mean = [0,0,0,0,0];
  var preds = new Array(n);
  for (var i = 0; i < n; i++) {
    var r = rows[i];
    var u = { qsat: r[0], b0c: r[1], b0n: r[2], dUc: r[3], dUn: r[4], rho: r[5],
              y: r[6], tads: r[7], PH: r[8], PI: r[9], PL: r[10], vF: r[11] };
    var k = MAPLE.evalKPI(cyc, u);
    preds[i] = [k.Pu, k.Re, k.En, k.Prod, k.En100];
    for (var j = 0; j < 5; j++) {
      var L = r[li[j]];
      if (isFinite(L)) { sum[j] += L; cnt[j]++; }
    }
  }
  for (j = 0; j < 5; j++) mean[j] = sum[j] / cnt[j];
  var st = [0,0,0,0,0];
  for (i = 0; i < n; i++) {
    var r2 = rows[i];
    for (j = 0; j < 5; j++) {
      var L2 = r2[li[j]];
      if (isFinite(L2)) {
        se[j] += (L2 - preds[i][j]) * (L2 - preds[i][j]);
        st[j] += (L2 - mean[j]) * (L2 - mean[j]);
      }
    }
  }
  for (j = 0; j < 5; j++) {
    var R2 = 1 - se[j] / st[j];
    check(cyc + ".R2." + names[j], R2, MAPLE_ANCHORS[cyc].R2[j], 1e-6);
    ok(cyc + ".R2claim." + names[j] + " (" + R2.toFixed(6) + ")", R2 > 0.99);
    print("  " + cyc + " R2(" + names[j] + ") = " + R2.toFixed(6) + " over " + cnt[j] + " rows");
  }
});

/* ---------- Section C: guards, isotherms, optimizer smoke ---------- */
var base = Object.assign({ y: 0.15, tads: 60, PH: 1, PI: 0.15, PL: 0.03, vF: 0.8 },
                         MAPLE.PRESETS.Z13X);
ok("guard.base-inside", MAPLE.guards(base).length === 0);
ok("guard.y-out", MAPLE.guards(Object.assign({}, base, { y: 0.7 })).length === 1);
ok("guard.PL-ratio", MAPLE.guards(Object.assign({}, base, { PL: 0.001 })).length === 2); // ratio + abs
ok("guard.PIPL-order", MAPLE.guards(Object.assign({}, base, { PI: 0.02 })).length >= 1);
// absolute sampled-pressure windows: the LHS sampled ABSOLUTE PI/PL, so a
// point can be ratio-legal yet in a joint data void at low P_H
ok("guard.PI-absfloor", MAPLE.guards(Object.assign({}, base, { PI: 0.05 })).length === 1);
ok("guard.PL-absfloor", MAPLE.guards(Object.assign({}, base, { PL: 0.008 })).length === 1);
ok("guard.PL-abscap", MAPLE.guards(Object.assign({}, base, { PH: 5, PI: 3, PL: 1.5 })).length === 1);
ok("guard.abs-inside", MAPLE.guards(Object.assign({}, base, { PL: 0.011 })).length === 0);
ok("guard.selectivity", MAPLE.guards(Object.assign({}, base, { b0n: base.b0c, dUn: base.dUc }))
    .length >= 1);
// Henry limit: q*/P -> qsat*b*C/P as P->0; also q(b298 identity)
var qlo = MAPLE.qStar(base, 1e-6, 298, "c");
var H = base.qsat * MAPLE.bT(base.b0c, base.dUc, 298) * (1e-6 * 1e5 / (8.314 * 298));
check("iso.henry", qlo, H, 1e-4);
var bb = MAPLE.bT(base.b0c, base.dUc, 298);
check("iso.b298", bb, base.b0c * Math.exp(-base.dUc * 1000 / (8.314 * 298)), 1e-14);
// binary reduces to pure when y=1
var qb = MAPLE.qStarBinary(base, 0.5, 1.0, 303.15);
check("iso.binary-pure", qb.qc, MAPLE.qStar(base, 0.5, 303.15, "c"), 1e-12);
// monotone saturation
ok("iso.saturating", MAPLE.qStar(base, 5, 303.15, "c") < base.qsat);

// NSGA-II smoke: Pu-Re Pareto for 13X at y=0.15, small budget
MAPLE.setSeed(20260815);
var prob = MAPLE.optProblem("LPP", MAPLE.PRESETS.Z13X, 0.15, "PuRe", true);
var front = MAPLE.nsga2({ nvar: prob.nvar, lo: prob.lo, hi: prob.hi,
                          fobj: prob.fobj, pop: 60, gens: 30 });
ok("nsga.front-nonempty", front.length >= 10);
var maxPu = 0, maxRe = 0;
front.forEach(function (p) {
  maxPu = Math.max(maxPu, -p.f[0]); maxRe = Math.max(maxRe, -p.f[1]);
});
print("  Pu-Re front: " + front.length + " pts, max Pu " + maxPu.toFixed(1) +
      "%, max Re " + maxRe.toFixed(1) + "%");
ok("nsga.attains-highPu", maxPu > 90);
ok("nsga.attains-highRe", maxRe > 90);
// front sanity: no member dominated by another
var dom = 0;
for (var a = 0; a < front.length; a++) for (var b = 0; b < front.length; b++) {
  if (a !== b) {
    var fa = front[a].f, fb = front[b].f;
    if (fa[0] <= fb[0] && fa[1] <= fb[1] && (fa[0] < fb[0] || fa[1] < fb[1])) dom++;
  }
}
ok("nsga.front-nondominated", dom === 0);

// stepwise API must reproduce the monolithic run exactly (same rnd sequence)
MAPLE.setSeed(20260815);
var S = MAPLE.nsga2init({ nvar: prob.nvar, lo: prob.lo, hi: prob.hi,
                          fobj: prob.fobj, pop: 60, gens: 30 });
var ngen = 0;
while (MAPLE.nsga2step(S)) ngen++;
ngen++;                                  // last step ran a generation too
var front2 = MAPLE.nsga2front(S);
ok("nsga.stepwise-gens", ngen === 30);
ok("nsga.stepwise-size", front2.length === front.length);
var sdev = 0;
for (var q = 0; q < Math.min(front.length, front2.length); q++) {
  sdev = Math.max(sdev, Math.abs(front[q].f[0] - front2[q].f[0]),
                        Math.abs(front[q].f[1] - front2[q].f[1]));
}
ok("nsga.stepwise-identical", sdev === 0);

// EnPr mode: penalized bi-objective (min En, max Prod) with Pu/Re constraints
MAPLE.setSeed(777);
var probE = MAPLE.optProblem("LPP", MAPLE.PRESETS.Z13X, 0.15, "EnPr", true,
                             { puMin: 95, reMin: 90 });
// budget = the app's defaults: with the penalty added to BOTH objectives,
// infeasible points form a domination chain, so front size is budget- and
// feasibility-sensitive (worst on VSA 13X, whose feasible set is marginal).
var frE = MAPLE.nsga2({ nvar: probE.nvar, lo: probE.lo, hi: probE.hi,
                        fobj: probE.fobj, pop: 120, gens: 70 });
ok("nsga.enpr-front-nonempty", frE.length >= 5);
var nFeas = 0, minEn = Infinity, maxPr = 0;
frE.forEach(function (p) {
  var k = MAPLE.evalKPI("LPP", probE.toU(p.x));
  if (Math.min(k.Pu, 100) >= 94.5 && Math.min(k.Re, 100) >= 89.5) {
    nFeas++; minEn = Math.min(minEn, k.En); maxPr = Math.max(maxPr, k.Prod);
  }
});
print("  En-Prod front: " + frE.length + " pts, " + nFeas + " near-feasible, min En " +
      (isFinite(minEn) ? minEn.toFixed(0) : "-") + " kWh/t, max Prod " + maxPr.toFixed(2));
ok("nsga.enpr-feasible-points", nFeas >= 5);
ok("nsga.enpr-energy-sane", isFinite(minEn) && minEn > 50 && minEn < 2000);
// constraint defaults: 6-arg call without cons behaves like puMin=95
var pd = MAPLE.optProblem("LPP", MAPLE.PRESETS.Z13X, 0.15, "EnRe", true);
var pc = MAPLE.optProblem("LPP", MAPLE.PRESETS.Z13X, 0.15, "EnRe", true, { puMin: 95, reMin: 90 });
var xprobe = [60, 1.0, -0.9, -1.5, 0.8];
check("nsga.cons-default-f0", pd.fobj(xprobe)[0], pc.fobj(xprobe)[0], 1e-15);
check("nsga.cons-default-f1", pd.fobj(xprobe)[1], pc.fobj(xprobe)[1], 1e-15);

// penalty must carry a LINEAR term: with a pure quadratic the GA parks the
// EnPr front a hair below the constraint and VSA runs tabulate as infeasible.
// Regression: VSA EnPr on 13X must deliver actually-feasible front points.
MAPLE.setSeed(4242);
var probV = MAPLE.optProblem("LPP", MAPLE.PRESETS.Z13X, 0.15, "EnPr", true,
                             { puMin: 95, reMin: 90 });
var frV = MAPLE.nsga2({ nvar: probV.nvar, lo: probV.lo, hi: probV.hi,
                        fobj: probV.fobj, pop: 120, gens: 70 });
var nV = 0;
frV.forEach(function (p) {
  var k = MAPLE.evalKPI("LPP", probV.toU(p.x));
  if (k.Pu > 0 && k.Pu < 102 && k.Re > 0 && k.Re < 102 && k.En > 0 && k.Prod > 0 &&
      Math.min(k.Pu, 100) >= 95 && Math.min(k.Re, 100) >= 90) nV++;
});
print("  EnPr VSA 13X: " + nV + "/" + frV.length + " strictly feasible front points");
ok("nsga.enpr-vsa-feasible", nV >= 5);

// trained-support regression: every optimizer front point must respect the
// ABSOLUTE sampled-pressure windows (PI >= 0.07, PL >= 0.01 bar). Without
// them, VSA fronts park in a joint data void (e.g. P_L ~ 2 mbar at P_H = 1,
// zero training rows) and produce artifact fronts (flat Pu ~ 99.5-100).
MAPLE.setSeed(97532);
var probU = MAPLE.optProblem("LPP", MAPLE.PRESETS.UTSA16, 0.15, "PuRe", true);
var frU = MAPLE.nsga2({ nvar: probU.nvar, lo: probU.lo, hi: probU.hi,
                        fobj: probU.fobj, pop: 120, gens: 70 });
var minPI = 1e9, minPL = 1e9, puLo = 1e9, puHi = -1e9;
frU.forEach(function (p) {
  var u = probU.toU(p.x);
  minPI = Math.min(minPI, u.PI); minPL = Math.min(minPL, u.PL);
  puLo = Math.min(puLo, -p.f[0]); puHi = Math.max(puHi, -p.f[0]);
});
print("  PuRe VSA UTSA-16: front PI >= " + minPI.toFixed(4) + ", PL >= " + minPL.toFixed(4) +
      ", Pu " + puLo.toFixed(1) + "-" + puHi.toFixed(1));
ok("nsga.support-PI", minPI >= MAPLE.PABS.PImin * 0.999);
ok("nsga.support-PL", minPL >= MAPLE.PABS.PLmin * 0.999);

/* ---- user-narrowable DV ranges (absolute units, narrow-only) ---- */
var rNarrow = MAPLE.dvRanges({ tads: [40, 60], PL: [0.03, 0.2], vF: [0.5, 0.5] }, false);
check("dvr.tads-lo", rNarrow.tads[0], 40, 1e-12);
check("dvr.tads-hi", rNarrow.tads[1], 60, 1e-12);
check("dvr.PL-lo", rNarrow.PL[0], 0.03, 1e-12);
check("dvr.vF-degenerate", rNarrow.vF[1], 0.5, 1e-12);
check("dvr.untouched-PH", rNarrow.PH[1], 5, 1e-12);
// widening is refused (clamped back to the trained support)
var rWide = MAPLE.dvRanges({ PL: [0.0001, 9], tads: [0, 500], PH: [0.2, 40] }, false);
check("dvr.noWiden-PL-lo", rWide.PL[0], MAPLE.PABS.PLmin, 1e-12);
check("dvr.noWiden-PL-hi", rWide.PL[1], MAPLE.PABS.PLmax, 1e-12);
check("dvr.noWiden-tads", rWide.tads[1], 110, 1e-12);
check("dvr.noWiden-PH-lo", rWide.PH[0], 1, 1e-12);
// inverted input collapses to a point rather than an empty interval
var rInv = MAPLE.dvRanges({ tads: [90, 30] }, false);
ok("dvr.inverted-safe", rInv.tads[1] >= rInv.tads[0]);
// fixPH pins P_H regardless of the user's range
check("dvr.fixPH", MAPLE.dvRanges({ PH: [2, 4] }, true).PH[1], 1, 1e-12);
// the optimizer must actually respect a narrowed range
MAPLE.setSeed(31415);
var probN = MAPLE.optProblem("LPP", MAPLE.PRESETS.Z13X, 0.15, "PuRe", false,
                             null, { tads: [40, 60], PL: [0.03, 0.2], PH: [1, 2] });
var frN = MAPLE.nsga2({ nvar: probN.nvar, lo: probN.lo, hi: probN.hi,
                        fobj: probN.fobj, pop: 80, gens: 40 });
var bad = 0;
frN.forEach(function (p) {
  var u = probN.toU(p.x);
  if (u.tads < 39.99 || u.tads > 60.01 || u.PL < 0.0299 || u.PL > 0.2001 ||
      u.PH < 0.999 || u.PH > 2.001) bad++;
});
print("  narrowed-range front: " + frN.length + " pts, " + bad + " outside the requested box");
ok("dvr.optimizer-respects", bad === 0 && frN.length >= 5);

/* ---------- Section D: vs published optimum points (external) ---------- */
(function () {
  var devs = { Pu: [], Re: [], En: [], Prod: [] };
  MAPLE_PAPER_PTS.forEach(function (r) {
    var u = { qsat: r.qsat, b0c: r.b0c, b0n: r.b0n, dUc: r.dUc, dUn: r.dUn, rho: r.rho,
              y: r.y, tads: r.tads, PH: r.PH, PI: r.PI, PL: r.PL, vF: r.vF };
    var k = MAPLE.evalKPI("LPP", u);            // E5-E7/P5-P7 are 4-step LPP
    devs.Pu.push(Math.abs(Math.min(k.Pu, 100) - r.Pu) / r.Pu);
    devs.Re.push(Math.abs(Math.min(k.Re, 100) - r.Re) / r.Re);
    devs.En.push(Math.abs(k.En - r.En) / r.En);
    devs.Prod.push(Math.abs(k.Prod - r.Prod) / r.Prod);
  });
  ok("paper.count (" + MAPLE_PAPER_PTS.length + ")", MAPLE_PAPER_PTS.length >= 200);
  var LIM = { Pu: [0.02, 0.05], Re: [0.02, 0.05], En: [0.04, 0.12], Prod: [0.03, 0.10] };
  Object.keys(devs).forEach(function (kk) {
    var v = devs[kk].slice().sort(function (a, b) { return a - b; });
    var med = v[(v.length / 2) | 0], p90 = v[(0.9 * v.length) | 0];
    print("  paper-pts " + kk + ": median " + (100 * med).toFixed(2) + "%  p90 " +
      (100 * p90).toFixed(2) + "%  max " + (100 * v[v.length - 1]).toFixed(2) + "%");
    ok("paper." + kk + ".median<" + LIM[kk][0], med < LIM[kk][0]);
    ok("paper." + kk + ".p90<" + LIM[kk][1], p90 < LIM[kk][1]);
  });
})();

print("=== " + NPASS + " passed, " + NFAIL + " failed; worst rel dev " +
      WORST.v.toExponential(2) + " (" + WORST.id + ")");
