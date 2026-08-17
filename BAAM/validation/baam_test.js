/* BAAM engine test harness (jsc). Load after baam_engine.js:
   jsc baam_engine.js baam_test.js
   Anchors: Balashankar's BAAM.m run in MATLAB R2026a at ode15s tolerances
   tightened to 1e-10/1e-12 (converged), cross-checked by baam_ref.py (scipy
   BDF 1e-9/1e-11, agreement <= 4.4e-6). Liske worked-example values from the
   e-BAAM paper text (reproduce at T = 298.15 K). */
"use strict";
const E = globalThis.BAAMENG;
let nPass = 0, nFail = 0;
function ck(name, got, ref, tol) {
  const e = Math.abs(got - ref) / Math.max(Math.abs(ref), 1e-3);
  if (e <= tol) { nPass++; }
  else { nFail++; print(`FAIL ${name}: got ${got} ref ${ref} rel ${e.toExponential(2)}`); }
}

const ADS2019 = {
 MgMOF74: { rho: 588.25,
   co2: { qsb: 6.80, qsd: 9.90, b0: 1.81e-7, d0: 1.06e-6, dUb: -39.30e3, dUd: -21.20e3 },
   n2:  { qsb: 14.00, qsd: 0, b0: 3.45e-6, d0: 0, dUb: -15.50e3, dUd: 0 } },
 Z13X: { rho: 1130.00,
   co2: { qsb: 3.09, qsd: 2.54, b0: 8.65e-7, d0: 2.63e-8, dUb: -36.60e3, dUd: -35.70e3 },
   n2:  { qsb: 5.84, qsd: 0, b0: 2.50e-6, d0: 0, dUb: -15.80e3, dUd: 0 } },
 UTSA16: { rho: 1092.00,
   co2: { qsb: 5.00, qsd: 3.00, b0: 6.24e-7, d0: 1.87e-23, dUb: -30.60e3, dUd: -44.70e3 },
   n2:  { qsb: 12.70, qsd: 0, b0: 2.96e-6, d0: 0, dUb: -9.77e3, dUd: 0 } },
 CSAC: { rho: 799.50,
   co2: { qsb: 0.59, qsd: 7.51, b0: 9.40e-6, d0: 1.04e-5, dUb: -25.61e3, dUd: -17.55e3 },
   n2:  { qsb: 0.16, qsd: 41.30, b0: 1.81e-3, d0: 1.72e-12, dUb: -8.67e3, dUd: -44.90e3 } },
};
const Z13Xssl = { rho: 1130.0,
  co2: { qsb: 4.390, qsd: 0, b0: 2.50e-6, d0: 0, dUb: -31.19e3, dUd: 0 },
  n2:  { qsb: 4.390, qsd: 0, b0: 2.70e-6, d0: 0, dUb: -16.38e3, dUd: 0 } };

/* MATLAB tight-tolerance anchors: P: [y,qA,qB,Ncum,NAcum,W] */
const ML_PATH = {
 MgMOF74: { 0.50: [0.2876043823, 6.3383679264, 0.0396635443, 0.0894555817, 0.0184100700, 104.8157029653],
            0.15: [0.8835919471, 6.2751732940, 0.0020956179, 0.2043150025, 0.0820583554, 754.5232051312],
            0.06: [1.0, 5.3721158107, 0.0, 1.1130930110, 0.9880374654, 11802.4904151252],
            0.03: [1.0, 4.3305637932, 0.0, 2.1558533310, 2.0307977854, 30404.5081886472] },
 Z13X:    { 0.50: [0.2946411412, 3.4266930010, 0.0085133640, 0.0278515688, 0.0058008436, 32.5118342658],
            0.15: [0.9399060447, 3.4028814541, 0.0002269780, 0.0672879787, 0.0297452098, 264.6437144563],
            0.06: [1.0, 2.8688867949, 0.0, 0.6033966529, 0.5654379069, 6611.5649427036],
            0.03: [1.0, 2.3936397656, 0.0, 1.0792726945, 1.0413139484, 15083.4355067695] },
 UTSA16:  { 0.50: [0.2970936618, 2.3076889628, 0.0147455512, 0.0398508441, 0.0083354252, 46.4731757634],
            0.15: [0.9587463112, 2.2701179983, 0.0002638075, 0.0994973976, 0.0460091211, 404.9578363313],
            0.06: [1.0, 1.2879327817, 0.0, 1.0838991248, 1.0300127804, 11866.8648673425],
            0.03: [1.0, 0.7391660036, 0.0, 1.6333168039, 1.5794304596, 21496.6658651577] },
 CSAC:    { 0.50: [0.2538324510, 0.7040245098, 0.1164685187, 0.1637831932, 0.0316877768, 192.7810982330],
            0.15: [0.6363238539, 0.6331997389, 0.0203630713, 0.3410854929, 0.1034450763, 1130.1177505578],
            0.06: [0.9962494436, 0.4592213069, 0.0000972681, 0.5379968347, 0.2784806785, 3540.7479980509],
            0.03: [0.9999994154, 0.2633040184, 0.0000000087, 0.7349004183, 0.4752803344, 6989.7709482427] },
};
/* (Plow,Pint): [Pu,Re,EnBLO,EnEVAC,EnT,WC] (EnBLO/EnEVAC null where unrecorded) */
const ML_KPI = {
 MgMOF74: { "0.03,0.15": [99.856580, 77.995657, 2.444353, 96.054093, 98.498446, 1146.345970],
            "0.05,0.1":  [99.99999981, 47.69462764, null, null, 114.56186475, 509.75552449],
            "0.03,0.3":  [98.77610577, 79.86407374, null, null, 96.19407867, 1173.80714689],
            "0.06,0.5":  [94.72371444, 60.36603298, null, null, 76.84470527, 570.38331536] },
 Z13X:    { "0.03,0.15": [99.958895, 79.305127, 1.651625, 92.483140, 94.134765, 1143.072675],
            "0.05,0.1":  [100.0, 44.27210790, null, null, 121.66474956, 496.29396545],
            "0.03,0.3":  [99.32383846, 80.75682881, null, null, 92.44257775, 1163.99693601],
            "0.06,0.5":  [97.23600786, 61.85209698, null, null, 74.58353601, 632.38988148] },
 UTSA16:  { "0.03,0.15": [99.974047, 79.977967, 1.667222, 86.835122, 88.502344, 1674.496102],
            "0.05,0.1":  [100.0, 41.63041760, null, null, 128.32050829, 777.09535057],
            "0.03,0.3":  [99.38094577, 81.53369868, null, null, 86.81364468, 1707.06840793],
            "0.06,0.5":  [97.85729014, 63.11301939, null, null, 73.32758795, 1115.67167187] },
 CSAC:    { "0.03,0.15": [94.418783, 55.639333, 19.187496, 99.487041, 118.674537, 297.282289],
            "0.05,0.1":  [96.63457041, 30.18082691, null, null, 148.12687854, 152.57646164],
            "0.03,0.3":  [85.79980864, 62.48748814, null, null, 105.66870698, 333.87214661],
            "0.06,0.5":  [65.94973416, 40.01795593, null, null, 90.57482927, 197.31092495] },
};

/* A: paths vs converged MATLAB */
print("== A: BLO/EVAC paths vs MATLAB (tight) ==");
const paths = {};
for (const name in ADS2019) {
  const s = ADS2019[name];
  const mix = E.edslMix(s.co2, s.n2, 298.15);
  const t0 = Date.now();
  const path = E.blowPath(mix, { T: 298.15, rho: s.rho, PH: 1.0, yF: 0.15, Plow: 0.03 });
  paths[name] = { mix, path, spec: s };
  for (const Pq in ML_PATH[name]) {
    const i = E.baIdx(path, +Pq);
    const r = ML_PATH[name][Pq];
    ck(`${name} y@${Pq}`,    path.y[i],  r[0], 2e-5);
    ck(`${name} qA@${Pq}`,   path.q1[i], r[1], 2e-5);
    ck(`${name} qB@${Pq}`,   path.q2[i], r[2], 2e-4);
    ck(`${name} N@${Pq}`,    path.N[i],  r[3], 5e-5);
    ck(`${name} NA@${Pq}`,   path.N1[i], r[4], 5e-4);
    ck(`${name} W@${Pq}`,    path.W[i],  r[5], 5e-5);
  }
  print(`  ${name} path ${Date.now() - t0} ms`);
}

/* B: KPIs vs converged MATLAB */
print("== B: cycle KPIs vs MATLAB (tight) ==");
for (const name in ADS2019) {
  const { mix, path } = paths[name];
  for (const key in ML_KPI[name]) {
    const [Plow, Pint] = key.split(",").map(Number);
    const k = E.cycleKPIs(mix, path, { Pint, Plow, yF: 0.15, mode: "LPP" });
    const r = ML_KPI[name][key];
    ck(`${name} Pu@${key}`, k.Pu, r[0], 5e-6);
    ck(`${name} Re@${key}`, k.Re, r[1], 5e-5);
    if (r[2] !== null) ck(`${name} EnBLO@${key}`, k.EnBLO, r[2], 2e-4);
    if (r[3] !== null) ck(`${name} EnEVAC@${key}`, k.EnEVAC, r[3], 5e-5);
    ck(`${name} En@${key}`, k.En, r[4], 5e-5);
    ck(`${name} WC@${key}`, k.WC, r[5], 5e-5);
    if (Math.abs(k.mbal) > 1e-8) { nFail++; print(`FAIL ${name} mbal@${key}: ${k.mbal}`); } else nPass++;
  }
}

/* C: Python-ref cross anchors (baam_ref.py check E, BDF 1e-9) */
print("== C: cross anchors vs Python reference ==");
{
  const { mix, path } = paths.Z13X;
  const kL = E.cycleKPIs(mix, path, { Pint: 0.15, Plow: 0.03, yF: 0.15, mode: "LPP" });
  ck("py Z13X LPP Pu", kL.Pu, 99.95889492, 1e-6);
  ck("py Z13X LPP Re", kL.Re, 79.30513286, 1e-5);
  ck("py Z13X LPP En", kL.En, 94.13468934, 1e-5);
  ck("py Z13X LPP WC", kL.WC, 1143.07269, 1e-5);
  ck("py Z13X LPP ydel", kL.ydel, 0.03154846254, 1e-5);
  ck("py Z13X LPP Npr", kL.Npr, 0.109465564, 1e-4);
  ck("py Z13X LPP Nfeed", kL.Nfeed, 8.503600509, 1e-5);
  ck("py Z13X LPP Nraff", kL.Nraff, 7.533793378, 1e-5);
  const kF = E.cycleKPIs(mix, path, { Pint: 0.15, Plow: 0.03, yF: 0.15, mode: "FP" });
  ck("py Z13X FP Re", kF.Re, 78.93000395, 1e-5);
  ck("py Z13X FP ydel", kF.ydel, 0.03218977122, 1e-5);
  ck("py Z13X FP Npr", kF.Npr, 0.1231668739, 1e-4);
}

/* D: Liske worked example (13X SSL refit; PH=2, Pint=0.8, Plow=0.07,
   eta = Maruyama). Paper text: LPP 96.6/74.8, FP 96.6/74.0 at 298.15 K. */
print("== D: Liske e-BAAM worked example ==");
{
  const mix = E.edslMix(Z13Xssl.co2, Z13Xssl.n2, 298.15);
  const path = E.blowPath(mix, { T: 298.15, rho: Z13Xssl.rho, PH: 2.0, yF: 0.15,
                                 Plow: 0.07, eta: { type: "maruyama" } });
  const kL = E.cycleKPIs(mix, path, { Pint: 0.8, Plow: 0.07, yF: 0.15, mode: "LPP" });
  const kF = E.cycleKPIs(mix, path, { Pint: 0.8, Plow: 0.07, yF: 0.15, mode: "FP" });
  ck("liske LPP Pu", kL.Pu, 96.57, 6e-4);
  ck("liske LPP Re", kL.Re, 74.77, 6e-4);
  ck("liske FP Pu", kF.Pu, 96.57, 6e-4);
  ck("liske FP Re", kF.Re, 73.99, 6e-4);
  print(`  LPP ${kL.Pu.toFixed(2)}/${kL.Re.toFixed(2)} (paper 96.6/74.8), FP ${kF.Pu.toFixed(2)}/${kF.Re.toFixed(2)} (paper 96.6/74.0)`);
  // PVSA energy parts present
  if (kF.EnPR > 0 && kF.EnADS > 0 && kL.EnADS > 0 && kL.EnPR === 0) nPass++;
  else { nFail++; print("FAIL PVSA compression-work flags"); }
}

/* E: IAST(equal-qs Langmuir) === extended Langmuir (route equivalence).
   Python: Pu=99.879583 Re=78.509721/30 En=97.521436/2 WC=1618.2224/6 */
print("== E: IAST route equivalence ==");
{
  const T = 298.15, RE_ = 8.314;
  const qs = 4.39;
  const bC = 2.50e-6 * Math.exp(31.19e3 / (RE_ * T));  // m3/mol
  const bN = 2.70e-6 * Math.exp(16.38e3 / (RE_ * T));
  const specE = {
    co2: { qsb: qs, qsd: 0, b0: bC, d0: 0, dUb: 0, dUd: 0 },
    n2:  { qsb: qs, qsd: 0, b0: bN, d0: 0, dUb: 0, dUd: 0 } };
  const mixE = E.edslMix(specE.co2, specE.n2, T);
  // p-basis (1/bar): b_p = b_c / (RB*T)
  const RB = 8.314e-5;
  const c1 = E.makeComponent({ id: "langmuir", pvals: { qs, b: bC / (RB * T) } }, 1e4);
  const c2 = E.makeComponent({ id: "langmuir", pvals: { qs, b: bN / (RB * T) } }, 1e4);
  const mixI = E.iastMix(c1, c2, T);
  const t0 = Date.now();
  const pE = E.blowPath(mixE, { T, rho: 1130, PH: 1.0, yF: 0.15, Plow: 0.03, dP: 5e-4 });
  const pI = E.blowPath(mixI, { T, rho: 1130, PH: 1.0, yF: 0.15, Plow: 0.03, dP: 5e-4 });
  print(`  IAST path ${Date.now() - t0} ms`);
  const kE = E.cycleKPIs(mixE, pE, { Pint: 0.15, Plow: 0.03, yF: 0.15, mode: "LPP" });
  const kI = E.cycleKPIs(mixI, pI, { Pint: 0.15, Plow: 0.03, yF: 0.15, mode: "LPP" });
  ck("iast==edsl Pu", kI.Pu, kE.Pu, 1e-5);
  ck("iast==edsl Re", kI.Re, kE.Re, 1e-4);
  ck("iast==edsl En", kI.En, kE.En, 1e-4);
  ck("iast==edsl WC", kI.WC, kE.WC, 1e-4);
  ck("edsl vs python Pu", kE.Pu, 99.879583, 1e-5);
  ck("edsl vs python Re", kE.Re, 78.509721, 1e-5);
  ck("edsl vs python En", kE.En, 97.521436, 1e-5);
  ck("edsl vs python WC", kE.WC, 1618.2224, 1e-5);
}

/* F: grid scan + r_max sanity on 13X */
print("== F: gridScan / rMax ==");
{
  const { mix, path } = paths.Z13X;
  const plows = [], pints = [];
  for (let p = 0.03; p <= 0.1 + 1e-9; p += 0.01) plows.push(+p.toFixed(4));
  for (let p = 0.05; p <= 0.95 + 1e-9; p += 0.05) pints.push(+p.toFixed(4));
  const t0 = Date.now();
  const g = E.gridScan(mix, path, { plows, pints, yF: 0.15, mode: "LPP" });
  print(`  grid ${plows.length}x${pints.length} in ${Date.now() - t0} ms`);
  const { rmax } = E.rMax(g);
  // r at (0.03,0.15) = hypot(99.9589,79.3051) = 127.60; rmax >= that
  if (rmax >= 127.5 && rmax < 141.5) nPass++;
  else { nFail++; print(`FAIL rmax ${rmax}`); }
  // grid cell (Pint=0.15,Plow=0.03) must equal the direct call
  const i = pints.indexOf(0.15), j = plows.indexOf(0.03);
  ck("grid Pu cell", g.Pu[i * plows.length + j], 99.958895, 1e-5);
}

/* G: PVSA (PH > 1 bar) vs Python reference on the same grid the app uses
   (dP is NOT scaled with PH: the vacuum work is a first-order sum over the
   pressure slices, so its value is grid-dependent at O(dP) and a scaled grid
   would silently shift the energy). Anchors: baam_ref.py check F. */
print("== G: PVSA, PH > 1 bar ==");
{
  const CALF20 = { rho: 570.0,
    co2: { qsb: 2.387, qsd: 3.271, b0: 5.52e-7, d0: 5.19e-8, dUb: -35.06e3, dUd: -28.95e3 },
    n2:  { qsb: 2.387, qsd: 3.271, b0: 8.14e-7, d0: 0, dUb: -17.96e3, dUd: 0 } };
  const SPEC = { Z13Xssl, CALF20 };
  const CASES = [
    {ads:"Z13Xssl", T:298.15, PH:2.0, yF:0.15, Pint:0.8, Plow:0.07, mode:"LPP", eta:"mar", ref:{Pu: 96.57387619, Re: 74.76742953, En: 224.4741691, EnADS: 140.6307167, EnPR: 0, WC: 1041.081318, ydel: 0.04082364828}},
    {ads:"Z13Xssl", T:298.15, PH:3.0, yF:0.15, Pint:1.5, Plow:0.08, mode:"LPP", eta:"mar", ref:{Pu: 94.49486599, Re: 79.56487444, En: 292.3836744, EnADS: 222.4938896, EnPR: 0, WC: 1057.59531, ydel: 0.03327541514}},
    {ads:"Z13Xssl", T:298.15, PH:3.0, yF:0.15, Pint:1.5, Plow:0.08, mode:"FP", eta:"mar", ref:{Pu: 94.49486599, Re: 78.63813968, En: 295.0057204, EnADS: 215.1377459, EnPR: 9.978189721, WC: 1057.59531, ydel: 0.03484618718}},
    {ads:"Z13Xssl", T:298.15, PH:5.0, yF:0.15, Pint:2.5, Plow:0.1, mode:"LPP", eta:"mar", ref:{Pu: 91.94383292, Re: 82.2355829, En: 395.1815346, EnADS: 340.8319907, EnPR: 0, WC: 969.3440939, ydel: 0.02821946117}},
    {ads:"Z13Xssl", T:298.15, PH:5.0, yF:0.15, Pint:0.6, Plow:0.05, mode:"LPP", eta:"mar", ref:{Pu: 99.51680572, Re: 87.8731664, En: 413.438878, EnADS: 318.9656021, EnPR: 0, WC: 1661.6954, ydel: 0.01408813131}},
    {ads:"CALF20", T:298.15, PH:5.0, yF:0.15, Pint:2.0, Plow:0.06, mode:"FP", eta:"const", ref:{Pu: 95.038738, Re: 88.12553473, En: 383.6032678, EnADS: 320.3650516, EnPR: 15.6129618, WC: 668.5636075, ydel: 0.01601590363}},
  ];
  let worst = 0, worstName = "";
  for (const c of CASES) {
    const s = SPEC[c.ads];
    const mix = E.edslMix(s.co2, s.n2, c.T);
    const eta = c.eta === "mar" ? { type: "maruyama" } : { val: 0.72 };
    const path = E.blowPath(mix, { T: c.T, rho: s.rho, eps: 0.37, PH: c.PH,
                                   yF: c.yF, Plow: c.Plow, dP: 1e-4, eta });
    // no vacuum work may be charged above 1 bar
    let wAbove = 0;
    for (let i = 0; i < path.P.length; i++) if (path.P[i] >= 1) wAbove = Math.max(wAbove, Math.abs(path.W[i]));
    if (wAbove === 0) nPass++;
    else { nFail++; print(`FAIL ${c.ads} PH=${c.PH}: W=${wAbove} above 1 bar`); }
    const k = E.cycleKPIs(mix, path, { Pint: c.Pint, Plow: c.Plow, yF: c.yF, mode: c.mode });
    const tag = `${c.ads} PH=${c.PH} ${c.mode}`;
    for (const key of ["Pu", "Re", "En", "EnADS", "EnPR", "WC", "ydel"]) {
      if (c.ref[key] === 0) { if (k[key] === 0) nPass++; else { nFail++; print(`FAIL ${tag} ${key} should be 0, got ${k[key]}`); } continue; }
      ck(`${tag} ${key}`, k[key], c.ref[key], 2e-5);
      const e = Math.abs(k[key] - c.ref[key]) / Math.abs(c.ref[key]);
      if (e > worst) { worst = e; worstName = `${tag} ${key}`; }
    }
    // energy split must close exactly
    ck(`${tag} split`, k.EnBLO + k.EnEVAC + k.EnADS + k.EnPR, k.En, 1e-12);
  }
  print(`  worst rel err vs Python = ${worst.toExponential(2)} (${worstName})`);
}

print(`\n${nPass} passed, ${nFail} failed`);
if (nFail > 0) throw new Error("BAAM engine tests FAILED");
