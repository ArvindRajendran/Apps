// DAC engine test harness (jsc).  Load order:
//   jsc dac_engine.js dac_anchors.js dac_test.js
var npass = 0, nfail = 0, worst = 0, worstName = "";

function check(name, got, ref, tol) {
  var denom = Math.max(Math.abs(ref), 1e-10);
  var rel = Math.abs(got - ref) / denom;
  if (rel > worst) { worst = rel; worstName = name; }
  if (rel <= tol) { npass++; }
  else { nfail++; print("FAIL " + name + ": got " + got + " ref " + ref + " rel " + rel.toExponential(2)); }
}

// ---------- A. anchor cases: KPIs + CSS state vs Python ----------
var TOL_KPI = 2e-4, TOL_STATE = 5e-4;
for (var name in DAC_ANCHORS) {
  var a = DAC_ANCHORS[name];
  var res = DacEngine.runToCSS(a.u, { record: false });
  var k = res.kpi;
  check(name + ".purity", k.purityDry, a.purityDry, TOL_KPI);
  check(name + ".recovery", k.recovery, a.recovery, TOL_KPI);
  check(name + ".Eth", k.Eth, a.Eth, TOL_KPI);
  check(name + ".Efan", k.Efan, a.Efan, TOL_KPI);
  check(name + ".Evac", k.Evac, a.Evac, TOL_KPI);
  check(name + ".Pr", k.prContact, a.prContact, TOL_KPI);
  check(name + ".w/c", k.waterPerCO2, a.waterPerCO2, TOL_KPI);
  check(name + ".mCO2", k.mCO2cycle, a.mCO2cycle, TOL_KPI);
  check(name + ".t_evac", k.times.evac, a.t_evac, 2e-3);
  // simplified-model comparison values (working capacity, swing, energies)
  if (a.simp) {
    ["dqC", "dqW", "dT", "etaCap", "Eth", "Efan", "Evac"].forEach(function (key) {
      check(name + ".simp." + key, k.simp[key], a.simp[key], TOL_KPI);
    });
  }
  for (var i = 0; i < 7; i++)
    check(name + ".state[" + i + "]", res.state[i], a.state[i], TOL_STATE);

  // mole balance on a CSS cycle (trapezoid quadrature => looser than Python)
  var cyc2 = DacEngine.runCycle(res.M, res.state, null);
  var mb = DacEngine.moleBalance(res.M, res.state, cyc2);
  if (Math.abs(mb.C) < 2e-5 && Math.abs(mb.N) < 2e-5 &&
      (isNaN(mb.W) || Math.abs(mb.W) < 2e-5)) npass++;
  else { nfail++; print("FAIL " + name + " mole balance " + mb.C + " " + mb.W + " " + mb.N); }
}

// ---------- B. vacuum-work identity on the heat step ----------
// At constant P_regen with the isothermal pump, W_vac(heat) must equal
// N_out,heat * R * Tamb * ln(Pamb/Pregen) / etaVac exactly.
(function () {
  var u = JSON.parse(JSON.stringify(DAC_ANCHORS["Lewatit-base"].u));
  var res = DacEngine.runToCSS(u, { record: false });
  var M = res.M;
  // isolate the heat step: run one more cycle, capture Q per step
  var Q1 = DacEngine.newQuad();
  var s1 = DacEngine.integrateStep(M, "ads", res.state, 0, u.t_ads, Q1, null, null);
  var Q2 = DacEngine.newQuad();
  var s2 = DacEngine.integrateStep(M, "evac", s1.x, 0, 4 * 3600, Q2, null, function (xx) {
    return xx[0] * DacEngine.R * xx[5] / M.VC - u.Pregen;
  });
  var x2 = s2.x.slice();
  x2[0] = u.Pregen * M.VC / (DacEngine.R * x2[5]);
  var Q3 = DacEngine.newQuad();
  DacEngine.integrateStep(M, "heat", x2, 0, u.t_heat, Q3, null, null);
  var Nout = Q3.prodC + Q3.prodW + Q3.prodN;
  var Wident = Nout * DacEngine.R * u.Tamb * Math.log(u.Pamb / u.Pregen) / u.etaVac;
  check("vac-identity(heat)", Q3.Wvac, Wident, 1e-6);
})();

// ---------- C. isotherm sanity ----------
// dry Toth reduction: Lewatit at 400 ppm, 293.15 K, qw=0
(function () {
  var S = { toth: DacEngine.PRESETS["Lewatit"].toth, gab: null, wet: null };
  var Rg = 8.314, T = 293.15, p = 400e-6 * 101325;
  var t = S.toth;
  var b = t.b0 * Math.exp((t.dH0 / (Rg * t.T0)) * (t.T0 / T - 1));
  var w = t.w0 + t.alpha * (1 - t.T0 / T);
  var f = b * p;
  var qref = t.qs0 * f / Math.pow(1 + Math.pow(f, w), 1 / w);
  check("toth-dry-reduction", DacEngine.qstarC(p, T, 0, S), qref, 1e-14);
  // GAB at 50% RH
  var g = DacEngine.PRESETS["Lewatit"].gab;
  var Ka = g.K * 0.5;
  var qwref = g.Cm * g.CG * Ka / ((1 - Ka) * (1 + (g.CG - 1) * Ka));
  check("gab-value", DacEngine.qstarW(0.5, { gab: g }), qwref, 1e-14);
})();

// ---------- D. energy-split closure ----------
(function () {
  var res = DacEngine.runToCSS(DAC_ANCHORS["Lewatit-base"].u, { record: false });
  var k = res.kpi;
  check("split-closure", k.Etot, k.Eth + k.Efan + k.Evac, 1e-14);
  check("elec-closure", k.Eel, k.Efan + k.Evac, 1e-14);
  var s = k.simp;
  check("simp-th-closure", s.Eth, s.Esens + s.EdesC + s.EdesW, 1e-14);
  check("simp-tot-closure", s.Etot, s.Eth + s.Efan + s.Evac, 1e-14);
  // the CO2 desorption term of the simplified model is dH_CO2 per mol CO2
  check("simp-dHc", s.EdesC * 1e6 * 0.04401, DacEngine.PRESETS["Lewatit"].dHc, 1e-12);
})();

// ---------- E. sorbent-property overrides ----------
(function () {
  var P = DacEngine.PRESETS["Lewatit"];
  var dH0_before = P.toth.dH0;
  var base = JSON.parse(JSON.stringify(DAC_ANCHORS["Lewatit-base"].u));
  // geometry overrides
  var u1 = JSON.parse(JSON.stringify(base));
  u1.eb = 0.30; u1.ep = 0.45; u1.rhop = 900; u1.cps = 1200;
  var M1 = DacEngine.buildModel(u1);
  check("prop.ms", M1.ms, 900 * 0.70, 1e-14);
  check("prop.epsT", M1.epsT, 0.30 + 0.70 * 0.45, 1e-14);
  check("prop.VC", M1.VC, M1.epsT, 1e-14);
  check("prop.cps", M1.cps, 1200, 1e-14);
  // linked dH: the Toth temperature coefficient follows the heat of adsorption
  var u2 = JSON.parse(JSON.stringify(base)); u2.dHc = 80000;
  var M2 = DacEngine.buildModel(u2);
  check("prop.dHc.linked.energy", M2.dHc, 80000, 1e-14);
  check("prop.dHc.linked.toth", M2.S.toth.dH0, 80000, 1e-14);
  // unlinked: energy balance only, isotherm untouched
  var u3 = JSON.parse(JSON.stringify(base)); u3.dHc = 80000; u3.unlinkDH = true;
  var M3 = DacEngine.buildModel(u3);
  check("prop.dHc.unlinked.energy", M3.dHc, 80000, 1e-14);
  check("prop.dHc.unlinked.toth", M3.S.toth.dH0, dH0_before, 1e-14);
  // the preset object itself must never be mutated
  check("prop.preset-intact", DacEngine.PRESETS["Lewatit"].toth.dH0, dH0_before, 1e-14);
  // defaults reproduce the preset exactly
  var M0 = DacEngine.buildModel(base);
  check("prop.default.ms", M0.ms, P.rho_p * (1 - P.eb), 1e-14);
  check("prop.default.dHc", M0.dHc, P.dHc, 1e-14);
  // water-inert sorbent ignores an edited dHw
  var u4 = JSON.parse(JSON.stringify(base)); u4.sorbent = "TMCM-41"; u4.dHw = 50000;
  check("prop.inert-dHw", DacEngine.buildModel(u4).dHw, 0, 1e-14);
})();

/* ---------------- Arrhenius LDF temperature dependence -------------- */
(function () {
  var R = DacEngine.R, TREF = 293.15;
  // Ea = 0 (or absent) is exactly the constant-k model
  check("ldf.Ea0", DacEngine.kLDF(4e-3, 0, 393.15), 4e-3, 1e-15);
  check("ldf.Ea-undef", DacEngine.kLDF(4e-3, undefined, 393.15), 4e-3, 1e-15);
  // at the reference temperature the slider value is returned unchanged
  check("ldf.at-Tref", DacEngine.kLDF(4e-3, 30000, TREF), 4e-3, 1e-14);
  // analytic Arrhenius ratio between two temperatures
  var T1 = 293.15, T2 = 393.15, Ea = 30000;
  var ratio = DacEngine.kLDF(1, Ea, T2) / DacEngine.kLDF(1, Ea, T1);
  check("ldf.ratio", ratio, Math.exp(-(Ea / R) * (1 / T2 - 1 / T1)), 1e-14);
  // k(T_des) > k(T_amb) for Ea > 0 (desorption faster than capture)
  check("ldf.monotone", DacEngine.kLDF(4e-3, 30000, 393.15) > 4e-3 ? 1 : 0, 1, 0);
  // full-cycle Ea = 0 identity: adding EaC/EaW = 0 to a converged anchor
  // case must not change a single KPI bit
  var u0 = JSON.parse(JSON.stringify(DAC_ANCHORS["Lewatit-base"].u));
  u0.EaC = 0; u0.EaW = 0;
  var k0 = DacEngine.runToCSS(u0, { record: false }).kpi;
  var kr = DacEngine.runToCSS(DAC_ANCHORS["Lewatit-base"].u, { record: false }).kpi;
  check("ldf.Ea0-cycle-identity.Pu", k0.purityDry, kr.purityDry, 0);
  check("ldf.Ea0-cycle-identity.Re", k0.recovery, kr.recovery, 0);
  check("ldf.Ea0-cycle-identity.Etot", k0.Etot, kr.Etot, 0);
})();

/* --------------------- Custom isotherm (v1.7) ----------------------- */
(function () {
  function isoFrom(name, over) {
    var S = JSON.parse(JSON.stringify(DacEngine.PRESETS[name]));
    for (var k in over) S[k] = over[k];
    return S;
  }
  // 1. Identity: Custom seeded with a preset's exact numbers must give a
  // bit-identical cycle (wet/modulated preset and dry/inert preset)
  ["Lewatit-base", "TMCM-dry"].forEach(function (anchor) {
    var u1 = JSON.parse(JSON.stringify(DAC_ANCHORS[anchor].u));
    var u2 = JSON.parse(JSON.stringify(u1));
    u2.iso = isoFrom(u1.sorbent, {});
    u2.sorbent = "Custom";
    var k1 = DacEngine.runToCSS(u1, { record: false }).kpi;
    var k2 = DacEngine.runToCSS(u2, { record: false }).kpi;
    check("cus.identity." + anchor + ".Pu", k2.purityDry, k1.purityDry, 0);
    check("cus.identity." + anchor + ".Re", k2.recovery, k1.recovery, 0);
    check("cus.identity." + anchor + ".Etot", k2.Etot, k1.Etot, 0);
  });
  // 2. bw0 = 0 (or cmp absent) reduces exactly to the plain wet-less Toth
  var Sd = { toth: DacEngine.PRESETS["Lewatit"].toth, gab: null, wet: null };
  var Sc0 = { toth: Sd.toth, gab: null, wet: null, cmp: { bw0: 0, dHw: 43800 } };
  var Sc = { toth: Sd.toth, gab: null, wet: null, cmp: { bw0: 2e-3, dHw: 43800 } };
  var p = 40, T = 298.15, pw = 1200;
  check("cus.cmp-bw0-reduction", DacEngine.qstarC(p, T, 0, Sc0, pw),
        DacEngine.qstarC(p, T, 0, Sd), 0);
  check("cus.cmp-pw0-reduction", DacEngine.qstarC(p, T, 0, Sc, 0),
        DacEngine.qstarC(p, T, 0, Sd), 0);
  // 3. competition only ever DECREASES q*, monotonically in pw
  var prev = DacEngine.qstarC(p, T, 0, Sd), mono = 1;
  [200, 600, 1200, 2400, 4800].forEach(function (pwv) {
    var q = DacEngine.qstarC(p, T, 0, Sc, pwv);
    if (!(q < prev)) mono = 0;
    prev = q;
  });
  check("cus.cmp-monotone-suppression", mono, 1, 0);
})();

print("");
print(npass + " passed, " + nfail + " failed; worst rel dev " +
      worst.toExponential(2) + " (" + worstName + ")");
