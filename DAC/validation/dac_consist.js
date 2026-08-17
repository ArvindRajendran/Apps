// Instant-kinetics consistency case, run on the JS engine (validated
// against the Python reference to <1e-5).  Prints a JSON blob with the
// P-struct for MATLAB dac.energyModel plus the dynamic energy terms.
//   jsc dac_engine.js dac_consist.js > dac_consist.json
var R = DacEngine.R, MW_C = 0.04401;
var u = {
  sorbent: "Lewatit",
  Tamb: 293.15, RH: 50, Tdes: 393.15,
  Pregen: 0.20e5, Pamb: 101325,
  t_ads: 6 * 3600, t_heat: 4 * 3600, t_cool: 2 * 3600,
  Vfan: 20000, Spump: 500,
  dPfan: 1013, etaFan: 0.70, etaVac: 0.50, pumpModel: "iso",
  tauHeat: 20, tauCool: 60, Ua: 1e7, wallCp: 0,
  kC: 0.5, kW: 0.5,
  cpaC: 0, cpaW: 0            // zero adsorbed-phase cp in BOTH models
};
var res = DacEngine.runToCSS(u, { record: false, tolCSS: 1e-7, maxCycles: 80 });
var M = res.M;

// per-step pass at CSS to isolate the heat step
var s1 = DacEngine.integrateStep(M, "ads", res.state, 0, u.t_ads, DacEngine.newQuad(), null, null);
var s2 = DacEngine.integrateStep(M, "evac", s1.x, 0, 4 * 3600, DacEngine.newQuad(), null,
  function (xx) { return xx[0] * R * xx[5] / M.VC - u.Pregen; });
var x2 = s2.x.slice();
x2[0] = u.Pregen * M.VC / (R * x2[5]);
var Q3 = DacEngine.newQuad();
var s3 = DacEngine.integrateStep(M, "heat", x2, 0, u.t_heat, Q3, null, null);

var xA = s1.x, xB = s3.x;
var dqc = xA[3] - xB[3], dqw = xA[4] - xB[4], dT = u.Tdes - u.Tamb;

// full-cycle quantities at CSS
var cyc = DacEngine.runCycle(M, res.state, null);
var Q = cyc.Q;
var mCO2 = Q.prodC * MW_C;

var CCO2air = M.yAir.c * u.Pamb / (R * u.Tamb);
var Vair = (u.Vfan / 3600) * u.t_ads;
var etaCap = Q.prodC / (CCO2air * Vair);

var H_CO2 = M.dHc * dqc * M.ms / mCO2 / 1e6;    // MJ/kg
var H_H2O = M.dHw * dqw * M.ms / mCO2 / 1e6;
var E_th = Q.Qheat / mCO2 / 1e6;
var NoutHeat = Q3.prodC + Q3.prodW + Q3.prodN;

var out = {
  P: {
    Cps: M.cps / 1000, mstr: 0, Cpstr: 0, CpCO2ads: 0, CpH2Oads: 0,
    dT: dT, dqCO2: dqc, dHCO2: M.dHc / 1000, dqH2O: dqw, dHH2O: M.dHw / 1000,
    Patm: u.Pamb, dPair: u.dPfan, etaBlow: u.etaFan,
    CCO2air: CCO2air, etaCap: etaCap, COP: 1,
    Tamb: u.Tamb, etaVac: u.etaVac, Pregen: u.Pregen
  },
  dyn: {
    E_th: E_th, E_sens: E_th - H_CO2 - H_H2O, E_CO2: H_CO2, E_H2O: H_H2O,
    E_fan: Q.Wfan / mCO2 / 1e6,
    E_vac_heat: Q3.Wvac / (Q3.prodC * MW_C) / 1e6,
    E_vac_total: Q.Wvac / mCO2 / 1e6,
    dqc: dqc, dqw: dqw, mCO2: mCO2, prodC: Q.prodC,
    NoutHeat: NoutHeat, msdq: M.ms * (dqc + dqw),
    purity: res.kpi.purityDry, recovery: res.kpi.recovery,
    cycles: res.cycles
  }
};
print(JSON.stringify(out, null, 1));
