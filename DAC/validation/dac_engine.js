"use strict";
/* ===================================================================== *
 *  DAC 0-D cycle engine  (Glaser-style dynamics + Rajendran-group
 *  isotherms and Simplified-DAC energy forms)
 *
 *  Basis: 1 m3 of packed contactor (bed).  Solid volume = (1-eb) m3.
 *  4-step TVSA cycle:
 *    1. adsorption      fan-driven CSTR at P_amb, duration t_ads
 *    2. evacuation      closed inlet, pump-speed input; ends at P_regen
 *    3. heat + vacuum   P held at P_regen, wall -> T_des (1st order),
 *                       outflow = product; duration t_heat
 *    4. cooling/repress instantaneous repressurisation with ambient air,
 *                       then constant-P cooling with inflow only, t_cool
 *
 *  State vector x = [N, yc, yw, qc, qw, T, Tw]
 *    N  total gas moles in the bed voids      [mol]
 *    yc, yw  CO2 / H2O gas mole fractions     [-]   (yn = 1-yc-yw)
 *    qc, qw  CO2 / H2O loadings               [mol/kg sorbent]
 *    T   bed temperature (gas+solid, lumped)  [K]
 *    Tw  wall temperature                     [K]
 *
 *  Isotherms:
 *    CO2: Toth, modified for water via q_H2O (Stampi-Bombelli):
 *         qs_wet = qs/(1 - psi*qw),  b_wet = b*(1 + beta*qw)
 *    H2O: GAB on relative-humidity fraction x (T enters via psat only)
 *  Both reduce to the published dry forms at qw = 0 / x = 0.
 *
 *  Energy forms follow the Simplified DAC model (dac.energyModel):
 *    fan: isentropic blower work over dP; vacuum: isothermal ln(P) work
 *    at T_amb (default) or Glaser's adiabatic form (toggle).
 * ===================================================================== */
var DacEngine = (function () {

  var R = 8.314;                    // J/mol/K
  var TREF_LDF = 293.15;            // [K] reference T for the LDF sliders (20 C)
  var GAM = 1.4;                    // adiabatic constant, air
  var MW_C = 0.04401;               // kg/mol CO2
  var MW_W = 0.018015;              // kg/mol H2O
  var CPG = { c: 37.1, w: 33.6, n: 29.1 };  // gas molar cp [J/mol/K]
  var CPA_C = 50.0, CPA_W = 75.0;   // adsorbed-phase cp [J/mol/K]
  var Y_CO2_AIR = 400e-6;           // 400 ppm CO2 in air

  // Antoine (paper SI): log10 P[bar] = A - B/(T+C)
  function psat(T) { return Math.pow(10, 4.6543 - 1435.264 / (T - 64.848)) * 1e5; } // Pa

  /* ------------------------------------------------------------------ *
   *  Sorbent presets.
   *  Toth: qs0 [mol/kg], chi [-], T0 [K], b0 [1/Pa], dH0 [J/mol],
   *        w0 [-], alpha [-]
   *  GAB:  Cm [mol/kg], CG [-], K [-]   (null -> H2O inert)
   *  wet:  psi, beta [kg/mol == g/mmol],  xmax = RH-fraction clamp for
   *        the coupling (SI: NbOFFIVE 0.75, SIFSIX 0.90)
   *  dHc, dHw [J/mol] heats used in the energy balance & thermal KPI
   *  rho_p particle density [kg/m3], eb bed voidage, ep particle voidage
   *  cps sorbent heat capacity [J/kg/K], kC/kW default LDF [1/s]
   * ------------------------------------------------------------------ */
  var PRESETS = {
    "TMCM-41": {
      toth: { qs0: 2.7, chi: 0, T0: 298, b0: 0.38776, dH0: 90000, w0: 0.3123, alpha: 0.3472 },
      gab: null, wet: null, dHc: 90000, dHw: 0,
      rho_p: 930, eb: 0.40, ep: 0.40, cps: 1000, kC: 2e-4, kW: 0.2,
      label: "TRI-PE-MCM-41", cls: "chemisorbent (amine silica)"
    },
    "TPMS": {
      toth: { qs0: 1.0, chi: 0, T0: 298, b0: 2.0123, dH0: 90000, w0: 0.2696, alpha: 0.08069 },
      gab: null, wet: null, dHc: 90000, dHw: 0,
      rho_p: 1250, eb: 0.40, ep: 0.40, cps: 900, kC: 2e-4, kW: 0.2,
      label: "SI-AEATPMS", cls: "chemisorbent (amine silica)"
    },
    "NbOFFIVE": {
      toth: { qs0: 2.22, chi: 0, T0: 273, b0: 0.17567, dH0: 50000, w0: 1.166, alpha: -0.4937 },
      gab: { Cm: 6.159, CG: 36.81, K: 0.4063 },
      wet: { psi: 0, beta: -0.11707, xmax: 0.75 }, dHc: 50000, dHw: 45036,
      rho_p: 1173.6, eb: 0.40, ep: 0.40, cps: 1000, kC: 2e-4, kW: 0.2,
      label: "NbOFFIVE-1-Ni", cls: "physisorbent (MOF)"
    },
    "SIFSIX": {
      toth: { qs0: 3.0, chi: 0, T0: 273, b0: 0.04572, dH0: 52000, w0: 0.5455, alpha: 1.592 },
      gab: { Cm: 0.4067, CG: 93.13, K: 0.7859 },
      wet: { psi: 0, beta: -0.6772, xmax: 0.90 }, dHc: 52000, dHw: 45036,
      rho_p: 786, eb: 0.40, ep: 0.40, cps: 1000, kC: 2e-4, kW: 0.2,
      label: "SIFSIX-18-Ni-β", cls: "physisorbent (MOF)"
    },
    "APDES": {
      toth: { qs0: 2.2, chi: 0, T0: 296, b0: 0.373, dH0: 60000, w0: 0.4247, alpha: -0.4921 },
      gab: { Cm: 36.48, CG: 0.1489, K: 0.5751 },
      wet: { psi: 0.00958, beta: 3.448, xmax: 1.0 }, dHc: 60000, dHw: 43800,
      rho_p: 61, eb: 0.092, ep: 0.9614, cps: 2070, kC: 2e-4, kW: 0.2,
      label: "APDES-NFC-FD-S", cls: "chemisorbent (amine cellulose)"
    },
    "Lewatit": {
      toth: { qs0: 3.93873034, chi: 0, T0: 303.15, b0: 20112.4573 / 1e5, dH0: 70999, w0: 0.26906407, alpha: 0.04501535 },
      gab: { Cm: 3.8640, CG: 1.7950, K: 0.76830 },
      wet: { psi: 4.70e-3, beta: 0.569, xmax: 1.0 }, dHc: 70999, dHw: 43800,
      rho_p: 760, eb: 0.37, ep: 0.35, cps: 1580, kC: 4e-3, kW: 0.2,
      label: "Lewatit VP OC 1065", cls: "chemisorbent (amine resin)"
    }
  };

  /* --------------------------- isotherms ---------------------------- */

  // GAB water loading [mol/kg] from RH fraction x (clamped to 0.999)
  function qstarW(x, S) {
    if (!S.gab) return 0;
    var g = S.gab, a = Math.min(Math.max(x, 0), 0.999);
    var Ka = g.K * a;
    return g.Cm * g.CG * Ka / ((1 - Ka) * (1 + (g.CG - 1) * Ka));
  }

  // Arrhenius LDF coefficient: k(T) = k_ref exp[-Ea/R (1/T - 1/Tref)].
  // k_ref is the slider value at TREF_LDF (20 C); Ea in J/mol, 0/absent
  // reproduces the constant-k behaviour bit-for-bit.
  function kLDF(kref, Ea, T) {
    return Ea ? kref * Math.exp(-(Ea / R) * (1 / T - 1 / TREF_LDF)) : kref;
  }

  // Modified-Toth CO2 loading [mol/kg]; p in Pa, qw in mol/kg.
  // Water coupling is EITHER S.wet (Stampi-Bombelli modulation of qs and
  // b by the water LOADING qw) OR S.cmp (extended/competitive Toth: the
  // water PARTIAL PRESSURE pw [Pa] joins the shared denominator,
  // q* = qs b p / [1 + (b p + b_w pw)^t]^(1/t), b_w Arrhenius on the
  // Toth T0). pw is only read by the cmp branch and may be omitted.
  function qstarC(p, T, qw, S, pw) {
    var t = S.toth;
    var qs = t.qs0 * Math.exp(t.chi * (1 - T / t.T0));
    var b = t.b0 * Math.exp((t.dH0 / (R * t.T0)) * (t.T0 / T - 1));
    var w = t.w0 + t.alpha * (1 - t.T0 / T);
    if (S.wet && qw > 0) {
      var qwEff = S.wet.qwClamp != null ? Math.min(qw, S.wet.qwClamp) : qw;
      var den = 1 - S.wet.psi * qwEff;
      qs = den > 0.05 ? qs / den : qs / 0.05;      // guard qs_wet > 0
      b = Math.max(b * (1 + S.wet.beta * qwEff), 0);
    }
    var f = b * Math.max(p, 0);
    if (f <= 0) return 0;
    var g = f;
    if (S.cmp && S.cmp.bw0 > 0) {
      var bw = S.cmp.bw0 * Math.exp((S.cmp.dHw / (R * t.T0)) * (t.T0 / T - 1));
      g = f + bw * Math.max(pw || 0, 0);
    }
    return qs * f / Math.pow(1 + Math.pow(g, w), 1 / w);
  }

  /* ------------------------ parameter assembly ----------------------- */
  // User params (all per m3 of bed):
  //   sorbent, Tamb [K], RH [%], Tdes [K], Pregen [Pa], Pamb [Pa],
  //   t_ads, t_heat, t_cool [s], Vfan [m3/h], Spump [m3/h],
  //   dPfan [Pa], etaFan, etaVac, pumpModel 'iso'|'adia',
  //   tauHeat, tauCool [s], Ua [W/m3/K], wallCp [J/K], kC, kW [1/s at 20 C],
  //   EaC, EaW [J/mol] optional Arrhenius activation energies for the LDF ks
  function buildModel(u) {
    // "Custom" carries its full isotherm/property spec in u.iso — the
    // same shape as a PRESETS entry (toth, gab, wet, cmp, eb, ep, rho_p,
    // cps, dHc, dHw), seeded by the UI from a source preset.
    var S = (u.sorbent === "Custom" && u.iso) ? u.iso : PRESETS[u.sorbent];
    if (!S) throw new Error("unknown sorbent " + u.sorbent);
    var pick = function (v, d) { return (v != null && isFinite(v)) ? v : d; };
    // sorbent properties may be overridden by the user (defaults = preset)
    var eb = pick(u.eb, S.eb), ep = pick(u.ep, S.ep);
    var rho_p = pick(u.rhop, S.rho_p), cps = pick(u.cps, S.cps);
    var dHc = pick(u.dHc, S.dHc);
    var dHw = S.gab ? pick(u.dHw, S.dHw) : 0;
    var epsT = eb + (1 - eb) * ep;              // total voidage
    var ms = rho_p * (1 - eb);                  // kg sorbent per m3 bed
    var VC = epsT * 1.0;                        // gas volume [m3]
    var Vsolid = (1 - eb) * 1.0;                // sorbent volume [m3]
    // The Toth temperature coefficient IS the heat of adsorption
    // (Clausius-Clapeyron), so an edited dHc moves the isotherm too unless
    // the user explicitly unlinks the two to test sensitivity.
    var toth = S.toth;
    if (!u.unlinkDH && dHc !== S.dHc) {
      toth = { qs0: toth.qs0, chi: toth.chi, T0: toth.T0, b0: toth.b0,
               dH0: dHc, w0: toth.w0, alpha: toth.alpha };
    }
    // clamp loading used in the CO2-H2O coupling (SI RH clamps)
    var wet = S.wet ? {
      psi: S.wet.psi, beta: S.wet.beta,
      qwClamp: (S.wet.xmax < 1 && S.gab) ? qstarW(S.wet.xmax, S) : null
    } : null;
    // competitive (shared-denominator) coupling: b_w0 from the spec, its
    // Arrhenius heat from the ΔH(H2O) box — works with inert water too,
    // since the competition acts through the gas-phase pw
    var cmp = (S.cmp && S.cmp.bw0 > 0) ? {
      bw0: S.cmp.bw0, dHw: pick(u.dHw, S.dHw || 0)
    } : null;
    var Seff = { toth: toth, gab: S.gab, wet: wet, cmp: cmp };
    // ambient air
    var pw = (u.RH / 100) * psat(u.Tamb);
    var ywA = pw / u.Pamb;
    var yA = { c: Y_CO2_AIR * (1 - ywA), w: ywA, n: 0 };
    yA.n = 1 - yA.c - yA.w;
    return {
      u: u, S: Seff, dHc: dHc, dHw: dHw,
      eb: eb, epsT: epsT, ms: ms, VC: VC, Vsolid: Vsolid,
      cps: cps, yAir: yA,
      cpaC: (u.cpaC != null) ? u.cpaC : CPA_C,
      cpaW: (u.cpaW != null) ? u.cpaW : CPA_W,
      // structure carried along with the sorbent: mass ratio m_str/m_sorb
      // times its cp [J/kg/K] (same term as dac.energyModel's mstr*Cpstr)
      cpStr: (u.mstr || 0) * (u.cpstr || 0),
      ndotFan: u.Pamb * (u.Vfan / 3600) / (R * u.Tamb),  // mol/s at ambient
      inertW: !S.gab
    };
  }

  /* ---------------------- step right-hand sides ---------------------- *
   * x = [N, yc, yw, qc, qw, T, Tw];  returns dx/dt.
   * mode: "ads" | "evac" | "heat" | "cool"
   * Also fills out.rates = {ndotO, ndotI, Qext, Wvac, Wfan} for
   * quadrature (evaluated at accepted solution points).
   * ------------------------------------------------------------------ */
  function rhs(M, mode, x, out) {
    var u = M.u;
    var N = x[0], yc = x[1], yw = x[2], qc = x[3], qw = x[4], T = x[5], Tw = x[6];
    var yn = 1 - yc - yw;
    var P = N * R * T / M.VC;

    // equilibria & LDF rates (mol/s per m3 bed)
    var xRH = Math.max(yw, 0) * P / psat(T);
    var qwS = M.inertW ? 0 : qstarW(xRH, M.S);
    var qcS = qstarC(Math.max(yc, 0) * P, T, Math.max(qw, 0), M.S,
                     Math.max(yw, 0) * P);
    var rc = kLDF(u.kC, u.EaC, T) * (qcS - qc) * M.ms;
    var rw = M.inertW ? 0 : kLDF(u.kW, u.EaW, T) * (qwS - qw) * M.ms;
    var dqc = rc / M.ms, dqw = M.inertW ? 0 : rw / M.ms;
    var rSum = rc + rw;

    // heat capacity of the lumped bed [J/K per m3 bed]
    var Cgas = N * (yc * CPG.c + yw * CPG.w + yn * CPG.n);
    var Ctot = M.ms * (M.cps + M.cpStr) + Cgas + M.ms * (qc * M.cpaC + qw * M.cpaW);

    var Qext = u.Ua * (Tw - T);       // W per m3 bed (wall -> bed)
    var Hads = M.dHc * rc + M.dHw * rw;  // W released by adsorption

    var ndotI = 0, ndotO = 0, dN, dT, dyc, dyw, dTw;

    if (mode === "heat") {
      dTw = (u.Tdes - Tw) / u.tauHeat;
    } else {
      dTw = (u.Tamb - Tw) / u.tauCool;
    }

    if (mode === "ads" || mode === "cool") {
      // constant pressure: dN/dt = -(N/T) dT/dt
      if (mode === "ads") {
        ndotI = M.ndotFan;
        // energy: Ctot*dT = ndotI*cpg*(Tamb-T) + Hads + Qext
        var cpgI = CPG.c * M.yAir.c + CPG.w * M.yAir.w + CPG.n * M.yAir.n;
        dT = (ndotI * cpgI * (u.Tamb - T) + Hads + Qext) / Ctot;
        dN = -(N / T) * dT;
        ndotO = ndotI - rSum - dN;
        dyc = (ndotI * (M.yAir.c - yc) - rc + yc * rSum) / N;
        dyw = (ndotI * (M.yAir.w - yw) - rw + yw * rSum) / N;
      } else {
        // cooling: inflow only (ndotO = 0), inflow at ambient composition
        // dT appears in ndotI = rSum + dN/dt = rSum - (N/T) dT  -> linear solve
        var cpgI2 = CPG.c * M.yAir.c + CPG.w * M.yAir.w + CPG.n * M.yAir.n;
        // Ctot*dT = ndotI*cpgI*(Tamb-T) + Hads + Qext
        //         = (rSum - (N/T)dT)*cpgI*(Tamb-T) + Hads + Qext
        var A = Ctot + (N / T) * cpgI2 * (u.Tamb - T);
        dT = (rSum * cpgI2 * (u.Tamb - T) + Hads + Qext) / A;
        dN = -(N / T) * dT;
        ndotI = rSum + dN;
        if (ndotI < 0) {  // gas would need to leave: vent excess instead
          ndotI = 0;
          dT = (Hads + Qext) / Ctot;
          dN = -(N / T) * dT;
          ndotO = -(dN + rSum) > 0 ? -(dN + rSum) : 0;
        }
        dyc = (ndotI * (M.yAir.c - yc) - rc + yc * rSum) / N;
        dyw = (ndotI * (M.yAir.w - yw) - rw + yw * rSum) / N;
      }
    } else if (mode === "evac") {
      // pump withdraws Spump m3/h at column conditions; V*dP/dt work term
      ndotO = (u.Spump / 3600) * N / M.VC;
      dN = -ndotO - rSum;
      // energy with gas-expansion term V dP/dt = R(T dN + N dT)
      // Ctot*dT = Hads + Qext + R*(T*dN + N*dT)  ->  (Ctot - R N) dT = ...
      dT = (Hads + Qext + R * T * dN) / (Ctot - R * N);
      // species: d(N yi)/dt = -ndotO yi - ri  ->  N dyi = -ri + yi*rSum
      dyc = (-rc + yc * rSum) / N;
      dyw = (-rw + yw * rSum) / N;
    } else { // heat: constant P = Pregen, outflow only
      // dN = -(N/T) dT ; ndotO = -rSum - dN
      // Ctot*dT = Hads + Qext  (no inflow; P const -> no pdV term at const P
      //           beyond enthalpy of leaving gas at T which cancels for CSTR)
      dT = (Hads + Qext) / Ctot;
      dN = -(N / T) * dT;
      ndotO = -rSum - dN;
      if (ndotO < 0) ndotO = 0;   // no reverse flow through the pump
      dyc = (-rc + yc * rSum) / N;
      dyw = (-rw + yw * rSum) / N;
    }

    if (out) {
      out.ndotI = ndotI; out.ndotO = ndotO; out.P = P;
      out.Qext = Qext; out.rc = rc; out.rw = rw;
      out.yc = yc; out.yw = yw; out.T = T;
    }
    return [dN, dyc, dyw, dqc, dqw, dT, dTw];
  }

  /* ------------------- TR-BDF2 implicit integrator ------------------- */
  var TRGAM = 2 - Math.SQRT2, TRD = TRGAM / 2;

  function solve7(Aa, bb) {  // Gaussian elim with partial pivoting, n=7
    var n = 7, i, j, k;
    var A = [], b = bb.slice();
    for (i = 0; i < n; i++) A.push(Aa[i].slice());
    for (k = 0; k < n; k++) {
      var p = k, mx = Math.abs(A[k][k]);
      for (i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > mx) { mx = Math.abs(A[i][k]); p = i; }
      if (p !== k) { var t = A[p]; A[p] = A[k]; A[k] = t; var tb = b[p]; b[p] = b[k]; b[k] = tb; }
      if (A[k][k] === 0) return null;
      for (i = k + 1; i < n; i++) {
        var m = A[i][k] / A[k][k];
        if (m === 0) continue;
        for (j = k; j < n; j++) A[i][j] -= m * A[k][j];
        b[i] -= m * b[k];
      }
    }
    var xx = new Array(n);
    for (i = n - 1; i >= 0; i--) {
      var s = b[i];
      for (j = i + 1; j < n; j++) s -= A[i][j] * xx[j];
      xx[i] = s / A[i][i];
    }
    return xx;
  }

  function jacFD(M, mode, x, f0) {
    var n = 7, J = [];
    for (var i = 0; i < n; i++) J.push(new Array(n));
    for (var j = 0; j < n; j++) {
      var h = 1e-7 * Math.max(Math.abs(x[j]), 1e-6);
      var xp = x.slice(); xp[j] += h;
      var fp = rhs(M, mode, xp, null);
      for (var i2 = 0; i2 < n; i2++) J[i2][j] = (fp[i2] - f0[i2]) / h;
    }
    return J;
  }

  // solve z - c*h*f(z) = rvec  by Newton, initial guess zg
  function newtonStage(M, mode, ch, rvec, zg) {
    var z = zg.slice();
    for (var it = 0; it < 12; it++) {
      var f = rhs(M, mode, z, null);
      var res = new Array(7), nrm = 0;
      for (var i = 0; i < 7; i++) {
        res[i] = z[i] - ch * f[i] - rvec[i];
        var sc = Math.max(Math.abs(z[i]), 1e-8);
        nrm = Math.max(nrm, Math.abs(res[i]) / sc);
      }
      if (nrm < 1e-12) return z;
      var J = jacFD(M, mode, z, f);
      var A = [];
      for (var r = 0; r < 7; r++) {
        A.push(new Array(7));
        for (var c = 0; c < 7; c++) A[r][c] = (r === c ? 1 : 0) - ch * J[r][c];
      }
      var dz = solve7(A, res);
      if (!dz) return null;
      var ok = true;
      for (var i3 = 0; i3 < 7; i3++) { z[i3] -= dz[i3]; if (!isFinite(z[i3])) ok = false; }
      if (!ok) return null;
      if (nrm < 1e-10 && it > 1) return z;
    }
    return z; // accept last iterate (residual checked via error control)
  }

  /* Integrate one step-mode from t0 to tend.
   * opts: { event: function(x, P) -> signed value (stop at zero cross,
   *         from positive to negative), record: fn(t, x, aux) }
   * Accumulates quadratures in Q = {feedC, feedW, feedN, outC, outW, outN,
   *   Wfan, Wvac, Qheat, QheatPos}
   */
  function integrateStep(M, mode, x0, t0, tend, Q, rec, eventFn) {
    var x = x0.slice(), t = t0;
    var h = Math.min(5, (tend - t0) / 20);
    var hmin = 1e-6, hmax = (tend - t0) / 4;
    var aux = {};
    var u = M.u;
    var quadLast = null;

    function auxRates(xx) {
      var o = {};
      rhs(M, mode, xx, o);
      // instantaneous powers
      // isentropic blower over dPfan (same form as dac.energyModel)  [W]
      o.Wfan = (mode === "ads")
        ? (GAM / (GAM - 1)) * (u.Pamb / u.etaFan) *
          (Math.pow((u.Pamb + u.dPfan) / u.Pamb, (GAM - 1) / GAM) - 1) * (u.Vfan / 3600)
        : 0;
      var Wv = 0;
      if ((mode === "evac" || mode === "heat") && o.ndotO > 0 && o.P < u.Pamb) {
        if (u.pumpModel === "adia") {
          Wv = (1 / u.etaVac) * (GAM / (GAM - 1)) * o.ndotO * R * xx[5] *
               (Math.pow(u.Pamb / o.P, (GAM - 1) / GAM) - 1);
        } else {
          Wv = o.ndotO * R * u.Tamb * Math.log(u.Pamb / o.P) / u.etaVac;
        }
      }
      o.Wvac = Wv;
      return o;
    }

    var steps = 0;
    while (t < tend - 1e-9) {
      if (t + h > tend) h = tend - t;
      var f0 = rhs(M, mode, x, null);
      // TR stage to t + gamma*h
      var r1 = new Array(7);
      for (var i = 0; i < 7; i++) r1[i] = x[i] + TRD * h * f0[i];
      var zg = new Array(7);
      for (var i1 = 0; i1 < 7; i1++) zg[i1] = x[i1] + TRGAM * h * f0[i1];
      var z = newtonStage(M, mode, TRD * h, r1, zg);
      if (!z) { h /= 4; if (h < hmin) throw new Error("step failure (TR)"); continue; }
      // BDF2 stage to t + h
      var c1 = 1 / (TRGAM * (2 - TRGAM)), c2 = (1 - TRGAM) * (1 - TRGAM) / (TRGAM * (2 - TRGAM));
      var r2 = new Array(7);
      for (var i4 = 0; i4 < 7; i4++) r2[i4] = c1 * z[i4] - c2 * x[i4];
      var x1 = newtonStage(M, mode, TRD * h, r2, z);
      if (!x1) { h /= 4; if (h < hmin) throw new Error("step failure (BDF2)"); continue; }

      // error estimate from 3-point divided difference of f
      var fz = rhs(M, mode, z, null), f1 = rhs(M, mode, x1, null);
      var err = 0;
      for (var i5 = 0; i5 < 7; i5++) {
        var dd = ((f1[i5] - fz[i5]) / ((1 - TRGAM) * h) - (fz[i5] - f0[i5]) / (TRGAM * h));
        var lte = 0.06 * h * h * Math.abs(dd);   // ~C*h^3*|x'''|
        var sc = 1e-10 + 1e-8 * Math.abs(x1[i5]) +
                 (i5 === 1 ? 1e-12 : i5 === 2 ? 1e-12 : i5 >= 5 ? 1e-7 : 1e-10);
        err = Math.max(err, lte / sc);
      }
      if (err > 1 && h > hmin * 2) { h = Math.max(h * Math.max(0.2, 0.9 * Math.pow(err, -1 / 3)), hmin); continue; }

      // event check (linear in x over the step)
      if (eventFn) {
        var e0 = eventFn(x), e1 = eventFn(x1);
        if (e0 > 0 && e1 <= 0) {
          // bisect on h for the crossing
          var lo = 0, hi = h;
          for (var bit = 0; bit < 40; bit++) {
            var mid = 0.5 * (lo + hi);
            var rm = new Array(7);
            for (var i6 = 0; i6 < 7; i6++) rm[i6] = x[i6] + TRD * mid * f0[i6];
            var zm = newtonStage(M, mode, TRD * mid, rm,
              x.map(function (v, k) { return v + TRGAM * mid * f0[k]; }));
            var r2m = new Array(7);
            for (var i7 = 0; i7 < 7; i7++) r2m[i7] = c1 * zm[i7] - c2 * x[i7];
            var xm = newtonStage(M, mode, TRD * mid, r2m, zm);
            if (eventFn(xm) <= 0) hi = mid; else { lo = mid; }
            if (hi - lo < 1e-8 * h + 1e-9) { x1 = xm; h = mid; break; }
            x1 = xm;
          }
          // accept the crossing step then stop
          acceptStep();
          return { x: x, t: t, stopped: true };
        }
      }

      acceptStep();

      function acceptStep() {
        // trapezoidal quadrature of flows/powers over [t, t+h]
        var a0 = quadLast || auxRates(x);
        var a1 = auxRates(x1);
        var w = 0.5 * h;
        if (Q) {
          Q.feedC += w * (a0.ndotI * M.yAir.c + a1.ndotI * M.yAir.c);
          Q.feedW += w * (a0.ndotI * M.yAir.w + a1.ndotI * M.yAir.w);
          Q.feedN += w * (a0.ndotI * M.yAir.n + a1.ndotI * M.yAir.n);
          var key = (mode === "heat") ? "prod" : "vent";
          Q[key + "C"] += w * (a0.ndotO * a0.yc + a1.ndotO * a1.yc);
          Q[key + "W"] += w * (a0.ndotO * a0.yw + a1.ndotO * a1.yw);
          Q[key + "N"] += w * (a0.ndotO * (1 - a0.yc - a0.yw) + a1.ndotO * (1 - a1.yc - a1.yw));
          Q.Wfan += w * (a0.Wfan + a1.Wfan);
          Q.Wvac += w * (a0.Wvac + a1.Wvac);
          if (mode === "heat") {
            Q.Qheat += w * (Math.max(a0.Qext, 0) + Math.max(a1.Qext, 0));
            // wall sensible while heating (prescribed Tw ODE)
            Q.Qwall += u.wallCp * Math.max(x1[6] - x[6], 0);
          }
        }
        quadLast = a1;
        x = x1; t += h;
        if (rec) rec(t, x, a1, mode);
        steps++;
        if (err > 0) {
          var fac = Math.min(3, Math.max(0.3, 0.85 * Math.pow(Math.max(err, 1e-8), -1 / 3)));
          h = Math.min(h * fac, hmax);
        }
      }
      if (steps > 200000) throw new Error("too many steps in " + mode);
    }
    return { x: x, t: t, stopped: false };
  }

  /* --------------------------- one cycle ----------------------------- */
  function newQuad() {
    return { feedC: 0, feedW: 0, feedN: 0, prodC: 0, prodW: 0, prodN: 0,
             ventC: 0, ventW: 0, ventN: 0, Wfan: 0, Wvac: 0, Qheat: 0, Qwall: 0 };
  }

  function runCycle(M, x0, rec) {
    var u = M.u, Q = newQuad();
    var t = 0;

    // 1) adsorption
    var s1 = integrateStep(M, "ads", x0, t, t + u.t_ads, Q, rec, null);
    t = s1.t;
    var tAdsEnd = t;

    // 2) evacuation: until P <= Pregen (or guard time)
    var evGuard = t + Math.max(4 * 3600, u.t_ads);
    var s2 = integrateStep(M, "evac", s1.x, t, evGuard, Q, rec, function (xx) {
      return xx[0] * R * xx[5] / M.VC - u.Pregen;   // + above Pregen
    });
    t = s2.t;
    var tEvacEnd = t, evacReached = s2.stopped;

    // clamp state exactly onto P = Pregen (event lands within tolerance)
    var x2 = s2.x.slice();
    x2[0] = u.Pregen * M.VC / (R * x2[5]);

    // 3) heat + vacuum at P = Pregen
    var s3 = integrateStep(M, "heat", x2, t, t + u.t_heat, Q, rec, null);
    t = s3.t;
    var tHeatEnd = t;

    // 4) instantaneous repressurisation with ambient air, then cooling
    var x3 = s3.x.slice();
    var Nadd = Math.max((u.Pamb - u.Pregen) * M.VC / (R * x3[5]), 0);
    var N3 = x3[0];
    var Nn = N3 + Nadd;
    var yc4 = (N3 * x3[1] + Nadd * M.yAir.c) / Nn;
    var yw4 = (N3 * x3[2] + Nadd * M.yAir.w) / Nn;
    Q.feedC += Nadd * M.yAir.c; Q.feedW += Nadd * M.yAir.w; Q.feedN += Nadd * M.yAir.n;
    var x4 = [Nn, yc4, yw4, x3[3], x3[4], x3[5], x3[6]];
    var s4 = integrateStep(M, "cool", x4, t, t + u.t_cool, Q, rec, null);
    t = s4.t;

    return { x: s4.x, Q: Q, t_cycle: t,
             times: { ads: tAdsEnd, evac: tEvacEnd - tAdsEnd, heat: tHeatEnd - tEvacEnd,
                      cool: t - tHeatEnd },
             // step-boundary states: richest (end of adsorption) and leanest
             // (end of the heating step) — the working capacity and the
             // temperature swing the simplified energy model asks for
             marks: { qcAds: s1.x[3], qwAds: s1.x[4], Tads: s1.x[5],
                      qcDes: s3.x[3], qwDes: s3.x[4], Tdes: s3.x[5] },
             evacReached: evacReached };
  }

  /* ------------------------ cycles to CSS ---------------------------- */
  function initialState(M) {
    var u = M.u;
    var N0 = u.Pamb * M.VC / (R * u.Tamb);
    var x = (u.RH / 100);
    var qw0 = M.inertW ? 0 : qstarW(x, M.S);
    var qc0 = qstarC(M.yAir.c * u.Pamb, u.Tamb, qw0, M.S, M.yAir.w * u.Pamb);
    return [N0, M.yAir.c, M.yAir.w, qc0, qw0, u.Tamb, u.Tamb];
  }

  function runToCSS(u, opts) {
    opts = opts || {};
    var M = buildModel(u);
    var x = opts.x0 || initialState(M);
    var tolCSS = opts.tolCSS || 1e-5, maxCyc = opts.maxCycles || 60;
    var cyc = null, nCyc = 0, delta = Infinity;
    var scale = [1, 1e-4, 1e-2, 1, 1, 100, 100];
    for (var k = 0; k < maxCyc; k++) {
      var rec = (opts.recordLast && k === maxCyc - 1) ? opts.recordLast : null;
      cyc = runCycle(M, x, rec);
      nCyc++;
      delta = 0;
      for (var i = 0; i < 7; i++) {
        delta = Math.max(delta, Math.abs(cyc.x[i] - x[i]) /
          (Math.abs(x[i]) + scale[i] * 1e-3));
      }
      x = cyc.x;
      if (delta < tolCSS) break;
    }
    // one more recorded cycle at CSS for profiles/KPIs
    var profile = [];
    var recFn = opts.record !== false ? function (t, xx, a, mode) {
      profile.push({ t: t, P: a.P, T: xx[5], Tw: xx[6], yc: a.yc, yw: a.yw,
                     qc: xx[3], qw: xx[4], mode: mode });
    } : null;
    cyc = runCycle(M, x, recFn);
    var kpi = computeKPIs(M, cyc);
    return { M: M, cycles: nCyc, cssDelta: delta, state: cyc.x, cyc: cyc,
             kpi: kpi, profile: profile };
  }

  /* ------------- Simplified-DAC energy model (comparison) ------------ *
   * The forms of dac.energyModel, evaluated with THIS cycle's own working
   * capacity, temperature swing and capture fraction, so the comparison
   * isolates the model-form differences only:
   *   - the void-gas inventory (the simplified model neglects it),
   *   - sensible heat swept out with the desorbing gas,
   *   - CO2 released during evacuation instead of in the product step.
   * The wall thermal mass has no analogue in the simplified model; it is
   * added to the sensible term so a non-zero wall does not by itself make
   * the two disagree.
   * ------------------------------------------------------------------- */
  function simplifiedEnergy(M, cyc) {
    var u = M.u, mk = cyc.marks;
    var dqc = mk.qcAds - mk.qcDes;              // working capacity [mol/kg]
    var dqw = mk.qwAds - mk.qwDes;
    var dT = mk.Tdes - mk.Tads;                 // actual bed swing [K]
    if (!(dqc > 1e-9) || !(dT > 1e-9)) return null;
    // sensible: sorbent + structure + adsorbed phases  [J/mol CO2]
    var thermalMass = M.cps + M.cpStr + dqc * M.cpaC + dqw * M.cpaW;  // J/kg/K
    var Esens = thermalMass * dT / dqc;
    if (u.wallCp > 0) Esens += u.wallCp * dT / (M.ms * dqc);
    var Ec = M.dHc, Ew = (dqw / dqc) * M.dHw;   // desorption [J/mol CO2]
    // blower: same isentropic work per m3 air, over the CO2 actually captured
    var wAir = (GAM / (GAM - 1)) * (u.Pamb / u.etaFan) *
      (Math.pow((u.Pamb + u.dPfan) / u.Pamb, (GAM - 1) / GAM) - 1);   // J/m3
    var cCO2air = M.yAir.c * u.Pamb / (R * u.Tamb);                   // mol/m3
    var nBlown = M.ndotFan * M.yAir.c * u.t_ads;   // mol CO2 offered on ads
    var etaCap = cyc.Q.prodC / Math.max(nBlown, 1e-30);
    var Eblow = wAir / (cCO2air * Math.max(etaCap, 1e-12));
    // vacuum: isothermal, CO2 + co-desorbed water (void gas neglected)
    var Evac = (R * u.Tamb / u.etaVac) * Math.log(u.Pamb / u.Pregen) * (1 + dqw / dqc);
    var conv = 1 / (1e6 * MW_C);                // J/mol CO2 -> MJ/kg CO2
    return {
      Esens: Esens * conv, EdesC: Ec * conv, EdesW: Ew * conv,
      Eth: (Esens + Ec + Ew) * conv,
      Efan: Eblow * conv, Evac: Evac * conv, Eel: (Eblow + Evac) * conv,
      Etot: (Esens + Ec + Ew + Eblow + Evac) * conv,
      dqC: dqc, dqW: dqw, dT: dT, etaCap: etaCap * 100
    };
  }

  /* ----------------------------- KPIs -------------------------------- */
  function computeKPIs(M, cyc) {
    var Q = cyc.Q, u = M.u;
    var mCO2 = Q.prodC * MW_C;                     // kg per cycle per m3 bed
    var tc = cyc.t_cycle;
    var puDry = Q.prodC / Math.max(Q.prodC + Q.prodN, 1e-30);
    var re = Q.prodC / Math.max(Q.feedC, 1e-30);
    var Eth = (Q.Qheat + Q.Qwall) / Math.max(mCO2, 1e-30) / 1e6;   // MJ/kg
    var Efan = Q.Wfan / Math.max(mCO2, 1e-30) / 1e6;
    var Evac = Q.Wvac / Math.max(mCO2, 1e-30) / 1e6;
    var prContact = mCO2 / tc * 3600;              // kg/m3 bed/h
    var prSorb = prContact / (1 - M.eb);           // kg/m3 sorbent/h
    return {
      purityDry: puDry * 100, recovery: re * 100,
      Eth: Eth, Efan: Efan, Evac: Evac, Eel: Efan + Evac, Etot: Eth + Efan + Evac,
      prContact: prContact, prSorb: prSorb,
      prContactTPD: prContact * 24 / 1000, prSorbTPD: prSorb * 24 / 1000,
      waterPerCO2: Q.prodW / Math.max(Q.prodC, 1e-30),
      simp: simplifiedEnergy(M, cyc),
      mCO2cycle: mCO2, t_cycle: tc, times: cyc.times,
      evacReached: cyc.evacReached,
      prodC: Q.prodC, prodW: Q.prodW, prodN: Q.prodN,
      feedC: Q.feedC, ventC: Q.ventC,
      Qheat: Q.Qheat, Qwall: Q.Qwall, Wfan: Q.Wfan, Wvac: Q.Wvac
    };
  }

  /* ------------------- mole-balance closure check -------------------- */
  // over a CSS cycle: feed - vent - prod - d(inventory) = 0 per species
  function moleBalance(M, x0, cyc) {
    var Q = cyc.Q, x1 = cyc.x;
    var invC0 = x0[0] * x0[1] + x0[3] * M.ms, invC1 = x1[0] * x1[1] + x1[3] * M.ms;
    var invW0 = x0[0] * x0[2] + x0[4] * M.ms, invW1 = x1[0] * x1[2] + x1[4] * M.ms;
    var nb0 = x0[0] * (1 - x0[1] - x0[2]), nb1 = x1[0] * (1 - x1[1] - x1[2]);
    return {
      C: (Q.feedC - Q.ventC - Q.prodC - (invC1 - invC0)) / Math.max(Q.feedC, 1e-30),
      W: (Q.feedW - Q.ventW - Q.prodW - (invW1 - invW0)) / Math.max(Q.feedW, 1e-30),
      N: (Q.feedN - Q.ventN - Q.prodN - (nb1 - nb0)) / Math.max(Q.feedN, 1e-30)
    };
  }

  return {
    R: R, PRESETS: PRESETS, psat: psat,
    qstarW: qstarW, qstarC: qstarC, kLDF: kLDF,
    buildModel: buildModel, initialState: initialState,
    runCycle: runCycle, runToCSS: runToCSS,
    computeKPIs: computeKPIs, simplifiedEnergy: simplifiedEnergy,
    moleBalance: moleBalance,
    integrateStep: integrateStep, rhs: rhs, newQuad: newQuad
  };
})();
if (typeof module !== "undefined") module.exports = DacEngine;
