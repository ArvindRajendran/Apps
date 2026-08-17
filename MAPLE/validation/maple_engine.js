/* MAPLE engine — ANN surrogate of a 4-step PVSA process (LPP / FP cycles)
 * Ports the trained networks of github.com/ArvindRajendran/MAPLE
 * (Pai, Prasad & Rajendran, I&EC Res 2020; ACS Sustain. Chem. Eng. 2021)
 * to JavaScript. Pure JS, no DOM — runs in both the page and a Worker.
 *
 * Net roles (as in the repo's MAPLE.m wrapper, and CONFIRMED against the
 * published optimum points of the Limits paper SI, where the purity
 * column is pinned at the active >=95% constraint — unambiguous):
 *   net1 = CO2 recovery [%]        net2 = CO2 purity [%]
 *   net3 = log10(energy [kWh/t])   net4 = productivity [mol/m3 ads/s]
 *   net5 = log10(energy at 100% pump efficiency)
 * NOTE the repo's LHS-samples xlsx has its Pu/Re column HEADERS swapped
 * (its "Pu" column is recovery and vice versa) — see validation README.
 *
 * Physical inputs (object u):
 *   qsat [mol/kg], b0c, b0n [m3/mol], dUc, dUn [kJ/mol, NEGATIVE],
 *   rho [kg/m3], y [-], tads [s], PH [bar], PI [bar], PL [bar], vF [m/s]
 * Net input vector (before z-score):
 *   [qsat*rho, log10 b_CO2(298K), log10 b_N2(298K), dUc, dUn, rho, y,
 *    tads, PH, log10(PI/PH), log10(PL/PH), vF]
 */
"use strict";
var MAPLE = (function () {
  var Rg = 8.314, TREF = 298.0;

  var W = (typeof MAPLE_WEIGHTS !== "undefined") ? MAPLE_WEIGHTS : null;
  function setWeights(w) { W = w; }

  /* ---------- network forward pass ---------- */
  function mmSettings(p) {
    var s = p.settings;
    return (s && s.values) ? s.values : s;
  }
  function applyProc(p, x) {
    var s = mmSettings(p);
    if (p.fcn === "mapminmax") {
      var y = new Array(x.length);
      var g = s.gain, o = s.xoffset, ymin = s.ymin;
      var gv = Array.isArray(g), ov = Array.isArray(o);
      for (var i = 0; i < x.length; i++)
        y[i] = (x[i] - (ov ? o[i] : o)) * (gv ? g[i] : g) + ymin;
      return y;
    }
    if (p.fcn === "removeconstantrows") {
      if (s && s.keep && s.keep.length && s.keep.length < x.length) {
        var k = new Array(s.keep.length);
        for (var j = 0; j < s.keep.length; j++) k[j] = x[s.keep[j] - 1];
        return k;
      }
      return x;
    }
    throw new Error("unknown processFcn " + p.fcn);
  }
  function revProc(p, v) {
    var s = mmSettings(p);
    if (p.fcn === "mapminmax")
      return (v - s.ymin) / s.gain + s.xoffset;
    if (p.fcn === "removeconstantrows") return v;
    throw new Error("unknown processFcn " + p.fcn);
  }
  function fwd(net, z) {                       // z: 12-vector (z-scored)
    var x = z;
    for (var m = 0; m < net.inProc.length; m++) x = applyProc(net.inProc[m], x);
    var a = x;
    for (var l = 0; l < net.layers.length; l++) {
      var L = net.layers[l], Wl = L.W, bl = L.b;
      if (!Array.isArray(Wl[0])) Wl = [Wl];    // 1xN row collapsed by jsonencode
      if (!Array.isArray(bl)) bl = [bl];       // 1x1 bias collapsed to scalar
      var n = Wl.length;
      var out = new Array(n);
      for (var i = 0; i < n; i++) {
        var row = Wl[i], sum = Array.isArray(bl[i]) ? bl[i][0] : bl[i];
        for (var j = 0; j < a.length; j++) sum += row[j] * a[j];
        out[i] = (L.transfer === "tansig") ? (2 / (1 + Math.exp(-2 * sum)) - 1)
               : (L.transfer === "purelin") ? sum
               : (L.transfer === "logsig") ? (1 / (1 + Math.exp(-sum)))
               : NaN;
      }
      a = out;
    }
    var v = a[0];
    for (var q = net.outProc.length - 1; q >= 0; q--) v = revProc(net.outProc[q], v);
    return v;
  }

  /* ---------- physical -> net inputs ---------- */
  function b298(b0, dU) { return b0 * Math.exp(-dU * 1000 / (Rg * TREF)); }
  function pack(u) {
    return [
      u.qsat * u.rho,
      Math.log10(b298(u.b0c, u.dUc)),
      Math.log10(b298(u.b0n, u.dUn)),
      u.dUc, u.dUn, u.rho, u.y, u.tads, u.PH,
      Math.log10(u.PI / u.PH), Math.log10(u.PL / u.PH), u.vF,
    ];
  }
  function evalKPI(cyc, u) {
    var C = W[cyc], c = pack(u), z = new Array(12);
    for (var i = 0; i < 12; i++) z[i] = (c[i] - C.mue[i]) / C.sig[i];
    return {
      Pu: fwd(C.net2, z), Re: fwd(C.net1, z),
      En: Math.pow(10, fwd(C.net3, z)), Prod: fwd(C.net4, z),
      En100: Math.pow(10, fwd(C.net5, z)),
    };
  }

  /* ---------- training-domain guards (data-driven box) ----------
   * Bounds are the min/max of the transformed inputs over the ~48k
   * labelled LHS rows shipped with the repo (union of both cycles),
   * plus the sampling constraints of the Limits paper. Outside this box
   * the surrogate extrapolates and its predictions are unreliable. */
  var BOX = [
    { k: "qrho",  lo: 404,     hi: 11382,   lbl: "q_sat·ρ [mol/m³]" },
    { k: "lbc",   lo: -4.998,  hi: 1.004,   lbl: "log₁₀ b_CO2(298K)" },
    { k: "lbn",   lo: -7.002,  hi: -0.0002, lbl: "log₁₀ b_N2(298K)" },
    { k: "dUc",   lo: -48.0,   hi: -7.0,    lbl: "ΔU_CO2 [kJ/mol]" },
    { k: "dUn",   lo: -25.0,   hi: -3.0,    lbl: "ΔU_N2 [kJ/mol]" },
    { k: "rho",   lo: 800,     hi: 1200,    lbl: "ρ [kg/m³]" },
    { k: "y",     lo: 0.05,    hi: 0.65,    lbl: "y_F" },
    { k: "tads",  lo: 10,      hi: 110,     lbl: "t_ADS [s]" },
    { k: "PH",    lo: 1.0,     hi: 5.0,     lbl: "P_H [bar]" },
    { k: "lPI",   lo: -1.846,  hi: -0.0519, lbl: "P_I/P_H" },
    { k: "lPL",   lo: -2.695,  hi: -0.2029, lbl: "P_L/P_H" },
    { k: "vF",    lo: 0.1,     hi: 1.5,     lbl: "v_F [m/s]" },
  ];
  /* The LHS sampled ABSOLUTE pressures: P_I ∈ [0.07, 4.0], P_L ∈
   * [0.01, 1.0] bar (identical in both cycles). The ratio bounds in BOX
   * are only their envelope over P_H ∈ [1, 5]: e.g. log10(P_L/P_H) =
   * -2.695 exists ONLY at P_H = 5. A point can sit inside every ratio
   * bound yet in a joint data void (P_L = 2 mbar at P_H = 1 — zero
   * training rows). The absolute windows below close that hole; every
   * consumer (guards, sweeps, optimizer) must apply BOTH. */
  var PABS = { PImin: 0.07, PImax: 4.0, PLmin: 0.01, PLmax: 1.0 };
  function guards(u) {
    var c = pack(u), out = [];
    for (var i = 0; i < 12; i++)
      if (!(c[i] >= BOX[i].lo && c[i] <= BOX[i].hi)) out.push(BOX[i].lbl);
    var d = c[1] - c[2];                       // log10 selectivity of b(298K)
    if (!(d >= 0.50 && d <= 7.00)) out.push("b_CO2/b_N2 selectivity");
    if (!((-u.dUc) - (-u.dUn) >= 3.0)) out.push("ΔU ordering (|ΔU_CO2| ≥ |ΔU_N2|+3)");
    if (!(u.PI > u.PL)) out.push("P_I > P_L");
    if (!(u.PI >= PABS.PImin && u.PI <= PABS.PImax)) out.push("P_I ∈ [0.07, 4.0] bar (absolute)");
    if (!(u.PL >= PABS.PLmin && u.PL <= PABS.PLmax)) out.push("P_L ∈ [0.01, 1.0] bar (absolute)");
    return out;
  }

  /* ---------- SSL isotherm (for plots; pure-component) ---------- */
  function bT(b0, dU, T) { return b0 * Math.exp(-dU * 1000 / (Rg * T)); }
  function qStar(u, P, T, comp) {              // P [bar] -> q [mol/kg]
    var C = P * 1e5 / (Rg * T);
    var b = (comp === "c") ? bT(u.b0c, u.dUc, T) : bT(u.b0n, u.dUn, T);
    return u.qsat * b * C / (1 + b * C);
  }
  function qStarBinary(u, P, y, T) {           // competitive at total P [bar]
    var Cc = y * P * 1e5 / (Rg * T), Cn = (1 - y) * P * 1e5 / (Rg * T);
    var bc = bT(u.b0c, u.dUc, T), bn = bT(u.b0n, u.dUn, T);
    var den = 1 + bc * Cc + bn * Cn;
    return { qc: u.qsat * bc * Cc / den, qn: u.qsat * bn * Cn / den };
  }

  /* ---------- presets ----------
   * SSL parameters as reported by Pai et al.:
   *  - first four: Limits paper SI Table S2 (qsat mol/kg, rho per material)
   *    except Mg-MOF-74, from the MAPLE paper SI (mol/m3 at the fixed
   *    rho=1130 of that study; converted).
   *  - library: MAPLE paper SI screening table (36 SSL-representable
   *    materials of Khurana & Farooq + 3 named), same mol/m3 basis. */
  var PRESETS = {
    Z13X:   { label: "Zeolite 13X",  qsat: 4.390, b0c: 2.50e-6, b0n: 2.70e-6, dUc: -31.19, dUn: -16.38, rho: 1130 },
    UTSA16: { label: "UTSA-16",      qsat: 4.478, b0c: 4.70e-7, b0n: 1.40e-6, dUc: -30.57, dUn: -9.91,  rho: 1000 },
    IISERP: { label: "IISERP-MOF2",  qsat: 5.000, b0c: 2.02e-7, b0n: 2.64e-7, dUc: -31.13, dUn: -11.89, rho: 1000 },
    MgMOF74:{ label: "Mg-MOF-74",    qsat: 4.7180, b0c: 6.38e-7, b0n: 2.06e-6, dUc: -33.73, dUn: -18.32, rho: 1130 },
  };
  var LIBRARY = [
    ["NAB", 7477.28, 2.53e-7, 5.94e-6, 25.71, 8.55],
    ["h8291835", 5278.94, 5.98e-7, 1.46e-5, 24.03, 8.00],
    ["h8155527", 4300.50, 2.13e-7, 1.17e-5, 30.67, 8.00],
    ["CaX", 7745.77, 6.52e-7, 8.50e-5, 32.52, 8.00],
    ["MgX", 8026.19, 2.15e-7, 3.27e-7, 33.80, 20.00],
    ["NaA", 5112.23, 3.09e-6, 2.83e-5, 27.59, 8.00],
    ["NaX", 8125.71, 2.82e-7, 1.04e-5, 35.88, 12.65],
    ["PS-MFI", 6322.72, 2.08e-6, 5.12e-7, 23.71, 19.74],
    ["Zn-MOF-74", 12000.00, 4.07e-6, 1.00e-5, 23.02, 10.94],
    ["Co-MOF-74", 10223.22, 1.19e-7, 3.82e-5, 34.02, 9.88],
    ["Ni-MOF-74", 8510.18, 4.63e-6, 1.54e-6, 25.20, 20.00],
    ["MOF-177 (a)", 800.00, 4.66e-6, 1.27e-4, 20.00, 8.00],
    ["MOF-177 (b)", 880.55, 5.07e-6, 7.22e-5, 20.00, 8.74],
    ["CuBTC", 11342.63, 2.87e-7, 3.65e-6, 27.27, 12.60],
    ["mmen-CuBTTri", 5242.98, 1.18e-8, 4.28e-7, 40.39, 17.72],
    ["ZIF-68", 5126.00, 1.79e-6, 3.16e-6, 21.58, 12.86],
    ["ZIF-69", 5850.63, 3.76e-6, 8.48e-6, 20.01, 10.68],
    ["ZIF-70", 2452.20, 4.90e-6, 1.03e-5, 20.00, 11.45],
    ["ZIF-78", 3654.30, 3.04e-6, 1.78e-5, 24.11, 10.75],
    ["ZIF-79", 3339.58, 1.49e-6, 6.78e-6, 23.46, 12.06],
    ["ZIF-81", 4359.92, 1.88e-6, 7.90e-6, 23.00, 11.50],
    ["ZIF-82", 4033.62, 2.87e-6, 9.17e-6, 22.40, 11.15],
    ["PPN-4", 800.00, 4.38e-6, 1.24e-4, 20.00, 8.00],
    ["PPN-6-SO3H", 3852.00, 1.39e-6, 3.01e-5, 25.30, 8.00],
    ["ZIF-36-CAG", 2970.81, 3.21e-7, 3.84e-6, 45.00, 20.00],
    ["ZIF-39-DIA", 9811.80, 1.00e-7, 5.80e-7, 26.51, 13.82],
    ["ZIF-39-ZNI", 2111.13, 8.01e-8, 6.80e-7, 27.14, 19.76],
    ["ZIF-116-MER", 10718.71, 1.89e-6, 1.33e-5, 20.16, 9.02],
    ["HMOF-MOF-5", 875.51, 7.86e-6, 1.52e-4, 20.02, 8.00],
    ["HMOF-16", 11762.36, 1.13e-6, 9.68e-6, 21.35, 8.00],
    ["HMOF-27", 845.91, 5.57e-6, 1.73e-4, 20.00, 8.00],
    ["HMOF-96", 906.68, 6.17e-6, 1.67e-3, 20.00, 8.00],
    ["HMOF-602", 2605.99, 3.50e-6, 5.93e-5, 20.05, 8.00],
    ["HMOF-972", 848.34, 7.02e-6, 1.64e-4, 20.00, 8.00],
    ["HMOF-992", 6894.91, 1.62e-6, 1.52e-5, 20.00, 8.00],
    ["HMOF-1041", 855.48, 7.34e-6, 1.38e-4, 20.00, 8.00],
  ].map(function (r) {
    return { label: r[0], qsat: r[1] / 1130, b0c: r[2], b0n: r[3],
             dUc: -r[4], dUn: -r[5], rho: 1130, fromM3: true };
  });

  /* ---------- NSGA-II (real-coded) ----------
   * opts: { nvar, lo[], hi[], pop, gens, fobj(x)->[f1,f2] (minimize),
   *         seedRows (optional array of x), onGen (optional cb) }     */
  var seed = 123456789;
  function rnd() {                            // deterministic LCG
    seed = (1103515245 * seed + 12345) % 2147483648;
    return seed / 2147483648;
  }
  function setSeed(s) { seed = s >>> 0 || 1; }
  function dominates(a, b) {
    var be = false;
    for (var i = 0; i < a.length; i++) {
      if (a[i] > b[i] + 1e-14) return false;
      if (a[i] < b[i] - 1e-14) be = true;
    }
    return be;
  }
  function ndSort(F) {
    var n = F.length, S = [], cnt = new Array(n).fill(0), fronts = [[]], rank = new Array(n);
    for (var i = 0; i < n; i++) S.push([]);
    for (i = 0; i < n; i++) for (var j = 0; j < n; j++) {
      if (i === j) continue;
      if (dominates(F[i], F[j])) S[i].push(j);
      else if (dominates(F[j], F[i])) cnt[i]++;
    }
    for (i = 0; i < n; i++) if (cnt[i] === 0) { rank[i] = 0; fronts[0].push(i); }
    var f = 0;
    while (fronts[f].length) {
      var nxt = [];
      fronts[f].forEach(function (p) {
        S[p].forEach(function (q) { if (--cnt[q] === 0) { rank[q] = f + 1; nxt.push(q); } });
      });
      f++; fronts.push(nxt);
    }
    fronts.pop();
    return { fronts: fronts, rank: rank };
  }
  function crowding(F, idx) {
    var d = {}, m = F[0].length;
    idx.forEach(function (i) { d[i] = 0; });
    for (var k = 0; k < m; k++) {
      var s = idx.slice().sort(function (a, b) { return F[a][k] - F[b][k]; });
      var lo = F[s[0]][k], hi = F[s[s.length - 1]][k], span = hi - lo || 1;
      d[s[0]] = d[s[s.length - 1]] = Infinity;
      for (var i = 1; i < s.length - 1; i++)
        d[s[i]] += (F[s[i + 1]][k] - F[s[i - 1]][k]) / span;
    }
    return d;
  }
  /* Stepwise NSGA-II: init → step (one generation each) → front.
   * nsga2(opts) below wraps them and is call-for-call identical to the
   * original monolithic version (same rnd() sequence, same results). */
  function nsga2init(opts) {
    var nv = opts.nvar, lo = opts.lo, hi = opts.hi;
    var NP = opts.pop || 120, NG = opts.gens || 80;
    function newX() {
      var x = new Array(nv);
      for (var i = 0; i < nv; i++) x[i] = lo[i] + rnd() * (hi[i] - lo[i]);
      return x;
    }
    var X = [], Fv = [];
    (opts.seedRows || []).slice(0, NP).forEach(function (r) { X.push(r.slice()); });
    while (X.length < NP) X.push(newX());
    X.forEach(function (x) { Fv.push(opts.fobj(x)); });
    var nd0 = ndSort(Fv);
    return { opts: opts, nv: nv, lo: lo, hi: hi, NP: NP, NG: NG,
             etaC: 15, etaM: 20, pm: 1 / nv,
             g: 0, X: X, Fv: Fv, nd: nd0, cd: crowdAll(Fv, nd0) };
  }
  function crowdAll(F, nd) {          // crowding distance for every member
    var d = {};
    nd.fronts.forEach(function (fr) {
      var c = crowding(F, fr);
      fr.forEach(function (i) { d[i] = c[i]; });
    });
    return d;
  }
  function nsga2step(S) {              // one generation; false when finished
    if (S.g >= S.NG) return false;
    var nv = S.nv, lo = S.lo, hi = S.hi, NP = S.NP;
    var etaC = S.etaC, etaM = S.etaM, pm = S.pm;
    var X = S.X, Fv = S.Fv, opts = S.opts;
    function clamp(v, i) { return v < lo[i] ? lo[i] : v > hi[i] ? hi[i] : v; }
    // crowded-comparison operator (Deb et al.): rank first, then the less
    // crowded of two equals. Rank-only selection degenerates to a random
    // walk once the whole population reaches rank 0, which loses spread
    // along near-degenerate (flat) fronts.
    function tourn() {
      var a = (rnd() * X.length) | 0, b = (rnd() * X.length) | 0;
      if (S.nd.rank[a] !== S.nd.rank[b]) return S.nd.rank[a] < S.nd.rank[b] ? a : b;
      var da = S.cd[a], db = S.cd[b];
      if (da === undefined) da = -1;
      if (db === undefined) db = -1;
      return da >= db ? a : b;
    }
    function sbx(p1, p2) {
      var c1 = p1.slice(), c2 = p2.slice();
      if (rnd() < 0.9) {
        for (var i = 0; i < nv; i++) if (rnd() < 0.5) {
          var u = rnd(), bq = (u <= 0.5)
            ? Math.pow(2 * u, 1 / (etaC + 1))
            : Math.pow(1 / (2 * (1 - u)), 1 / (etaC + 1));
          var x1 = 0.5 * ((1 + bq) * p1[i] + (1 - bq) * p2[i]);
          var x2 = 0.5 * ((1 - bq) * p1[i] + (1 + bq) * p2[i]);
          c1[i] = clamp(x1, i); c2[i] = clamp(x2, i);
        }
      }
      return [c1, c2];
    }
    function mut(x) {
      for (var i = 0; i < nv; i++) if (rnd() < pm) {
        var u = rnd(), d = (u < 0.5)
          ? Math.pow(2 * u, 1 / (etaM + 1)) - 1
          : 1 - Math.pow(2 * (1 - u), 1 / (etaM + 1));
        x[i] = clamp(x[i] + d * (hi[i] - lo[i]), i);
      }
    }
    var kids = [];
    while (kids.length < NP) {
      // binary tournament by rank then crowding is approximated by
      // rank-only here (crowding used at survival); adequate in practice
      var a = tourn(), b = tourn();
      var c = sbx(X[a], X[b]);
      mut(c[0]); mut(c[1]);
      kids.push(c[0]); if (kids.length < NP) kids.push(c[1]);
    }
    var kidF = kids.map(opts.fobj);
    var allX = X.concat(kids), allF = Fv.concat(kidF);
    var nd = ndSort(allF), nX = [], nF = [];
    for (var fi = 0; fi < nd.fronts.length && nX.length < NP; fi++) {
      var fr = nd.fronts[fi];
      if (nX.length + fr.length <= NP) {
        fr.forEach(function (i) { nX.push(allX[i]); nF.push(allF[i]); });
      } else {
        var cd = crowding(allF, fr);
        fr.slice().sort(function (p, q) { return cd[q] - cd[p]; })
          .slice(0, NP - nX.length)
          .forEach(function (i) { nX.push(allX[i]); nF.push(allF[i]); });
      }
    }
    S.X = nX; S.Fv = nF; S.nd = ndSort(nF); S.cd = crowdAll(nF, S.nd);
    if (opts.onGen) opts.onGen(S.g, S.X, S.Fv);
    S.g++;
    return S.g < S.NG;
  }
  function nsga2front(S) {             // current first front
    return S.nd.fronts[0].map(function (i) { return { x: S.X[i], f: S.Fv[i] }; });
  }
  function nsga2(opts) {
    var S = nsga2init(opts);
    while (nsga2step(S)) { }
    return nsga2front(S);
  }

  /* ---------- canned optimization problems ----------
   * DVs (5): tads, PH, lPI = log10(PI/PH), lPL = log10(PL/PH), vF
   * (log-ratio DVs keep the search inside the trained manifold).
   * mode: "PuRe" | "EnRe" | "PrRe" | "EnPr"; fixPH pins P_H = 1 bar (VSA).
   * cons: {puMin, reMin} in % — penalty constraints (defaults 95/90).
   * PuRe ignores them (the whole trade-off is the answer); EnRe/PrRe use
   * puMin only (recovery is an objective); EnPr uses both.              */
  /* DV_FULL: the full user-space search box = the trained support. P_I and
   * P_L are ABSOLUTE bar (log-sampled), matching how the LHS was built;
   * the network's ratio inputs are derived per candidate and rejected if
   * they leave BOX. `dvr` may only NARROW these (never widen). */
  var DV_FULL = { tads: [10, 110], PH: [1, 5],
                  PI: [PABS.PImin, PABS.PImax], PL: [PABS.PLmin, PABS.PLmax],
                  vF: [0.1, 1.5] };
  var DV_KEYS = ["tads", "PH", "PI", "PL", "vF"];
  function dvRanges(dvr, fixPH) {         // clamp user ranges into DV_FULL
    dvr = dvr || {};
    var out = {};
    DV_KEYS.forEach(function (k) {
      var f = DV_FULL[k], r = dvr[k];
      var a = (r && isFinite(+r[0])) ? +r[0] : f[0];
      var b = (r && isFinite(+r[1])) ? +r[1] : f[1];
      a = Math.min(Math.max(a, f[0]), f[1]);
      b = Math.min(Math.max(b, f[0]), f[1]);
      if (b < a) b = a;                    // inverted input → degenerate, not empty
      out[k] = [a, b];
    });
    if (fixPH) out.PH = [1, 1];
    return out;
  }
  function optProblem(cyc, mat, y, mode, fixPH, cons, dvr) {
    cons = cons || {};
    var puMin = (cons.puMin != null ? cons.puMin : 95) / 100;
    var reMin = (cons.reMin != null ? cons.reMin : 90) / 100;
    var R = dvRanges(dvr, fixPH);
    var lo = [R.tads[0], R.PH[0], Math.log10(R.PI[0]), Math.log10(R.PL[0]), R.vF[0]];
    var hi = [R.tads[1], R.PH[1], Math.log10(R.PI[1]), Math.log10(R.PL[1]), R.vF[1]];
    function toU(x) {
      return { qsat: mat.qsat, b0c: mat.b0c, b0n: mat.b0n, dUc: mat.dUc,
               dUn: mat.dUn, rho: mat.rho, y: y, tads: x[0], PH: x[1],
               PI: Math.pow(10, x[2]), PL: Math.pow(10, x[3]), vF: x[4] };
    }
    function fobj(x) {
      var u = toU(x);
      if (!(u.PI > u.PL)) return [1e6, 1e6];
      // absolute sampled-pressure windows (see PABS)
      if (!(u.PI >= PABS.PImin && u.PI <= PABS.PImax &&
            u.PL >= PABS.PLmin && u.PL <= PABS.PLmax)) return [1e6, 1e6];
      // the networks see RATIOS: reject candidates whose ratios leave BOX
      // (absolute DVs can produce out-of-support ratios at some P_H)
      var qI = Math.log10(u.PI / u.PH), qL = Math.log10(u.PL / u.PH);
      if (!(qI >= BOX[9].lo && qI <= BOX[9].hi)) return [1e6, 1e6];
      if (!(qL >= BOX[10].lo && qL <= BOX[10].hi)) return [1e6, 1e6];
      var k = evalKPI(cyc, u);
      // repo wrapper's validity window: reject nonphysical net output
      // (the surrogate can stray a little past 100% in DV corners)
      if (!(k.Pu > 0 && k.Pu < 102 && k.Re > 0 && k.Re < 102 &&
            k.Prod > 0 && k.En > 0)) return [1e6, 1e6];
      var Pu = Math.min(k.Pu, 100), Re = Math.min(k.Re, 100);
      if (mode === "PuRe") return [-Pu, -Re];
      // penalty P = sum over active constraints of 20 v + 5000 v^2.
      // The linear term matters: a pure quadratic has ~zero gradient at
      // tiny violations, so the GA parks the front a hair OUTSIDE the
      // constraint (e.g. Pu = 94.9) and the results table rejects it all.
      // Keep this form in sync with the objective display in the UI.
      var pv = Math.max(0, puMin - Pu / 100);
      var rv = mode === "EnPr" ? Math.max(0, reMin - Re / 100) : 0;
      var pen = 20 * (pv + rv) + 5000 * (pv * pv + rv * rv);
      if (mode === "EnRe") return [k.En / 100 + pen, -Re / 100 + pen];
      if (mode === "PrRe") return [10 / Math.max(k.Prod, 1e-6) + pen, -Re / 100 + pen];
      /* EnPr */           return [k.En / 100 + pen, -k.Prod + pen];
    }
    return { nvar: 5, lo: lo, hi: hi, fobj: fobj, toU: toU, ranges: R };
  }

  return {
    setWeights: setWeights, setSeed: setSeed,
    pack: pack, evalKPI: evalKPI, guards: guards, BOX: BOX, PABS: PABS,
    qStar: qStar, qStarBinary: qStarBinary, bT: bT,
    PRESETS: PRESETS, LIBRARY: LIBRARY,
    nsga2: nsga2, nsga2init: nsga2init, nsga2step: nsga2step,
    nsga2front: nsga2front, optProblem: optProblem, fwd: fwd,
    DV_FULL: DV_FULL, DV_KEYS: DV_KEYS, dvRanges: dvRanges,
  };
})();
if (typeof module !== "undefined") module.exports = MAPLE;
