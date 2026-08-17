/* UptakeCurves engine — shared by main thread and worker.
   Models: micro | macro | bidisperse | ldf. Spheres, step y0->y1 at constant P.
   FV cell-centered on normalized radius; CN time stepping on sqrt-spaced grid
   with quasi-Newton (colored FD Jacobian, reuse + rebuild policy, sub-stepping).
   Mirrors the Python reference uptake_ref.py exactly (same discretization). */
"use strict";
const R_GAS = 8.314462618;

/* ---------- isotherm models: q [mol/kg] = f(y), y = mole fraction ---------- */
const MODELS = {
  langmuir: { label: "Langmuir", params: [["qs", "q<sub>s</sub> [mol/kg]", 5], ["b", "b [-]", 10]],
    make: p => y => p.qs*p.b*y/(1 + p.b*y) },
  antilang: { label: "Anti-Langmuir", params: [["qs", "q<sub>s</sub> [mol/kg]", 1], ["b", "b [-]", 0.9]],
    make: p => y => p.qs*p.b*y/(1 - p.b*y) },
  linear: { label: "Linear", params: [["H", "H [mol/kg]", 3]], make: p => y => p.H*y },
  sips: { label: "Sips / Hill", params: [["qs", "q<sub>s</sub> [mol/kg]", 5], ["b", "b [-]", 10], ["n", "n [-]", 0.7]],
    make: p => y => { const t = Math.pow(Math.max(p.b*y, 0), p.n); return p.qs*t/(1 + t); } },
  quad: { label: "Quadratic", params: [["qs", "q<sub>s</sub> [mol/kg]", 3], ["b1", "b<sub>1</sub> [-]", 5], ["b2", "b<sub>2</sub> [-]", 10]],
    make: p => y => p.qs*y*(p.b1 + 2*p.b2*y)/(1 + p.b1*y + p.b2*y*y) },
  typev: { label: "Type V (Langmuir + Hill step)",
    params: [["qs1", "q<sub>s1</sub> [mol/kg]", 1], ["b1", "b<sub>1</sub> [-]", 5],
             ["qs2", "q<sub>s2</sub> [mol/kg]", 8], ["b2", "b<sub>2</sub> [-]", 2.5], ["n", "n [-]", 6]],
    make: p => y => { const t = Math.pow(Math.max(p.b2*y, 0), p.n);
      return p.qs1*p.b1*y/(1 + p.b1*y) + p.qs2*t/(1 + t); } },
  bet: { label: "BET (Type II/III)",
    params: [["qs", "q<sub>s</sub> [mol/kg]", 2], ["bs", "b<sub>s</sub> [-]", 30], ["bl", "b<sub>l</sub> [-]", 0.8]],
    make: p => y => p.qs*p.bs*y/((1 - p.bl*y)*(1 - p.bl*y + p.bs*y)) },
  dsl: { label: "Dual-site Langmuir",
    params: [["qs1", "q<sub>s1</sub> [mol/kg]", 3], ["b1", "b<sub>1</sub> [-]", 50],
             ["qs2", "q<sub>s2</sub> [mol/kg]", 4], ["b2", "b<sub>2</sub> [-]", 1]],
    make: p => y => p.qs1*p.b1*y/(1 + p.b1*y) + p.qs2*p.b2*y/(1 + p.b2*y) },
  dodoA: { label: "Do & Do water/AC — ads. branch (Hefti 2015)",
    params: [["ns", "n<sub>s</sub><sup>&infin;</sup> [mol/kg]", 10.791], ["Ks", "K<sub>s</sub>", 0.258],
             ["Kc", "K<sub>c</sub>", 0.222], ["nma", "n<sub>&mu;,a</sub><sup>&infin;</sup> [mol/kg]", 2.227],
             ["Kma", "K<sub>&mu;,a</sub>", 28.579], ["ma", "m<sub>a</sub>", 8.8],
             ["ys", "y<sub>sat</sub>", 0.0968]],
    make: p => y => { const x = y/p.ys, t = Math.pow(Math.max(x, 0), p.ma);
      return p.ns*p.Ks*x/((1 - p.Kc*x)*(1 + (p.Ks - p.Kc)*x)) +
             p.nma*p.ma*p.Kma*t/(1 + p.Kma*t); } },
  dodoD: { label: "Do & Do water/AC — des. branch with hysteresis (Hefti 2015)",
    params: [["ns", "n<sub>s</sub><sup>&infin;</sup> [mol/kg]", 10.791], ["Ks", "K<sub>s</sub>", 0.258],
             ["Kc", "K<sub>c</sub>", 0.222], ["nma", "n<sub>&mu;,a</sub><sup>&infin;</sup> [mol/kg]", 2.227],
             ["Kma", "K<sub>&mu;,a</sub>", 28.579], ["ma", "m<sub>a</sub>", 8.8],
             ["nmd", "n<sub>&mu;,d</sub><sup>&infin;</sup> [mol/kg]", 1.02], ["Kmd", "K<sub>&mu;,d</sub>", 24715],
             ["md", "m<sub>d</sub>", 18.243], ["ys", "y<sub>sat</sub>", 0.0968]],
    make: p => y => { const x = y/p.ys, xc = Math.max(x, 0);
      const ta = Math.pow(xc, p.ma), td = Math.pow(xc, p.md);
      return p.ns*p.Ks*x/((1 - p.Kc*x)*(1 + (p.Ks - p.Kc)*x)) +
             Math.max(p.nma*p.ma*p.Kma*ta/(1 + p.Kma*ta),
                      p.nmd*p.md*p.Kmd*td/(1 + p.Kmd*td)); } },
  custom: { label: "Custom q = f(y) …", params: [], make: null }
};

/* ---------- expression parser for custom isotherms (no eval) ---------- */
function parseExpr(src) {
  const toks = [];
  const re = /\s*(\d+\.?\d*(?:[eE][+-]?\d+)?|[A-Za-z_][A-Za-z_0-9]*|\*\*|[-+*/^(),])/y;
  let i = 0;
  while (i < src.length) {
    re.lastIndex = i;
    const m = re.exec(src);
    if (!m) throw "Bad token near '" + src.slice(i, i + 8) + "'";
    toks.push(m[1]); i = re.lastIndex;
  }
  const FUN = { exp: Math.exp, log: Math.log, ln: Math.log, sqrt: Math.sqrt,
                abs: Math.abs, tanh: Math.tanh, pow: Math.pow, min: Math.min, max: Math.max };
  const CONST = { pi: Math.PI, e: Math.E };
  const out = [], ops = [];
  const prec = { "+": 1, "-": 1, "*": 2, "/": 2, "u-": 3, "^": 4, "**": 4 };
  const rightAssoc = { "^": 1, "**": 1, "u-": 1 };
  let prev = null;
  for (const t0 of toks) {
    let t = t0;
    if (t === "-" && (prev === null || (prec[prev] || prev === "(" || prev === ","))) t = "u-";
    if (/^[\d.]/.test(t)) { out.push(parseFloat(t)); }
    else if (/^[A-Za-z_]/.test(t)) {
      if (FUN[t]) ops.push(t);
      else if (t in CONST) out.push(CONST[t]);
      else if (t === "y" || t === "c" || t === "x") out.push("VAR");
      else throw "Unknown symbol '" + t + "'";
    }
    else if (t === "(") ops.push(t);
    else if (t === ",") { while (ops.length && ops[ops.length-1] !== "(") out.push(ops.pop()); }
    else if (t === ")") {
      while (ops.length && ops[ops.length-1] !== "(") out.push(ops.pop());
      ops.pop();
      if (ops.length && FUN[ops[ops.length-1]]) out.push(ops.pop());
    } else {
      while (ops.length) {
        const o = ops[ops.length-1];
        if (o === "(" || !(prec[o] >= (rightAssoc[t] ? prec[t] + 1 : prec[t]))) break;
        out.push(ops.pop());
      }
      ops.push(t);
    }
    prev = t;
  }
  while (ops.length) out.push(ops.pop());
  const f = y => {
    const st = [];
    for (const tk of out) {
      if (typeof tk === "number") st.push(tk);
      else if (tk === "VAR") st.push(y);
      else if (FUN[tk]) {
        if (tk === "pow" || tk === "min" || tk === "max") { const b = st.pop(), a = st.pop(); st.push(FUN[tk](a, b)); }
        else st.push(FUN[tk](st.pop()));
      }
      else if (tk === "u-") st.push(-st.pop());
      else { const b = st.pop(), a = st.pop();
        st.push(tk === "+" ? a + b : tk === "-" ? a - b : tk === "*" ? a * b :
                tk === "/" ? a / b : Math.pow(a, b)); }
    }
    if (st.length !== 1) throw "Malformed expression";
    return st[0];
  };
  f(0.123);
  return f;
}

function makeIsotherm(iso) {
  if (iso.id === "custom") return parseExpr(iso.expr);
  const M = MODELS[iso.id];
  const p = {};
  M.params.forEach(([k, , d]) => { p[k] = (iso.pvals && iso.pvals[k] !== undefined) ? +iso.pvals[k] : d; });
  return M.make(p);
}

/* ---------- Darken factor table: Gamma(q) = q/(y f'(y)) ---------- */
function gammaTable(f, ylo, yhi, npts) {
  npts = npts || 3000;
  const y0 = Math.max(yhi * 1e-6, ylo <= 0 ? yhi * 1e-6 : ylo * 0.5);
  const qA = new Float64Array(npts), gA = new Float64Array(npts);
  let mono = true;
  for (let k = 0; k < npts; k++) {
    const y = y0 + (yhi - y0) * k / (npts - 1);
    const h = Math.max(1e-9, y * 1e-6);
    const fp = (f(y + h) - f(y - h)) / (2 * h);
    qA[k] = f(y);
    gA[k] = qA[k] / (y * fp);
    if (k > 0 && qA[k] <= qA[k-1]) mono = false;
    if (!isFinite(qA[k]) || !isFinite(gA[k])) mono = false;
  }
  const G = q => {
    if (q <= qA[0]) return gA[0];
    if (q >= qA[npts-1]) return gA[npts-1];
    let lo = 0, hi = npts - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (qA[m] <= q) lo = m; else hi = m; }
    const w = (q - qA[lo]) / (qA[hi] - qA[lo]);
    return gA[lo] + w * (gA[hi] - gA[lo]);
  };
  return { G, mono, qA, gA };
}

/* ---------- FV geometry ---------- */
function fvGeom(N) {
  const h = 1 / N, rf2 = new Float64Array(N + 1), Vi = new Float64Array(N);
  for (let i = 0; i <= N; i++) rf2[i] = (i * h) * (i * h);
  for (let i = 0; i < N; i++) {
    const a = i * h, b = (i + 1) * h;
    Vi[i] = (b * b * b - a * a * a) / 3;
  }
  return { h, rf2, Vi };
}

/* crystal RHS: writes dq into out[off..off+Nc), returns surface flux Fs */
function crystalRhs(x, off, Nc, qs, amu, geo, G, out) {
  const { h, rf2, Vi } = geo;
  let Fprev = 0;
  for (let i = 0; i < Nc; i++) {
    let Fi;
    if (i < Nc - 1) {
      const g = G ? G(0.5 * (x[off + i] + x[off + i + 1])) : 1;
      Fi = g * (x[off + i + 1] - x[off + i]) / h;
    } else {
      const g = G ? G(0.5 * (qs + x[off + i])) : 1;
      Fi = g * (qs - x[off + i]) / (h / 2);
    }
    out[off + i] = amu * (rf2[i + 1] * Fi - rf2[i] * Fprev) / Vi[i];
    Fprev = Fi;
  }
  return Fprev; // flux at rho=1
}

/* ---------- model assembly ---------- */
function buildModel(cfg) {
  const fRaw = makeIsotherm(cfg.iso);
  const CT = cfg.P / (R_GAS * cfg.T);
  const ylo = Math.min(cfg.y0, cfg.y1), yhi = Math.max(cfg.y0, cfg.y1);
  // max principle: y stays in [ylo,yhi]; clamp so Newton wander can't hit
  // isotherm poles (e.g. Langmuir at y=-1/b, BET at y=1/bl)
  const f = y => fRaw(y < ylo ? ylo : y > yhi ? yhi : y);
  let G = null, gmono = true;
  if (cfg.darken && cfg.model !== "macro" && cfg.model !== "ldf") {
    const gt = gammaTable(f, ylo, yhi);
    G = gt.G; gmono = gt.mono;
    if (!gmono) G = null; // fall back to constant D; UI warns
  }
  const amu = cfg.Dmu0 / (cfg.Rc * cfg.Rc);
  const De = cfg.epsp * (1 / (1 / cfg.Dm + 1 / cfg.DK)) / cfg.taup;
  const aDe = De / (cfg.Rp * cfg.Rp);
  const q0 = f(cfg.y0), q1 = f(cfg.y1);
  const hh = Math.max(1e-9, yhi * 1e-6);
  // derivative of the RAW isotherm at a point nudged just inside [ylo,yhi],
  // so the clamp above never truncates the central-difference stencil
  const fp = y => {
    const yc = Math.min(Math.max(y, ylo + 2 * hh), Math.max(ylo + 2 * hh, yhi - 2 * hh));
    return (fRaw(yc + hh) - fRaw(yc - hh)) / (2 * hh);
  };

  const m = { f, CT, G, gmono, amu, De, aDe, q0, q1, cfg };

  if (cfg.model === "micro") {
    const Nc = cfg.Nc, geo = fvGeom(Nc);
    m.n = Nc;
    m.x0 = new Float64Array(Nc).fill(q0);
    m.rhs = (x, out) => { crystalRhs(x, 0, Nc, q1, amu, geo, G, out); };
    m.sparsity = triSparsity(Nc);
    m.uptake = x => {
      let qb = 0;
      for (let i = 0; i < Nc; i++) qb += 3 * geo.Vi[i] * x[i];
      return (qb - q0) / (q1 - q0);
    };
    m.profiles = x => ({ crystal: [Array.from(x)], geoN: Nc });
  }
  else if (cfg.model === "macro") {
    const Np = cfg.Np, geo = fvGeom(Np);
    m.n = Np;
    m.x0 = new Float64Array(Np).fill(cfg.y0);
    m.rhs = (x, out) => {
      let Fprev = 0;
      for (let j = 0; j < Np; j++) {
        const Fj = (j < Np - 1) ? (x[j + 1] - x[j]) / geo.h
                                : (cfg.y1 - x[j]) / (geo.h / 2);
        const lap = (geo.rf2[j + 1] * Fj - geo.rf2[j] * Fprev) / geo.Vi[j];
        const beta = cfg.epsp + (cfg.rhop / CT) * fp(x[j]);
        out[j] = aDe * lap / beta;
        Fprev = Fj;
      }
    };
    m.sparsity = triSparsity(Np);
    const m0 = cfg.epsp * CT * cfg.y0 + cfg.rhop * q0;
    const minf = cfg.epsp * CT * cfg.y1 + cfg.rhop * q1;
    m.uptake = x => {
      let mm = 0;
      for (let j = 0; j < Np; j++) mm += 3 * geo.Vi[j] * (cfg.epsp * CT * x[j] + cfg.rhop * f(x[j]));
      return (mm - m0) / (minf - m0);
    };
    m.profiles = x => ({ particleY: Array.from(x), qbar: Array.from(x, v => f(v)), geoN: Np });
  }
  else if (cfg.model === "bidisperse") {
    const Np = cfg.Np, Nc = cfg.Nc;
    const geoC = fvGeom(Nc), geoP = fvGeom(Np);
    const blk = Nc + 1;
    m.n = Np * blk;
    m.x0 = new Float64Array(m.n);
    for (let j = 0; j < Np; j++) {
      for (let i = 0; i < Nc; i++) m.x0[j * blk + i] = q0;
      m.x0[j * blk + Nc] = cfg.y0;
    }
    m.rhs = (x, out) => {
      let Fprev = 0;
      for (let j = 0; j < Np; j++) {
        const yj = x[j * blk + Nc];
        const qs = f(yj);
        const Fs = crystalRhs(x, j * blk, Nc, qs, amu, geoC, G, out);
        const dqbar = 3 * amu * Fs;
        const yR = (j < Np - 1) ? x[(j + 1) * blk + Nc] : 0;
        const Fj = (j < Np - 1) ? (yR - yj) / geoP.h : (cfg.y1 - yj) / (geoP.h / 2);
        const lap = (geoP.rf2[j + 1] * Fj - geoP.rf2[j] * Fprev) / geoP.Vi[j];
        out[j * blk + Nc] = (aDe * lap - (cfg.rhop / CT) * dqbar) / cfg.epsp;
        Fprev = Fj;
      }
    };
    m.sparsity = bidiSparsity(Np, Nc);
    const m0 = cfg.epsp * CT * cfg.y0 + cfg.rhop * q0;
    const minf = cfg.epsp * CT * cfg.y1 + cfg.rhop * q1;
    m.uptake = x => {
      let mm = 0;
      for (let j = 0; j < Np; j++) {
        let qb = 0;
        for (let i = 0; i < Nc; i++) qb += 3 * geoC.Vi[i] * x[j * blk + i];
        mm += 3 * geoP.Vi[j] * (cfg.epsp * CT * x[j * blk + Nc] + cfg.rhop * qb);
      }
      return (mm - m0) / (minf - m0);
    };
    m.profiles = x => {
      const pY = [], qbar = [];
      for (let j = 0; j < Np; j++) {
        pY.push(x[j * blk + Nc]);
        let qb = 0;
        for (let i = 0; i < Nc; i++) qb += 3 * geoC.Vi[i] * x[j * blk + i];
        qbar.push(qb);
      }
      // crystal profiles at particle center, mid, surface shells
      const shells = [0, Math.floor(Np / 2), Np - 1];
      const crystal = shells.map(j => Array.from(x.subarray(j * blk, j * blk + Nc)));
      return { particleY: pY, qbar, crystal, shells, geoN: Np, geoNc: Nc };
    };
  }
  return m;
}

function triSparsity(n) {
  const rows = [];
  for (let i = 0; i < n; i++) {
    const c = [i];
    if (i > 0) c.push(i - 1);
    if (i < n - 1) c.push(i + 1);
    rows.push(c);
  }
  return rows;
}
function bidiSparsity(Np, Nc) {
  const blk = Nc + 1, rows = [];
  for (let j = 0; j < Np; j++) {
    for (let i = 0; i < Nc; i++) {
      const c = [j * blk + i];
      if (i > 0) c.push(j * blk + i - 1);
      if (i < Nc - 1) c.push(j * blk + i + 1);
      else c.push(j * blk + Nc);           // surface cell sees y_j
      rows.push(c);
    }
    const c = [j * blk + Nc, j * blk + Nc - 1];
    if (j > 0) c.push((j - 1) * blk + Nc);
    if (j < Np - 1) c.push((j + 1) * blk + Nc);
    rows.push(c);
  }
  return rows;
}

/* ---------- greedy column coloring for FD Jacobian ---------- */
function colorColumns(rows, n) {
  // columns conflict if they appear in the same row
  const colRows = Array.from({ length: n }, () => []);
  rows.forEach((cols, r) => cols.forEach(c => colRows[c].push(r)));
  const color = new Int32Array(n).fill(-1);
  let ncol = 0;
  const rowMark = new Int32Array(rows.length).fill(-1);
  for (let c = 0; c < n; c++) {
    const used = new Set();
    for (const r of colRows[c]) for (const c2 of rows[r]) if (color[c2] >= 0) used.add(color[c2]);
    let k = 0;
    while (used.has(k)) k++;
    color[c] = k;
    ncol = Math.max(ncol, k + 1);
  }
  return { color, ncol, colRows };
}

/* ---------- dense LU ---------- */
function luFactor(A, n, piv) {
  for (let k = 0; k < n; k++) {
    let p = k, mx = Math.abs(A[k * n + k]);
    for (let i = k + 1; i < n; i++) { const v = Math.abs(A[i * n + k]); if (v > mx) { mx = v; p = i; } }
    if (mx === 0) return false;
    piv[k] = p;
    if (p !== k) for (let j = 0; j < n; j++) { const t = A[k * n + j]; A[k * n + j] = A[p * n + j]; A[p * n + j] = t; }
    const d = A[k * n + k];
    for (let i = k + 1; i < n; i++) {
      const l = A[i * n + k] / d;
      A[i * n + k] = l;
      if (l !== 0) for (let j = k + 1; j < n; j++) A[i * n + j] -= l * A[k * n + j];
    }
  }
  return true;
}
function luSolve(A, n, piv, b) {
  for (let k = 0; k < n; k++) { const p = piv[k]; if (p !== k) { const t = b[k]; b[k] = b[p]; b[p] = t; } }
  for (let i = 1; i < n; i++) { let s = b[i]; for (let j = 0; j < i; j++) s -= A[i * n + j] * b[j]; b[i] = s; }
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i];
    for (let j = i + 1; j < n; j++) s -= A[i * n + j] * b[j];
    b[i] = s / A[i * n + i];
  }
}

/* ---------- CN integrator with quasi-Newton + colored FD Jacobian ---------- */
function integrateCN(model, times, opts) {
  opts = opts || {};
  const n = model.n, tol = opts.tol || 1e-10;
  const rhs = model.rhs;
  const { color, ncol, colRows } = colorColumns(model.sparsity, n);
  const x = Float64Array.from(model.x0);
  const f0 = new Float64Array(n), fx = new Float64Array(n), r = new Float64Array(n);
  const xn = new Float64Array(n), xp = new Float64Array(n), fb = new Float64Array(n);
  const Jraw = new Float64Array(n * n);   // dF/dx (sparse in dense storage)
  const M = new Float64Array(n * n);      // I - dt/2 J, LU-factored
  const piv = new Int32Array(n);
  let haveJ = false, dtJ = 0;
  const out = { states: [], stateAt: [] };
  const keep = opts.keepStates || null;    // array of step indices to snapshot
  const uptake = new Float64Array(times.length);
  uptake[0] = model.uptake(x);
  if (keep && keep.includes(0)) { out.states.push(Float64Array.from(x)); out.stateAt.push(0); }

  const buildJ = (xa) => {
    Jraw.fill(0);
    rhs(xa, fb);
    for (let cgrp = 0; cgrp < ncol; cgrp++) {
      xp.set(xa);
      const eps = [];
      for (let c = 0; c < n; c++) if (color[c] === cgrp) {
        const e = Math.max(1e-8, 1e-8 * Math.abs(xa[c]));
        xp[c] += e; eps[c] = e;
      }
      rhs(xp, fx);
      for (let c = 0; c < n; c++) if (color[c] === cgrp) {
        for (const rr of colRows[c]) Jraw[rr * n + c] = (fx[rr] - fb[rr]) / eps[c];
      }
    }
    haveJ = true;
  };
  const factorM = (wdt) => {          // M = I - wdt*J  (wdt = stage weight * dt)
    for (let i = 0; i < n * n; i++) M[i] = -wdt * Jraw[i];
    for (let i = 0; i < n; i++) M[i * n + i] += 1;
    dtJ = wdt;
    return luFactor(M, n, piv);
  };

  /* TR-BDF2 (L-stable), gamma = 2-sqrt(2): both stages share M = I - 0.29289*dt*J */
  const GAM = 2 - Math.SQRT2;
  const A1 = 1 / (GAM * (2 - GAM)), A2 = (1 - GAM) * (1 - GAM) / (GAM * (2 - GAM));
  const W = GAM / 2;                             // = (1-GAM)/(2-GAM)
  const xg = new Float64Array(n), rhsFix = new Float64Array(n);

  // solve xs - w*dt*f(xs) = rhsFix, starting from guess xs (in xn); true on success
  const solveStage = (w, dt, guess) => {
    xn.set(guess);
    let rn0 = Infinity, restarted = false;
    for (let it = 0; it < 40; it++) {
      rhs(xn, fx);
      let rn = 0;
      for (let i = 0; i < n; i++) { r[i] = xn[i] - w * dt * fx[i] - rhsFix[i]; rn = Math.max(rn, Math.abs(r[i])); }
      if (!isFinite(rn)) {
        if (restarted) return false;
        restarted = true;
        xn.set(guess); buildJ(x);
        if (!factorM(w * dt)) return false;
        rn0 = Infinity;
        continue;
      }
      if (rn < tol * Math.max(1, maxAbs(xn))) return true;
      if (it > 2 && rn > 0.5 * rn0) {
        buildJ(xn);
        if (!factorM(w * dt)) return false;
      }
      rn0 = rn;
      luSolve(M, n, piv, r);
      for (let i = 0; i < n; i++) xn[i] -= r[i];
    }
    return false;
  };

  const stepCN = (dt) => {
    rhs(x, f0);
    if (!haveJ) buildJ(x);
    if (dtJ === 0 || Math.abs (W * dt / dtJ - 1) > 0.3) { if (!factorM(W * dt)) return false; }
    // stage 1: TR over gamma*dt
    for (let i = 0; i < n; i++) rhsFix[i] = x[i] + W * dt * f0[i];
    if (!solveStage(W, dt, x)) return false;
    xg.set(xn);
    // stage 2: BDF2 to t+dt
    for (let i = 0; i < n; i++) rhsFix[i] = A1 * xg[i] - A2 * x[i];
    if (!solveStage(W, dt, xg)) return false;
    return true;
  };

  for (let k = 1; k < times.length; k++) {
    const dt = times[k] - times[k - 1];
    let nsub = 1, done = false;
    while (!done && nsub <= 64) {
      const xsave = Float64Array.from(x);
      let okAll = true;
      for (let s = 0; s < nsub; s++) {
        if (stepCN(dt / nsub)) { x.set(xn); }
        else { okAll = false; break; }
      }
      if (okAll) done = true;
      else { x.set(xsave); haveJ = false; nsub *= 2; }
    }
    if (!done) throw new Error("CN step failed at t=" + times[k]);
    uptake[k] = model.uptake(x);
    if (keep && keep.includes(k)) { out.states.push(Float64Array.from(x)); out.stateAt.push(k); }
    if (opts.onProgress && (k & 15) === 0) opts.onProgress(k / (times.length - 1));
  }
  out.uptake = uptake;
  return out;
}
function maxAbs(v) { let m = 0; for (let i = 0; i < v.length; i++) m = Math.max(m, Math.abs(v[i])); return m; }

/* ---------- top-level solve ---------- */
function timesSqrt(T, K) {
  const t = new Float64Array(K + 1);
  for (let k = 0; k <= K; k++) t[k] = (k / K) * (k / K) * T;
  return t;
}

function autoTend(cfg) {
  const f = makeIsotherm(cfg.iso);
  const CT = cfg.P / (R_GAS * cfg.T);
  const q0 = f(cfg.y0), q1 = f(cfg.y1);
  const dq = Math.abs(q1 - q0), dy = Math.abs(cfg.y1 - cfg.y0);
  const Kchord = dy > 0 ? (cfg.rhop * dq / (dy * CT)) : 0;
  const De = cfg.epsp * (1 / (1 / cfg.Dm + 1 / cfg.DK)) / cfg.taup;
  const tauMac = cfg.Rp * cfg.Rp * (cfg.epsp + Kchord) / De;
  let gbar = 1;
  if (cfg.darken) {
    const gt = gammaTable(f, Math.min(cfg.y0, cfg.y1), Math.max(cfg.y0, cfg.y1), 400);
    if (gt.mono) { let s = 0; for (let i = 0; i < gt.gA.length; i++) s += gt.gA[i]; gbar = s / gt.gA.length; }
  }
  const tauMic = cfg.Rc * cfg.Rc / (cfg.Dmu0 * Math.max(gbar, 1e-6));
  if (cfg.model === "ldf") return 5 / cfg.kldf;
  if (cfg.model === "micro") return 0.5 * tauMic;
  if (cfg.model === "macro") return 0.5 * tauMac;
  return 0.5 * tauMic + 0.5 * tauMac;
}

function solveUptake(cfg, onProgress) {
  const T = cfg.tEnd || autoTend(cfg);
  const K = cfg.K || 400;
  if (cfg.model === "ldf") {
    const f = makeIsotherm(cfg.iso);
    const t = timesSqrt(T, K);
    const U = Float64Array.from(t, tv => 1 - Math.exp(-cfg.kldf * tv));
    return { t: Array.from(t), U: Array.from(U), snaps: null, T,
             meta: { q0: f(cfg.y0), q1: f(cfg.y1) } };
  }
  const model = buildModel(cfg);
  const t = timesSqrt(T, K);
  const nSnap = Math.min(60, K);
  const keep = [];
  for (let s = 0; s <= nSnap; s++) keep.push(Math.round(s * K / nSnap));
  const res = integrateCN(model, t, { keepStates: keep, onProgress });
  const snaps = res.stateAt.map((ki, idx) => ({
    t: t[ki], ...model.profiles(res.states[idx])
  }));
  return {
    t: Array.from(t), U: Array.from(res.uptake), snaps, T,
    meta: { q0: model.q0, q1: model.q1, De: model.De, amu: model.amu,
            CT: model.CT, gmono: model.gmono }
  };
}

/* exports for jsc / worker */
if (typeof globalThis !== "undefined") {
  globalThis.UPTAKE = { MODELS, makeIsotherm, gammaTable, buildModel,
                        integrateCN, timesSqrt, solveUptake, autoTend, parseExpr, fvGeom };
}
