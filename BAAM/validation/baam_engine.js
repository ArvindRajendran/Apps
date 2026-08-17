/* BAAM engine — Batch Adsorber Analogue Model (Balashankar 2019 / e-BAAM
   Liske 2026). 4-step VSA/PVSA cycles with LPP or FP pressurization.
   Mixture equilibria: extended dual-site Langmuir (analytic partials,
   concentration basis) or binary IAST over pure-component isotherms
   (partial-pressure basis, FD partials). No DOM access.

   Sections 1-2 (parseExprP, IAST_MODELS, buildPsiTable, makeComponent,
   iastSolve) are copied VERBATIM from the IAST app engine (ia-engine block in
   Apps/IAST/IAST.html) — keep in sync with that app. */
"use strict";

/* ================= 1. expression parser (from IAST app) ================= */
function parseExprP(src) {
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
    if (/^[\d.]/.test(t)) out.push(parseFloat(t));
    else if (/^[A-Za-z_]/.test(t)) {
      if (FUN[t]) ops.push(t);
      else if (t in CONST) out.push(CONST[t]);
      else if (t === "p" || t === "P") out.push("VAR");
      else throw "Unknown symbol '" + t + "' — use p, numbers, exp/log/sqrt/…";
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
  const f = p => {
    const st = [];
    for (const tk of out) {
      if (typeof tk === "number") st.push(tk);
      else if (tk === "VAR") st.push(p);
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

/* ================= 2. IAST pieces (from IAST app) ======================= */
const IAST_MODELS = {
  linear: { label: "Linear", params: [["H", "H [mol/(kg·bar)]", 3]],
    make: c => ({ f: p => c.H * p, psi: p => c.H * p, henry: c.H, pole: null, psiKind: "exact" }) },
  langmuir: { label: "Langmuir", params: [["qs", "q<sub>s</sub> [mol/kg]", 4], ["b", "b [1/bar]", 8]],
    make: c => ({ f: p => c.qs * c.b * p / (1 + c.b * p),
                  psi: p => c.qs * Math.log(1 + c.b * p),
                  henry: c.qs * c.b, pole: null, psiKind: "exact" }) },
  dsl: { label: "Dual-site Langmuir",
    params: [["qs1", "q<sub>s1</sub> [mol/kg]", 3], ["b1", "b<sub>1</sub> [1/bar]", 50],
             ["qs2", "q<sub>s2</sub> [mol/kg]", 4], ["b2", "b<sub>2</sub> [1/bar]", 1]],
    make: c => ({ f: p => c.qs1 * c.b1 * p / (1 + c.b1 * p) + c.qs2 * c.b2 * p / (1 + c.b2 * p),
                  psi: p => c.qs1 * Math.log(1 + c.b1 * p) + c.qs2 * Math.log(1 + c.b2 * p),
                  henry: c.qs1 * c.b1 + c.qs2 * c.b2, pole: null, psiKind: "exact" }) },
  antilang: { label: "Anti-Langmuir", params: [["qs", "q<sub>s</sub> [mol/kg]", 1], ["b", "b [1/bar]", 0.05]],
    make: c => ({ f: p => c.qs * c.b * p / (1 - c.b * p),
                  psi: p => -c.qs * Math.log(1 - c.b * p),
                  henry: c.qs * c.b, pole: 1 / c.b, psiKind: "exact" }) },
  quad: { label: "Quadratic", params: [["qs", "q<sub>s</sub> [mol/kg]", 3], ["b1", "b<sub>1</sub> [1/bar]", 5], ["b2", "b<sub>2</sub> [1/bar²]", 10]],
    make: c => ({ f: p => c.qs * p * (c.b1 + 2 * c.b2 * p) / (1 + c.b1 * p + c.b2 * p * p),
                  psi: p => c.qs * Math.log(1 + c.b1 * p + c.b2 * p * p),
                  henry: c.qs * c.b1, pole: null, psiKind: "exact" }) },
  sips: { label: "Sips / Hill", params: [["qs", "q<sub>s</sub> [mol/kg]", 5], ["b", "b [1/bar]", 2], ["n", "n [-]", 0.85]],
    make: c => ({ f: p => { const t = Math.pow(Math.max(c.b * p, 0), c.n); return c.qs * t / (1 + t); },
                  psi: p => c.qs / c.n * Math.log(1 + Math.pow(Math.max(c.b * p, 0), c.n)),
                  henry: c.n === 1 ? c.qs * c.b : (c.n < 1 ? Infinity : 0),
                  pole: null, psiKind: "exact" }) },
  toth: { label: "Toth", params: [["qs", "q<sub>s</sub> [mol/kg]", 4], ["b", "b [1/bar]", 8], ["t", "t [-]", 0.6]],
    make: c => ({ f: p => c.qs * c.b * p / Math.pow(1 + Math.pow(Math.max(c.b * p, 0), c.t), 1 / c.t),
                  psi: null, henry: c.qs * c.b, pole: null, psiKind: "table" }) },
  freundlich: { label: "Freundlich", params: [["k", "k [mol/(kg·bar<sup>1/n</sup>)]", 2], ["n", "n [-]", 2]],
    make: c => ({ f: p => c.k * Math.pow(Math.max(p, 0), 1 / c.n),
                  psi: p => c.n * c.k * Math.pow(Math.max(p, 0), 1 / c.n),
                  henry: c.n === 1 ? c.k : (c.n > 1 ? Infinity : 0),
                  pole: null, psiKind: "exact" }) },
  bet: { label: "BET (Type II/III)", params: [["qs", "q<sub>s</sub> [mol/kg]", 2], ["bs", "b<sub>s</sub> [1/bar]", 30], ["bl", "b<sub>l</sub> [1/bar]", 0.8]],
    make: c => ({ f: p => c.qs * c.bs * p / ((1 - c.bl * p) * (1 - c.bl * p + c.bs * p)),
                  psi: p => c.qs * Math.log((1 + (c.bs - c.bl) * p) / (1 - c.bl * p)),
                  henry: c.qs * c.bs, pole: 1 / c.bl, psiKind: "exact" }) },
  typev: { label: "Type V (Langmuir + Hill step)",
    params: [["qs1", "q<sub>s1</sub> [mol/kg]", 1], ["b1", "b<sub>1</sub> [1/bar]", 5],
             ["qs2", "q<sub>s2</sub> [mol/kg]", 8], ["b2", "b<sub>2</sub> [1/bar]", 2.5], ["n", "n [-]", 6]],
    make: c => ({ f: p => { const t = Math.pow(Math.max(c.b2 * p, 0), c.n);
                    return c.qs1 * c.b1 * p / (1 + c.b1 * p) + c.qs2 * t / (1 + t); },
                  psi: p => c.qs1 * Math.log(1 + c.b1 * p) +
                            c.qs2 / c.n * Math.log(1 + Math.pow(Math.max(c.b2 * p, 0), c.n)),
                  henry: c.n < 1 ? Infinity : c.qs1 * c.b1,
                  pole: null, psiKind: "exact" }) },
  custom: { label: "Custom q(p) …", params: [], make: null }
};

function buildPsiTable(f, pMax, N) {
  N = N || 20000;
  const uLo = Math.log(1e-10), uHi = Math.log(pMax);
  const h = (uHi - uLo) / N;
  const F = new Float64Array(N + 1), Psi = new Float64Array(N + 1);
  for (let i = 0; i <= N; i++) F[i] = f(Math.exp(uLo + i * h));
  const p0 = Math.exp(uLo);
  const q0 = f(p0), qh = f(p0 / 2);
  const alpha = (q0 > 0 && qh > 0) ? Math.log(q0 / qh) / Math.LN2 : 1;
  Psi[0] = alpha > 0.01 ? q0 / alpha : 0;
  for (let i = 0; i + 2 <= N; i += 2) {
    Psi[i + 1] = Psi[i] + h / 12 * (5 * F[i] + 8 * F[i + 1] - F[i + 2]);
    Psi[i + 2] = Psi[i + 1] + h / 12 * (-F[i] + 8 * F[i + 1] + 5 * F[i + 2]);
  }
  return { uLo, h, N, Psi,
    eval(p) {
      if (p <= 0) return 0;
      const u = Math.log(p);
      if (u <= this.uLo) return this.Psi[0] * Math.exp(alpha * (u - this.uLo));
      const s = (u - this.uLo) / this.h;
      const i = Math.min(this.N - 1, Math.floor(s));
      const w = s - i;
      return this.Psi[i] + w * (this.Psi[i + 1] - this.Psi[i]);
    } };
}

function makeComponent(cfg, pMaxNeeded) {
  let inst;
  if (cfg.id === "custom") {
    const f = parseExprP(cfg.expr);
    inst = { f, psi: null, henry: NaN, pole: null, psiKind: "table" };
    const h1 = f(1e-8) / 1e-8, h2 = f(1e-9) / 1e-9;
    inst.henry = Math.abs(h1 / h2 - 1) < 0.02 ? h2 : (h2 > h1 ? Infinity : 0);
  } else {
    const M = IAST_MODELS[cfg.id];
    const par = {};
    M.params.forEach(([k, , d]) => {
      const v = cfg.pvals ? cfg.pvals[k] : undefined;
      par[k] = (v === undefined || !isFinite(v)) ? d : +v;
    });
    inst = M.make(par);
    inst.par = par;
  }
  const cap = inst.pole ? inst.pole * (1 - 1e-9) : (pMaxNeeded || 1e10);
  inst.pCap = cap;
  if (inst.psiKind === "table" || !inst.psi) {
    const tab = buildPsiTable(inst.f, cap);
    inst.psi = p => tab.eval(Math.min(p, cap));
    inst.psiTab = tab;
  } else {
    const psiRaw = inst.psi;
    inst.psi = p => psiRaw(Math.min(p, cap));
  }
  inst.mono = true;
  let prevQ = 0;
  for (let k = 0; k <= 240; k++) {
    const p = Math.pow(10, -8 + 10.5 * k / 240);
    if (p >= cap) break;
    const q = inst.f(p);
    if (!isFinite(q) || q < prevQ - 1e-12 * Math.abs(prevQ)) { inst.mono = false; break; }
    prevQ = q;
  }
  return inst;
}

function iastSolve(P, y1, c1, c2) {
  if (y1 <= 0) { const q = c2.f(P); return { x1: 0, q1: 0, q2: q, qT: q, p1o: 0, p2o: P, psi: c2.psi(P), edge: true }; }
  if (y1 >= 1) { const q = c1.f(P); return { x1: 1, q1: q, q2: 0, qT: q, p1o: P, p2o: 0, psi: c1.psi(P), edge: true }; }
  const y2 = 1 - y1;
  const G = x1 => c1.psi(P * y1 / x1) - c2.psi(P * y2 / (1 - x1));
  let lo = 1e-15, hi = 1 - 1e-15;
  for (let it = 0; it < 50; it++) {
    const mid = 0.5 * (lo + hi);
    if (G(mid) > 0) lo = mid; else hi = mid;
  }
  let x1 = 0.5 * (lo + hi);
  for (let it = 0; it < 6; it++) {
    const p1o = P * y1 / x1, p2o = P * y2 / (1 - x1);
    const g = c1.psi(p1o) - c2.psi(p2o);
    const dg = -c1.f(Math.min(p1o, c1.pCap)) / x1 - c2.f(Math.min(p2o, c2.pCap)) / (1 - x1);
    const step = g / dg;
    const x1n = x1 - step;
    if (!(x1n > 0 && x1n < 1)) break;
    x1 = x1n;
    if (Math.abs(step) < 1e-16) break;
  }
  const p1o = P * y1 / x1, p2o = P * y2 / (1 - x1);
  const q1o = c1.f(Math.min(p1o, c1.pCap)), q2o = c2.f(Math.min(p2o, c2.pCap));
  const qT = 1 / (x1 / q1o + (1 - x1) / q2o);
  return { x1, q1: x1 * qT, q2: (1 - x1) * qT, qT,
           p1o, p2o, psi: c1.psi(p1o),
           capped: (p1o >= c1.pCap * 0.999) || (p2o >= c2.pCap * 0.999) };
}

/* ================= 3. BAAM core ========================================= */
const BA_RB = 8.314e-5;   // m3 bar / (mol K)
const BA_RE = 8.314;      // J / (mol K)
const BA_K  = 1.4;        // adiabatic constant
const BA_J2KWH = 1 / 3.6e6;

/* mixture: extended DSL, concentration basis (paper Eq. 20-21).
   co2/n2: {qsb,qsd,b0,d0,dUb,dUd}; dU in J/mol, NEGATIVE for exothermic. */
function edslMix(co2, n2, T) {
  const act = g => [(g.b0 > 0) ? g.b0 * Math.exp(-g.dUb / (BA_RE * T)) : 0,
                    (g.d0 > 0) ? g.d0 * Math.exp(-g.dUd / (BA_RE * T)) : 0];
  const [bA, dA] = act(co2), [bB, dB] = act(n2);
  const sites = [[bA, bB, co2.qsb, n2.qsb], [dA, dB, co2.qsd, n2.qsd]];
  return { T, kind: "edsl", sites,
    loadings(P, y) {
      const RT = BA_RB * T;
      const c1 = P * y / RT, c2 = P * (1 - y) / RT;
      const c1P = y / RT, c2P = (1 - y) / RT, c1y = P / RT, c2y = -P / RT;
      let q1 = 0, q2 = 0, q1c1 = 0, q1c2 = 0, q2c1 = 0, q2c2 = 0;
      for (const [kA, kB, qsA, qsB] of sites) {
        const D = 1 + kA * c1 + kB * c2, D2 = D * D;
        q1 += qsA * kA * c1 / D;
        q2 += qsB * kB * c2 / D;
        q1c1 += qsA * kA * (1 + kB * c2) / D2;
        q1c2 += -qsA * kA * c1 * kB / D2;
        q2c2 += qsB * kB * (1 + kA * c1) / D2;
        q2c1 += -qsB * kB * c2 * kA / D2;
      }
      return [q1, q2,
              q1c1 * c1P + q1c2 * c2P, q1c1 * c1y + q1c2 * c2y,
              q2c1 * c1P + q2c2 * c2P, q2c1 * c1y + q2c2 * c2y];
    } };
}

/* mixture: binary IAST over two makeComponent instances (p-basis, bar);
   partials by central finite differences (same deltas as the Python ref). */
function iastMix(comp1, comp2, T) {
  function qq(P, y) {
    if (P <= 0) return [0, 0];
    const r = iastSolve(P, Math.min(Math.max(y, 0), 1), comp1, comp2);
    return [r.q1, r.q2];
  }
  return { T, kind: "iast", comp1, comp2,
    q(P, y) { return qq(P, y); },     // fast path: loadings only, no partials
    loadings(P, y) {
      const dP = Math.max(1e-6, 1e-6 * P);
      const dy = 1e-7;
      const yl = Math.max(0, y - dy), yh = Math.min(1, y + dy);
      const c0 = qq(P, y);
      const cp = qq(P + dP, y), cm = qq(Math.max(P - dP, 1e-12), y);
      const cyh = qq(P, yh), cyl = qq(P, yl);
      return [c0[0], c0[1],
              (cp[0] - cm[0]) / (2 * dP), (cyh[0] - cyl[0]) / (yh - yl),
              (cp[1] - cm[1]) / (2 * dP), (cyh[1] - cyl[1]) / (yh - yl)];
    } };
}

function etaOf(mode, Pbar) {
  if (mode && mode.type === "maruyama") {
    const Patm = Pbar / 1.01325;
    return 0.8 * (19.55 * Patm) / (1 + 19.55 * Patm);
  }
  return (mode && mode.val) || 0.72;
}

/* BLO/EVAC path: RK4 fixed step on the uniform dP grid from PH down to Plow.
   Returns arrays + cumulative removed moles and vacuum work (MATLAB
   convention: slice dN charged at the slice-END pressure; only P < 1 bar). */
function blowPath(mix, opt) {
  const T = opt.T, w = opt.w ?? 1, eps = opt.eps ?? 0.37;
  const V = w / (opt.rho * (1 - eps));
  const RT = BA_RB * T;
  const rhs = (P, yv) => {
    const y = Math.min(Math.max(yv, 0), 1);
    const L = mix.loadings(P, y);
    const a1 = V * eps / RT + w * (L[2] + L[4]);
    const a2 = y * V * eps / RT + w * L[2];
    const f1 = w * (L[3] + L[5]);
    const f2 = P * V * eps / RT + w * L[3];
    return (a1 * y - a2) / (f2 - f1 * y);
  };
  const dP = opt.dP ?? 1e-4;
  const n = Math.round((opt.PH - opt.Plow) / dP);
  const P = new Float64Array(n + 1), y = new Float64Array(n + 1);
  P[0] = opt.PH; y[0] = opt.yF;
  const h = -dP;
  for (let i = 0; i < n; i++) {
    const Pi = opt.PH - i * dP;
    const k1 = rhs(Pi, y[i]);
    const k2 = rhs(Pi + h / 2, y[i] + h / 2 * k1);
    const k3 = rhs(Pi + h / 2, y[i] + h / 2 * k2);
    const k4 = rhs(Pi + h, y[i] + h * k3);
    P[i + 1] = opt.PH - (i + 1) * dP;
    y[i + 1] = Math.min(Math.max(y[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), 0), 1);
  }
  const q1 = new Float64Array(n + 1), q2 = new Float64Array(n + 1);
  const N1s = new Float64Array(n + 1), N2s = new Float64Array(n + 1);
  for (let i = 0; i <= n; i++) {
    const L = mix.loadings(P[i], y[i]);
    q1[i] = L[0]; q2[i] = L[1];
    N1s[i] = P[i] * y[i] * V * eps / RT + w * L[0];
    N2s[i] = P[i] * (1 - y[i]) * V * eps / RT + w * L[1];
  }
  const N = new Float64Array(n + 1), N1 = new Float64Array(n + 1), W = new Float64Array(n + 1);
  const kk = (BA_K - 1) / BA_K, fac = BA_K / (BA_K - 1) * BA_RE * T;
  for (let i = 1; i <= n; i++) {
    const dN = (N1s[i - 1] + N2s[i - 1]) - (N1s[i] + N2s[i]);
    N[i] = N[i - 1] + dN;
    N1[i] = N1[i - 1] + (N1s[i - 1] - N1s[i]);
    W[i] = W[i - 1] + (P[i] < 1
      ? fac / etaOf(opt.eta, P[i]) * dN * (Math.pow(1 / P[i], kk) - 1) : 0);
  }
  return { P, y, q1, q2, N, N1, W, N1s, N2s, V, w, eps, T, dP, PH: opt.PH, eta: opt.eta };
}

function baIdx(path, Pq) {
  const i = Math.round((path.PH - Pq) / path.dP);
  return Math.min(Math.max(i, 0), path.P.length - 1);
}

function baInv(mix, T, V, eps, w, P, y) {
  const L = mix.q ? mix.q(P, y) : mix.loadings(P, y);
  return [P * y * V * eps / (BA_RB * T) + w * L[0],
          P * (1 - y) * V * eps / (BA_RB * T) + w * L[1]];
}

/* bisection helper on [lo,hi] for scalar g with bracket scan fallback */
function baRoot(g, lo, hi) {
  let glo = g(lo), ghi = g(hi);
  if (glo * ghi > 0) {
    let found = false;
    let pl = lo, gl = glo;
    for (let k = 1; k <= 400; k++) {
      const t = lo + (hi - lo) * k / 400;
      const gt = g(t);
      if (gl * gt <= 0) { lo = pl; hi = t; found = true; break; }
      pl = t; gl = gt;
    }
    if (!found) return null;
  }
  glo = g(lo);
  for (let it = 0; it < 200; it++) {
    const mid = 0.5 * (lo + hi);
    if (glo * g(mid) <= 0) hi = mid; else { lo = mid; glo = g(lo); }
    if (hi - lo < 1e-16) break;
  }
  return 0.5 * (lo + hi);
}

/* LPP: pressurize (Plow,ygam)->(PH,ydel) with gas of composition ydel */
function solveLPP(mix, T, V, eps, w, PH, Plow, ygam) {
  const [N1i, N2i] = baInv(mix, T, V, eps, w, Plow, ygam);
  const Ni = N1i + N2i;
  const g = yd => {
    const [N1f, N2f] = baInv(mix, T, V, eps, w, PH, yd);
    return N1i - N1f + ((N1f + N2f) - Ni) * yd;
  };
  const yd = baRoot(g, 1e-12, ygam);
  if (yd === null) return null;
  const [N1f, N2f] = baInv(mix, T, V, eps, w, PH, yd);
  return { ydel: yd, Npr: (N1f + N2f) - Ni };
}

/* FP: pressurize with FEED gas */
function solveFP(mix, T, V, eps, w, PH, Plow, ygam, yF) {
  const [N1i, N2i] = baInv(mix, T, V, eps, w, Plow, ygam);
  const Ni = N1i + N2i;
  const g = ye => {
    const [N1f, N2f] = baInv(mix, T, V, eps, w, PH, ye);
    return N1i - N1f + ((N1f + N2f) - Ni) * yF;
  };
  const ye = baRoot(g, 1e-12, Math.max(ygam, yF));
  if (ye === null) return null;
  const [N1f, N2f] = baInv(mix, T, V, eps, w, PH, ye);
  return { ydel: ye, Npr: (N1f + N2f) - Ni };
}

/* ADS: linear solve; raffinate leaves at y_delta (BAAM.m convention) */
function solveADS(mix, T, V, eps, w, PH, ydel, yF) {
  const [N1i, N2i] = baInv(mix, T, V, eps, w, PH, ydel);
  const [N1f, N2f] = baInv(mix, T, V, eps, w, PH, yF);
  const A = (N1f + N2f) - (N1i + N2i);
  const B = N1f - N1i;
  const Nfeed = (B - A * ydel) / (yF - ydel);
  return { Nfeed, Nraff: Nfeed - A };
}

/* full-cycle KPIs at (Pint, Plow) reusing a precomputed path.
   mode: "LPP" | "FP". Compression work for PH>1 bar per e-BAAM. */
function cycleKPIs(mix, path, opt) {
  const { Pint, Plow, yF } = opt;
  const mode = opt.mode || "LPP";
  const T = path.T, V = path.V, eps = path.eps, w = path.w, PH = path.PH;
  const ib = baIdx(path, Pint), ig = baIdx(path, Plow);
  const n1e = path.N1[ig] - path.N1[ib];
  const ne = path.N[ig] - path.N[ib];
  const Wblo = path.W[ib], Wevac = path.W[ig] - path.W[ib];
  const ygam = path.y[ig];
  // pressurization + ADS depend on Plow only — memoizable across a Pint sweep
  const ck = opt.cache ? `${ig}|${mode}|${yF}` : null;
  let pa = ck && opt.cache[ck];
  if (!pa) {
    const pr = (mode === "LPP")
      ? solveLPP(mix, T, V, eps, w, PH, Plow, ygam)
      : solveFP(mix, T, V, eps, w, PH, Plow, ygam, yF);
    if (!pr) return null;
    pa = { pr, ads: solveADS(mix, T, V, eps, w, PH, pr.ydel, yF) };
    if (ck) opt.cache[ck] = pa;
  }
  const { ydel, Npr } = pa.pr;
  const ads = pa.ads;
  const kk = (BA_K - 1) / BA_K;
  const comp = N => (PH > 1)
    ? N * BA_K / (BA_K - 1) * BA_RE * T / etaOf(path.eta, 1.0) * (Math.pow(PH, kk) - 1) : 0;
  const Wads = comp(ads.Nfeed);
  const Wpr = (mode === "FP") ? comp(Npr) : 0;
  const mT = n1e * 44e-6;                      // tonnes CO2 (44 g/mol, as BAAM.m)
  const fed = ads.Nfeed * yF + (mode === "FP" ? Npr * yF : 0);
  const nin = ads.Nfeed + (mode === "FP" ? Npr : 0);
  const nout = (ads.Nraff - (mode === "LPP" ? Npr : 0)) + path.N[ib] + ne;
  return {
    Pu: 100 * n1e / ne,
    Re: 100 * n1e / fed,
    EnBLO: Wblo * BA_J2KWH / mT,
    EnEVAC: Wevac * BA_J2KWH / mT,
    EnPR: Wpr * BA_J2KWH / mT,
    EnADS: Wads * BA_J2KWH / mT,
    En: (Wblo + Wevac + Wads + Wpr) * BA_J2KWH / mT,
    WC: n1e / (V * (1 - eps)),
    ydel, ygam, Npr, Nfeed: ads.Nfeed, Nraff: ads.Nraff,
    n1evac: n1e, nevac: ne,
    mbal: (nin - nout) / Math.max(Math.abs(nout), 1e-300) * 100
  };
}

/* design-space scan: KPI matrices over (Plow columns, Pint rows).
   One path + one pressurization/ADS solve per Plow; buffer diffs per Pint. */
function gridScan(mix, path, opt) {
  const { plows, pints, yF } = opt;
  const nr = pints.length, nc = plows.length;
  const mk = () => { const a = new Float64Array(nr * nc); a.fill(NaN); return a; };
  const Pu = mk(), Re = mk(), En = mk(), WC = mk();
  const cache = {};
  for (let j = 0; j < nc; j++) {
    const Plow = plows[j];
    for (let i = 0; i < nr; i++) {
      const Pint = pints[i];
      if (Pint <= Plow + 1e-12 || Pint >= path.PH - 1e-12) continue;
      const k = cycleKPIs(mix, path, { Pint, Plow, yF, mode: opt.mode, cache });
      if (!k) continue;
      const m = i * nc + j;
      Pu[m] = k.Pu; Re[m] = k.Re; En[m] = k.En; WC[m] = k.WC;
    }
  }
  return { plows, pints, Pu, Re, En, WC };
}

/* r_max over a scanned grid + e-BAAM DOE classification (r_cut = 124.5) */
function rMax(grid) {
  let rm = 0, im = -1;
  for (let m = 0; m < grid.Pu.length; m++) {
    if (!isFinite(grid.Pu[m])) continue;
    const r = Math.hypot(grid.Pu[m], grid.Re[m]);
    if (r > rm) { rm = r; im = m; }
  }
  return { rmax: rm, idx: im };
}

const BA_RCUT_EBAAM = 124.5;      // Liske 2026 linear-SVM PRT cutoff
const BA_EN_SCALE_EBAAM = 1.58;   // En_full ~= 1.58 x En_BAAM (Liske Fig. 9a)
const BA_RCUT_2019 = 110.25;      // Balashankar 2019 r(95,90)
const BA_EN_SCALE_2019 = [1.1446, 66.528]; // En_full = 1.1446 En + 66.528

if (typeof globalThis !== "undefined") {
  globalThis.BAAMENG = { parseExprP, IAST_MODELS, buildPsiTable, makeComponent,
    iastSolve, edslMix, iastMix, blowPath, cycleKPIs, gridScan, rMax, baIdx,
    solveLPP, solveFP, solveADS, etaOf,
    BA_RB, BA_RE, BA_K, BA_RCUT_EBAAM, BA_EN_SCALE_EBAAM, BA_RCUT_2019,
    BA_EN_SCALE_2019 };
}
