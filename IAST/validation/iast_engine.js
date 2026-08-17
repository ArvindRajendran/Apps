/* IAST engine — binary Ideal Adsorbed Solution Theory.
   Pure isotherms q_i(p) [mol/kg], p = partial pressure [bar].
   psi(p) = int_0^p q/p' dp' (reduced grand potential), analytic for most
   models, cumulative Newton-Cotes table in ln p for Toth/custom.
   Solver: bisection + Newton on x1 in (0,1):
     G(x1) = psi1(P y1/x1) - psi2(P y2/(1-x1)), strictly decreasing,
     G(0+)=+inf, G(1-)=-inf -> unique root, forward psi evals only.
   Shared by main thread; no DOM access. */
"use strict";

/* ---------- expression parser (variable p) ---------- */
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

/* ---------- isotherm menu ----------
   make(par) -> { f, psi, henry, pole, psiKind } ; pole = pressure where the
   model diverges (null if none); henry: finite>0, 0, or Infinity. */
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

/* ---------- numerical psi: cumulative Newton-Cotes table in ln p ----------
   Per interval pair (u0,u1,u2), f0,f1,f2 known:
     int_u0^u1 = h/12 (5 f0 + 8 f1 - f2),  int_u1^u2 = h/12 (-f0 + 8 f1 + 5 f2)
   4th-order accurate; tail below u_min via local power-law exponent. */
function buildPsiTable(f, pMax, N) {
  N = N || 20000;                           // intervals (even)
  const uLo = Math.log(1e-10), uHi = Math.log(pMax);
  const h = (uHi - uLo) / N;
  const F = new Float64Array(N + 1), Psi = new Float64Array(N + 1);
  for (let i = 0; i <= N; i++) F[i] = f(Math.exp(uLo + i * h));
  const p0 = Math.exp(uLo);
  const q0 = f(p0), qh = f(p0 / 2);
  const alpha = (q0 > 0 && qh > 0) ? Math.log(q0 / qh) / Math.LN2 : 1;
  Psi[0] = alpha > 0.01 ? q0 / alpha : 0;   // int_0^p0 q/p dp for q ~ C p^alpha
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

/* ---------- component assembly ---------- */
function makeComponent(cfg, pMaxNeeded) {
  // cfg: {id, pvals, expr}; pMaxNeeded: largest hypothetical pressure the
  // solver may request (cap; poles shrink it)
  let inst;
  if (cfg.id === "custom") {
    const f = parseExprP(cfg.expr);
    inst = { f, psi: null, henry: NaN, pole: null, psiKind: "table" };
    // numeric Henry estimate + monotonicity scan
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
  // monotonicity scan (for badges; over 12 decades)
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

/* ---------- binary IAST solve ----------
   Bisection (50) + Newton polish on x1; dG/dx1 = -q1/x1 - q2/(1-x1). */
function iastSolve(P, y1, c1, c2) {
  if (y1 <= 0) { const q = c2.f(P); return { x1: 0, q1: 0, q2: q, qT: q, p1o: 0, p2o: P, psi: c2.psi(P), edge: true }; }
  if (y1 >= 1) { const q = c1.f(P); return { x1: 1, q1: q, q2: 0, qT: q, p1o: P, p2o: 0, psi: c1.psi(P), edge: true }; }
  const y2 = 1 - y1;
  const G = x1 => c1.psi(P * y1 / x1) - c2.psi(P * y2 / (1 - x1));
  let lo = 1e-15, hi = 1 - 1e-15;
  // capped-p sanity: G(lo) uses p1o = huge -> clipped at c1.pCap; still fine
  for (let it = 0; it < 50; it++) {
    const mid = 0.5 * (lo + hi);
    if (G(mid) > 0) lo = mid; else hi = mid;
  }
  let x1 = 0.5 * (lo + hi);
  for (let it = 0; it < 6; it++) {          // Newton polish
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

function iastSelectivity(r, y1) {
  if (y1 <= 0 || y1 >= 1) return NaN;
  return (r.x1 / y1) / ((1 - r.x1) / (1 - y1));
}

/* ---------- sweep over (log P, y1) ---------- */
function iastSweep(c1, c2, Pmin, Pmax, NP, NY) {
  const q1 = new Float64Array(NP * NY), q2 = new Float64Array(NP * NY);
  const qT = new Float64Array(NP * NY), x1A = new Float64Array(NP * NY);
  const S = new Float64Array(NP * NY);
  const ext = new Float64Array(NP * NY);   // max extrapolation factor p_io/P
  const Ps = new Float64Array(NP), Ys = new Float64Array(NY);
  const lp0 = Math.log10(Pmin), lp1 = Math.log10(Pmax);
  for (let i = 0; i < NP; i++) Ps[i] = Math.pow(10, lp0 + (lp1 - lp0) * i / (NP - 1));
  for (let j = 0; j < NY; j++) Ys[j] = j / (NY - 1);
  let extMax = 0, cappedAny = false;
  for (let i = 0; i < NP; i++) {
    for (let j = 0; j < NY; j++) {
      const r = iastSolve(Ps[i], Ys[j], c1, c2);
      const k = i * NY + j;
      q1[k] = r.q1; q2[k] = r.q2; qT[k] = r.qT; x1A[k] = r.x1;
      S[k] = iastSelectivity(r, Ys[j]);
      ext[k] = Math.max(r.p1o, r.p2o) / Ps[i];
      if (isFinite(ext[k])) extMax = Math.max(extMax, ext[k]);
      if (r.capped) cappedAny = true;
    }
  }
  return { Ps, Ys, NP, NY, q1, q2, qT, x1: x1A, S, ext, extMax, cappedAny };
}

/* x-y diagram at fixed P + azeotrope detection */
function xyCurve(c1, c2, P, N) {
  N = N || 121;
  const y = new Float64Array(N), x = new Float64Array(N);
  for (let j = 0; j < N; j++) {
    y[j] = j / (N - 1);
    x[j] = iastSolve(P, y[j], c1, c2).x1;
  }
  const az = [];
  for (let j = 1; j < N - 1; j++) {
    const d0 = x[j] - y[j], d1 = x[j + 1] - y[j + 1];
    if (d0 === 0 && j > 0) az.push(y[j]);
    else if (d0 * d1 < 0) az.push(y[j] + (y[j + 1] - y[j]) * d0 / (d0 - d1));
  }
  return { y, x, az };
}

/* Selectivity-reversal pressures: p* where psi1(p) = psi2(p).
   Binary IAST cannot produce composition azeotropes — S = p2o/p1o depends on
   psi alone, so S=1 forces p1o=p2o=P: the x-y curve touches the diagonal only
   at P = p*, where it degenerates onto it. S>1 on one side of p*, <1 on the
   other; the S=1 contour in (P,y1) is the horizontal line P=p*. */
function findSReversal(c1, c2, pLo, pHi) {
  pLo = pLo || 1e-8; pHi = Math.min(pHi || 1e8, c1.pCap, c2.pCap);
  const g = p => c1.psi(p) - c2.psi(p);
  const scale = p => Math.max(c1.psi(p), c2.psi(p), 1e-300);
  const N = 400;
  // degeneracy first: psi1 == psi2 (identical isotherms, e.g. quad-quad with
  // the same parameters) makes every sample a "crossing" — detect and report
  // it as one condition, not 400 bogus p* values
  let maxRel = 0;
  const us = [], gs = [];
  for (let k = 0; k <= N; k++) {
    const u = Math.log(pLo) + (Math.log(pHi) - Math.log(pLo)) * k / N;
    const p = Math.exp(u);
    us.push(u); gs.push(g(p));
    maxRel = Math.max(maxRel, Math.abs(gs[k]) / scale(p));
  }
  if (maxRel < 1e-7) return { ps: [], degenerate: true };
  const ps = [];
  for (let k = 1; k <= N && ps.length < 10; k++) {
    if (!(gs[k - 1] * gs[k] < 0)) continue;
    let a = us[k - 1], b = us[k];
    const sPrev = gs[k - 1];
    for (let it = 0; it < 60; it++) {
      const m = 0.5 * (a + b);
      if (g(Math.exp(m)) * sPrev > 0) a = m; else b = m;
    }
    const pStar = Math.exp(0.5 * (a + b));
    // transversality: flanks 5% away must have opposite signs and real size
    const tol = 1e-9 + 1e-7 * scale(pStar);
    const gl = g(pStar * 0.95), gr = g(pStar * 1.05);
    if (gl * gr < 0 && Math.abs(gl) > tol && Math.abs(gr) > tol) ps.push(pStar);
  }
  return { ps, degenerate: false };
}

/* extended (competitive) Langmuir for the overlay: both comps Langmuir */
function extendedLangmuir(P, y1, par1, par2) {
  const p1 = P * y1, p2 = P * (1 - y1);
  const d = 1 + par1.b * p1 + par2.b * p2;
  return { q1: par1.qs * par1.b * p1 / d, q2: par2.qs * par2.b * p2 / d };
}

if (typeof globalThis !== "undefined") {
  globalThis.IAST = { IAST_MODELS, parseExprP, makeComponent, buildPsiTable,
                      iastSolve, iastSelectivity, iastSweep, xyCurve,
                      findSReversal, extendedLangmuir };
}
