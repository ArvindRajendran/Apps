/* Isotherm Fitter engine — 16 models + custom, Levenberg-Marquardt with
   logit-bounded parameters and multi-start, weighting choice, AICc ranking,
   confidence intervals, per-T and global van 't Hoff multi-temperature modes,
   model-based isosteric heat. DOM-free; p in bar, q in mol/kg, T in K. */
"use strict";
const R_GAS = 8.314462618;

/* ---------- model library ----------
   f(th, p, T): loading; init(stats): heuristic start; bounds(stats): [lo, hi]
   arrays; vh: indices of affinity params that become b0*exp(U/RT) in global
   mode ("intrinsic" = model already uses T explicitly). */
const FIT_MODELS = {
  linear: { label: "Linear (Henry)", pnames: ["H"],
    f: (th, p) => th[0] * p,
    init: s => [Math.max(s.henry, 1e-6)],
    bounds: s => [[1e-9], [1e6]], vh: [0] },
  langmuir: { label: "Langmuir", pnames: ["qs", "b"],
    f: (th, p) => th[0] * th[1] * p / (1 + th[1] * p),
    init: s => [1.05 * s.qmax, 1 / s.pmid],
    bounds: s => [[1e-9, 1e-9], [1e6, 1e6]], vh: [1] },
  dsl: { label: "Dual-site Langmuir", pnames: ["qs1", "b1", "qs2", "b2"],
    f: (th, p) => th[0] * th[1] * p / (1 + th[1] * p) + th[2] * th[3] * p / (1 + th[3] * p),
    init: s => [0.6 * s.qmax, 5 / s.pmid, 0.6 * s.qmax, 0.2 / s.pmid],
    bounds: s => [[1e-9, 1e-9, 1e-9, 1e-9], [1e6, 1e6, 1e6, 1e6]], vh: [1, 3] },
  tsl: { label: "Triple-site Langmuir", pnames: ["qs1", "b1", "qs2", "b2", "qs3", "b3"],
    f: (th, p) => th[0] * th[1] * p / (1 + th[1] * p) + th[2] * th[3] * p / (1 + th[3] * p) +
                  th[4] * th[5] * p / (1 + th[5] * p),
    init: s => [0.4 * s.qmax, 25 / s.pmid, 0.4 * s.qmax, 1 / s.pmid, 0.4 * s.qmax, 0.04 / s.pmid],
    bounds: s => [Array(6).fill(1e-9), Array(6).fill(1e6)], vh: [1, 3, 5] },
  antilang: { label: "Anti-Langmuir", pnames: ["qs", "b"],
    f: (th, p) => th[0] * th[1] * p / (1 - th[1] * p),
    init: s => [0.5 * s.qmax, 0.5 / s.pmax],
    bounds: s => [[1e-9, 1e-9], [1e6, 0.95 / s.pmax]], vh: [1] },
  quad: { label: "Quadratic", pnames: ["qs", "b1", "b2"],
    f: (th, p) => th[0] * p * (th[1] + 2 * th[2] * p) / (1 + th[1] * p + th[2] * p * p),
    init: s => [0.55 * s.qmax, 1 / s.pmid, 0.1 / (s.pmid * s.pmid)],
    bounds: s => [[1e-9, 1e-9, 1e-12], [1e6, 1e6, 1e6]], vh: [1, 2] },
  sips: { label: "Sips / Langmuir-Freundlich", pnames: ["qs", "b", "n"],
    f: (th, p) => { const t = Math.pow(Math.max(th[1] * p, 0), th[2]); return th[0] * t / (1 + t); },
    init: s => [1.05 * s.qmax, 1 / s.pmid, 1],
    bounds: s => [[1e-9, 1e-9, 0.1], [1e6, 1e6, 3]], vh: [1] },
  dss: { label: "Dual-site Sips", pnames: ["qs1", "b1", "n1", "qs2", "b2", "n2"],
    f: (th, p) => {
      const t1 = Math.pow(Math.max(th[1] * p, 0), th[2]);
      const t2 = Math.pow(Math.max(th[4] * p, 0), th[5]);
      return th[0] * t1 / (1 + t1) + th[3] * t2 / (1 + t2);
    },
    init: s => [0.6 * s.qmax, 5 / s.pmid, 1, 0.6 * s.qmax, 0.2 / s.pmid, 1],
    bounds: s => [[1e-9, 1e-9, 0.1, 1e-9, 1e-9, 0.1], [1e6, 1e6, 3, 1e6, 1e6, 3]], vh: [1, 4] },
  toth: { label: "Toth", pnames: ["qs", "b", "t"],
    f: (th, p) => th[0] * th[1] * p / Math.pow(1 + Math.pow(Math.max(th[1] * p, 0), th[2]), 1 / th[2]),
    init: s => [1.1 * s.qmax, 1 / s.pmid, 0.7],
    bounds: s => [[1e-9, 1e-9, 0.05], [1e6, 1e6, 5]], vh: [1] },
  freundlich: { label: "Freundlich", pnames: ["k", "n"],
    f: (th, p) => th[0] * Math.pow(Math.max(p, 0), 1 / th[1]),
    init: s => { const n = Math.min(10, Math.max(0.2, s.fexp)); return [s.qmax / Math.pow(s.pmax, 1 / n), n]; },
    bounds: s => [[1e-9, 0.2], [1e6, 10]], vh: [0] },
  bet: { label: "BET / GAB (Type II/III)", pnames: ["qs", "bs", "bl"],
    f: (th, p) => th[0] * th[1] * p / ((1 - th[2] * p) * (1 - th[2] * p + th[1] * p)),
    init: s => [0.5 * s.qmax, 10 / s.pmid, 0.3 / s.pmax],
    bounds: s => [[1e-9, 1e-9, 1e-9], [1e6, 1e6, 0.95 / s.pmax]], vh: [1, 2] },
  typev: { label: "Type V (Langmuir + Hill step)", pnames: ["qs1", "b1", "qs2", "b2", "n"],
    f: (th, p) => {
      const t = Math.pow(Math.max(th[3] * p, 0), th[4]);
      return th[0] * th[1] * p / (1 + th[1] * p) + th[2] * t / (1 + t);
    },
    init: s => [0.25 * s.qmax, 1 / s.pmid, s.qmax, 1 / s.pmid, 4],
    bounds: s => [[1e-9, 1e-9, 1e-9, 1e-9, 1], [1e6, 1e6, 1e6, 1e6, 20]], vh: [1, 3] },
  unilan: { label: "UNILAN", pnames: ["qs", "b", "s"],
    f: (th, p) => th[0] / (2 * th[2]) *
      Math.log((1 + th[1] * Math.exp(th[2]) * p) / (1 + th[1] * Math.exp(-th[2]) * p)),
    init: s => [1.1 * s.qmax, 1 / s.pmid, 1],
    bounds: s => [[1e-9, 1e-9, 0.01], [1e6, 1e6, 10]], vh: [1] },
  temkin: { label: "Temkin-type (B ln(1+Ap))", pnames: ["B", "A"],
    f: (th, p) => th[0] * Math.log(1 + th[1] * p),
    init: s => [s.qmax / Math.log(1 + s.pmax / s.pmid), 1 / s.pmid],
    bounds: s => [[1e-9, 1e-9], [1e6, 1e6]], vh: [1] },
  da: { label: "Dubinin–Astakhov", pnames: ["qs", "E", "n", "p0"],
    // E in kJ/mol; A = RT ln(p0/p); q = qs exp[-(A/1000E)^n]; T explicit
    f: (th, p, T) => {
      if (p <= 0) return 0;
      const A = R_GAS * (T || 298.15) * Math.log(th[3] / p);
      if (A <= 0) return th[0];
      return th[0] * Math.exp(-Math.pow(A / (1000 * th[1]), th[2]));
    },
    init: s => [1.05 * s.qmax, 10, 2, 3 * s.pmax],
    bounds: s => [[1e-9, 0.5, 0.5, s.pmax * 1.0001], [1e6, 200, 6, 1e7]], vh: "intrinsic" },
  custom: { label: "Custom q(p) …", pnames: null }
};

/* ---------- custom expression: p + free parameters a,b,c,... ---------- */
function compileCustom(src) {
  const FUN = { exp: Math.exp, log: Math.log, ln: Math.log, sqrt: Math.sqrt,
                abs: Math.abs, tanh: Math.tanh, pow: Math.pow, min: Math.min, max: Math.max };
  const CONST = { pi: Math.PI, e: Math.E };
  const ids = [];
  const re = /\s*(\d+\.?\d*(?:[eE][+-]?\d+)?|[A-Za-z_][A-Za-z_0-9]*|\*\*|[-+*/^(),])/y;
  const toks = [];
  let i = 0;
  while (i < src.length) {
    re.lastIndex = i;
    const m = re.exec(src);
    if (!m) throw "Bad token near '" + src.slice(i, i + 8) + "'";
    toks.push(m[1]); i = re.lastIndex;
  }
  for (const t of toks)
    if (/^[A-Za-z_]/.test(t) && !FUN[t] && !(t in CONST) && t !== "p" && t !== "T" && !ids.includes(t))
      ids.push(t);
  if (ids.length === 0) throw "No free parameters found (use names like a, b, qs …)";
  if (ids.length > 9) throw "Too many parameters (max 9)";
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
      else if (t === "p") out.push("VP");
      else if (t === "T") out.push("VT");
      else out.push("A" + ids.indexOf(t));
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
  const f = (th, p, T) => {
    const st = [];
    for (const tk of out) {
      if (typeof tk === "number") st.push(tk);
      else if (tk === "VP") st.push(p);
      else if (tk === "VT") st.push(T || 298.15);
      else if (tk[0] === "A" && tk.length <= 2) st.push(th[+tk.slice(1)]);
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
  f(ids.map(() => 1), 0.123, 298.15);
  return { f, pnames: ids };
}

/* ---------- data statistics for heuristics ---------- */
function dataStats(p, q) {
  let qmax = 0, pmax = 0, pmin = Infinity;
  for (let i = 0; i < p.length; i++) {
    qmax = Math.max(qmax, q[i]); pmax = Math.max(pmax, p[i]);
    if (p[i] > 0) pmin = Math.min(pmin, p[i]);
  }
  // p at half loading (first crossing, sorted by p)
  const idx = Array.from(p.keys()).sort((a, b) => p[a] - p[b]);
  let pmid = Math.sqrt(pmin * pmax);
  for (let k = 1; k < idx.length; k++) {
    const a = idx[k - 1], b = idx[k];
    if (q[a] <= qmax / 2 && q[b] >= qmax / 2) {
      const w = (qmax / 2 - q[a]) / (q[b] - q[a] || 1);
      pmid = p[a] + w * (p[b] - p[a]);
      break;
    }
  }
  // Henry estimate: slope through lowest-p points
  const lo = idx.slice(0, Math.max(2, Math.floor(idx.length / 4)));
  let sxy = 0, sxx = 0;
  lo.forEach(i => { sxy += p[i] * q[i]; sxx += p[i] * p[i]; });
  const henry = sxx > 0 ? sxy / sxx : 1;
  // Freundlich exponent estimate from end-to-end log slope
  const a0 = idx[0], a1 = idx[idx.length - 1];
  let fexp = 2;
  if (q[a0] > 0 && q[a1] > 0 && p[a1] > p[a0])
    fexp = Math.log(p[a1] / p[a0]) / Math.log(q[a1] / q[a0]);
  return { qmax: qmax || 1, pmax: pmax || 1, pmin, pmid: Math.max(pmid, pmin), henry, fexp };
}

/* ---------- weighting ---------- */
function makeWeights(q, mode) {
  let qmax = 0;
  for (const v of q) qmax = Math.max(qmax, v);
  const qf = 1e-3 * (qmax || 1);
  return q.map(v => mode === "abs" ? 1 :
                    mode === "rel" ? 1 / Math.max(v, qf) :
                    1 / Math.sqrt(Math.max(v, qf)));
}

/* ---------- Levenberg-Marquardt on logit-bounded parameters ---------- */
function logit01(x) { return Math.log(x / (1 - x)); }
function sig(u) { return 1 / (1 + Math.exp(-u)); }

function lmFit(resFun, nRes, th0, lo, hi, opts) {
  // resFun(theta, out[nRes]); bounds strict; returns {theta, ssew, iters, ok}
  opts = opts || {};
  const k = th0.length;
  const u = new Float64Array(k);
  for (let j = 0; j < k; j++) {
    const x = Math.min(0.999999, Math.max(1e-6, (th0[j] - lo[j]) / (hi[j] - lo[j])));
    u[j] = logit01(x);
  }
  const theta = uu => Array.from(uu, (v, j) => lo[j] + (hi[j] - lo[j]) * sig(v));
  const r = new Float64Array(nRes), r2 = new Float64Array(nRes);
  const J = new Float64Array(nRes * k);
  const A = new Float64Array(k * k), g = new Float64Array(k), d = new Float64Array(k);
  const sse = uu => { resFun(theta(uu), r); let s = 0; for (let i = 0; i < nRes; i++) s += r[i] * r[i]; return s; };
  let S = sse(u);
  if (!isFinite(S)) return { theta: theta(u), ssew: Infinity, iters: 0, ok: false };
  let lam = 1e-3;
  const uTry = new Float64Array(k);
  let it = 0;
  for (; it < (opts.maxIter || 120); it++) {
    resFun(theta(u), r);
    for (let j = 0; j < k; j++) {
      const du = 1e-6 * Math.max(1, Math.abs(u[j]));
      const us = u[j];
      u[j] = us + du;
      resFun(theta(u), r2);
      u[j] = us;
      for (let i = 0; i < nRes; i++) J[i * k + j] = (r2[i] - r[i]) / du;
    }
    A.fill(0); g.fill(0);
    for (let i = 0; i < nRes; i++)
      for (let a = 0; a < k; a++) {
        g[a] += J[i * k + a] * r[i];
        for (let b = a; b < k; b++) A[a * k + b] += J[i * k + a] * J[i * k + b];
      }
    for (let a = 0; a < k; a++) for (let b = 0; b < a; b++) A[a * k + b] = A[b * k + a];
    let improved = false;
    for (let tries = 0; tries < 12; tries++) {
      // solve (A + lam*diag(A)) d = -g  (Cholesky-free: Gaussian elim, small k)
      const M = Float64Array.from(A);
      for (let a = 0; a < k; a++) M[a * k + a] += lam * Math.max(A[a * k + a], 1e-12);
      for (let a = 0; a < k; a++) d[a] = -g[a];
      let singular = false;
      for (let c = 0; c < k; c++) {
        let piv = c;
        for (let rr = c + 1; rr < k; rr++) if (Math.abs(M[rr * k + c]) > Math.abs(M[piv * k + c])) piv = rr;
        if (Math.abs(M[piv * k + c]) < 1e-300) { singular = true; break; }
        if (piv !== c) {
          for (let cc = 0; cc < k; cc++) { const t = M[c * k + cc]; M[c * k + cc] = M[piv * k + cc]; M[piv * k + cc] = t; }
          const t = d[c]; d[c] = d[piv]; d[piv] = t;
        }
        for (let rr = c + 1; rr < k; rr++) {
          const fmul = M[rr * k + c] / M[c * k + c];
          if (fmul !== 0) {
            for (let cc = c; cc < k; cc++) M[rr * k + cc] -= fmul * M[c * k + cc];
            d[rr] -= fmul * d[c];
          }
        }
      }
      if (!singular) {
        for (let c = k - 1; c >= 0; c--) {
          let s2 = d[c];
          for (let cc = c + 1; cc < k; cc++) s2 -= M[c * k + cc] * d[cc];
          d[c] = s2 / M[c * k + c];
        }
        for (let j = 0; j < k; j++) uTry[j] = u[j] + Math.max(-20, Math.min(20, d[j]));
        const Snew = sse(uTry);
        if (isFinite(Snew) && Snew < S) {
          u.set(uTry);
          const rel = (S - Snew) / Math.max(S, 1e-300);
          S = Snew;
          lam = Math.max(lam * 0.3, 1e-12);
          improved = true;
          if (rel < 1e-12) it = 1e9;   // converged
          break;
        }
      }
      lam *= 2.5;
      if (lam > 1e10) break;
    }
    if (!improved) break;
  }
  return { theta: theta(u), ssew: S, iters: Math.min(it, 120), ok: true };
}

/* ---------- assemble residual function for a dataset/mode ---------- */
function buildProblem(modelId, model, data, mode, weighting) {
  // data: {p:[], q:[], T:[]}; mode: "single" | "global"
  const N = data.p.length;
  const w = makeWeights(data.q, weighting);
  const stats = dataStats(data.p, data.q);
  let f = model.f, pn = model.pnames, init, lo, hi, vhInfo = null;
  const [bLo, bHi] = model.bounds ? model.bounds(stats) : [null, null];
  if (mode === "global" && model.vh !== "intrinsic" && model.vh) {
    // replace each vh param b -> b0*exp(U/RT): vector = base params with each
    // vh slot expanded to [b0, U]; map at eval time
    const baseInit = model.init(stats);
    const Tmid = data.T.reduce((a, b) => a + b, 0) / N;
    pn = []; init = []; lo = []; hi = [];
    const slots = [];   // for each base param: {idx in new vec, isVH}
    model.pnames.forEach((nm, j) => {
      if (model.vh.includes(j)) {
        slots.push({ j0: init.length, vh: true });
        pn.push(nm + "0", "U_" + nm + " (J/mol)");
        const U0 = 25000;
        init.push(baseInit[j] / Math.exp(U0 / (R_GAS * Tmid)), U0);
        lo.push(1e-15, 1e3); hi.push(1e6, 1.5e5);
      } else {
        slots.push({ j0: init.length, vh: false });
        pn.push(nm);
        init.push(baseInit[j]); lo.push(bLo[j]); hi.push(bHi[j]);
      }
    });
    vhInfo = { slots, base: model };
    f = (th, p, T) => {
      const bth = slots.map(sl => sl.vh ? th[sl.j0] * Math.exp(th[sl.j0 + 1] / (R_GAS * T)) : th[sl.j0]);
      return model.f(bth, p, T);
    };
  } else {
    init = model.init(stats); lo = bLo; hi = bHi;
  }
  const resFun = (th, out) => {
    for (let i = 0; i < N; i++)
      out[i] = w[i] * (f(th, data.p[i], data.T[i]) - data.q[i]);
  };
  return { f, pn, init, lo, hi, resFun, N, w, stats, vhInfo };
}

/* ---------- multi-start driver + statistics ---------- */
function fitModel(modelId, data, opts) {
  opts = opts || {};
  const weighting = opts.weighting || "rel";
  const mode = opts.mode || "single";
  const model = modelId === "custom"
    ? Object.assign({}, FIT_MODELS.custom, opts.custom)
    : FIT_MODELS[modelId];
  if (!model.f) return null;
  let prob;
  if (modelId === "custom") {
    const k = model.pnames.length;
    const w = makeWeights(data.q, weighting);
    prob = {
      f: model.f, pn: model.pnames,
      init: model.pnames.map(() => 1),
      lo: model.pnames.map(() => -1e6), hi: model.pnames.map(() => 1e6),
      resFun: (th, out) => {
        for (let i = 0; i < data.p.length; i++)
          out[i] = w[i] * (model.f(th, data.p[i], data.T[i]) - data.q[i]);
      },
      N: data.p.length, w, stats: dataStats(data.p, data.q), vhInfo: null
    };
  } else prob = buildProblem(modelId, model, data, mode, weighting);
  const k = prob.init.length;
  if (prob.N < k + 2) return { modelId, error: "too few points (" + prob.N + ") for " + k + " parameters" };
  // multi-start
  const starts = [prob.init];
  let seed = 12345;
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
  for (let s = 0; s < (opts.nStarts || 7); s++)
    starts.push(prob.init.map((v, j) => {
      const span = prob.hi[j] - prob.lo[j];
      if (prob.lo[j] > 0 || modelId !== "custom") {
        const vv = v * Math.pow(10, 2 * rnd() - 1);
        return Math.min(prob.hi[j] - 1e-12 * span, Math.max(prob.lo[j] + 1e-12 * span, vv));
      }
      return v + (2 * rnd() - 1) * Math.max(1, Math.abs(v));
    }));
  let best = null;
  for (const th0 of starts) {
    const r = lmFit(prob.resFun, prob.N, th0, prob.lo, prob.hi, opts);
    if (r.ok && (!best || r.ssew < best.ssew)) best = r;
  }
  if (!best) return { modelId, error: "fit failed" };
  // statistics
  const th = best.theta;
  const N = prob.N;
  const rr = new Float64Array(N);
  prob.resFun(th, rr);
  let sseU = 0, sst = 0, qbar = 0;
  data.q.forEach(v => qbar += v / N);
  for (let i = 0; i < N; i++) {
    const qm = prob.f(th, data.p[i], data.T[i]);
    sseU += (qm - data.q[i]) ** 2;
    sst += (data.q[i] - qbar) ** 2;
  }
  const rmse = Math.sqrt(sseU / N);
  const r2 = 1 - sseU / Math.max(sst, 1e-300);
  const aicc = (N - k - 1 > 0)
    ? N * Math.log(best.ssew / N) + 2 * k + 2 * k * (k + 1) / (N - k - 1)
    : Infinity;
  // CI via weighted Jacobian in theta space
  const Jt = new Float64Array(N * k), rB = new Float64Array(N);
  for (let j = 0; j < k; j++) {
    const dth = 1e-6 * Math.max(Math.abs(th[j]), 1e-12);
    const ths = th.slice(); ths[j] += dth;
    prob.resFun(ths, rB);
    for (let i = 0; i < N; i++) Jt[i * k + j] = (rB[i] - rr[i]) / dth;
  }
  const AtA = [];
  for (let a = 0; a < k; a++) {
    AtA.push([]);
    for (let b = 0; b < k; b++) {
      let s = 0;
      for (let i = 0; i < N; i++) s += Jt[i * k + a] * Jt[i * k + b];
      AtA[a].push(s);
    }
  }
  const cov = invSym(AtA);
  let ci = null, corrWarn = [];
  if (cov) {
    const s2 = best.ssew / Math.max(N - k, 1);
    ci = th.map((v, j) => 1.96 * Math.sqrt(Math.max(cov[j][j] * s2, 0)));
    for (let a = 0; a < k; a++)
      for (let b = a + 1; b < k; b++) {
        const rho = cov[a][b] / Math.sqrt(Math.max(cov[a][a] * cov[b][b], 1e-300));
        if (Math.abs(rho) > 0.98)
          corrWarn.push(prob.pn[a] + "–" + prob.pn[b] + " " + (100 * Math.abs(rho)).toFixed(1) + "% correlated");
      }
  }
  return { modelId, label: model.label, pn: prob.pn, theta: th, ci, corrWarn,
           ssew: best.ssew, sseU, rmse, r2, aicc, k, N, iters: best.iters,
           f: prob.f, mode, weighting, stats: prob.stats, vhInfo: prob.vhInfo };
}

function invSym(A) {
  const k = A.length;
  const M = A.map((row, i) => row.concat(Array.from({ length: k }, (_, j) => i === j ? 1 : 0)));
  for (let c = 0; c < k; c++) {
    let piv = c;
    for (let r = c + 1; r < k; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    if (Math.abs(M[piv][c]) < 1e-300) return null;
    [M[c], M[piv]] = [M[piv], M[c]];
    const dv = M[c][c];
    for (let cc = 0; cc < 2 * k; cc++) M[c][cc] /= dv;
    for (let r = 0; r < k; r++) {
      if (r === c) continue;
      const fmul = M[r][c];
      if (fmul !== 0) for (let cc = 0; cc < 2 * k; cc++) M[r][cc] -= fmul * M[c][cc];
    }
  }
  return M.map(row => row.slice(k));
}

/* ---------- isosteric heat from a global fit: qst(q) = -R dlnp/d(1/T) ---- */
function isostericHeat(fitRes, data, nq) {
  if (!fitRes || fitRes.mode !== "global") return null;
  const Tmid = data.T.reduce((a, b) => a + b, 0) / data.T.length;
  const T1 = Tmid - 10, T2 = Tmid + 10;
  const pmax = Math.max(...data.p);
  const f = fitRes.f, th = fitRes.theta;
  const qcap = Math.min(f(th, pmax * 10, Tmid), f(th, pmax * 10, T2));
  const out = { q: [], qst: [] };
  const solveP = (q, T) => {
    let lo = 1e-12, hi = pmax * 10;
    if (f(th, hi, T) < q) return NaN;
    for (let it = 0; it < 80; it++) {
      const m = Math.sqrt(lo * hi);
      if (f(th, m, T) < q) lo = m; else hi = m;
    }
    return Math.sqrt(lo * hi);
  };
  for (let i = 1; i <= (nq || 40); i++) {
    const q = qcap * i / ((nq || 40) + 1) * 0.95;
    const p1 = solveP(q, T1), p2 = solveP(q, T2);
    if (!isFinite(p1) || !isFinite(p2) || p1 <= 0 || p2 <= 0) continue;
    const qst = R_GAS * Math.log(p2 / p1) / (1 / T1 - 1 / T2);
    out.q.push(q); out.qst.push(qst / 1000);   // kJ/mol
  }
  return out;
}

/* ---------- model-free isosteric heat from the raw data ----------
   For each loading level q in the range covered by >=2 temperature groups,
   interpolate ln p LINEARLY in q within each group (standard practice —
   far better behaved than interpolating p on the steep branch), then
   qst = -R * slope of ln p vs 1/T (least squares across groups). kJ/mol. */
function expIsostericHeat(data, nq) {
  const Tg = [...new Set(data.T.map(v => Math.round(v * 100) / 100))].sort((a, b) => a - b);
  if (Tg.length < 2) return null;
  const groups = Tg.map(Tv => {
    const idx = [];
    for (let i = 0; i < data.T.length; i++)
      if (Math.abs(data.T[i] - Tv) < 0.05 && data.p[i] > 0 && data.q[i] > 0) idx.push(i);
    idx.sort((a, b) => data.p[a] - data.p[b]);
    return { T: Tv, p: idx.map(i => data.p[i]), q: idx.map(i => data.q[i]) };
  }).filter(gr => gr.p.length >= 2);
  if (groups.length < 2) return null;
  let qlo = 0, qhi = Infinity;
  groups.forEach(gr => {
    qlo = Math.max(qlo, Math.min.apply(null, gr.q));
    qhi = Math.min(qhi, Math.max.apply(null, gr.q));
  });
  if (!(qhi > qlo)) return null;
  const lnPat = (gr, qv) => {
    for (let i = 1; i < gr.q.length; i++) {
      const a = gr.q[i - 1], b = gr.q[i];
      if ((a - qv) * (b - qv) <= 0 && a !== b) {
        const w = (qv - a) / (b - a);
        return Math.log(gr.p[i - 1]) + w * (Math.log(gr.p[i]) - Math.log(gr.p[i - 1]));
      }
    }
    return NaN;
  };
  const out = { q: [], qst: [] };
  const NQ = nq || 25;
  for (let k = 1; k <= NQ; k++) {
    const qv = qlo + (qhi - qlo) * k / (NQ + 1);
    const xs = [], ys = [];
    for (const gr of groups) {
      const lp = lnPat(gr, qv);
      if (isFinite(lp)) { xs.push(1 / gr.T); ys.push(lp); }
    }
    if (xs.length < 2) continue;
    const n = xs.length;
    let sx = 0, sy = 0, sxy = 0, sxx = 0;
    for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxy += xs[i] * ys[i]; sxx += xs[i] * xs[i]; }
    const den = n * sxx - sx * sx;
    if (Math.abs(den) < 1e-300) continue;
    const slope = (n * sxy - sx * sy) / den;
    out.q.push(qv);
    out.qst.push(-R_GAS * slope / 1000);
  }
  return out.q.length ? out : null;
}

if (typeof globalThis !== "undefined") {
  globalThis.IFIT = { FIT_MODELS, compileCustom, dataStats, makeWeights,
                      lmFit, fitModel, isostericHeat, expIsostericHeat, R_GAS };
}
