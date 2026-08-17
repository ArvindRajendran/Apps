# Isotherm Fitter — technical handoff

Implementation notes for future sessions. User-facing docs + validation
table in `README.md`.

## File layout

- `IsothermFitter.html` — whole app; `<script id="if-engine">` = the
  fitting engine (`globalThis.IFIT`, DOM-free), main `<script>` = UI.
  No worker: 15 models on ~60 points fit in 30–300 ms, debounced 450 ms.
- `validation/` — `fitter_ref.py` (scipy reference), `fitter_engine.js`
  (snapshot of the embedded engine — keep in sync), `fitter_test.js`
  (jsc harness; the seeded noisy datasets are PRINTED BY the Python run
  and pasted as JS literals, so both engines fit byte-identical data).

## Engine (`IFIT.*`)

- `FIT_MODELS` — 15 models: `{label, pnames, f(th,p,T), init(stats),
  bounds(stats), vh}`. `vh` lists the affinity-parameter indices that
  become b₀·e^(U/RT) in global mode; `"intrinsic"` (Dubinin–Astakhov)
  means the model uses T explicitly and is fitted globally unchanged.
  Bounds are data-dependent where physics demands it (anti-Langmuir/BET
  pole parameters capped at 0.95/p_max).
- `compileCustom(src)` — expression → `{f, pnames}`; free identifiers
  (not p, T, functions, pi, e) become parameters (max 9), fitted in
  *linear* space with wide bounds (they may legitimately be negative).
- `dataStats` — q_max, p_max, p_mid (half-loading crossing), Henry slope
  (regression through origin on the lowest-p quartile), Freundlich
  exponent estimate; feeds all `init` heuristics.
- `lmFit(resFun, nRes, th0, lo, hi)` — LM with the **logit-bounded
  transform** θ = lo + (hi−lo)·σ(u): strict interior feasibility, no
  active-set logic. FD Jacobian in u; damped normal equations via
  Gaussian elimination with pivoting; λ ∈ [1e-12, 1e10], ×0.3 / ×2.5;
  step clamp |δu| ≤ 20; converges on relative SSE change < 1e-12.
- `fitModel(modelId, data, opts)` — multi-start driver: heuristic init +
  `nStarts` (default 7) randomized starts (positive params × 10^U(−1,1),
  deterministic LCG seed → reproducible). Statistics: weighted SSE
  (objective), unweighted SSE/RMSE/R², AICc **on the weighted SSE**
  (models compared under the chosen objective), CI = ±1.96σ from the
  weighted Jacobian in θ-space with s² = SSE_w/(N−k), pairwise
  correlation warnings at |ρ| > 0.98. Returns Infinity AICc when
  N ≤ k+1.
- Global mode (`opts.mode = "global"`): `buildProblem` expands each vh
  parameter slot into (b₀, U) with bounds b₀ ∈ [1e-15, 1e6],
  U ∈ [1e3, 1.5e5] J/mol, init U = 25 kJ/mol and b₀ back-computed from
  the heuristic b at the mean temperature. **b₀ and U are always
  strongly correlated** (the van 't Hoff lever arm is short) — the >98%
  correlation warning on b₀–U pairs is expected and honest, not a bug.
- `isostericHeat(fitRes, data, nq)` — q_st(q) from the fitted global
  model: solve p(q,T) by log-bisection at T_mid ± 10 K, Clausius–
  Clapeyron; returns kJ/mol.

## Basis handling and hand-off (added after v1)

- `st.pu` may be a pressure unit (converted to bar at parse) or a
  concentration unit (`molm3`, `molL`, `gL` — values fitted AS ENTERED,
  no conversion; `isConc()/xUnit()/xSym()` drive labels/CSV/status). In
  c-basis a global-fit U is the c-basis internal energy (the p-basis
  fitted U differs by ≈ +RT̄ — demonstrated by the CO2/13X example).
- `expIsostericHeat(data, nq)` in the engine: model-free q_st — ln p
  linearly interpolated in q per temperature group, LSQ slope over 1/T;
  restricted-range check in the tests because the plateau region is
  ill-conditioned by nature.
- Hand-off: `handOff(R, key)` builds `{v:1, model, pvals, slot, note}`
  and opens `<target>.html#import=<encodeURIComponent(JSON)>`. Targets in
  `HO_URL`; p-basis → IAST slots + UptakeCurves; c-basis → step/pulse
  (+ Systems slots for Langmuir). Multi-T: `baseThetaAt(R, Tsel)`
  evaluates vh params b₀e^(U/RT) at the chosen T (global) or picks the
  nearest per-T fit. Receivers live in the five target apps right after
  their localStorage restore: parse hash, `history.replaceState` to clear
  it, apply to their native state (`st.pvals["model.key"]` maps in the
  chromatography apps; `st.p1/p2[model]` in IAST; y-basis conversion at
  `st.Pbar` in UptakeCurves — b×P, quad b2×P²), then a 6-s toast. If a
  receiver's state schema changes, its import block must follow.
- BAAM targets (`baamA`/`baamB` → `handOffBaam`): langmuir/dsl build an
  extended-DSL payload `{app:"baam", mode:"edsl", slot:"co2"|"n2", pvals}`
  with affinities in m³/mol — global van 't Hoff fits pass (b₀, U) exactly
  in c-basis, p-basis converted via U_c = U_p − R·T̄ and
  b₀_c = b₀_p·R_B·T̄·e; non-global fits pass b(T_sel) with dU = 0 and pin
  BAAM's T (payload.T). g/L basis is refused (needs molar mass). Other
  p-basis suite models reuse the IAST-style payload with slot comp1/comp2
  aimed at BAAM.html.

## Fit-quality philosophy (baked in, keep it)

SSE alone always ranks the biggest model first; the app therefore ranks
by AICc and displays ΔAICc. Per-T aggregation charges k_total = k·n_T so
per-temperature fitting also pays its parameter bill vs global van 't
Hoff. The B-check in the test suite pins this: DSL achieves lower SSE
than Langmuir on Langmuir data yet loses on AICc — if that ever flips,
something broke.

## UI notes

- Parse: lines → numeric tokens; 3 columns = (T, p, q), 2 = (p, q) with
  the T field; non-numeric lines skipped (headers/comments). Units
  converted at parse time (wt%/mg/g need molar mass; cm³STP/g uses
  22.414 L/mol). `st.excluded` holds row indices unticked in the data
  table; loading new data resets exclusions, selection, and overlay set.
- Per-T mode stores `perT: [{T, r}]` per model; `evalModel` picks the
  nearest-T fit when drawing. Global stores one `global: r`.
- Colors: model = MCOL by overlay order; temperature = TCOL by group.
- Suite hand-off block: for the 10 models shared with IAST/UptakeCurves
  the parameter object uses those apps' key names verbatim.
- Debug handle `globalThis.__if` (incl. `views`/`geo`); localStorage `isofit_v1`.
- Plot interactivity: `wireZoom(cvId, id, redraw, onClick?)` — box zoom
  (pointer drag > 7 px), wheel zoom about the cursor (log-aware on x),
  dbl-click reset; per-plot view overrides in `views` (session-only,
  cleared on new data), current mappings in `geo` (each draw stores its
  final ranges there — keep that invariant when touching draws). Legend
  hit-zones from `legendHits`; click toggles `st.hidden` (models) /
  `st.hiddenT` (temperature groups) — separate from the overlay checkbox
  set, both persisted. Curves/points are clipped to the axes rect so
  zooming doesn't spill ink outside the frame. PNG export =
  `canvas.toDataURL`; per-figure CSVs mirror exactly what is plotted.

## Deliberate scope / extension points

- Virial isotherm excluded (p(q)-implicit — doesn't fit the q(p)
  framework); Jensen–Seaton excluded (niche; expressible via custom).
- No uncertainty propagation to derived quantities (Henry constants,
  q_st bands) — a natural v2 along with bootstrap CIs.
- No isotherm-file import from the NIST ISODB (parked idea; see the
  Apps-suite memory note: the API is CORS-blocked, the GitHub mirror
  isn't — an "import from mirror" button is feasible if ever wanted).

## Environment reminders

- jsc for headless JS tests (path in validation/README.md); no node.
- scipy/pandas/pyiast installed user-level (numpy 2.0.2); the reference
  needs only scipy.
- Preview server: serve a scratchpad copy (launch.json "apps" →
  serve_apps.py, port 8642; recreate the script if the scratchpad was
  cleaned). Browser-pane compositor can wedge on scroll — use a tall
  resize_window viewport for full-page screenshots.
