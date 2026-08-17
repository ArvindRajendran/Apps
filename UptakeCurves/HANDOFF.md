# UptakeCurves — technical handoff

Implementation notes for future development sessions. The user-facing
description and validation table live in `README.md`; this file records the
architecture, the numerical decisions and their reasons, and the extension
points.

## File layout

- `UptakeCurves.html` — the entire app. Three parts:
  1. CSS + markup (suite conventions: CSS variables, light/dark via
     `prefers-color-scheme`, `.panel`/`.pcard` cards, footer credit block).
  2. `<script id="uc-engine">` — the solver engine. **This same text is
     used twice**: parsed normally for the main thread (exposes
     `globalThis.UPTAKE`), and read back via `textContent` to build a Blob
     Worker. Keep it dependency-free and DOM-free.
  3. Main `<script>` — UI, plotting, worker protocol, localStorage.
- `validation/` — Python reference + jsc harness + engine snapshot
  (see its README). `uptake_engine.js` there is a *copy*; sync manually.

## Engine (`UPTAKE.*`)

State conventions: loading q in mol/kg, gas as mole fraction y; c = y·c_T
with c_T = P/RT. Radii normalized to [0,1]; cell-centred FV with N cells,
faces at i/N; cell "volumes" V_i = (ρ_{i+1}³−ρ_i³)/3 (4π dropped
throughout, it cancels). Surface Dirichlet closure: half-cell face flux
F_s = Γ_face·(q_s − q_N)/(h/2). Symmetry at centre = zero flux (face 0
never written).

- `MODELS` — isotherm menu; same functional forms as the sorption
  chromatography app but in mol/kg (Do & Do without that app's `S`
  scaling). `makeIsotherm({id, pvals, expr})` → f(y).
- `gammaTable(f, ylo, yhi)` — Darken factor Γ(q) = q/(y·f′(y)) tabulated on
  3000 points of [max(ylo/2, yhi·10⁻⁶), yhi], binary-search interpolation,
  clamped at the ends; also returns `mono` (false → app disables Darken
  with a warning rather than feeding garbage into the PDE).
- `buildModel(cfg)` — returns `{n, x0, rhs(x,out), sparsity, uptake(x),
  profiles(x)}` for micro / macro / bidisperse. Bidisperse ordering:
  block j holds [q_{j,0..Nc−1}, y_j] (block size Nc+1) — keeps the
  Jacobian banded-ish and the coloring small. The macropore ODE is written
  explicitly (no mass matrix) because the crystal flux telescopes:
  dq̄_j/dt = 3·a_μ·F_s,j.
- Two protective details, both load-bearing:
  - **f is clamped to [min(y0,y1), max(y0,y1)]** (max principle) so Newton
    wander can never hit isotherm poles (Langmuir at y = −1/b, BET at
    y = 1/b_l). Safe *only* because the integrator is L-stable (see below).
  - **fp (macro-only β) differentiates the RAW isotherm at a point nudged
    2h inside the interval** — differentiating the clamped f would halve
    the centred difference at the interval edges (the same halving bug once
    hit in the sorption app's Do & Do derivative).
- `integrateCN(model, times, opts)` — despite the legacy name this is
  **TR-BDF2** (γ = 2−√2; both stages share M = I − 0.29289·Δt·J).
  Quasi-Newton with colored-FD Jacobian (greedy column coloring on the
  model sparsity, ~4–6 colors), dense LU, Jacobian reuse with rebuild on
  stall (residual not halving after it > 2) or NaN (restart stage from the
  step start), sub-step halving up to 64× as last resort. `opts.keepStates`
  = list of step indices to snapshot.
- `solveUptake(cfg, onProgress)` — √t-spaced grid t_k = (k/K)²·T
  (resolves the √t-linear early regime without adaptivity), K default 280;
  61 profile snapshots; LDF is closed-form (no ODE). `autoTend(cfg)` sums
  half of each relevant time constant, with a chord-averaged Γ for Darken.

**Why TR-BDF2 and not Crank–Nicolson** (the most important lesson here):
CN is A-stable but not L-stable. When Δt ≫ τ_macro (easy in the
micropore-controlled bidisperse limit) the stiff macropore mode rings
undamped (R(∞) = −1), and the y-clamp *rectifies* that oscillation —
qs = f(clamp(y)) is biased low — which shifted uptake curves by O(0.1).
An earlier "fix" attempt (clamping alone, keeping CN) silently produced
those wrong curves while every non-stiff check still passed. TR-BDF2
(R(∞) = 0) kills the mode in one step and makes the clamp inert.

## Worker protocol / UI

- Main thread posts the full cfg (isotherm params included — the custom
  expression is compiled inside the worker too; parser is in the engine).
- Worker replies `{type:'prog'|done|err}`. If a run is in flight when
  parameters change, the worker is terminated and respawned (cheap) rather
  than queued.
- Sync fallback path (`setTimeout` + direct `solveUptake`) exists for
  contexts where Blob workers fail (some `file://` configurations).
- Auto-extension: if t_end was auto and U(T) < 0.97, re-run once or twice
  with T×3 (S-shaped isotherms have slow tails because Γ collapses near
  the plateau).
- Uptake normalization is (m_t − m₀)/(m_∞ − m₀) — rises 0→1 for both
  adsorption and desorption; the y-axis label spells this out. m includes
  the macropore gas inventory ε_p·c as well as ρ_p·q̄ (what a balance sees).
- localStorage key `uptake_v1`. Debug handle: `globalThis.__uc`
  ({st, sol, cfg}).
- Grids: N_c, N_p, K from the numerics panel apply to **every** model
  (micro uses N_c, macro uses N_p, bidisperse both; defaults 20/20/280,
  clamped ≥ 8 in `buildCfg()`). The active grid is echoed in the summary
  line.
- CSV export (`downloadCSV` + `csvMeta`): uptake curve (t, √t, U) and both
  profiles at the current slider snapshot, each prefixed with `#` metadata
  lines (model, isotherm + params, T/P/y-step, transport parameters,
  q₀/q₁). Profile buttons are disabled when the model has no such field.

## Known limitations / deliberate scope

- Spheres only (user chose to start there). Slab/cylinder = change the
  geometry factors in `fvGeom` (face "areas" ρ^s, volumes ∫ρ^s dρ for
  s = 0,1,2) and the uptake weights — the rest is geometry-agnostic.
- Equilibrium at the outer surface: no external film. Adding film
  resistance = replace the particle-surface Dirichlet half-cell flux with
  a Robin closure k_f(c_b − c_s); one line in the macro/bidisperse rhs
  plus an input field.
- Isothermal. Heat effects (the usual cause of long gravimetric tails)
  would need a particle energy balance coupled to q̄ — a real extension,
  not a tweak.
- Single step only; no multi-step staircase or frequency response.
- Darken applies to the micropore diffusivity only (macropore is a
  gas-phase mechanism); macro-only + Darken checkbox is deliberately a
  no-op (checkbox hidden for that model).
- The Do & Do desorption branch uses max(ads, des) hysteresis composite;
  for an uptake step you pick the branch matching the step direction —
  the app does not switch branches mid-run.

## Editing workflow (from this project's sessions)

- No node on this machine — headless JS via JavaScriptCore `jsc` (path in
  `validation/README.md`).
- The preview server sandbox cannot read this Google Drive folder — serve
  a *copy* from the session scratchpad (`launch.json` "apps" config →
  `serve_apps.py` pattern) and re-copy after each edit.
- If you change the engine block in the HTML, re-copy it into
  `validation/uptake_engine.js` and re-run both validation suites before
  shipping.

## Parameter import

The app accepts `#import=<urlencoded JSON>` from the Isotherm Fitter
(p-basis parameters). Affinities are converted to the app's
mole-fraction basis at the CURRENT `st.Pbar` (b×P; quadratic b2×P²) —
re-import after changing P if exactness matters. Applied after the
localStorage restore; hash cleared; toast shown.
