# IAST — technical handoff

Implementation notes for future sessions. User-facing description and the
validation table are in `README.md`.

## File layout

- `IAST.html` — whole app. `<script id="ia-engine">` holds the solver
  engine (`globalThis.IAST`), DOM-free; main `<script>` holds UI/plots.
  No worker: a full 61×61 sweep is ~15–60 ms, computed synchronously with
  a 300 ms debounce.
- `validation/` — `iast_ref.py` (Python reference), `iast_engine.js`
  (snapshot of the embedded engine — re-copy when the engine changes),
  `iast_test.js` (jsc harness with baked anchors).

## Engine

- `IAST_MODELS` — 11 models, `make(par) → {f, psi, henry, pole, psiKind}`.
  ψ analytic for all but Toth/custom. Henry constants encode the
  consistency flags: finite (ok), `Infinity` (Freundlich n>1, Sips n<1 —
  flagged in UI), `0` (super-linear onset).
- `buildPsiTable(f, pMax)` — numeric ψ: cumulative Newton–Cotes on pairs
  of ln-p intervals (h/12·(5f₀+8f₁−f₂) then h/12·(−f₀+8f₁+5f₂); 4th
  order), 20 000 intervals from p = 10⁻¹⁰ to pMax; below the table the
  tail ∫₀^{p₀} q/p dp = q(p₀)/α with α the local power-law exponent
  estimated from q(p₀)/q(p₀/2). Accuracy ~10⁻⁷ rel except near poles
  (irrelevant: pole models all have analytic ψ).
- `makeComponent(cfg, pMaxNeeded)` — assembles f/ψ, clips at
  `pCap = pole·(1−10⁻⁹)` or 10¹⁰ bar, runs a 240-point monotonicity scan,
  estimates Henry numerically for custom expressions.
- `iastSolve(P, y1, c1, c2)` — the core: G(x₁) = ψ₁(Py₁/x₁) −
  ψ₂(Py₂/(1−x₁)) is strictly decreasing with G(0⁺)=+∞, G(1⁻)=−∞
  (fixed bracket!). 50 bisections + ≤6 Newton steps using
  dG/dx₁ = −q₁/x₁ − q₂/(1−x₁). Only forward ψ evaluations — no inversion
  tables, no initial guess. Edge cases y₁∈{0,1} short-circuit to pure.
  Returns p_i⁰ and a `capped` flag when a hypothetical pressure hits pCap.
- `findSReversal(c1, c2)` — returns `{ps, degenerate}`. Sign-scan +
  bisection on ψ₁(p) − ψ₂(p) over a 400-point log grid, but with two
  guards added after a real user bug: (i) **degeneracy detection** —
  if max|ψ₁−ψ₂|/max(ψ₁,ψ₂) < 10⁻⁷ over the grid, the pair is
  IAST-indistinguishable (identical model+parameters gives ψ₁ ≡ ψ₂, and
  the naive scan then flags all 400 samples as "crossings", which the UI
  printed as a wall of numbers); (ii) **transversality validation** —
  each bisected root must have opposite-signed flanks at p*·(1±0.05) with
  magnitude above 10⁻⁹+10⁻⁷ψ, killing noise wiggles; list capped at 10,
  UI shows at most 3. Regression: test J (identical quad–quad →
  degenerate, x=y to machine precision; quad(3,5,10) vs quad(3,8,3) →
  p* = 3/7 bar exactly, from 5p+10p² = 8p+3p²).
- `iastSweep`, `xyCurve` (with x=y crossing scan), `extendedLangmuir`.

**The key theoretical structure baked into the app** (derive before
touching the azeotrope/reversal code): from Py_i = x_i p_i⁰,
S₁₂ = (x₁/y₁)/(x₂/y₂) = p₂⁰/p₁⁰ — a function of ψ (and T) only.
Hence S = 1 ⇔ p₁⁰ = p₂⁰, and closure then gives P = p₁⁰ = p₂⁰ = p*.
Consequences: (i) binary IAST cannot produce composition azeotropes;
(ii) at P = p* the entire x–y curve lies on the diagonal; (iii) the S = 1
contour in (P, y₁) is exactly the horizontal line P = p*. The original
test suite briefly asserted an azeotrope for a crossing Langmuir pair —
wrong physics, caught by the code; check I now asserts the correct
statements.

## Python reference (`validation/iast_ref.py`)

Deliberately different algorithm: bisection on p₁⁰ with upward bracket
expansion (x₁ = Py₁/p₁⁰, p₂⁰ = Py₂/(1−x₁)); adaptive Simpson in ln p for
numeric ψ (tol 10⁻¹²). Check H installs nothing at runtime but expects
`pyiast`+`pandas` (`pip install --user pyiast scipy pandas` — worked on
this machine with numpy 2.0.2/pandas 2.3.3); pyIAST ModelIsotherms are
built by fitting exact synthetic Langmuir data (constructor requires a
DataFrame; bypassing `__init__` breaks its internals). Skips gracefully
if the import fails. Run with `MPLCONFIGDIR=/tmp/... python3 iast_ref.py`
(pyiast imports matplotlib; the default cache dir is not writable here).

## UI notes

- Inspection point (P, y₁): heatmap hover (transient) / click (pin).
  `drawInspection()` re-solves ~400 IAST points per frame — cheap.
- Heatmap: S uses log₁₀ with a diverging palette **only when the range
  straddles S=1** (`useDiv`), sequential viridis otherwise; colorbar
  labels show actual S values. S=1 contour drawn as the horizontal p*
  line(s).
- Extended-Langmuir overlay only offered when both models are Langmuir
  (`ia-elwrap` visibility).
- 3D view: hand-rolled orthographic projection, quad mesh decimated to
  ~40×40, painter-sorted; yaw/pitch by pointer drag; drawn only when the
  `<details>` is open.
- CSV: `meta()` writes model+params (or the custom expression) into `#`
  headers.
- localStorage `iast_v1`; debug handle `globalThis.__ia`.

## Deliberate scope / extension points

- Binary only. The solver generalizes to n components (ψ-bisection with
  Σ Py_i/p_i⁰(ψ) = 1 stays monotone) but every 2D view is binary-specific.
- Ideal adsorbed phase. **RAST** is the natural phase 2: route activity
  coefficients γ_i(x, ψ) into Py_i = γ_i x_i p_i⁰ — the G(x₁) solve stays
  1-D for spreading-pressure-independent γ models (Margules/Wilson), but
  ψ-dependent γ (the thermodynamically complete form) couples the loop.
- Ideal gas (p not fugacity); single T (parameters are at fixed T).
- No isotherm fitting — pairs naturally with a future Isotherm Fitter app.

## Environment reminders

- No node: jsc at /System/Library/Frameworks/JavaScriptCore.framework/
  Versions/Current/Helpers/jsc.
- Preview server can't read this Google Drive folder: serve a scratchpad
  copy (launch.json "apps" → serve_apps.py, port 8642); the scratchpad is
  session-scoped, so `serve_apps.py` may need recreating.
- The Browser-pane compositor occasionally wedges on scroll; the page JS
  keeps running (javascript_tool still works). Use a tall viewport
  (resize_window) for full-page screenshots instead of scrolling.

## Parameter import

The app accepts `#import=<urlencoded JSON>` (`{v:1, model, pvals, slot}`)
from the Isotherm Fitter: applied right after the localStorage restore,
hash cleared via replaceState, 6-s toast shown. Slot 1/2 → component
1/2. Keep the import block in sync if `st.m1/p1` naming changes.
