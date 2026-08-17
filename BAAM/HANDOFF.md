# BAAM — technical handoff

Implementation notes for future sessions. User-facing docs + validation table
in `README.md`.

## File layout

- `BAAM.html` — whole app. `<script id="ba-engine">` = the engine
  (`globalThis.BAAMENG`, DOM-free), main `<script>` = UI. No worker: the
  extended-DSL path integrates in ~20 ms; the IAST route in ~0.5 s
  (ΔP floored at 2.5·10⁻⁴ bar in that mode).
- `validation/` — `baam_ref.py` (scipy reference), `baam_engine.js`
  (snapshot of the embedded block — **keep in sync**), `baam_test.js`
  (jsc harness, 210 checks).
- Sections 1–2 of the engine (parseExprP, IAST_MODELS, buildPsiTable,
  makeComponent, iastSolve) are **verbatim copies from the IAST app's
  ia-engine block** — if that app's engine changes, re-copy.

## Engine (`BAAMENG.*`)

- `edslMix(co2, n2, T)` — extended DSL, **concentration basis**
  (b₀/d₀ in m³/mol, ΔU in J/mol stored NEGATIVE for exothermic; UI fields
  show −ΔU in kJ/mol, positive). `loadings(P, y)` returns
  `[q1, q2, ∂q1/∂P, ∂q1/∂y, ∂q2/∂P, ∂q2/∂y]` with analytic partials
  (chain rule through c₁ = Py/RT, c₂ = P(1−y)/RT).
- `iastMix(comp1, comp2, T)` — same interface over two `makeComponent`
  instances (partial-pressure basis, bar); partials by central FD
  (δP = max(1e-6, 1e-6·P), δy = 1e-7 — same as the Python ref). Also
  exposes `q(P, y)` (single IAST solve, no partials) — the fast path used
  by `baInv` and the family-curve plotting.
- `blowPath(mix, {T, rho, eps, PH, yF, Plow, dP, eta})` — BLO/EVAC ODE
  dy/dP = (a₁y−a₂)/(f₂−f₁y), **fixed-step RK4** on the uniform ΔP grid
  (default 1e-4 bar), y clamped to [0,1] (y = 1 is an invariant manifold —
  the clamp is benign, unlike the UptakeCurves CN story). Returns the path
  plus cumulative removed moles (telescoping inventory differences — exact
  for any grid) and vacuum work W (MATLAB convention: each slice's dN
  charged at the slice-END pressure, only where P < 1 bar).
- `cycleKPIs(mix, path, {Pint, Plow, yF, mode, cache})` — LPP or FP.
  LPP/FP reduce to a 1-D root find in the end-of-pressurization composition
  (`baRoot`: bisection with a 400-point bracket-scan fallback); ADS is a
  **linear** 2×2 solve (the MATLAB used fsolve, but it's linear once y_δ is
  known). Compression work (P_H > 1 bar) uses η at 1 bar. `opt.cache`
  (object) memoizes the pressurization+ADS solve per (Plow-index, mode, yF) —
  essential for grid scans in IAST mode (they are P_INT-independent).
- `gridScan`, `rMax`; constants `BA_RCUT_EBAAM = 124.5`,
  `BA_EN_SCALE_EBAAM = 1.58`, `BA_RCUT_2019 = 110.25`,
  `BA_EN_SCALE_2019 = [1.1446, 66.528]`.

## Pressure levels and the ODE grid (v2 — read before touching ΔP)

`pullLight()` is the **only** place P_H > P_INT > P_LOW is enforced: it
rewrites the slider `min`/`max` attributes on every pull (P_LOW ≤
min(0.9, P_H − 0.03), P_INT ∈ [P_LOW + 0.01, P_H − 0.01]), so an out-of-order
value is clamped by the range element itself. Verified against six adversarial
sequences including P_H collapsing below P_INT. P_H has both a slider
(`ba-PHs`, 0.2–5 bar, step 0.05) and a number box (`ba-PH`, canonical — the
slider listener writes into it); typing an exact P_H is what the validation
cases need.

**Do not scale ΔP with P_H.** It was tried (to hold the step count constant)
and reverted: the vacuum work is a first-order Riemann sum over the pressure
slices, so its value is grid-dependent at O(ΔP) — scaling shifted the energy
by up to 1.2e-3 relative, and worse, design pressures stopped landing on the
grid (`baIdx` rounds to the nearest point silently, giving a 1.2e-3 error in
WC/Re at P_H = 3). With ΔP fixed at 1e-4, every slider pressure (a multiple of
0.005 bar) is exactly on the grid for any P_H, and the JS engine matches the
Python reference to 2e-8 at P_H = 2/3/5. Cost at P_H = 5: 50k RK4 steps =
32 ms (dual-site) or 0.9 s (IAST). `dPeff()` only raises ΔP if the user types
one so fine that the path would exceed `BA_NMAX` = 2e5 points; `snapP()`
quantizes the scan grid so a heatmap cell's axis label is the pressure that
was actually evaluated.

Compression work uses η at **1 bar** (not η(P_H)) — the convention that
reproduces the P_H = 2 worked example exactly. LPP gas is raffinate already at
P_H, so only FP charges pressurization compression.

## Conventions nailed down from the sources

- Basis: w = 1 kg adsorbent, V = w/(ρ_s(1−ε)), ε = 0.37 default.
  R = 8.314e-5 m³·bar/(mol·K) for inventories, 8.314 J/mol·K for van 't Hoff
  and work; k (adiabatic) = 1.4; WC divides by V(1−ε); energy uses
  M_CO₂ = 44 g/mol (as `BAAM.m`).
- **Raffinate leaves at y_δ** (the pre-adsorption composition) — from
  `BAAM.m` `feed_VSB`, matching e-BAAM Table 2 (T11).
- LPP pressurizes with gas at the (unknown) end-of-LPP composition;
  FP with feed. Recovery denominator: CO₂ fed in ADS (+ FP when used).
- e-BAAM η(P) = 0.8·19.55·P_atm/(1+19.55·P_atm), P in **atm**.
- States α (end ADS = start BLO), β (end BLO), γ (end EVAC), δ (end LPP) —
  2019 lettering; the e-BAAM paper uses a,b,c,… (a = end LPP).
- The 2019 classifier picked grid points with r ∈ (110.2, 110.3) and took
  min energy there; this app instead reports min energy among points with
  Pu ≥ 95 ∧ Re ≥ 90 plus the r_max/r_cut verdict.

## Validation ground truth (how it was made)

`Papers/Adsorption/SimplifiedModel/MgMOF-74_Check/BAAM.m` (Balashankar's
own code) run via `/Applications/MATLAB_R2026a.app/bin/matlab -batch`
(**sandbox must be disabled** for MATLAB — it needs system temp access).
At published tolerances (ode15s 1e-6) it reproduces paper Table 2 exactly;
anchors embedded in the tests come from a run with tolerances tightened to
1e-10/1e-12 (converged — the 1e-6 path values are off by up to ~3e-3 in
places, so DON'T re-anchor against a default-tolerance run). The Liske
worked example (Pu/Re quoted in her §Visualization) reproduces at
**T = 298.15 K**, not the 30 °C in the figure caption.

e-BAAM sources: `Papers/Adsorption/e-BAAM/AfterReviews/eBAAMAfterReviewClean.tex`
(authoritative Table 1 — the PDF text extraction shuffles columns!), energy
scaling 1.58 read from Fig9a.pdf text. Her MATLAB was not found locally
(GitHub only) — PEQ cycles (6-step/8-step) are NOT yet implemented; her
paper's 6/8-step example values (99.9/73.0 with P_eq = 0.23; 100.0/55.5 with
P_eq,2 = 0.12, P_eq,1 = 0.37) are the validation anchors to use when adding
them. PE steps need the donor buffer from the path (already stored) plus an
iterative P_eq match (<100 Pa in the paper).

## UI notes

- State `st` (localStorage `baam_v1`); debug handle `globalThis.__ba`
  ({st, PATH, KPI, GRID, MIX}).
- Recompute split: `onHeavy` (mixture/T/PH/yF/η/ΔP/scan ranges; 250 ms
  debounce) rebuilds path + KPIs + async grid (chunked, cancellation token);
  `onLight` (P_INT/P_LOW sliders) reuses the cached path. The path is always
  integrated down to min(P_LOW, scan-min, 0.02) so slider moves stay light.
- Plots: CO₂/N₂ competitive isotherms with constant-y family + transitions
  (BLO solid --blo, EVAC solid --evac, LPP dashed --lpp, ADS dotted --ads,
  states α β γ δ), y(P) path, design-space heatmap (viridis-ish `vir()`,
  2–98 % percentile color scale, white DOE-region outline, 1 bar line when
  P_H > 1, click-to-set, hover readout), Pareto with r_cut arc. Canvas `_map`
  fields carry the mapping for hit-testing.
- Composition curves are coloured by `ycol(y)` — a single-hue (rose, 330°)
  lightness+saturation ramp, deliberately away from the four step colours
  (blue/green/orange/violet); the feed curve is drawn bold and labelled
  "(feed)". End labels are de-collided by an 11 px minimum spacing pass
  (curves bunch up near saturation).
- **Crosshair** (replaced the old click-to-pin lines): `drawIsoPlot` stores
  `c._snap = ctx.getImageData(...)` and `c._fam`; `drawCross(which, P, q)`
  blits the snapshot back and draws the overlay, so a mousemove costs no
  isotherm evaluations (important in IAST mode). Loadings come from `qOn()`,
  which indexes the log-uniform 161-point family grid directly (interpolation
  error ≤ 7e-5 rel, checked in-browser). `hovA`/`hovB` persist the position so
  a redraw (slider move, resize) restores the crosshair. `putImageData`
  ignores the canvas transform — reset to identity before it, back to dpr
  after (see the function).
- Schematic: SVG built per cycle mode in `drawScheme` — a **step-list
  loop**, so adding steps (PEQ) means extending the `steps` array and the
  letters row.
- Import receiver in `applyImport()` (runs after restore, before UI build):
  accepts `{app:"baam", mode:"edsl", slot:"co2"|"n2", pvals, T?, note}` and
  IAST-style `{model, pvals, slot:"comp1"|"comp2"}`. Fitter side:
  `handOffBaam` in IsothermFitter.html (targets `baamA`/`baamB`) — c-basis
  van 't Hoff DSL passes b₀/U exactly; p-basis global fits converted with
  U_c = U_p − R·T̄ and b₀_c = b₀_p·R_B·T̄·e; single-T fits pin BAAM's T.

## Extension points (agreed with Arvind)

- **Liske cycles next**: 4-step FP is done; 6-step 1×PEQ and 8-step 2×PEQ
  remain (donor path from `blowPath` buffers + iterative receiver
  pressurization; validation anchors above). The schematic/step-list and
  the KPI plumbing were built with this in mind.
- More steps beyond the papers are anticipated ("we will introduce other
  steps") — keep the step-list architecture.
- Possible later: wireZoom-style box zoom on the isotherm plots (pattern in
  IsothermFitter), overlay of a second adsorbent's Pareto for comparison.
- **Parked by Arvind (2026-08-15)**: a "minimum energy vs P_H" panel with the
  per-step breakdown — the analogue of e-BAAM Fig. 10 (one grid scan per P_H
  value, ~10 scans, async like the existing map). It would give a second
  published figure to validate against. Offered, deferred, not commissioned.

## Figure exports & report (v2.1)

The shared module is `Apps/_shared/report_core.js`, pasted into `BAAM.html`
between the sentinels `/* == report-core v1 == */ … /* == /report-core == */`
as its own `<script id="ba-report-core">` between the engine block and the
UI block. Never edit the embedded copy — edit the canonical file, then:

```
sh _shared/sync_report_core.sh BAAM/BAAM.html
sh _shared/check_report_core.sh
```

**Read `_shared/README.md` first** — the module's contract, the
forced-light-palette export, the automatic equation numbering and the
integration checklist live there. BAAM-specific points:

- **`E` is taken** — it is the engine handle (`globalThis.BAAMENG`). The
  MathML builders are therefore named `ML`, not `E` as in DAC. Anything
  copied from another app's report block must be renamed.
- The old `dl()`/`csvPath()`/`csvIso()`/`csvMap()`/`csvPar()` functions are
  gone. Their bodies survive as `tblPath`/`tblIso`/`tblMap`/`tblPar`, which
  **return `{head, rows, notes}`** instead of writing a file; the core owns
  the download. They still export more than the plotted series (the path
  carries cumulative moles and vacuum work), which is why BAAM does not use
  the spec-derived CSV that DAC and MAPLE get for free.
- The `data-png` attribute handler was removed; the `.pfoot` rows now hold
  only their live readout spans (`ba-infoA`, `ba-infoB`, `ba-mapinfo`) and
  the two cards that had nothing else lost their `.pfoot` entirely.
- `Report.refresh()` is called at the end of `drawAll()` and in the
  `computeGridAsync` callback — the map and Pareto figures only become
  exportable once the async scan finishes.
- The report's classification block recomputes `r_max` and the feasible-point
  count from `GRID` directly. If the classifier constants
  (`BA_RCUT_EBAAM`, `BA_EN_SCALE_EBAAM`) ever change, the report follows
  automatically, but the *wording* of the verdict is duplicated from
  `renderClass()` — keep the two in step.

## Environment reminders

- jsc for headless JS tests; no node. MATLAB headless needs
  `dangerouslySkipSandbox` (writes to /tmp and its own caches).
- Preview server: launch.json "apps" → serve_apps.py port 8642 serving
  scratchpad `appserve/`; copy BAAM/ + IsothermFitter/ there to test the
  hand-off chain (relative `../BAAM/BAAM.html` needs sibling folders).
