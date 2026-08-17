# MAPLE app — handoff notes

Read this before editing `MAPLE.html`.

## File assembly

```
{ cat maple_head.html
  echo '<script id="maple-weights">'; cat maple_weights.js; echo '</script>'
  echo '<script id="maple-engine">';  cat maple_engine.js;  echo '</script>'
  cat maple_ui.html; } > MAPLE.html
```

- `maple_head.html` — CSS + all DOM (source parts live in the build
  scratchpad; the deployed app is the single file)
- `maple_weights.js` — `var MAPLE_WEIGHTS = {...}` (from
  `validation/export_maple_weights.m`; snapshot in `validation/`)
- `maple_engine.js` — pure JS, no DOM (snapshot in `validation/`,
  **keep in sync**)
- `maple_ui.html` — one `<script>` with the UI, closes `</body></html>`

No Web Worker: a full network evaluation is ~5 µs, so sweeps, maps and
even the NSGA-II comparisons run on the main thread (chunked per
sorbent with setTimeout).

After any engine change, re-run the harness in `validation/` (2,189
checks) and re-sync the `maple_engine.js` / `maple_test.js` snapshots
there — the harness runs against the snapshot, not the deployed HTML,
so a drifted snapshot means you are testing the wrong code.

## Rules that must not be broken

1. **Net roles: net1 = recovery, net2 = purity** (then log₁₀En, Prod,
   log₁₀En at η=100%). The repo's LHS-samples xlsx has its Pu/Re column
   HEADERS swapped, so any "obvious" re-derivation against that file
   will tell you to swap them back — it is wrong. The proof of the
   correct mapping is the pinned ≥95 % purity constraint at the
   Limits-paper SI optimum points (validation README, Findings).
2. **The input convention** (order, q_sat·ρ, log₁₀ b(298 K) — T=298.0
   not 298.15, log₁₀ P-*ratios*, z-score) is pinned by 2,120 anchor
   checks at 10⁻¹². Any change must fail the harness loudly. Note the
   networks take ratios while the optimizer's decision variables are
   absolute pressures (rule 3 and the Optimizer section) — these are
   deliberately different spaces; `optProblem.toU` is the only bridge.
3. **Guards are data-driven** — the BOX in the engine is the min/max of
   the transformed inputs over the ~96k labelled rows (NOT the rounded
   ranges in the repo README, which are wrong for y_F (0.65 not 0.5)
   and q_sat·ρ). **And the BOX alone is NOT the trained support**: the
   LHS sampled ABSOLUTE pressures (P_I ∈ [0.07, 4], P_L ∈ [0.01, 1]
   bar), so `PABS` must be enforced alongside the ratio bounds
   everywhere (`guards()`, optimizer validity, any new DV space).
   Ratio-only bounds admit joint data voids at low P_H (P_L ~ 2 mbar at
   P_H = 1 — zero training rows) and the GA will park entire fronts in
   them, producing convincing-looking artifacts (flat Pu ≈ 99.5–100
   VSA fronts; FP "beating" LPP on energy). See validation README,
   "absolute-pressure support lesson".
4. **Nonphysical-output window:** raw Pu/Re outside (0, 102) →
   invalid (repo wrapper convention); above 100 → clamp for display
   (asterisk note) and for objectives. Without the clamp the Pu–Re
   optimizer rides a Pu ≈ 101.8 artifact pocket at P_H ≈ 4 bar and the
   front collapses to a fake (100, 100) point.

## Optimizer

**Decision variables are absolute pressures** (since v1.4):
x = [t_ADS, P_H, log₁₀ P_I, log₁₀ P_L, v_F], matching how the LHS was
sampled. The networks' ratio inputs are derived per candidate in
`toU()` and rejected in `fobj` when they leave BOX. Do not "simplify"
this back to ratio DVs: ratio-space search was half of why high-
selectivity fronts came out gappy (validation README, "Fragmented
Pu–Re fronts").

**Selection is the crowded-comparison operator** (rank, then crowding
distance). `S.cd` must be recomputed alongside `S.nd` every generation
— a stale `cd` silently degrades diversity and was the other half of
the gappy-front problem. Rank-only selection degenerates to a random
walk once the population is all rank 0.

**`dvRanges(dvr, fixPH)`** clamps user ranges into `DV_FULL` and can
only NARROW. Inverted input collapses to a point rather than an empty
interval. The UI panel writes the clamped values back into the boxes,
so a rejected widening is visible immediately.

**The penalty must keep its linear term**: P = 20 v + 5000 v² per
active constraint. A pure quadratic has ~zero gradient at small
violations, so the GA parks the front a hair OUTSIDE the constraint
(Pu ≈ 94.9) and the table then rejects everything as infeasible — a
real v1.1 bug, worst in VSA EnPr, now guarded by harness check
`nsga.enpr-vsa-feasible`. The `mp-objbox` display must be kept in sync
with `optProblem`: it claims to show the exact objective functions.

**Modes:** PuRe (unconstrained), EnRe/PrRe (Pu_min penalty), EnPr
(Pu_min + Re_min). GA pop/gens/seed are user-editable state
(gaPop/gaGen/gaSeed); the operator η's stay fixed in the engine.

**Stepwise API:** the UI drives `nsga2init` / `nsga2step` /
`nsga2front`, two generations per setTimeout tick, so the front
animates live and per-generation fronts feed the hypervolume plot
(`hv2`, normalized to the final front). `nsga2(opts)` is a thin wrapper
over the same functions and is harness-checked to be *identical* to the
stepwise path — don't let them drift apart.

**Known limitation:** penalties are added to both objectives, so
infeasible points form a domination chain and front size drops when the
feasible set is marginal (VSA 13X in EnPr) at small GA budgets. At the
app defaults (120×70) fronts are 32–120 points. The textbook fix is
Deb's constrained-domination instead of penalties — but that would
change the objective display Arvind explicitly asked for, so it needs
his go-ahead.

**Front annotations:** `decodeFront` flags DVs within 0.5 % of a bound
(skipping the VSA-pinned P_H, whose lo == hi) plus the absolute
P_I/P_L window edges; the table shows ⚠ with a footnote. `resLimited`
compares a front's purity span across Re 80–95 % against
`PU_RESOLUTION` (0.15 pp, the surrogate's median purity error) and
prints the "resolution-limited" note — rarely triggered since v1.4,
which is the intended outcome.

## Layout / UI notes

- Plot/crosshair machinery is the DAC pattern (baseH cache, snapshot
  blit; see DAC/HANDOFF for the canvas-height bug it prevents). The 2-D
  map reuses the same canvas as the 1-D sweep by legitimately changing
  `dataset.baseH` through `setBaseH()`.
- Isotherm panel is pure-component SSL at 303.15 K, 4 reference
  sorbents toggleable; current material always in accent colour.
- The Pareto panel re-evaluates every front point through `evalKPI`
  before tabulating, and filters constraint feasibility (editable
  Pu_min/Re_min, defaults 95/90) on clamped values.
- Pareto hover: series carry `meta` (decoded points incl. `u` and the
  `pinned` bound list); `drawCross` snaps to the nearest scatter point
  within 20 px and calls `drawPointTip`, otherwise it draws both
  crosshair lines with an x/y readout.
- Both collapsible panels (`mp-dvwrap` DV trends, `mp-dvrwrap` DV
  ranges) are `<details>`. Canvases cannot size themselves while
  `display:none`, so never draw into a closed panel —
  `optState.dvDirty` defers the DV-trend render to the toggle event.
- The process schematic routes the LPP light-product line *above* the
  step boxes (branch dot on the N₂-rich stream → dashed line at y=22 →
  down into the LPP box top). Don't route it through the box band —
  that was a v1.0 visual bug.
- Isotherm Fitter hand-off: accepts the BAAM-style `#import=` payload
  (per-slot co2/n2; `dUb` arrives in J/mol, negative); DSL second sites
  are dropped with a note.

## Parked / possible next increments

- Downloadable HTML report (suite-wide feature, parked by Arvind).
- Nonlinearity-plot mode (Limits paper case 5: an optimization per
  (H_CO2, H_N2) pixel — feasible at ~0.2 s/pixel, needs a worker +
  progress UI).
- Material-DV optimization ("what is the best possible adsorbent",
  curves E3/P3) — same NSGA-II with 10 DVs; overlay the paper's
  power-law limit curves for reference.
- Adding MAPLE as an explicit target button in the Isotherm Fitter's
  hand-off row (currently the payload works but the Fitter's UI only
  offers BAAM).

## Figure exports & report (v1.5)

Shared module `Apps/_shared/report_core.js`, pasted into `MAPLE.html`
between the sentinels as `<script id="mp-report-core">` between the engine
block and the UI block. Re-paste with
`sh _shared/sync_report_core.sh MAPLE/MAPLE.html`, verify with
`sh _shared/check_report_core.sh`. Never edit the embedded copy.

**Read `_shared/README.md` first** — the module's contract, the
forced-light-palette export, the automatic equation numbering and the
integration checklist live there. MAPLE-specific points:

- **Eight canvases live in one `.pcard`** (Pareto, convergence, five DV
  trends) and two more share the sweep card. The core detects this: a card
  holding one registered figure gets its buttons in the card header, a card
  holding several gets a right-aligned `.figrow` under each canvas. This
  generalization was added for MAPLE — do not assume one header per card.
- `plot()` ends with `Report.spec(cv.id, spec)`, so every CSV here is
  spec-derived. `renderSweep2()` must therefore call
  `Report.spec("mp-sw-a", null)` before handing the canvas to `drawMap`,
  or the 2-D map would export the previous 1-D sweep's series. It stashes
  the grid in `lastMap` for the map CSV instead.
- A figure's `csv()` may return `null` to defer to the spec — that is how
  `mp-sw-a` serves a grid in 2-D mode and the plotted series in 1-D mode.
  The core falls back automatically.
- **DV-trend `ready()` checks `$("mp-dvwrap").open`.** Those canvases are
  rendered lazily on the `<details>` toggle, so before the panel is opened
  they have never been drawn; without the check the report would embed five
  blank frames. The toggle handler calls `Report.refresh()` so the buttons
  enable the moment the panel opens.
- The report clamps purity and recovery to 100 % exactly as the KPI cards
  do — the surrogate can return ~101.8 %, and an unclamped value in a
  report would read as a result rather than as the artifact it is.
