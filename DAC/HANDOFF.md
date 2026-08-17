# DAC app — handoff notes

Read this before editing `DAC.html`.

## File assembly

`DAC.html` is assembled by concatenation (BAAM pattern):

```
cat dac_head.html dac_engine.js dac_ui.html > DAC.html
```

- `dac_head.html` — CSS + all DOM + opens `<script id="dac-engine">`
- `dac_engine.js` — pure JS, no DOM (identical copy in `validation/`)
- `dac_ui.html` — closes the engine script, then the UI `<script>`

The UI builds its Web Worker from the *text content* of the
`dac-engine` script tag (Blob URL), so the engine block must stay pure
JS with no DOM references — it runs in both the page and the worker.

## Model/implementation decisions (agreed with Arvind)

- 4 steps: ads / evac / heat+vac / cool. Evacuation duration is an
  outcome (pump speed input, event at P = P_regen). Heating is a
  first-order wall approach to T_des (τ_heat input); the wall has no
  thermal resistance (couples via volumetric U_a); wall thermal mass is
  an input (default 0) added to the thermal duty.
- Repressurisation at the start of cooling is an instantaneous
  isothermal air addition (documented in README limits).
- Energy forms follow the Simplified DAC model: isentropic blower over
  ΔP; vacuum pump isothermal at T_amb by default with Glaser's
  adiabatic form as a toggle. Thermal and electrical reported
  separately (no COP).
- Product = heat-step outflow only; evacuation stream is vented; purity
  is dry-basis. Recovery counts all CO₂ entering over the cycle
  (ads + repressurisation + cooling inflow).
- Productivity on both bases: per m³ contactor and per m³ sorbent
  (= contactor / (1−ε_b)).
- Basis 1 m³ bed; all flows are per m³ of contactor.
- Custom sorbent (v1.7): "Custom" in the preset menu is a VIRTUAL
  preset. `CUSBASE` (UI) holds a deep copy of the source preset and is
  what `preset()` returns for Custom — it feeds the property boxes and
  their defaults, NEVER the isotherm evaluation. Evaluation always goes
  through `u.iso` (built by `readIso()` each collectU), which has the
  exact shape of a PRESETS entry; both `buildModel` (JS) and `Model`
  (Python) just swap `u.iso` in for the preset lookup. Keep that shape
  in lockstep with PRESETS. Coupling is EITHER `wet` (ψ/β modulation)
  OR `cmp` ({bw0}); `qstarC` gained a 5th arg `pw` [Pa] read only by
  the cmp branch — every call site must pass the local water partial
  pressure (rhs, initialState, the UI isotherm curves). cmp's b_w heat
  comes from the ΔH(H₂O) box (`pick(u.dHw, S.dHw||0)`), which is why
  that box is enabled in cmp mode even for inert water. Identity
  guarantee: Custom seeded from a preset must run bit-identical cycles
  (`cus.identity.*`, tolerance 0) — don't add custom-only behaviour to
  shared code paths.
- LDF coefficients are Arrhenius in bed temperature (v1.6):
  k(T) = k_ref·exp[−E_a/R (1/T − 1/293.15)] via `kLDF()` in the engine
  and `k_ldf()` in `dac_ref.py` — keep the two in lockstep. The k
  sliders are the 20 °C values (`TREF_LDF` is fixed, NOT tied to
  T_amb, so moving the ambient slider never silently rescales
  kinetics). `u.EaC`/`u.EaW` are in J/mol (UI converts from kJ/mol);
  absent or 0 must stay bit-identical to constant k — harness checks
  `ldf.*` enforce this, and preset switches reset both E_a sliders to 0
  because the published preset models are constant-k. E_a(CO₂) is a
  sweep axis (`SWAXES.EaC`, tr ×1000).

## Numerical structure (engine)

- Single state vector x = [N, yc, yw, qc, qw, T, Tw]; per-step
  constraint swap is algebraic, not a DAE: constant-P steps eliminate
  the unknown flow via dN/dt = −(N/T)dT/dt and solve the energy
  equation (linear in dT/dt where the inflow depends on it); the evac
  step advances N directly and keeps the V_C dP/dt expansion term
  (the (C_tot − RN) denominator).
- TR-BDF2 (γ = 2−√2) with Newton + FD Jacobian (7×7 Gaussian solve),
  divided-difference error estimate, adaptive h, bisection event
  location on P − P_regen.
- Flow/power quadratures are trapezoid sums on accepted steps — good to
  ~1e-5. The Python reference integrates them as ODE states (exact);
  that asymmetry is deliberate and documented in validation/README.
- CSS: repeat cycles until the start-of-cycle state repeats to 1e-5
  (component-relative with floors); sweep uses 1e-4/40 cycles for speed.

## Parameter provenance & assumptions

- Five paper sorbents: SI Tables S1 (Toth), S2 (GAB), S3 (ΔH_H2O),
  S4 (ψ, β), S5 (ρ, ε_b, c_p) of the five-sorbent DAC screening paper
  (Papers/Adsorption/DAC-1/AfterReview). b₀ given in Pa⁻¹ there.
- TMCM-41 & TPMS: water inert (agreed base case). The paper also has an
  "APDES-water" variant for them (Sec. on TMCM/TPMS water effects) —
  a possible future toggle.
- NbOFFIVE/SIFSIX coupling clamps (75 %/90 % RH) implemented as a cap
  on the q_w used inside the modified Toth (q_w,eff = min(q_w,
  q_GAB(x_max))), plus b_wet ≥ 0 and qs_wet denominator ≥ 0.05 guards.
- Lewatit: `params_dac_binary.m` / `isotherm_co2_binary.m` /
  `isotherm_h2o_gab.m` from FullModelLDFComparison, b₀ = 20112.4573
  bar⁻¹ converted to Pa⁻¹. **Assumed** (no measured value on hand):
  ΔH_H2O = 43.8 kJ/mol (as APDES), c_p,s = 1580 J/kg/K, k_CO2 default
  4e-3 s⁻¹ (the FullModelLDF kamine fit at 30 °C).
- Gas cp (37.1/33.6/29.1 J/mol/K), adsorbed-phase cp (50/75 J/mol/K —
  overridable via u.cpaC/u.cpaW, used by the consistency check),
  Antoine constants: engine header.
- Air: 400 ppm CO₂ of the dry fraction, P_amb = 101325 Pa.

## Layout (v1.1, same day)

Compact layout per Arvind: panel 1 is three columns — sorbent |
ambient+regeneration stacked | 2×2 cycle-schematic SVG (viewBox-scaled,
`.schewrap svg {width:100%}`); the separate schematic panel is gone.
Column rows carry `.ctrls.cols` (align-items: flex-start) so groups are
top-aligned; plain `.ctrls` stays flex-end for slider rows.
Profile plots are one row of three `.third` cards: T+P combined
(dual axes: T left °C, P right bar), loadings, gas composition.

## Layout (v1.3) & the simplified-model line

- Panel 1 is now sorbent | horizontal 1×4 schematic (`renderScheme`,
  viewBox 944×214, dashed "next cycle" return underneath). Ambient air
  and regeneration moved into the second panel, stacked under "Step
  durations & flows" in a flex column, which fills the space that column
  used to leave empty.
- "CO₂ per cycle" was dropped (redundant with productivity) and replaced
  by the working capacity Δq_CO₂ = q(end of ads) − q(end of heat), with
  Δq_H₂O and ΔT in the detail line.
- Every energy card carries a small `<span class="simp">` line with the
  Simplified DAC energy model's value, its kWh/t and the ratio to the
  dynamic number. `simplifiedEnergy(M, cyc)` in the engine implements the
  `dac.energyModel` forms but is fed **this cycle's own** Δq, ΔT and
  η_cap, so the gap is model-form only (void gas, swept-out sensible
  heat, CO₂ lost in evacuation). Do not "fix" it to use equilibrium Δq or
  a fixed η_cap = 0.30 — that was considered and rejected: it turns the
  line into a statement about the CSTR idealisation (factors of 2–5),
  not about the energy model. Wall mass is added to the simplified
  sensible term because it has no analogue there.
- The footer now carries the standard `../index.html` back-link that all
  the other apps have.

## Layout (v1.5) — dead-space pass

Four spots of empty space removed; keep the pattern when adding content:

- The two long explanatory paragraphs (`#dac-sorbnote`, the "Structure…"
  note) live at PANEL level, full-width below the `.ctrls.cols` row —
  not inside a column group, where they made one column ~150–290 px
  taller than the other and left a blank field beside it.
- The KPI panel bottom row is `.kpibot` (flex): the Working-capacity
  card (`#dac-kpiwc`, basis 320 px) sits beside `#dac-kpinote`, which
  fills the rest — previously the card sat alone on a grid row with
  ~1,050 px of empty cells. `renderKPIs` writes the card into
  `#dac-kpiwc`, and the dim class is toggled on BOTH `#dac-kpis` and
  `#dac-kpibot`.
- The sweep plot card `#dac-swwrap` starts `display:none` and is revealed
  at the top of `runSweep()` — before any draw, so `prepCanvas` never
  measures a hidden canvas (that would give width 0; see THE canvas
  bug below). Don't draw into it while hidden; the resize/theme
  redraw paths are safe because they check `sweepState.grid` first.

## Editable sorbent properties (v1.4)

`buildModel` takes optional `u.eb / u.ep / u.rhop / u.cps / u.dHc /
u.dHw`, each defaulting to the preset. Two rules that must not be
broken if you touch this:

1. **Never mutate `PRESETS`.** An edited ΔH builds a *copy* of the Toth
   object; a test pins `PRESETS.Lewatit.toth.dH0` after an override.
2. **ΔH(CO₂) is the Toth ΔH₀.** They are the same physical quantity, so
   an edit drives both unless `u.unlinkDH` is set. The UI exposes the
   unlink as a checkbox (enabled only once ΔH is edited) and badges it
   red, because unlinked runs are thermodynamically inconsistent by
   construction. Base case, 71 → 80 kJ/mol: linked gives Re 34.5→44.0 %
   and E_th 11.54→9.61; unlinked leaves Re alone and moves E_th to
   11.74. Both are pinned as anchors (`Lewatit-props`, the linked case
   with geometry edits too, and `Lewatit-unlinkDH`).

ΔH(H₂O) is free (GAB is on RH, no T-dependence) and the input is
disabled for the inert-water sorbents; `buildModel` also forces dHw = 0
for them regardless of what is passed. UI side: `PROPS` map →
`loadProps` / `markProps` / `propsU`, preset switch and Reset both call
`loadProps`. ε_b, ρ_p and −ΔH(CO₂) were added to the sweep axes.

The isotherm parameters themselves are deliberately NOT editable — that
is the Isotherm Fitter's `#import=` hand-off job, not a 13-field panel.

## Crosshairs (v1.2)

`plot()` caches `cv._geo` (margins, dpr, the ym mapper, the spec) plus
`cv._snap` (ImageData) and calls `wireCross(cv)` once. `drawCross`
inverts the x mapping, blits the snapshot, draws the vertical line,
interpolates every series at that x (`interpAt`, binary search — series
must be ascending in x, which both time and isotherm series are), dots
them, and draws a readout box. Series may carry `label`/`unit`; the
spec may carry `x.sym`/`x.unit` and `stepAt(x)` (profile plots use it to
name and colour the active step in the header).

The 2-D sweep heatmap keeps its own cell-snapping hover, so `drawSweep`
sets `cv._geo = null` in that branch — otherwise both handlers would
draw on the same canvas. The 1-D sweep branch goes through `plot()` and
gets the generic crosshair (its `sweepGeom` is nulled instead).

## Structure thermal mass (v1.2)

`u.mstr` [kg structure / kg sorbent] × `u.cpstr` [J/kg/K] is folded into
`M.cpStr` at build time and enters C_tot as `m_s·(c_p,s + cpStr)` — the
same term as `mstr*Cpstr` in `dac.energyModel`. Default 0. It changes
only the sensible duty (KPI E_thermal); mole balances and equilibria are
untouched. Anchor case `Lewatit-struct` (m_str = 2, c_p,str = 710)
pins it: E_th 11.55 → 15.97 MJ/kg.

## THE canvas bug (v1.0 sweep crash — sad-face tab)

`prepCanvas` originally did `h = +cv.getAttribute("height");
cv.height = h*dpr` — assigning `cv.height` **mutates the height
attribute**, so each redraw compounded the height by dpr
(260→520→1040→…). The sweep redrew fully on every mousemove, so
hovering grew the canvas exponentially until Chrome killed the GPU
context (frowny-face placeholder). Fixes now in place — keep all
three if you touch plotting:

1. design height cached once in `cv.dataset.baseH`;
2. backing store reallocated only when the target size changes;
3. sweep hover never re-renders — full render caches an ImageData
   snapshot (`sweepSnap`), hover blits it (transform reset to identity
   for putImageData) and draws only the crosshair; `sweepLeave`
   restores the snapshot and base readout.

The other apps do not use the `getAttribute("height")` pattern
(checked 2026-08-15) — this was DAC-only.

## Pitfalls discovered while building

- The default fan flow matters enormously: DAC beds are shallow, so
  specific flows are huge (v/L ≈ 10⁴–10⁵ m³/h per m³ bed). 500 m³/h
  starves the bed (recovery ≈ 96 % but productivity ≈ nothing);
  20,000 m³/h is the sensible default (recovery ≈ 35 % ≈ the
  Simplified model's η_cap = 0.30).
- scipy BDF cannot handle the instant-kinetics consistency case at
  rtol 1e-10 in reasonable time; the case runs on the JS engine (see
  validation README).
- Sweep colour scale needs robust clipping (95th percentile) + log
  mapping: zero-working-capacity corners (low T_des, high P_regen)
  produce 10³–10⁴ MJ/kg outliers that flatten a linear scale.
- jsc path: `/System/Library/Frameworks/JavaScriptCore.framework/Versions/Current/Helpers/jsc`.

## Validation provenance (know what the numbers rest on)

The 373 checks are JS ↔ our own Python re-implementation, plus the
instant-kinetics collapse onto our own MATLAB `dac.energyModel`. **No
result has been crosschecked against Glaser et al.'s published
numbers.** Two independent implementations of the same reading of a
paper agree even if the reading is wrong, so do not describe the app as
"validated against Glaser" — say "validated against an independent
reference implementation". Section D of validation/README states what
an external anchor would take.

## Figure exports & report (v1.8) — DAC is the pilot app

The shared module lives in `Apps/_shared/report_core.js`. DAC carries a
**byte-identical copy** pasted into `dac_ui.html` between the sentinels
`/* == report-core v1 == */ … /* == /report-core == */`. Never edit the
embedded copy:

```
sh _shared/sync_report_core.sh <path-to-dac_ui.html>   # re-paste
sh _shared/check_report_core.sh                        # verify all apps
```

Because DAC is assembled by `cat`, sync the **source part**
(`dac_ui.html`), then rebuild — syncing `DAC.html` would be overwritten
by the next build.

**Read `_shared/README.md` first** — the module's contract, the
forced-light-palette mechanism, the automatic equation numbering and the
integration checklist live there, not here. DAC-specific points:

- **The core must load before the UI IIFE.** It sits in its own
  `<script id="dac-report-core">` between the engine script and the UI
  script. It must NOT go inside `dac-engine`: the worker is built from
  that tag's text content, and the core touches the DOM.
- **CSV comes from the plot spec.** `plot()` ends with
  `Report.spec(cv.id, spec)`, so exports are derived from exactly what
  was drawn. The corollary is that a canvas drawn by anything *other*
  than `plot()` must clear its spec — `isoPlots` does this
  (`Report.spec("dac-isoW", null)`) in the inert-water branch, or
  switching Lewatit → TMCM would export the previous sorbent's GAB
  curve under the new sorbent's name.
- **The `redraw` hook is `renderAll(lastRes)` + `drawSweep()`** and must
  keep repainting *every* canvas, or figures captured for the report
  keep the screen palette.
- The report button is mounted by the core, top right beside the page
  title. It used to sit in the KPI panel header; do not put it back
  there — placement is uniform across the suite on purpose.
- The CO₂-isotherm card's `pinfo` was shortened to
  `Δq(CO₂) = … mol/kg (A→B)` so the header still fits on one line beside
  the export buttons. Lengthen it and that header wraps.

## Parked / possible next increments

- Exports + report are now also in BAAM (v2.1) and MAPLE (v1.5). PSA
  Simulator was descoped by Arvind mid-rollout and is untouched. Still
  to do: SMB Triangle and the four Equilibrium Chromatography apps —
  all bespoke plotters, so one `csv()` per figure, and the four
  chromatography apps can share most of their theory prose.
- Weather-driven year loop (the driver signature (ambient T, RH,
  initial state) → cycle result was kept for it).
- Steam-assisted desorption (s-TVSA), paper eqn 11 of Glaser.
- APDES-water variant for TMCM-41/TPMS as a toggle.
- Purity/recovery vs 1-D model calibration factor (like BAAM's ×1.58).
