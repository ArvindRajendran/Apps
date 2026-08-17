# Isotherm Fitter

`IsothermFitter.html` — fit pure-component adsorption equilibrium data to a
library of isotherm models and rank them honestly. The front door to the
rest of the suite: fitted parameters export in the exact format used by the
IAST, Uptake Curves, and chromatography apps (p in bar, q in mol/kg).

## Data input

Paste-first: whitespace/comma/tab-separated columns, header lines skipped
automatically. Two columns = (p, q) at a single temperature (entered in a
field); three columns = (T, p, q) — multiple isotherms at once. Units are
converted on entry — pressure: bar, kPa, Pa, MPa, atm, mmHg; loading:
mol/kg ≡ mmol/g, mmol/kg, cm³(STP)/g, mg/g and wt% (the last two ask for
the molar mass). The x-axis can also be a **concentration** (mol/m³ ≡
mmol/L, mol/L, g/L) for liquid-phase/chromatography data — concentration
values are fitted exactly as entered (no conversion), so fitted affinities
carry the same units as your data and a global-fit U is the
concentration-basis internal energy. A parsed-data table lets you exclude
outlier points individually. Two built-in examples (single-T Langmuir-like; 3-temperature
dual-site Langmuir with van 't Hoff affinities).

## Model library (16 + custom)

Linear (Henry), Langmuir, Dual-site Langmuir, Triple-site Langmuir,
Anti-Langmuir, Quadratic, Sips/Langmuir–Freundlich, Dual-site Sips, Toth,
Freundlich, BET/GAB (identical functional family once GAB's activity basis
is absorbed into pressure units), Type V (Langmuir + Hill step), UNILAN,
Temkin-type B·ln(1+Ap), Dubinin–Astakhov (T-explicit, with fitted
saturation pressure). Plus **custom expressions**: any q(p) formula — free
identifiers become fitted parameters, and T is available, so e.g.
`qs*b*exp(U/(8.314*T))*p/(1+b*exp(U/(8.314*T))*p)` fits a van 't Hoff
Langmuir directly.

## Fitting

Levenberg–Marquardt on logit-bounded parameters (positivity and physical
bounds enforced by construction; pole models bounded by the data range),
with heuristic initialization from the data (plateau, half-loading
pressure, Henry slope) plus randomized multi-start. **Residual weighting is
a choice** — relative Σ(Δq/q)² (default, emphasizes the Henry region),
absolute ΣΔq², or √q-weighted — and the ranking re-sorts live when
changed.

**Ranking is by AICc**, not SSE:

AICc = N ln(SSE_w/N) + 2k + 2k(k+1)/(N−k−1)

so extra parameters must earn their keep — dual-site Langmuir does *not*
beat Langmuir on Langmuir data. Reported per model: RMSE (mol/kg), R²,
AICc, ΔAICc, parameters ±95% CI from the Jacobian at the optimum, and an
over-parametrization warning whenever two parameters are >98% correlated
(the honest "your data cannot distinguish these sites" diagnostic).

## Multi-temperature modes (pull-down)

- **Per-temperature**: independent fits per T; aggregate AICc uses
  k_total = k·n_T so per-T fitting also pays its parameter bill.
- **Global van 't Hoff**: one parameter set across all temperatures with
  temperature-independent q_s and b(T) = b₀e^(U/RT) (ΔH_ads = −U) for
  every affinity parameter; Dubinin–Astakhov is fitted globally as-is
  (its temperature dependence is intrinsic). The selected global fit also
  produces the **model-based isosteric heat** q_st(q) via
  Clausius–Clapeyron, plotted vs loading.

## Display & export

Ranked table (tick to overlay, click for details), data + fit overlay
(points colored by temperature, log/lin p), residual plot (Δq/q %), parity
plot, and numeric-formula rendering of the selected fit (global fits shown
evaluated at the mean temperature). **All plots are interactive,
plotly-style**: drag a box to zoom, mouse-wheel to zoom about the cursor,
double-click to reset autoscale, and click legend entries (models or
temperature groups) to toggle their visibility without touching the
overlay selection. Every figure has **PNG and CSV download buttons**
(fit curves, parity pairs, residuals, isosteric-heat series).

**Isosteric-heat plot** (any data with ≥2 temperatures): the model-based
q_st(q) from a selected global fit (line) is overlaid with the
**model-free experimental q_st** (points) — ln p interpolated linearly in
q within each temperature's isotherm, slope of ln p vs 1/T by least
squares. Expect the points to scatter around the line and to blow up near
the plateau (where dq/dp → 0 makes the inversion ill-conditioned) — the
contrast between the noisy model-free route and the fitted-model route is
the point.

**One-click hand-off** (beyond the Copy-parameters block): buttons in the
details panel open the sibling app with the fitted parameters passed in
the URL fragment — no server, no clipboard, works offline. Pressure-basis
fits → IAST (component 1 or 2, verbatim) and Uptake Curves (converted to
its mole-fraction basis at that app's current total pressure, noted in
the confirmation). Concentration-basis fits → Equilibrium Chromatography
and Pulse (parameters verbatim — those apps are unit-agnostic, the
provenance note records your units), and, for Langmuir fits, the 2×2
Systems app (component A or B of the competitive Langmuir). Multi-T fits
hand off at a temperature you pick (global fits: b(T) evaluated there).
**BAAM** (CO₂ or N₂ slot): Langmuir/DSL fits land in the extended-DSL
parameters — global van 't Hoff fits keep their full temperature dependence
(concentration-basis b₀/U passed exactly; pressure-basis converted with
U_c = U_p − RT̄), single-T fits pin BAAM's temperature; other pressure-basis
suite models land as a BAAM IAST component.

## Validation

Reference: `validation/fitter_ref.py` (scipy.optimize.least_squares, trf,
native bounds) generating seeded synthetic datasets that are embedded
verbatim in the jsc harness (`validation/fitter_test.js`):

| Check | Result |
|---|---|
| Zero-noise recovery (Langmuir, Sips, Toth, DSL) | ≤ 1.3·10⁻¹⁵ rel (JS), machine (scipy) |
| Noisy Langmuir (seed 7): params vs scipy | 2.2·10⁻⁸ rel; SSE_w to 2.5·10⁻¹⁰ |
| AICc prefers Langmuir over DSL on Langmuir data | confirmed both engines (−138.3 vs −135.2) |
| DSL on the degenerate overfit valley | JS multi-start found a *deeper* minimum than scipy (8.39·10⁻³ vs 8.71·10⁻³) |
| Weighting changes the optimum (Toth, seed 11) | rel/abs params match scipy to ≤ 6·10⁻⁷ |
| Global van 't Hoff DSL, 3 T (seed 3): q_s, U recovery | within 4% of truth; SSE_w matches scipy to 1.4·10⁻⁹ |
| Isosteric heat from the global fit | 29.6 → 12.9 kJ/mol across loading, spanning U₁ = 30 → U₂ = 12 kJ/mol |
| Custom expression recovers Langmuir | 2.9·10⁻¹¹ |
| All 15 built-in models on 19 points | 26–60 ms total |

Notes: CI are asymptotic (±1.96σ from J_w at the optimum); AICc uses the
*weighted* SSE, i.e. models are compared under the same objective you
chose. Settings persist in localStorage. Works offline — single file.
