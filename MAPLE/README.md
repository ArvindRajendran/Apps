# MAPLE — ANN PVSA Emulator

The MAPLE framework (Machine-assisted Adsorption Process Learning and
Emulation) as a self-contained offline web app: the trained neural
networks of [github.com/ArvindRajendran/MAPLE](https://github.com/ArvindRajendran/MAPLE)
ported to JavaScript. One file — open `MAPLE.html`.

Sources: Pai, Prasad & Rajendran, *Ind. Eng. Chem. Res.* **59**:16730
(2020) (framework) and *ACS Sustain. Chem. Eng.* **9**:3838 (2021)
(performance limits; the repo's 12-input networks are this generation).

## What it emulates

A 4-step PVSA cycle for CO₂/N₂ separation on a fixed rig (1 m column,
fixed pump curves): adsorption at P_H → blowdown to P_I → evacuation to
P_L (CO₂ product) → pressurization back to P_H, either with light
product (LPP) or feed (FP). Instead of integrating PDEs to cyclic
steady state (~minutes to hours), a dense feed-forward network
(12 → 30 → 30 → 30 → 1, tansig, Bayesian-regularized, trained on ~48,000
detailed-model CSS solutions per cycle) returns each KPI in
microseconds:

- CO₂ purity and recovery [%]
- specific energy [kWh/tonne CO₂, also shown in MJ/kg], plus the
  100 %-pump-efficiency variant as a small second line
- productivity [mol CO₂ / m³ adsorbent / s]

## Adsorbent description

Competitive single-site Langmuir:

$$q_i^* = \frac{q_\mathrm{sat}\, b_i C_i}{1 + b_\mathrm{CO_2} C_\mathrm{CO_2} + b_\mathrm{N_2} C_\mathrm{N_2}} \tag{1}$$

$$b_i = b_{0,i}\, e^{-\Delta U_i / (R T)} \tag{2}$$

Six editable properties: q_sat [mol/kg], b₀ for both components
[m³/mol], −ΔU for both [kJ/mol], particle density ρ [kg/m³]. The
network inputs derived from them are q_sat·ρ, log₁₀ b_i(298 K), ΔU_i, ρ
plus the six process conditions (y_F, t_ADS, P_H, log₁₀ P_I/P_H,
log₁₀ P_L/P_H, v_F), z-scored with the training moments.

Presets: Zeolite 13X, UTSA-16, IISERP-MOF2 (Limits paper SI Table S2)
and Mg-MOF-74 (MAPLE paper SI, converted from its mol/m³ basis at
ρ = 1130), plus the 36-material screening library of the MAPLE paper SI
(SSL refits of the Khurana & Farooq set, same basis). Custom SSL
parameters can be typed in or handed off from the Isotherm Fitter
(`#import=` payload; DSL second sites are ignored with a warning).

## Trained-range guards

A surrogate is only valid inside its training manifold. Every
evaluation is checked against the data-driven box of the ~96k labelled
training rows (e.g. q_sat·ρ 404–11,382 mol/m³, y_F 0.05–0.65, the
P_I/P_H and P_L/P_H ratio bounds, and the b and ΔU
ordering/selectivity constraints) **and** against the absolute pressure
windows the LHS actually sampled, P_I ∈ [0.07, 4] and P_L ∈ [0.01, 1]
bar. Both are needed: the ratio bounds are only an envelope over
P_H ∈ [1, 5], so ratios alone would admit joint data voids at low P_H
(P_L ≈ 2 mbar at P_H = 1 bar contains no training rows at all) where
the networks still return confident-looking numbers. Violations show a
red *extrapolating* badge and dim the KPIs. Raw network outputs outside
0–102 % are treated as invalid, as in the published wrapper; values
between 100 and 102 are shown clamped with an asterisk.

## Panels

1. **Cycle & sorbent** — LPP/FP, preset/library/custom, six property
   boxes with edit badges, 4-step schematic.
2. **KPIs** — the four indicators + US-DOE badge (purity ≥ 95 % &
   recovery ≥ 90 %).
3. **Isotherms at 30 °C** — your material against the four reference
   sorbents (toggle each; optional N₂ curves), crosshair readout.
4. **Parametric sweep** — any input 1-D (all four KPIs) or any pair as
   a full 2-D map with DOE-feasible cells dotted, hover readout.
5. **Optimizer** — NSGA-II run on the surrogate, in-browser. A
   four-sorbent comparison takes about 2 s, so no worker is needed.

   *Objectives* — purity–recovery, energy–recovery,
   productivity–recovery, or energy–productivity. The purity and
   recovery constraint levels are editable (defaults = US-DOE 95/90);
   they are enforced by penalty in every mode except purity–recovery
   (there the whole trade-off is the answer), and energy–productivity
   uses both. The objective functions actually being minimized are
   printed above the plot with the penalty embedded and the current
   constraint values substituted:

   $$J_1 = E/100 + P, \qquad J_2 = -Prod + P \tag{3}$$

   $$P = \sum_c \left(20\,v_c + 5000\,v_c^2\right), \qquad
     v_\mathrm{Pu} = \max\!\left(0,\ \tfrac{Pu_\mathrm{min}}{100} - \tfrac{Pu}{100}\right) \tag{4}$$

   The linear term in eq. (4) matters: a pure quadratic has almost no
   gradient at small violations, and the GA parks the front just
   outside the constraint where the results table then rejects it.

   *Decision variables* — t_ADS, P_H (unless VSA pins it at 1 bar),
   P_I, P_L and v_F. P_I and P_L are searched as **absolute** pressures,
   the way the training set sampled them, not as ratios of P_H. A
   collapsible range panel lets you narrow the search in physical units
   ("the pump cannot pull below 30 mbar", "t_ADS ≤ 60 s"); entries clamp
   to the trained support, so the search can be restricted but never
   pushed outside the data.

   *Algorithm* — real-coded NSGA-II with SBX crossover (η = 15),
   polynomial mutation (η = 20, rate 1/n_DV), crowded-comparison
   tournament and crowding-distance survival. Population, generations
   and seed are editable (defaults 120 × 70); the seed is deterministic,
   so identical settings reproduce an identical front.

   *Reading the result* — the front is redrawn live each generation,
   with a normalized hypervolume-vs-generation convergence plot after
   the run. Hovering a point gives its material, all four KPIs and all
   five decision variables; a collapsible panel plots each decision
   variable along the front. Optima sitting at a decision-variable bound
   are flagged ⚠, since bounds are the edge of the trained range. If a
   front's purity varies by less than the surrogate's own accuracy, a
   note says so rather than letting you read structure into noise.

## Exports and reports

Every figure carries **⬇ PNG** and **⬇ CSV** — the isotherms, both sweep
panels, the Pareto front, the convergence trace and each of the five
decision-variable trends. Because all of them are drawn through the
declarative plotter, their CSVs are derived from the very spec that was
plotted; the two-dimensional sweep map is the exception and exports its
grid in long format, one row per point, with a flag marking the cells
that meet the US-DOE targets. Figures are exported in the light palette
on an opaque white background whatever the screen theme.

Buttons stay disabled until the figure exists: the Pareto, convergence
and trend figures need an optimization run, and the decision-variable
trends additionally need their collapsible panel to have been opened,
since they are rendered lazily.

**Generate report** (top right) writes a self-contained HTML document:
introduction, objectives, the surrogate and optimizer formulation in
MathML with a nomenclature table, the numerical method including the
surrogate's accuracy against the detailed model, the adsorbent
description and operating point with any trained-range violations
flagged, the predicted indicators, the optimization outcome, and the
selected figures with numbered captions. Printable to PDF; no network
access. The report states plainly that its numbers are surrogate
predictions rather than detailed-model results.

## Validation (see `validation/README.md`)

The harness is 2,189 checks; all pass.

- JS forward pass vs the MATLAB network objects: ≤ 10⁻¹² (2,120 checks).
- R² vs the ~96k detailed-model labels shipped with the repo: 0.994–0.9997
  for all five outputs, matching MATLAB's values to 10⁻⁶.
- **External anchor:** the 216 published optimum points of the Limits
  paper SI (detailed-model KPIs at MAPLE-Opt conditions) are reproduced
  at the surrogate's own accuracy — purity median deviation 0.0 %,
  recovery 0.02 %, energy 1.8 %, productivity 0.9 %.
- Two findings about the public repo are documented there: the
  `MAPLE.m` wrapper's net1/net2 roles are correct but the LHS-samples
  spreadsheet has its Pu/Re column *headers* swapped; and the .mat files
  contain a sixth network that matches none of the shipped labels.

## Known limits

- Everything is rig-specific: 1 m column, fixed particle size, fixed
  pump curves, T_feed = 298 K. Only the 12 trained inputs vary.
- SSL equilibria only — materials needing DSL (strong dual-site
  character) are outside scope; the paper's own screening filter
  applied an R² ≥ 0.95 SSL-representability test.
- The optimizer optimizes the *surrogate*; near sparse corners of the
  training box (e.g. v_F at its lower bound) it can be a few per cent
  optimistic vs the detailed model. The published protocol re-ran every
  optimum through the detailed model — this app cannot, and says so.
  Treat a ⚠-flagged optimum as a candidate for a detailed-model run,
  not as an answer.
- Surrogate error is ~1–2 % typical (up to ~19 % on energy in the worst
  published-point case); do not read meaning into differences smaller
  than that. Purity is reproduced to ±0.15 pp (median) against the
  detailed-model labels.
- Pareto fronts can be genuinely disconnected — UTSA-16 has an isolated
  high-purity branch in purity–recovery, which a dense random sweep of
  the decision space confirms is a sparse region rather than a solver
  artifact. Gaps are not automatically a sign of a failed run.
