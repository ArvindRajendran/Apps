# Uptake Curves

`UptakeCurves.html` — predicts a gravimetric uptake experiment: step change
of sorbate mole fraction y₀ → y₁ in an inert carrier at constant total
pressure and temperature, single component, spherical particle and crystals.

## Models (pull-down)

- **Micropore only** — intracrystalline diffusion
  ∂q/∂t = r⁻²∂/∂r(r²D_μ(q)∂q/∂r) with equilibrium loading at the crystal
  surface and symmetry at the centre.
- **Macropore only** — gas diffusion in the particle voids with local
  crystal equilibrium: ε_p∂c/∂t + ρ_p∂q*/∂t = R⁻²∂/∂R(R²D_e∂c/∂R); the
  effective diffusivity of the equivalent nonlinear diffusion problem is
  D_e/[ε_p + ρ_p dq*/dc].
- **Macro + micro (bidisperse)** — the Ruckenstein-type two-scale model: a
  crystal diffusion problem lives at every particle radial location, its
  surface loading set by the local macropore composition; its mean-loading
  rate feeds back into the macropore balance.
- **LDF** — dq̄/dt = k(q* − q̄), k entered directly or via Glueckauf
  k = 15 D_μ0/r_c².

Macropore transport is assembled from parts the way you would in practice:
molecular D_m and Knudsen D_K combine by Bosanquet
(1/D_pore = 1/D_m + 1/D_K), then D_e = ε_p D_pore/τ_p; a helper estimates
D_K = (2r_pore/3)√(8RT/πM). The app displays D_pore, D_e and both diffusion
time constants (τ_micro = r_c²/D_μ0, τ_macro = R_p²(ε_p+K)/D_e with K the
dimensionless chord slope) plus which resistance controls.

**Darken correction** (checkbox): D_μ = D_μ0·Γ(q) with
Γ = d ln p/d ln q = q/(y f′(y)), evaluated from the *forward* isotherm on a
fine tabulated y↔q map — no symbolic inversion, so it works for every model
in the menu including custom expressions. For Langmuir this reduces to the
exact 1/(1−θ). If the isotherm is non-monotonic over the step interval the
correction is disabled with a warning.

Isotherm menu (same set as the chromatography apps, in mol/kg): Linear,
Langmuir, Anti-Langmuir, Sips/Hill, Quadratic, Type V (Langmuir + Hill
step), BET, Dual-site Langmuir, Do & Do water/AC (ads. and des. branches,
Hefti 2015 parameters), custom f(y). The selected isotherm is displayed as
a rendered mathematical expression with the current parameter values
substituted in, updating live as you type. States are draggable on the
isotherm panel; the Γ(q) curve is overlaid (dashed) when Darken is active.

## Outputs

- **Intracrystalline profiles** q(r) — for the bidisperse model, crystals at
  the particle centre / mid / surface shells are shown together.
- **Intraparticle profiles** — gas fraction y(R) and crystal-average loading
  q̄(R).
- **Uptake curve** plotted as (m_t − m₀)/(m_∞ − m₀) — the y-axis label
  spells out the normalization, so desorption runs (where both numerator
  and denominator are negative) unambiguously read as the fraction of the
  total mass change released, rising 0 → 1 like the standard gravimetric
  convention. m includes the macropore gas inventory as well as the
  adsorbed phase (what a balance sees). **t ↔ √t axis toggle**, hover
  crosshair, time slider + animation, ghost profiles at earlier times, PNG
  export. The summary reports t₅₀, t₉₀ and the apparent D/R² from the
  initial √t slope (π s²/36) — exactly the short-time analysis used on
  gravimetric data.
- **CSV export** on every result: the uptake curve (t, √t, uptake
  fraction) and both profiles at the slider time (crystal q(r) — three
  shell columns in bidisperse mode; particle y(R) and q̄(R)), each with a
  metadata header recording the model, isotherm, conditions and transport
  parameters, so the file is self-describing.

## Numerics

Conservative finite volumes on the normalized radius at both scales
(cell-centred, half-cell Dirichlet closure at the surface), method of
lines, L-stable TR-BDF2 with quasi-Newton (colored finite-difference
Jacobian, reuse + rebuild policy, sub-stepping fallback), √t-spaced time
grid. Runs in a Web Worker; a typical bidisperse solve (20×20 cells, 280
steps) takes well under a second. t_end is estimated automatically from
both time constants (chord-averaged Darken factor included) and
auto-extends when the curve has not plateaued; a manual override is in the
numerics panel, alongside the grid controls (crystal cells N_c, particle
cells N_p, time steps K — honored by every model, and the active grid is
echoed in the summary line so convergence can be checked by refining).

During development an early Crank–Nicolson stepper was rejected: it is
A-stable but not L-stable, and when Δt ≫ τ_macro its undamped ringing of
the stiff macropore mode biased bidisperse uptake curves by ~10⁻¹ in the
micropore-controlled limit. TR-BDF2 eliminated this.

## Validation

The engine was validated against an independent Python reference and exact
solutions; the shipped JS reproduces it. The reference scripts and a jsc
test harness ship in [`validation/`](validation/) with run instructions;
implementation notes for future development are in
[`HANDOFF.md`](HANDOFF.md).

| Check | Reference | App engine |
|---|---|---|
| Micropore, constant D vs exact series 1−(6/π²)Σn⁻²exp(−n²π²τ) | U(τ=0.01…0.4) = 0.308514/0.606940/0.770479/0.915496/0.988269 | max err 3.6·10⁻⁵ |
| Macropore, linear isotherm vs same series with D_eff = D_e/(R_p²β), β = ε_p+ρ_pH/c_T | β = 74.719, D_eff = 0.38664 s⁻¹ | max err 1.1·10⁻⁵ |
| LDF vs 1−e^(−kt) | U(1/k) = 1−1/e = 0.63212 | 0.63211 in-app |
| Bidisperse → macropore-controlled limit (D_μ → ∞) | series | 2.1·10⁻³ (grid-limited) |
| Bidisperse → micropore-controlled limit (D_e → ∞) | series | 1.6·10⁻³ (grid-limited) |
| Γ table vs exact Langmuir 1/(1−θ) | — | 1.2·10⁻⁶ |
| Darken Langmuir (q_s 5, b 10, y 0↔0.5, D_μ0/r_c² = 10⁻³ s⁻¹): ads t₅₀/t₉₀ | 11.24 / 59.69 s | 11.241 / 59.689 s |
| … desorption t₅₀/t₉₀ (slower, favorable isotherm) | 19.97 / 140.85 s | 19.965 / 140.854 s |
| Bidisperse linear anchor U(5,20,60,150,400 s) | 0.599232 / 0.912529 / 0.998233 / 1 / 1 | max diff 6·10⁻⁴ |
| Bidisperse Langmuir+Darken anchor U(5,20,60,150 s) | 0.841206 / 0.999919 / 1 / 1 | max diff 1.3·10⁻⁵ |

The desorption-slower-than-adsorption asymmetry under Darken (t₅₀ ratio
1.78, t₉₀ ratio 2.36 for the case above) is the classic signature of a
favorable isotherm with thermodynamically corrected diffusion.

Physics notes: equilibrium at the particle surface (no film resistance);
isothermal (no heat effects — real gravimetric data often show thermal
tails); constant total pressure with the sorbate in an inert carrier.

Settings persist in localStorage. Works offline — single file, share freely.
