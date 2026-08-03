# Equilibrium Chromatography — Sorption Effect

`EquilibriumChromatographySorption.html` — equilibrium theory for a
single adsorbable component in an inert carrier when the ad/de-sorption
of the solute changes the fluid velocity (the "sorption effect").
Isothermal, isobaric, fixed inlet velocity; the Riemann problem
(uniform initial state y_I, constant feed y_F) is solved **analytically**.

Follows Ortner, Joss & Mazzotti, *J. Chromatogr. A* 1373 (2014) 131
(liquid formulation, one adsorbing + one inert compound) and Hefti,
Joss, Marx & Mazzotti, *IECR* 54 (2015) 12165 (gas formulation, §4).

## Model

With y the solute mole (or volume) fraction, ω = v/v_F, ν = (1−ε)/ε,
and dimensionless loading f(y) = q*/c_T:

1. [1 + ν f′(y)(1−y)] ∂y/∂τ + ω ∂y/∂ξ = 0
2. ν f′(y) ∂y/∂τ + ∂ω/∂ξ = 0,  ω(0, τ) = 1

Eq. 2 is the sorption effect: adsorption removes moles from the fluid
and decelerates it; desorption accelerates it. The dilute limit y → 0
recovers the classic constant-velocity app exactly.

## Analytical solution

In inert-flow coordinates (stream function of the carrier) the system
reduces **exactly** to a scalar conservation law in the mole ratio
Φ = y/(1−y) with flux νf. Therefore:

- **Wave/shock selection** = Oleinik envelope construction on νq vs Φ
  (loading → lower convex envelope, desorption → upper concave), drawn
  in the app mapped back onto the isotherm (shock chords are straight
  in Φ, curved in y).
- **Simple waves:** slope η = ω/(1+νf′(1−y)) with velocity from the
  quadrature d ln ω/dy = νf′/(1+νf′(1−y)) (exact trapezoid on a
  4000-point grid).
- **Shocks:** solute + inert jump balances give both the slope
  η_s = ω_U Δy/[Δy + νΔf(1−y_D)] and the velocity jump
  ω_D = ω_U (Δy + νΔf(1−y_U))/(Δy + νΔf(1−y_D)) in closed form.

The velocity plateau ahead of the leading front (ω_I, the outlet
velocity before breakthrough) falls out of the construction and is
reported in the summary line.

Isotherms: Langmuir, anti-Langmuir, linear, Sips, quadratic, Type V,
BET, dual-site Langmuir, the Do & Do water/activated-carbon model
(both branches, with Hefti et al. 2015 Table 2 defaults — the
desorption branch implements the hysteresis switching via the
max(n_A, n_D) composite), and custom f(y). All evaluated with
y ∈ [0, 0.995].

## Display

- Twin colour-mapped physical planes: composition y (viridis) and
  velocity v/v_F (diverging, white at ω = 1 — blue decelerated, red
  accelerated), same fan/shock skeleton, **linked crosshair** (hover or
  pin in either plane).
- Column profile and breakthrough cards show y (solid, left axis) and
  v/v_F (dashed, right axis) together at the crosshair location.
- RH ↔ y helper for water-vapour cases (Buck psat, y = RH·psat/100P).

## Validation

The construction was validated against an independent Python reference
before porting, and the shipped JS engine reproduces all six reference
cases to 6 decimals (ν = 1.7027, i.e. ε = 0.37):

| Case | Setup | Result | Checked against |
|---|---|---|---|
| A | Langmuir (qs 5, b 10), 0 → 0.3 | single shock η = 0.044876, ω_I = 0.713463 | hand algebra (exact) |
| B | same, desorption 0.3 → 0 | wave η 0.011610→0.292147, ω_I = 1.380297 | independent PDE solver, rms 6×10⁻⁴ |
| C | dilute 0 → 0.001 | η = 0.011724 | classic constant-velocity formula (exact) |
| D | Type V (1, 5, 8, 2.5, 6), 0 → 0.8 | single shock η = 0.051367, ω_I = 0.241094 | Godunov solver in Φ-coordinates (provably entropic) |
| E | anti-Langmuir (1, 0.9), 0 → 0.6 | **shock** (classic theory: all-wave!), ω_I = 0.538524 | Godunov |
| F | Type V desorption 0.8 → 0 | wave (0.8→0.4484) + shock (0.4484→0), ω_I = 3.357746 | Godunov (wave positions + shock speed) |

Case E shows genuine sorption-effect physics: with y of order 1 the
velocity coupling can sharpen a front that classic dilute theory says
should spread. During development, the hull-variable question (y vs Φ)
was settled by case D: hulling in y gives a plausible-looking but
non-entropic solution; the Φ-hull matches the provably-convergent
Godunov reference.

## Reproducing Hefti et al. (2015) Figures 3–4

Water vapor + He on activated carbon at 45 °C, Table 2 isotherm
(Do-&-Do type: BET term + Hill term; desorption branch via the
max(n_A, n_D) composite, which reproduces the hysteresis switching at
the closure points including the derivative kink). Mapping y → RH uses
x = y·P/psat with psat(45 °C) = 9.59 kPa (Buck), P = 1 atm, so
x = 10.331·y; loading scale assumes ρ_b = 450 kg/m³, ε_t = 0.68
(Table S3 of their SI was not available — these affect only absolute
times, not the structure, H states, velocities, or time ratios).

Both branches are built into the isotherm menu ("Do & Do water/AC …",
adsorption and desorption-with-hysteresis) with the Table 2 values as
defaults, so the reproduction is: set ε = 0.68, pick the branch, and

- **Adsorption run**: ads. branch, y_I = 0, y_F = 0.092
- **Desorption run**: des. branch, y_I = 0.092, y_F = 0

(y_sat = psat/P and S = ρ_b/((1−ε)c) are editable model parameters,
so different T, P, or bed properties are one edit away.)

Comparison against their equilibrium-theory solution:

| Quantity | Hefti Fig 3/4 | This app |
|---|---|---|
| Ads. structure | wave + trailing shock (semishock) | same |
| Ads. H state | ≈ 5.0 % | 5.4 % |
| Ads. velocity plateaus | 0.91 → ~0.955 → 1.00 | 0.908 → 0.960 → 1.00 |
| Des. structure | shock 2 – wave – shock 1 | same |
| Des. H₁ | 6.2 % | 6.5 % |
| Des. H₂ | 7.0 % | 8.9 % (sensitive — see below) |
| Des. velocity plateaus | 1.10 → ~1.07 → ~1.065 → 1.00 | 1.101 → 1.097 → 1.069 → 1.00 |

The H₂ tangency sits on the nearly-flat top of the desorption branch
(Δn ≈ 0.3 mol/kg over Δy ≈ 2 %), so its location shifts by percent-
points for ~0.2 mol/kg differences in isotherm shape; their inset shows
n(I) ≈ 22.5 mol/kg where Table 2 parameters give 21.9. All velocity
plateaus, both front structures, and the other intermediate states
reproduce well.

The initial and feed states can be set by typing values or by dragging the state markers directly along the isotherm curve — the full solution updates live during the drag.

Settings persist in localStorage. Works offline — single file, share freely.
