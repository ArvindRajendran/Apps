# Equilibrium Theory of Chromatography — Pulse Injection

Single-file interactive solver for a **finite pulse** (injection of c_F for
a duration τ_inj, then purge back to c_I) under single-component
equilibrium theory: `EquilibriumChromatographyPulse.html` (~57 kB).
The initial and feed states can be set by typing values or by dragging the state markers directly along the isotherm curve — the full solution updates live during the drag.

Works offline; share by sending the file. Forked from
`EquilibriumChromatography/` (the single-step app), which stays as the
validated reference for pure Riemann problems.

## Method: wavefront tracking (Holden–Risebro)

The pulse launches two Riemann fans — loading at (ξ,τ) = (0,0) and purge
at (0, τ_inj) — that inevitably interact. The solver:

1. Resolves each Riemann problem by envelope construction (Oleinik):
   lower convex envelope for loading, upper concave for purge.
2. Discretizes rarefaction fans into ~140 small jumps, so **every wave
   is a straight-line front** separating two constant states, moving at
   the chord (Rankine–Hugoniot) speed
   σ = 1/(1 + ν·Δq/Δc),  ν = (1−ε)/ε.
3. Marches through **collision events**: when two fronts meet, the
   middle state vanishes and the local Riemann problem (left state,
   right state) is re-resolved by the same envelope construction; the
   new fronts replace the old ones.

This scheme converges to the entropy solution for any isotherm, is
**exactly conservative**, and is *exact* (no discretization error) for
piecewise-linear isotherms — so the "use data as-is" table mode works
natively. Curved decaying shocks emerge as polylines through the event
cascade. Interactions beyond ξ ≈ 1.06 are ignored (they cannot affect
the column), and an event cap of 30 000 guards against pathological
inputs.

## Features

Everything from the single-step app (isotherm library + live formula,
custom f(c) expressions, data table with as-is/smoothed interpolation,
zoomable isotherm plot, dimensionless ↔ dimensional axes, RH↔c helper,
pinnable crosshair, PNG export) plus:

- **Injection duration input** that follows the axis toggle
  (τ_inj dimensionless, or t_inj in seconds via t₀ = L/v_i).
- **Physical plane** shows both fans, all tracked fronts (waves thin,
  shocks bold with width ∝ strength) and the curved decaying shock
  paths over the concentration colour map.
- **Isotherm plot** overlays both envelopes: loading solid, purge dashed.
- **Peak metrics at ξ = 1**: elution window, peak concentration and its
  retention time, width at half height.
- **Mass-balance check**: eluted area / injected area (should read 1).
- **Interaction event table**: every collision with τ, ξ, the incident
  front types, and the local states.

## Validation

Against the closed-form Rhee–Aris–Amundson results for Langmuir
(q = 5·10·c/(1+10·c), ε = 0.37, c_I = 0 → c_F = 1, τ_inj = 1):

- First wave–shock collision: analytic (τ, ξ) = (1.2421, 0.1421);
  app (1.2429, 0.1422).
- Decayed-shock exit: integrating dξ/dτ = σ_sh along the rear fan gives
  ξ_s = K·(u/(u−1))², u = 1 + b·c_s → arrival τ = 38.77, peak
  c_s = 0.0522; app: τ = 38.7568, peak 0.050 (within the Δc = 1/140
  fan discretization).
- Mass balance = 1.000 in all tested cases (Langmuir, Type V composite
  interactions, negative/desorption pulse, tabular data as-is).
- τ_inj → large reproduces the single-step solution exactly
  (0 interactions, breakthrough at τ = 8.7396).

## Increasing accuracy

Fan resolution is set by `dc = (cmax−cmin)/140` in `buildPulse()`.
Raising 140 sharpens the tail quantization (peak height resolution)
at roughly linear cost in events.
