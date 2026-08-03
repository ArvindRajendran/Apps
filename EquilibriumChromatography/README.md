# Equilibrium Theory of Chromatography

Single-file interactive solver for the single-component equilibrium-theory
(Rhee–Aris–Amundson) Riemann problem: `EquilibriumChromatography.html`
(~40 kB, all closed-form — no embedded data). Works offline; share by
sending the file.

## Theory

Column mass balance ∂c/∂τ + σ(c)·∂c/∂ξ = 0 with

1. σ(c) = 1/(1 + ν·f′(c)),  ν = (1−ε)/ε

Riemann solution by envelope construction (Oleinik condition) on
[min(cI,cF), max(cI,cF)]: **loading** (cF > cI) uses the lower *convex*
envelope of q = f(c), **desorption** the upper *concave* envelope.
Envelope-on-curve arcs are simple waves; chords are shocks with

2. σ_sh = 1/(1 + ν·Δq/Δc)

Implementation: f sampled on a 4000-point grid, Andrew monotone-chain
hull → consecutive hull indices = wave arcs, index jumps = shocks;
f′(c) by central finite difference (h = max(1e-8, 1e-6·c)) — no
derivative input needed. Tiny σ inversions at wave↔shock tangent
points (FD vs chord slope) are clamped; genuine ordering violations
raise an error.

## Features

- **Inputs**: L, interstitial velocity v_i, voidage ε, initial state cI,
  feed state cF (mol/L), swap button (adsorption ↔ desorption). Axis
  toggle: dimensionless (ξ, τ) ↔ dimensional (z in m, t in s via
  t₀ = L/v_i).
- **Four linked plots**:
  1. *Physical plane* (ξ vs τ) — concentration colour map (viridis),
     characteristics thin white, shock paths bold red, region labels.
     Hover for a crosshair; **click to pin**.
  2. *Isotherm & envelope* — f(c) muted, wave arcs thick blue on the
     curve, shock chords red, cI/cF markers.
  3. *Column profile* c(ξ) at the crosshair time.
  4. *Breakthrough curve* c(τ) at the crosshair position.
  Profile and breakthrough export as PNG (self-captioned with the
  τ / ξ, or t / z, they were taken at).
- **Isotherms**: Langmuir, anti-Langmuir, linear, Sips/Hill, quadratic,
  Type V (Langmuir + Hill step — water-on-MOF shape), BET, dual-site
  Langmuir, each with editable parameters; plus a **custom f(c)**
  expression box (safe parser: numbers, `c`, `+ - * / ^`, parentheses,
  exp/log/log10/sqrt/abs/tanh/…, `pow/min/max`) — no raw JS eval.
- **Live formula display**: the selected isotherm is rendered in
  mathematical form (stacked fractions, superscripts) with the current
  parameter values substituted, updating as you edit them.
- **Data table mode**: paste two columns (c, q) straight from Excel
  (tab/comma/semicolon/space separated; header rows ignored) or load a
  .csv/.txt file. Two interpolation choices: **use data as-is**
  (exact piecewise-linear through the points) or **smoothen**
  (Gaussian kernel smoother with an adjustable bandwidth slider,
  1–25% of the c-range; note kernel smoothers have some edge bias at
  the data boundaries — reduce the bandwidth if the ends matter).
  Raw points are drawn on the isotherm plot over the active curve.
  The data must span the cI–cF range (checked, with a clear message).
- **Isotherm zoom**: default view anchored at (0, 0); scroll-wheel
  zoom in/out about the cursor, drag to pan, −/+/⟲ buttons,
  double-click to reset.
- **RH ↔ c helper** for water vapour (Buck psat correlation, ideal gas).
- Structure summary line (segments front→back with σ values, exit
  window at ξ = 1 in τ and seconds). The initial and feed states can be set by typing values or by dragging the state markers directly along the isotherm curve — the full solution updates live during the drag.

Settings persist (localStorage).

## Validation

Checked against the analytic limits from the source MATLAB toolbox spec:

- Langmuir, cI = 0 → cF: **single shock** at σ = 1/(1 + ν·f(cF)/cF)
  (defaults give σ = 0.11442, exit τ = 8.740, t = 0.874 s ✓).
- Langmuir desorption: **pure simple wave** spanning
  σ = 1/(1+ν·f′(0)) … 1/(1+ν·f′(cF)) ✓.
- Type V step from cI ≈ 0: shock → wave → shock composite
  (semi-shocks: shock speed = adjacent wave speed at tangency ✓);
  desorption gives the mirrored structure ✓.

## Extending

The solver is a pure function of (f, cI, cF, ν) returning a segment
list — a finite-pulse mode (two Riemann problems + wave–shock
interaction) can be layered on later without touching the envelope
machinery.
