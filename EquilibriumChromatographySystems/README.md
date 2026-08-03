# Equilibrium Theory — Coupled 2×2 Systems

Single-file interactive Riemann solver for **systems of two conservation
laws** in equilibrium chromatography:
`EquilibriumChromatographySystems.html` (~50 kB). The I and F states can be set by typing values or by dragging the I/F markers directly in the hodograph plane — the full solution updates live during the drag.

Works offline; share by
sending the file.

**v1 implements isothermal binary chromatography.** The adiabatic
single-component mode (c, T) is designed in (mode switch present,
disabled) — the engine is state-vector agnostic and the adiabatic flux
drops into the same machinery.

## Theory

System: (I + ν·∂q/∂c)·∂c/∂τ + ∂c/∂ξ = 0 with c = (c₁, c₂),
ν = (1−ε)/ε. Two characteristic families with

1. σ_k = 1/(1 + ν·λ_k),  λ_k = eigenvalues of ∂q/∂c

Riemann solution: feed F → intermediate state M via the **slow** family,
M → initial state I via the **fast** family; each transition is a simple
wave (σ increasing frontward) or a shock (Rankine–Hugoniot:
Δq₁/Δc₁ = Δq₂/Δc₂ = m, σ_sh = 1/(1+νm)).

## Two engines (hybrid architecture)

1. **Exact hodograph construction** — for the thermodynamically
   consistent constant-selectivity Langmuir
   (q_i = q_s·b_i·c_i/(1+b₁c₁+b₂c₂)): the hodograph characteristics are
   straight lines (h-transformation, Rhee–Aris–Amundson) and wave/shock
   images coincide, so M is the intersection of the slow line through F
   and the fast line through I. The construction self-verifies (RH
   consistency across both components, speed ordering, positivity) and
   falls back to the numerical engine if any check fails.
2. **Numerical reference solution** — for arbitrary isotherms
   (independent-q_s Langmuir, anti-Langmuir, custom q₁/q₂ expressions):
   conservative 2nd-order MUSCL scheme (minmod, local Rusanov flux,
   Newton inversion of w = c + νq) marched to τ = 1 and read off in the
   self-similar variable s = ξ/τ — the **numerically generated
   hodograph**. Wave structure detected from gradients; shocks appear
   slightly smeared but plateau values and speeds are accurate.

When both engines run (shared-q_s Langmuir), the cross-check panel
reports the intermediate-plateau deviation between them.

## Panels

Physical plane (colour map of c₁ / c₂ / c₁+c₂ with feature rays,
pinnable crosshair) · hodograph plane, plotted as (c₂, c₁), showing the
**full characteristic web** — integral curves of both families
numerically integrated through composition space (exactly straight for
constant-selectivity Langmuir, curved otherwise) — with the solution
path overlaid (wave/shock segments; numerical path dashed) ·
fluid-phase profiles c₁,c₂(ξ) · breakthrough curves c₁,c₂(τ) ·
adsorbed-phase profiles q₁,q₂(ξ) · engine cross-check.

## Validation record

Shared-q_s Langmuir (q_s = 5, b₁ = 1, b₂ = 2, ε = 0.37), saturation
(0,0) → (0.5,0.5):

- Exact: slow shock σ = 0.12803, fast shock σ = 0.17525, intermediate
  plateau **c_M = (0.809, 0)** — the classic component-1 **roll-up**
  (+61.8%) with complete displacement of component 2.
- Hand check: Δq₁/Δc₁ = Δq₂/Δc₂ = 4.000 exactly on the slow shock ✓.
- Independent numerical engine: plateau (0.8086, 0), speeds
  0.1278/0.1749 — 0.06% deviation.
- Desorption (swap): two rarefactions along the same hodograph lines
  (wave/shock image coincidence) ✓.
- Scalar limit (c₂ ≡ 0): reproduces the single-component shock
  σ = 0.11442 from the scalar app; represented as a degenerate
  zero-width-plateau double shock (flagged as such in the summary —
  eigenvalues cross on the composition-space boundary).
- Custom-expression path (same isotherm typed as text, numerical only):
  speeds within 0.2% of exact.

## Known limitations (v1)

- Numerical engine resolution: 900 cells over s ∈ [0, 1.06]; very slow
  features (σ < ~0.01) compress near the inlet and lose resolution.
- Non-genuinely-nonlinear isotherms (e.g. Type V generalizations) are
  handled by the numerical engine only; composite single-family waves
  appear correctly in profiles but are classified heuristically.
- Adiabatic (c, T) mode: next version.
