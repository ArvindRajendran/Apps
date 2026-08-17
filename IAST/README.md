# IAST — Binary Ideal Adsorbed Solution Theory

`IAST.html` — mixture adsorption equilibria predicted from two
pure-component isotherms (Myers & Prausnitz 1965). Partial-pressure basis:
q in mol/kg, p in bar, ideal gas phase.

## Model

Reduced grand potential ("spreading pressure")
ψ_i(p) = ∫₀^p q_i⁰/p′ dp′, with the IAST conditions
ψ₁(p₁⁰) = ψ₂(p₂⁰), Py_i = x_i p_i⁰, Σx_i = 1, and
1/q_T = Σ x_i/q_i⁰(p_i⁰), q_i* = x_i q_T.

**Solver.** G(x₁) = ψ₁(Py₁/x₁) − ψ₂(Py₂/(1−x₁)) is strictly decreasing
with G(0⁺) = +∞ and G(1⁻) = −∞, so the root is unique and bracketed in
(0,1) from the start: 50 bisections + Newton polish
(dG/dx₁ = −q₁/x₁ − q₂/(1−x₁)), forward ψ evaluations only. No initial
guess, no divergence, closure residual ~10⁻¹⁵.

**Isotherm menu** (independent model + parameters per component): Linear,
Langmuir, Dual-site Langmuir, Anti-Langmuir, Quadratic, Sips/Hill, Toth,
Freundlich, BET, Type V (Langmuir + Hill step), custom q(p) expression.
ψ is analytic for all but Toth and custom (e.g. Langmuir q_s ln(1+bp);
quadratic q_s ln(1+b₁p+b₂p²); BET q_s ln[(1+(b_s−b_l)p)/(1−b_l p)]); the
numeric path uses a cumulative Newton–Cotes table in ln p (4th order,
20 000 intervals) with the p→0 tail handled by a local power-law exponent.

**Thermodynamic consistency, made visible.**
- Badges per component: finite Henry constant (Freundlich and Sips n<1
  are flagged — no Henry limit, so their low-pressure selectivities are
  artifacts of the input, not physics), monotonicity, pole location
  (anti-Langmuir, BET), ψ analytic vs quadrature.
- The **hypothetical pressures p_i⁰ are always displayed** with their
  extrapolation factors p_i⁰/P: the weakly adsorbed component is evaluated
  far above P (exponentially far when the strong component sits on a
  plateau). A user-set trust limit turns the readout red and the summary
  reports the worst extrapolation over the whole sweep — nothing is hidden.
- Exact structural result displayed and exploited: **S = p₂⁰/p₁⁰ is a
  function of ψ alone**, so binary IAST admits **no composition
  azeotropes** — S = 1 forces p₁⁰ = p₂⁰ = P, which happens only at the
  pressure p* where the ψ curves cross. The app finds p* (selectivity
  reversal), draws the S = 1 line on the heatmap (a horizontal line at
  P = p*, exactly), and labels the x–y diagram when P ≈ p*, where the
  whole curve degenerates onto the diagonal. Crossings are validated for
  transversality, and the fully **degenerate case** — two isotherms with
  (near-)identical spreading pressures, e.g. the same model with the same
  parameters on both components — is detected and reported as a single
  condition ("IAST-indistinguishable: S ≈ 1 and x ≈ y everywhere") rather
  than as hundreds of spurious crossings.

## Panels

- **Pure isotherms** with markers at the current p₁⁰, p₂⁰ (log/lin).
- **ψ-construction**: both ψ_i(p) with the equal-ψ horizontal, the current
  operating pressure, p_i⁰ markers, and ψ-crossing points — the graphical
  heart of IAST, and the extrapolation problem made visible.
- **x–y diagram** at the inspection pressure, with diagonal and optional
  extended-Langmuir overlay.
- **Heatmap** over (y₁, log P) of q₁*, q₂*, q_T, S₁₂ or x₁ (diverging
  palette centred on S = 1 when a reversal is in range, sequential
  otherwise); hover to inspect any point, click to pin; linked **slices**
  q*(P) at fixed y₁ and q*(y₁) at fixed P.
- **Extended Langmuir overlay** (when both components are Langmuir):
  equal q_s → coincides with IAST exactly; unequal q_s → visibly deviates,
  a one-toggle demonstration of why thermodynamic consistency matters.
- **3D surface** (bonus view, drag to rotate) of the selected quantity.
- **CSV export**: full sweep, x–y curve, both slices, each with a
  self-describing metadata header.

## Validation

Independent Python reference (`validation/iast_ref.py`, different
algorithm: bisection on p₁⁰) plus jsc-tested shipped engine
(`validation/`): 

| Check | Result |
|---|---|
| Equal-q_s Langmuir pair ≡ extended Langmuir (exact IAST result) | 1.1·10⁻¹³ (py), 1.3·10⁻¹⁵ (JS) |
| **pyIAST cross-validation** (Langmuir pair, 3 states) | agrees to 1.7·10⁻¹³ |
| Pure-component limits y₁ → 0, 1 | 3·10⁻¹¹ |
| Henry-limit selectivity → H₁/H₂ = 15.4 | 6·10⁻⁹ |
| Closure Σx−1 at solution | 0 (machine) |
| Toth(t=1) numeric-ψ ≡ Langmuir analytic | 1·10⁻¹³ (py), 4·10⁻⁷ (JS table) |
| BET ψ closed form vs quadrature | 2·10⁻¹² (py) |
| Cross anchors py↔JS (DSL/Sips pair, 4 states) | 1.6·10⁻⁷ rel |
| Selectivity reversal: p* = 4.2805 bar for Langmuir (2,20)/(6,0.8); S(p*) = 1 ∀y | 7·10⁻¹⁶ |
| No azeotropes at P ≠ p* (theory requires none) | confirmed |
| Identical quad–quad → degenerate flag (no spurious p*), x₁ = y₁ | machine precision |
| Genuine quad–quad crossing (b: 5,10 vs 8,3) | p* = 3/7 bar exactly |

Assumptions: ideal gas phase (p, not fugacity), ideal adsorbed solution
(no activity coefficients — RAST is the natural extension), isothermal,
binary. Settings persist in localStorage. Works offline — single file.
