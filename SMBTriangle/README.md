# SMB Triangle Designer

`SMBTriangle.html` — equilibrium-theory ("triangle theory") design of
classical 4-section simulated moving bed units, entirely algebraic:
no numerical PDE solution anywhere.

Sources: Rajendran, Paredes & Mazzotti, *J. Chromatogr. A* 1216 (2009)
709 (§3, Table 2 — complete separation regions for the generalized
Langmuir isotherm); Rajendran, *J. Chromatogr. A* 1185 (2008) 216
(reduced-purity design, linear isotherms, restrictive case); triangle
theory per Storti/Mazzotti/Morbidelli.

## Model

Flow-rate ratios m_j = (Q_j t* − Vε)/(V(1−ε)). Isotherm menu:

- **Linear** n_i = H_i c_i — the exact triangle H_B < m₂ < m₃ < H_A,
  m₁ > H_A, m₄ < H_B, **plus** the reduced-purity pentagon PQRST:
  for specified extract/raffinate purities P_E, P_R the region is
  bounded by m₃ = H*_A, m₂ = H*_B and the two purity lines pivoting
  through (H_A, H_A) and (H_B, H_B). The exact-purity operating point
  P = (H*_B, H*_A) is marked (and the app notes its non-robustness).
- **Generalized Langmuir** n_i = H_i c_i/(1 + p_A K_A c_A + p_B K_B c_B),
  the four classes L (+,+), A (−,−), M1 (−,+), M2 (+,−). Boundaries
  from the ω-quadratic at feed composition and Table 2 of the review:
  curved sides ar/bs, straight sides rw/sw/bw/aw, and the
  (m₂,m₃)-dependent m₁,min / m₄,max square-root expressions where the
  class requires them.

A geometric point worth noting: the straight sides rw and sw are
**tangent** to the curved sides ar and bs at the junction points r and
s (the boundary is smooth there). The tangency points are computed in
closed form — m₂ʳ = H_A/[1 + K_A c_A^F(l₂+l₃)/l₃]² from
d/dm₂[rw(m₂, ar(m₂))] = 0, and the mirror expression for s — and both
are marked on the plot. (Naive intersection-finding fails here
precisely because tangency gives no sign change.)

## Panels

- **(m₂, m₃) plane**: complete-separation region (a-side chain red =
  crossing pollutes the raffinate; b-side chain blue = crossing
  pollutes the extract), dashed linear "ghost" triangle for reference,
  labelled points a, b, r, s, w, pentagon overlay (linear + reduced
  purity), zone labels in every region (pure extract / pure raffinate /
  both polluted / complete separation / flooded), and a draggable
  operating point W with live regime classification. A **hover
  crosshair** shows (m₂, m₃) at the cursor — and for linear isotherms
  also the **predicted purities P_E and P_R at that point**, from
  closed-form inversions of the 2008 design equations (region 2:
  P_R = c_B^F(m₃−m₂)/[c_B^F(m₃−m₂)+c_A^F(m₃−H_A)]; region 3 mirror;
  region 4 decoupled in m₂ and m₃; each → 100% on the triangle
  boundary, and the region-4 P_E reproduces the paper's Eq. 28 minimum
  at the zero-flow limit m₂ = −ε/(1−ε)). The **feed-dilution slider**
  morphs the region live — Langmuir shrinking down-left, anti-Langmuir
  up-right (review Fig. 9).
- **(m₁, m₄) plane**: feasible regeneration region for the current
  (m₂, m₃) point (m₁,min / m₄,max recompute as W moves for classes
  L, A, M1), corner q = minimum desorbent, zone labels, hover
  crosshair with (m₁, m₄) readout, draggable point with violation
  warnings.
- **Safety factor β**: β = 1 places the operating point at the vertex
  (max productivity, zero robustness); larger β retreats toward the
  region interior, and m₁ = β·m₁,min, m₄ = m₄,max/β.
- **Open-loop design card**: column L, d, ε, d_p, μ and ΔP_max
  (or a measured ΔP at reference flow) → Blake–Kozeny column
  resistance → closed-form minimum switch time
  t*_min = (R_col V/ΔP_max)·Σ S_j[(1−ε)m_j + ε], then all internal
  flows Q₁…Q₄, external streams (desorbent, extract, feed, raffinate,
  section-4 outlet), per-section ΔP, productivity, and solvent
  consumption. Open loop: sections 1→4 in hydraulic series, outlet at
  ambient, so the constraint is the desorbent-pump back pressure.

## Validation

The boundary engine was validated first in an independent Python
reference (scratchpad smb_ref.py, 19 checks), then the shipped JS was
verified to reproduce it exactly:

| Check | Reference | App |
|---|---|---|
| ω roots, case L (H 2/1, K 0.1, c 16 g/L) | 0.30724 / 1.54990 (hand) | identical |
| ω roots, case A (c 2) | 1.18350 / 2.81650 (hand) | identical |
| ω roots, case M1 (c 4) | 0.75660 / 2.64340 (hand) | identical |
| ω roots, case M2 (c 1.5) | 1.25 / 1.60 (hand, exact) | identical |
| Vertices w for the four review Fig. 8 cases | (0.775, 0.916), (1.449, 2.367), (1.182, 1.636), (1.000, 2.000) | identical |
| Dilute limit (all classes) | w → (H_B, H_A), linear triangle | pass |
| Region endpoints | boundaries meet diagonal at (H_B,H_B), (H_A,H_A) | pass |
| Pentagon H*_A, H*_B (2008 case: H 3/2, racemic, P 99%) | 3.0102 / 1.9898 | identical |
| Purity lines | Eq. 26 through Q and R; Eq. 27 through S and T | pass |
| t*_min | closed form reproduces ΔP_max exactly | pass |
| Open-loop node balance | Q_D + Q_F = Q_E + Q_R + Q_W | exact |

A notable exact result surfaced during validation: for class M2 the
vertex w sits exactly at the linear vertex (H_B, H_A) while the region
sides bulge — reproduced by both implementations.

Additional r/s tangency checks (Python + JS identical): case L at
16 g/L, r = (1.2011, 1.2644) lies on rw with vanishing derivative;
case A at 2 g/L, s = (1.2323, 1.4007) on sw likewise; case M2 both.
Purity-map anchors: along the H*_A locus (P = 99%) the inverted
formula returns P_R = 0.9900 exactly, and along the Eq.-26 line the
region-2 inversion returns 0.9900 exactly.

Equilibrium theory neglects axial dispersion and mass-transfer
resistance; treat the boundaries as the ideal limit and design with a
safety margin (β > 1) in practice.

Settings persist in localStorage. Works offline — single file, share freely.
