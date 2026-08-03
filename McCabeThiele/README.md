# McCabe–Thiele Distillation Designer

Single-file HTML app for teaching binary distillation design by the
McCabe–Thiele graphical method. Built as an educational tool for
undergraduate chemical engineering students.

## Model

- Constant pressure, constant relative volatility α:
  equilibrium y = αx / (1 + (α−1)x).
- Equilibrium stages by default; optional Murphree vapor efficiency
  E_MV (applied to plates via a pseudo-equilibrium curve drawn a
  fraction E_MV of the vertical distance from the operating line to
  the equilibrium curve; a partial condenser stays at true equilibrium).
- Feed at the optimal stage (operating-line switch at the q-line
  intersection); stages numbered from the top; partial reboiler counts
  as one equilibrium stage, and a partial condenser (selectable vs
  total) counts as stage 1.

## Inputs

α, x_D, z_F, x_B (mole fraction of the light key; x_B < z_F < x_D),
feed quality q (any value: subcooled to superheated, with a live
description of the feed state), reflux as R directly or as a multiple
of R_min, Murphree efficiency, condenser type. Every input has a
paired slider + number box with instant redraw.

## Outputs

- Diagram: equilibrium curve, y = x, q-line, both operating lines,
  R_min pinch construction (dashed), pseudo-equilibrium curve when
  E_MV < 100%, colour-coded staircase with stage numbers and feed
  stage highlighted, crosshair readout, PNG export.
- Numbers: total equilibrium stages (fractional + rounded up), feed
  stage, R_min (from the q-line pinch), N_min (Fenske), D/F from the
  overall mass balance, boil-up ratio, operating-line equations, and
  a stage-by-stage x/y table with section labels.
- Stage-by-stage reveal animation (step / play) for lecture use.
- Warns and suppresses the staircase when R ≤ R_min (infeasible).

## Validation

Cross-checked by hand for α = 2.5, z_F = 0.5, x_D = 0.95, x_B = 0.05,
q = 1, R = 1.5·R_min: R_min = 1.100, N_min (Fenske) = 6.43,
q-line/ROL intersection (0.500, 0.670), first stages x₁ = 0.8837,
x₂ = 0.7993 — all reproduced exactly. Pinch location and R_min also
verified for q = 0 and q = 0.5.

Self-contained, works offline, no external requests.
