# BAAM — Batch Adsorber Analogue Model

`BAAM.html` — rapid screening of 4-step VSA/PVSA cycles for binary CO₂/N₂
separation using the batch adsorber analogue model of Subramanian Balashankar,
Subraveti, Li & Rajendran (*Ind. Eng. Chem. Res.* **58** (2019) 3314–3328),
with extensions from the e-BAAM of Liske & Rajendran (*Can. J. Chem. Eng.*,
2026): feed-pressurization cycle, pressure-dependent vacuum-pump efficiency,
PVSA operation (P_H > 1 bar) with feed-compression work, and the e-BAAM
US-DOE classification and energy calibration.

## The model in one paragraph

A well-mixed, isothermal batch adsorber (1 kg adsorbent, bed voidage ε) at
instantaneous equilibrium. Each species' inventory is
N_i = Py_iVε/R_gT + wq_i*(P, y). Blowdown and evacuation reduce to a single
ODE for the gas composition along the pressure path,

dy/dP = (a₁y − a₂)/(f₂ − f₁y),

with coefficients built from the competitive-isotherm partial derivatives;
pressurization (with raffinate = LPP, or feed = FP) and adsorption are small
algebraic balances (the raffinate leaves at the pre-adsorption composition).
The cycle closes exactly in one pass — no cyclic-steady-state iteration — so
purity, recovery, energy (adiabatic vacuum work below 1 bar, η = 72 % constant
or the Maruyama η(P) = 0.8·19.55P/(1+19.55P)), and working capacity respond
instantly to the three design pressures P_H > P_INT > P_LOW.

## Features

- **Adsorbent presets**: the 2019 paper's four (Mg-MOF-74, zeolite 13X,
  UTSA-16, CS-AC; dual-site Langmuir, Table 1) and the e-BAAM four
  (13X and UTSA-16 single-site refits, IISERP MOF-2, CALF-20), all editable.
- **Mixture equilibrium, two routes**: extended dual-site Langmuir
  (concentration basis, analytic partials — the papers' route), or **binary
  IAST over any pair of pure isotherms** from the IAST app's menu (Langmuir,
  DSL, Sips, Toth, BET, quadratic, Type V, custom expressions, …) with
  finite-difference partials feeding the same ODE.
- **Cycle visualization**: the papers' signature plots — cycle transitions
  α→β→γ→δ drawn on the competitive CO₂ and N₂ isotherms, over a
  colour-coded family of constant-composition curves (log/lin pressure axis),
  plus the composition path y(P) and a live cycle schematic. Moving the
  cursor over either isotherm plot drops a **crosshair** that reads off the
  loading on every composition curve at that pressure, and on the cycle path
  itself.
- **Three pressure levels, all adjustable**: P_H (0.2–5 bar), P_INT and
  P_LOW on sliders that enforce P_H > P_INT > P_LOW. Above 1 bar the cycle is
  a PVSA: the feed (and, in the FP cycle, the pressurization gas) is charged
  adiabatic compression work from 1 bar to P_H at η(1 bar), while
  depressurization above 1 bar is free — only the sub-atmospheric part of
  blowdown and evacuation is charged to a vacuum pump.
- **Design-space scan**: purity/recovery/energy/working-capacity maps over
  (P_LOW, P_INT) at the chosen P_H, with the DOE-feasible region outlined
  (Pu ≥ 95 %, Re ≥ 90 %) and — for a PVSA — the 1 bar line above which no
  vacuum work is charged; click any point to load it into the sliders. Purity–recovery
  Pareto curves per P_LOW with the e-BAAM classification radius.
- **US-DOE verdict**: r_max = max √(Pu²+Re²) compared against the e-BAAM
  linear-SVM cutoff r_cut = 124.5 (82 % accuracy); minimum energy among
  DOE-feasible scanned points, with the e-BAAM full-model calibration
  E_full ≈ 1.58·E_BAAM. The 2019 VSA-specific calibration
  (r₉₅₋₉₀ = 110.25, E_full = 1.1446·E + 66.53) is quoted in the app note.
- **Feed presets** NGCC 4 % / coal 15 % / cement 25 % CO₂; free y_feed and T.
- **Hand-off**: receives fitted parameters from the Isotherm Fitter
  (`#import=` URL hash) — Langmuir/DSL fits land in the extended-DSL slots
  (global van 't Hoff fits keep their temperature dependence; pressure-basis
  fits are converted with U_c = U_p − RT̄), any other pressure-basis suite
  model lands as an IAST component.
- **Exports and reports**: ⬇ PNG / ⬇ CSV on every figure — the composition
  path (with cumulative inventories and vacuum work), the isotherm
  families, the design-space grid and the Pareto curves. Figures are
  exported in the light palette on an opaque white background whatever
  the screen theme, so they are publication-ready either way.
  **Generate report** (top right) writes a self-contained HTML document:
  introduction, objectives, the model equations in MathML with a
  nomenclature table, the numerical method, every input including the
  isotherm parameters actually used, the cycle results, the design-space
  summary with the classification verdict, and the selected figures with
  numbered captions. Printable to PDF; no network access.
  *Filenames changed slightly in v2.1* — one basename now serves both the
  PNG and the CSV of a figure (`baam_composition_path.*`,
  `baam_design_space.*`, `baam_isotherms_co2.*`, `baam_isotherms_n2.*`,
  `baam_pareto.*`).

## Validation

Ground truth is **the authors' original MATLAB `BAAM.m`** (Balashankar), run
headlessly in MATLAB R2026a with the Table-1 parameters. At its published
tolerances it reproduces Table 2 of the paper to the displayed digits; with
ode15s tolerances tightened to 1e-10/1e-12 it provides converged
machine-precision anchors. The Python reference (`validation/baam_ref.py`,
scipy BDF 1e-9) and the JS engine (`validation/baam_test.js`, jsc; fixed-step
RK4 on the ΔP = 10⁻⁴ bar grid) are checked against those anchors:

| Check | Result |
|---|---|
| Paper Table 2 (all 4 adsorbents; Pu/Re/En-BLO/En-EVAC/En/WC at P_LOW=0.03, P_INT=0.15) | reproduced to the displayed digits |
| BLO/EVAC path (y, q_CO₂, q_N₂, cumulative moles, vacuum work) vs converged MATLAB, 4 adsorbents × 4 pressures | ≤ 4.4·10⁻⁶ rel (Python), ≤ 5·10⁻⁴ rel (JS, noise-floor entries) |
| Cycle KPIs vs converged MATLAB, 4 adsorbents × 4 (P_LOW, P_INT) | ≤ 8.7·10⁻⁷ rel (Python), ≤ 5·10⁻⁵ (JS) |
| e-BAAM worked example (13X SSL refit, P_H=2, P_INT=0.8, P_LOW=0.07, η(P)): LPP 96.6/74.8, FP 96.6/74.0 | reproduced (96.57/74.77 and 96.57/73.99 at 298.15 K) |
| IAST route ≡ extended Langmuir (equal-q_s Langmuir pair) | KPIs agree to ≤ 1·10⁻⁴ (JS), 1·10⁻⁷ (Python) |
| PVSA, P_H = 2/3/5 bar (6 cases, LPP + FP, both η models): JS vs Python | ≤ 2.0·10⁻⁸ rel |
| Vacuum work charged above 1 bar | exactly zero, every grid point |
| Energy split En_BLO + En_EVAC + En_ADS + En_PR = En | 1.3·10⁻¹⁶ (Python) |
| Cycle mole-balance closure | ≤ 10⁻⁸ % (algebraic steps solved exactly) |

Notes: the e-BAAM paper's worked-example figure is captioned 30 °C but the
quoted Pu/Re values correspond to T = 298.15 K — both codes here reproduce
them at 298.15 K. Working capacity uses M = 44 g/mol as in `BAAM.m`.
Compression work is evaluated at η(1 bar), the convention that reproduces the
P_H = 2 bar worked example; the ODE grid ΔP is held fixed as P_H rises rather
than scaled, because the vacuum work is a first-order sum over pressure slices
and a scaled grid shifts it (see `validation/README.md`).

To re-run: `python3 validation/baam_ref.py` (needs scipy) and
`jsc validation/baam_engine.js validation/baam_test.js`.

## References

- V. Subramanian Balashankar, A. K. Rajagopalan, R. de Pauw, A. M. Avila,
  A. Rajendran, "Analysis of a Batch Adsorber Analogue for Rapid Screening of
  Adsorbents for Postcombustion CO₂ Capture", *Ind. Eng. Chem. Res.* 58
  (2019) 3314–3328.
- G. Liske, A. Rajendran, "Extended Batch Adsorber Analogue Model (e-BAAM)
  for Rapid Screening of PSA Cycles for CO₂ Capture", *Can. J. Chem. Eng.*
  (2026), doi:10.1002/cjce.70453.

Settings persist in localStorage. Works offline — single file, share freely.
