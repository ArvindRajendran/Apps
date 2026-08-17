# BAAM validation

Three-way cross-validation: the authors' original MATLAB `BAAM.m` →
Python reference → JS engine (the block embedded in `BAAM.html`).

## Files

- `baam_ref.py` — Python/scipy reference implementation (generic mixture
  interface: extended DSL analytic partials + binary IAST with FD partials;
  BDF integration at rtol 1e-9). Embeds the MATLAB anchors. Run:
  `python3 baam_ref.py` (~3 s; needs scipy).
- `baam_engine.js` — snapshot of the `ba-engine` block in `BAAM.html`.
  **Keep in sync with the app.**
- `baam_test.js` — jsc harness (264 checks). Run:
  `/System/Library/Frameworks/JavaScriptCore.framework/Versions/Current/Helpers/jsc baam_engine.js baam_test.js`

## Ground truth

`Papers/Adsorption/SimplifiedModel/MgMOF-74_Check/BAAM.m` (Subramanian
Balashankar's own code) driven headlessly in MATLAB R2026a with the 2019
Table-1 parameters:

- At the published ode15s tolerances (1e-6) it reproduces **Table 2 of the
  2019 paper to the displayed digits** for all four adsorbents.
- Anchors in the test files come from a run with tolerances tightened to
  RelTol 1e-10 / AbsTol 1e-12 (converged; the default-tolerance path values
  are up to ~3·10⁻³ off in places — do not re-anchor against them).

Agreement achieved: Python ≤ 4.4·10⁻⁶ (paths) / 8.7·10⁻⁷ (KPIs);
JS ≤ 5·10⁻⁴ (paths, noise-floor entries) / 5·10⁻⁵ (KPIs).

e-BAAM checks: Liske's worked example (13X SSL refit, P_H = 2 bar,
P_INT = 0.8, P_LOW = 0.07, η(P)) — paper values LPP 96.6/74.8 and
FP 96.6/74.0 reproduce at T = 298.15 K (the figure caption's 30 °C does not
match the quoted numbers; 298.15 K does, to all displayed digits).

Route equivalence: IAST with an equal-q_s Langmuir pair is mathematically
identical to the extended Langmuir — both engines confirm KPI agreement
(1·10⁻⁷ Python, ≤ 1·10⁻⁴ JS with the coarser IAST path grid).

## PVSA (P_H up to 5 bar) — checks A–G, added with the adjustable P_H

Six cases (13X SSL and CALF-20, P_H = 2 / 3 / 5 bar, LPP and FP, both
efficiency models) are checked in `baam_ref.py` (check F) and against the JS
engine (`baam_test.js` section G), plus two properties that only exist above
1 bar:

| Check | Result |
|---|---|
| JS vs Python across all six PVSA cases (Pu, Re, En, En_ADS, En_PR, WC, y_δ) | ≤ 2.0·10⁻⁸ |
| Vacuum work charged above 1 bar | exactly 0 (all grid points, both engines) |
| En_BLO + En_EVAC + En_ADS + En_PR = En | 1.3·10⁻¹⁶ (Python), machine (JS) |
| Live app at P_H = 5, P_INT = 2.5, P_LOW = 0.1 vs Python | 8.9·10⁻⁹, mole balance 2.6·10⁻¹⁴ |

**The ODE grid ΔP is deliberately not scaled with P_H.** The vacuum work is a
first-order sum over pressure slices, so its value depends on ΔP at O(ΔP): a
grid scaled to hold the step count constant shifted the energy by up to
1.2·10⁻³ relative and moved the design pressures off the grid (`baIdx` then
rounds them silently). Fixed ΔP keeps every slider pressure — a multiple of
0.005 bar — exactly on the grid at any P_H. Cost at P_H = 5 bar is 50k RK4
steps: 32 ms for the dual-site route, 0.9 s for the IAST route.

Compression work uses η evaluated at 1 bar (0.72 constant, or 0.761 with the
Maruyama correlation). That is the convention that reproduces the paper's
P_H = 2 bar worked example exactly — do not "improve" it to η(P_H).

## Not yet covered

- 6-step (1×PEQ) and 8-step (2×PEQ) e-BAAM cycles are not implemented.
  Anchors for the future implementation (from the e-BAAM paper's worked
  example, same conditions as above): 6-step 99.9/73.0 with P_eq = 0.23 bar;
  8-step 100.0/55.5 with P_eq,2 = 0.12, P_eq,1 = 0.37 bar.
