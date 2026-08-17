# DAC 0-D simulator — validation

Two independent implementations of the same model, plus a structural
consistency check against the group's Simplified DAC energy model in
MATLAB.

## Files

- `dac_ref.py` — independent Python/scipy reference (BDF, rtol 1e-10,
  flow/power quadratures carried as augmented ODE states, i.e. exact).
  Run `python3 dac_ref.py` (~1 min) → regenerates `dac_anchors.json`
  and the console table.
- `dac_anchors.json` / `dac_anchors.js` — 14 anchor cases (inputs, KPIs,
  CSS states). The `.js` file is the same data as a JS literal.
- `dac_engine.js` — snapshot of the `dac-engine` block in `DAC.html`.
  **Keep in sync with the app.**
- `dac_test.js` — jsc harness, 373 checks. Run:
  `/System/Library/Frameworks/JavaScriptCore.framework/Versions/Current/Helpers/jsc dac_engine.js dac_anchors.js dac_test.js`
- `dac_consist.js` — instant-kinetics consistency case (JS engine).
  `jsc dac_engine.js dac_consist.js > dac_consist.json`
- `dac_consist_check.m` — MATLAB comparison against `dac.energyModel`
  (needs the DAC EnergySimplifiedModel project on the path):
  `matlab -batch "addpath('…/DAC EnergySimplifiedModel'); run('dac_consist_check.m')"`

## A. JS engine vs Python reference

14 anchor cases: Lewatit base + adiabatic-pump variant + hot/humid/deep-
vacuum variant with wall mass + structure-mass variant (m_str = 2 kg/kg,
c_p,str = 710 J/kg/K) + edited-properties variant (ε_b 0.30, ε_p 0.45,
ρ_p 900, c_p,s 1200, −ΔH 80/48 kJ/mol, ΔH linked to the isotherm) +
unlinked-ΔH variant + Arrhenius-LDF variant (E_a 30/20 kJ/mol, both
coefficients T-dependent), APDES, APDES with E_a(CO₂) = 45 kJ/mol
(dry sorbent, water-inert branch), TMCM-41 (dry, inert water), SIFSIX
at 70 % RH, NbOFFIVE at 85 % RH (exercises the RH clamp), and two
custom site-competition cases (v1.7): Lewatit's numbers with the ψ/β
modulation replaced by b_w0 = 2·10⁻³ Pa⁻¹ (GAB water), and APDES's
numbers with the same b_w0 but inert water — competition through the
gas-phase p_w only, with ΔH(H₂O) supplied as a user input so b_w stays
T-dependent.
Per case: purity, recovery, E_th, E_fan, E_vac, productivity, water/CO₂,
CO₂/cycle, t_evac, all 7 CSS state components, and the simplified-model
comparison block (Δq_CO₂, Δq_H₂O, ΔT, η_cap, E_th/E_fan/E_vac).

**Result: 373 passed, 0 failed; worst relative deviation 6.4·10⁻⁶**
(a CSS state component; KPIs agree to better than 2·10⁻⁴ tolerance with
typical deviations ~10⁻⁶).

The Arrhenius LDF block (`ldf.*`) additionally checks the k(T) function
itself: E_a = 0 and E_a = undefined return the slider value unchanged,
k(T_ref) = k_ref, the two-temperature ratio matches the analytic
exp[−E_a/R (1/T₂ − 1/T₁)] to 10⁻¹⁴, and a full CSS run with
E_a explicitly set to 0 is bit-identical (tolerance 0) to one where the
fields are absent — the presets' published constant-k behaviour cannot
drift.

The custom-isotherm block (`cus.*`, v1.7) checks: a Custom sorbent
seeded with a preset's exact numbers runs a bit-identical cycle
(tolerance 0, checked for Lewatit — modulated/GAB — and TMCM-41 —
dry/inert); the site-competition form with b_w0 = 0 or p_w = 0 reduces
exactly to the plain Toth; and competition suppresses q*(CO₂)
monotonically in p_w — it can never increase uptake.

Additional property checks in the same harness:

| Check | Result |
|---|---|
| Mole balance over a CSS cycle (C, H₂O, N₂) | ≤ 2·10⁻⁵ (JS, trapezoid quadrature); ~10⁻¹³ in Python (exact quadrature) |
| Isothermal vacuum-work identity on the heat step: W = N_out·R·T_amb·ln(P_amb/P_regen)/η | 10⁻⁶ |
| Modified Toth reduces to dry Toth at q_w = 0; GAB closed form | machine precision |
| E_total = E_th + E_fan + E_vac closure | machine precision |
| Property overrides: m_s = ρ_p(1−ε_b), ε_T, c_p,s; ΔH linked → Toth ΔH₀ follows; unlinked → isotherm untouched; PRESETS never mutated; inert-water sorbent ignores an edited ΔH(H₂O) | machine precision |

## B. Instant-kinetics collapse onto `dac.energyModel` (MATLAB)

Lewatit, k_CO2 = k_H2O = 0.5 s⁻¹, U_a = 10⁷ W/m³/K, τ_heat = 20 s, long
steps, adsorbed-phase c_p zeroed in **both** models (the simplified
model heats Δq where the dynamic model heats the actual inventory — the
one acknowledged structural difference). The dynamic Δq's, ΔT and the
actual capture fraction feed the P-struct of `dac.energyModel`:

| Term | Simplified [MJ/kg] | Dynamic [MJ/kg] | rel. dev |
|---|---|---|---|
| E_sensible | 4.434 | 4.451 | 3.8·10⁻³ |
| E_CO2 des | 1.6132 | 1.6195 | 3.9·10⁻³ |
| E_H2O des | 3.666 | 3.680 | 3.9·10⁻³ |
| E_blower | 10.18228 | 10.18228 | **7·10⁻¹⁶** |
| E_vacuum | 0.8417 | 0.8453 | 4.3·10⁻³ |
| E_total | 20.738 | 20.778 | 2.0·10⁻³ |

The blower forms are identical, so that term is exact. Every Δq-based
term deviates by the same ~0.39 % factor = (moles desorbed)/(moles
collected as product) = 387.6/386.1: the simplified model neglects the
void-gas inventory (its own comment: "nvoid neglected"), while the
dynamic product excludes the CO₂ left in the voids at the end of
desorption. This is the expected, documented approximation of the
simplified model — not a discrepancy between implementations.

## C. In-app simplified-model comparison (the KPI "simplified" line)

`simplifiedEnergy()` (JS) / `simplified_energy()` (Python) evaluate the
same `dac.energyModel` forms with each cycle's own Δq, ΔT and capture
fraction — the quantity the app prints in small type under the energy
cards. Section B above is the same comparison driven all the way to the
instant-kinetics limit, where it must (and does) collapse to ≤0.43 %.
At the anchor conditions:

| case | Δq_CO₂ [mmol/g] | ΔT [K] | η_cap [%] | E_th simp/dyn | E_el simp/dyn |
|---|---|---|---|---|---|
| Lewatit-base | 0.728 | 100.0 | 35.1 | 1.002 | 0.999 |
| Lewatit-adia | 0.728 | 100.0 | 35.1 | 1.002 | 0.922 |
| APDES-base | 0.581 | 80.0 | 3.0 | 0.922 | 0.996 |
| TMCM-dry | 0.258 | 99.9 | 14.1 | 0.947 | 0.999 |
| SIFSIX-hum | 0.012 | 80.0 | 0.5 | 0.914 | 0.998 |
| NbOFFIVE-clamp | 0.015 | 100.0 | 1.0 | 0.948 | 0.998 |
| Lewatit-hotwet | 0.788 | 95.0 | 41.5 | 1.073 | 0.999 |
| Lewatit-struct | 0.730 | 100.0 | 35.2 | 1.000 | 0.999 |
| Lewatit-props | 0.890 | 100.0 | 56.6 | 1.008 | 0.999 |
| Lewatit-unlinkDH | 0.728 | 100.0 | 35.1 | 1.002 | 0.999 |

The electrical ratio is ~1 by construction (identical blower form, same
isothermal vacuum integral); `Lewatit-adia` is the exception because the
dynamic run uses the adiabatic pump while the simplified model is
isothermal — an intended structural difference, not an error. The
thermal ratio drifts a few per cent away from 1 as the CSTR's void gas
and swept-out sensible heat grow relative to a small working capacity,
and rises above 1 for `Lewatit-hotwet`, where the wall mass is charged
at the bed swing in the simplified term but at the actual wall swing in
the dynamic duty.

## D. Figure exports and the generated report (v1.8)

Checked in the browser, on the deployed build, in **both** colour
schemes. These are interface checks, not numerical ones — but the first
of them is the one that matters, because a report that quietly disagreed
with the screen would be worse than no report at all.

| Check | Result |
|---|---|
| Every KPI in the report equals the on-screen card | purity 98.70 %, recovery 34.46 %, E_th 11.537, E_el 6.644, E_tot 18.181 — identical |
| CO₂ balance within the report's own tables | product 346.96 + vented 659.92 = fed 1006.88 mol |
| Figure CSV derives from the drawn spec | 3-series isotherm CSV wide-format, 161 rows, markers A/B in the header |
| Sweep CSV | 9×9 = 81 rows, both axes + all 8 KPIs, blanks for failed points |
| Export forces the light palette | dark-mode screen → PNG corner pixel (255,255,255,255), 1.1 % dark ink in the x-axis band |
| Screen theme restored after export | `--bg` back to `#0f1420` |
| Report is self-contained | 0 external URLs (6 figures embedded as data URIs, 524 kB) |
| Equation/figure cross-references | 19 equations, 19 tags, 0 unresolved `{{eq:}}`/`{{fig:}}` |
| Uncomputed figures | named in place ("not computed in this session"), checkbox disabled, never a blank frame |
| No horizontal page scroll | body 420 vs 420 at a 420 px viewport; 14 tables and all equations scroll in their own boxes |

The 373-check numerical harness is unchanged by this work and was
re-run after it: **373 passed, 0 failed, worst 6.35·10⁻⁶** — identical to
the v1.7 figure, as it must be, since no engine file was touched.

## E. What is NOT validated

Nothing here is anchored to a **published** DAC cycle result. The chain
is (i) JS ↔ our own Python re-implementation of the same equations, and
(ii) the instant-kinetics collapse onto our own MATLAB simplified model.
Both catch implementation error; neither would catch a shared
misreading of Glaser et al. (2025) — a wrong constant, a mis-stated
step sequence, an omitted term. No number in this app has been compared
against a figure or table of that paper. If a reproducible case (given
sorbent, conditions and reported E_th / purity / recovery / productivity)
is ever digitised, add it as an eleventh anchor and record the deviation
here; that would be the first genuinely external check.

## Notes for future edits

- The Python quadratures are exact (ODE states); the JS engine uses
  trapezoid sums on accepted steps. If you tighten/loosen the JS error
  control, re-run the harness — KPI agreement should stay ≤ 1e-4.
- Anchor states are CSS states; they depend on the CSS tolerance
  (1e-5) and the initial state (equilibrated with feed). Change either
  and the anchors must be regenerated **in both implementations**.
- scipy struggles at the consistency-check settings (very stiff);
  that case is intentionally run on the JS engine (itself validated
  against Python at ordinary settings) — do not "fix" it back to
  Python without budgeting ~hours of BDF time.
