# Vapour Pressure Plotter

`VaporPressure.html` — saturation curves for **653 components**, drawn from
published correlation coefficients and overlaid on one plot. Single file
(~294 kB), works offline, share by sending the file.

## Features

- **Component picker** — search by name, formula, CAS number or common
  synonym (`n-hexane`, `r-134a`, `isobutane` all work). Up to ten
  components on one plot, each with its own colour and its own choice of
  correlation.
- **Coordinate modes**
  - *P vs T* — linear or logarithmic pressure axis.
  - *Cox chart* — ln P against 1/T, where the curves become nearly
    straight; the x-axis carries both 1000/T and the corresponding
    temperature.
  - *Reduced* — P/P<sub>c</sub> against T/T<sub>c</sub>, on which every
    component converges on the critical point (1, 1) and the
    corresponding-states collapse is visible.
- **Units** — T in K / °C / °F, P in Pa, kPa, MPa, bar, atm, mmHg or psia.
- **Validity ranges shown honestly** — each curve is solid inside the
  correlation's published range and dashed outside it. Antoine fits are
  narrow (the Landolt fit for water covers only 323–393 K), and
  extrapolating them is the classic way to get a badly wrong number.
- **Markers** for normal boiling points and critical points.
- **Isobar** — enter a pressure to draw a horizontal line and get each
  component's boiling temperature at that pressure (bisection on the
  correlation). Useful for vacuum-distillation reasoning.
- **Sublimation curves** for the 92 components that have them.
- **Crosshair** reporting every plotted curve's value at the cursor
  temperature, flagging extrapolated values.
- **Comparison points A and B** — click the chart to place two
  temperatures (drag to move, or type exact values), and every plotted
  curve is read at both, with the ratio B/A. On this chart a "point" is a
  temperature, because the x-axis is temperature in every mode, so one
  marker reads a value off every curve at once. In the *Compare methods*
  tab the rows become the correlations, so the markers give a direct
  read-off of how far the correlations disagree at a chosen temperature.
  With two points placed on a pressure curve the app also reports the
  **two-point Clausius–Clapeyron** mean ΔH<sub>vap</sub> between them,
  from ln(P_B/P_A) = −(ΔH/R)(1/T_B − 1/T_A).
  In reduced coordinates a single vertical line would be wrong — one
  temperature is a different T<sub>r</sub> for each fluid — so the markers
  are drawn there as one dot per curve instead.
- **Compare methods** — one component, every correlation it carries,
  overlaid with a deviation table relative to the highest-ranked fit. 610
  of the 653 components support this. The spread between independent fits
  is a fair estimate of how well that compound's vapour pressure is
  actually known.
- **Relative volatility α(T)** — α = P₁ˢᵃᵗ/P₂ˢᵃᵗ for any pair, with a
  geometric mean over the range to carry into a constant-α short-cut
  design (the [McCabe–Thiele designer](../McCabeThiele/McCabeThiele.html)
  or the Fenske equation).
- **ΔH<sub>vap</sub>** from the Clapeyron slope of the selected
  correlation.
- **Export** — chart as PNG (opaque background), data as CSV.
- Selection and settings persist in localStorage.

## Correlations

All return pressure in Pa. Each is a direct transcription of the
corresponding function in the [`chemicals`](https://github.com/CalebBell/chemicals)
package (MIT), so the JavaScript and Python engines agree to machine
precision.

| Method | Form | Components |
|---|---|---|
| Wagner (2.5,5) | ln(P/P<sub>c</sub>) = (aτ + bτ^1.5 + cτ^2.5 + dτ⁵)/T<sub>r</sub> | 104 |
| Wagner (3,6) | ln(P/P<sub>c</sub>) = (aτ + bτ^1.5 + cτ³ + dτ⁶)/T<sub>r</sub> | 245 |
| VDI PPDS | Wagner (2.5,5) form, VDI coefficients | 274 |
| DIPPR-101 | ln P = C₁ + C₂/T + C₃ ln T + C₄T^C₅ | 339 |
| Antoine extended (TRC) | log₁₀P = A − B/(T+C) + 0.43429xⁿ + Ex⁸ + Fx¹² | 97 |
| Antoine | log₁₀(P/Pa) = A − B/(T+C) | 325 |
| Antoine (Landolt) | **ln**(P/Pa) = A − B/(T+C) | 596 |
| Sublimation | ln(P/Pa) = A − B/(T+C) | 92 |
| Ambrose–Walton | corresponding states, needs T<sub>c</sub>, P<sub>c</sub>, ω | 648 |
| Lee–Kesler | corresponding states, needs T<sub>c</sub>, P<sub>c</sub>, ω | 648 |

`τ = 1 − T_r`. When "Auto" is in force the app picks the first available
method in the order listed above — fitted multi-parameter forms first,
narrow Antoine fits last, and the predictive corresponding-states methods
only when nothing was fitted.

**The conventions differ between tables and this matters.** Landolt's
Antoine coefficients are natural-log based; Poling's are base-10. Running
the Landolt coefficients for water through the log₁₀ form returns
3.4 × 10¹¹ Pa instead of 101 kPa — seven orders of magnitude out. Each
method therefore carries its own evaluator rather than sharing one
"Antoine" routine.

## Component set

The 653 components are the union of the six engineering coefficient
tables (Poling/Reid, Perry's, VDI, McGarry) — a compound is included if at
least one standard handbook carries it. The much larger
Landolt–Börnstein Antoine table (6346 compounds) is used only to add a
method to components already in the set, never to add new ones, since its
entries are the least curated.

- 610 of 653 carry two or more correlations
- 648 have T<sub>c</sub>, P<sub>c</sub> and ω, so the predictive methods work
- 92 have a sublimation curve

## Data sources

Coefficients are redistributed from the `chemicals` package (MIT licence),
which compiles them from:

- Poling, Prausnitz & O'Connell, *The Properties of Gases and Liquids*,
  5th ed., Appendix A — Antoine, extended Antoine, Wagner (2.5,5)
- McGarry, *Ind. Eng. Chem. Process Des. Dev.* **22** (1983) 313 — Wagner (3,6)
- Perry's *Chemical Engineers' Handbook*, 8th ed., Table 2-8 — DIPPR-101
- VDI *Heat Atlas*, 2nd ed. (2010) — PPDS
- Landolt–Börnstein, New Series IV/20 — Antoine and sublimation

Critical constants, acentric factors, boiling and melting points come from
the same package.

## Validation

**Engine.** `make_reference.py` evaluates every component/correlation pair
at nine points across its published range using the Python implementations
in `chemicals`, producing **29 776 reference cases**. Running the same
cases through the JavaScript engine shipped in the app:

| Correlation | Cases | Max relative error |
|---|---|---|
| DIPPR-101 | 2 933 | 5.7 × 10⁻¹⁴ |
| Antoine (Landolt) | 5 348 | 1.1 × 10⁻¹⁴ |
| Wagner (2.5,5) | 802 | 1.1 × 10⁻¹⁴ |
| VDI PPDS | 2 327 | 1.1 × 10⁻¹⁴ |
| Antoine | 2 921 | 6.2 × 10⁻¹⁵ |
| Antoine extended | 851 | 4.2 × 10⁻¹⁵ |
| Wagner (3,6) | 2 111 | 3.6 × 10⁻¹⁵ |
| Lee–Kesler | 5 832 | 3.6 × 10⁻¹⁵ |
| Sublimation | 819 | 1.3 × 10⁻¹⁵ |
| Ambrose–Walton | 5 832 | 3.5 × 10⁻¹⁶ |

All at the level of floating-point round-off; no non-finite results.

**Physical checks.**

| Check | Result |
|---|---|
| Benzene at 353.24 K (its NBP) | Wagner(3,6) 101 294 Pa, Wagner(2.5,5) 101 323, PPDS 101 331, Antoine 101 568, Landolt 101 431 — all within 0.25 % of 101 325 |
| Water at 373.15 K | Wagner(3,6) 101 285 Pa, DIPPR 101 261, Antoine 101 047 vs IAPWS-95 101 418 |
| T<sub>sat</sub>(P<sub>sat</sub>(T)) round trip | max error 1.1 × 10⁻¹³ K |
| Boiling point at 1 atm | benzene 353.24 K (tabulated 353.22), toluene 383.79 (383.75), water 373.16 (373.12) |
| α benzene/toluene at 373 K | 2.427 (textbook ≈ 2.4) |
| ΔH<sub>vap</sub> water at 373 K | 41.33 kJ/mol vs accepted 40.65 — the +1.7 % is the ideal-vapour assumption (ΔZ ≈ 0.984) |
| Reduced coordinates | all components terminate at (T<sub>r</sub>, P<sub>r</sub>) = (0.999, 0.99), i.e. on the critical point |
| Comparison points at 50 °C / 100 °C | water 12.337 / 101.28 kPa (accepted 12.35 / 101.32), benzene 36.188, toluene 12.286 kPa — and the marked temperatures survive a change of temperature unit or coordinate mode |
| Two-point Clausius–Clapeyron, water 323→373 K | 42.22 kJ/mol, matching the mean of the local values (43.0 at 323 K, 41.3 at 373 K) |
| Sublimation convention | P<sub>sub</sub>/P<sub>liq</sub> at the top of the sublimation range has median 0.985 over all 92 components — the solid and liquid branches meet at the triple point, confirming the natural-log reading |

## Known data limitations

These are properties of the published data, not of the code. They are
surfaced in the app rather than hidden.

- **Correlations and tabulated boiling points sometimes disagree.** Over
  the 613 components whose T<sub>b</sub> falls inside the auto-chosen
  correlation, the deviation from 1 atm at T<sub>b</sub> has a median of
  0.37 % and a 90th percentile of 3.8 % — but **51 exceed 5 %** and 19
  exceed 20 %. The coefficients and the boiling point generally come from
  different laboratories. The component table shows this check per
  component and flags deviations above 5 %; when it fires, the *Compare
  methods* tab usually shows the other correlations agreeing with each
  other against one outlier. By method, Wagner (3,6) is the most reliable
  (median 0.18 %, p90 1.5 %) and DIPPR-101 the least among auto-chosen
  fits (median 1.17 %, p90 7 %).
- **Landolt is the most frequent outlier**, which is why it ranks last
  among the fitted methods. Isobutyraldehyde's Landolt fit misses 1 atm at
  its boiling point by 654 % while its Wagner fit is exact.
- **Two records were internally impossible and have been repaired** by
  `generate_data.py`, which reports them on stderr: phenanthrene was
  listed with T<sub>c</sub> = 0.869 K (a factor-1000 slip — its own
  T<sub>b</sub> is 611 K), so its critical data is dropped and only its
  three fitted correlations remain; 1-buten-3-yne (vinyl acetylene) had
  T<sub>m</sub> = 476 K above its T<sub>c</sub> of 455 K, so the melting
  point is dropped.
- **Above the critical point there is no saturation state.** Several
  correlations silently clamp T<sub>r</sub> to 1 and return P<sub>c</sub>
  rather than failing, which would show up as a flat curve or a
  ΔH<sub>vap</sub> of exactly zero. The app refuses to evaluate there.
- **ΔH<sub>vap</sub> uses Clausius–Clapeyron**, ΔH = RT² d(ln P)/dT, which
  assumes an ideal vapour and neglects the liquid molar volume. Near the
  normal boiling point this costs a few per cent. Both assumptions
  collapse as T → T<sub>c</sub>, where the true ΔH<sub>vap</sub> must fall
  to zero but this estimate turns *upward*; curves are dashed above
  T<sub>r</sub> = 0.85 for that reason.
- **α(T) is the ideal-solution relative volatility**, assuming Raoult's
  law and an ideal vapour. It is a starting point for a non-ideal mixture,
  not an answer.
- **Sublimation ranges are narrow** — median span 36 K, typically ending a
  few kelvin below the melting point where the measurements stop. CO₂'s
  Landolt sublimation fit covers only 79–83 K and should not be
  extrapolated to its 194.7 K sublimation point.

## Rebuilding

```bash
pip install chemicals
python generate_data.py > vp_data.json     # curate + normalise the tables
python build.py VaporPressure.html         # inline the JSON into the template
python make_reference.py                   # 29 776 Python reference values
```

Then open the app and run the reference comparison against `reference.json`
in the browser console to re-verify the engine.
