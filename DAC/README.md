# DAC — 0-D Cycle Simulator

Dynamic 0-D (CSTR) simulation of a 4-step temperature–vacuum swing
direct-air-capture cycle, run to cyclic steady state (CSS), for six
sorbents. One self-contained HTML file — open `DAC.html`.

## The cycle

Basis: 1 m³ of packed contactor (bed). Four steps:

1. **Adsorption** — fan-driven ambient air (400 ppm CO₂, H₂O from RH,
   balance N₂) flows through the bed at P_amb; outlet is at the
   (well-mixed) column composition. Duration set by slider.
2. **Evacuation** — inlet closed; a vacuum pump with a set volumetric
   speed draws the column down. The step ends when P reaches P_regen,
   so t_evac is an *outcome*, not an input.
3. **Heating + vacuum desorption** — the wall approaches T_des as a
   first-order process (time constant τ_heat) while the pump holds
   P = P_regen; everything leaving the column in this step is the
   **product** stream.
4. **Cooling / repressurisation** — the bed is repressurised
   instantaneously with ambient air, then cools at constant P_amb with
   inflow only (the wall relaxes to T_amb with τ_cool). Passive — no
   energy charged.

## Model equations

State: x = [N, y_CO2, y_H2O, q_CO2, q_H2O, T, T_W] (total gas moles in
the voids, gas mole fractions, loadings, bed and wall temperature).

Molar balance (CSTR, per species i, Glaser et al. eqn 1):

$$\frac{d(N y_i)}{dt} = \dot n_I y_{I,i} - \dot n_O y_i - k_i (q_i^* - q_i)\, m_s \tag{1}$$

Energy balance (lumped bed = solid + structure + adsorbed phase + gas,
Glaser eqn 2):

$$C_{tot}\frac{dT}{dt} = \dot n_I c_{p,g}(T_{amb}-T) + V_C\frac{dP}{dt}
 + \sum_i \Delta H_i\, k_i (q_i^*-q_i)\, m_s + U_a (T_W - T) \tag{2}$$

$$C_{tot} = m_s\big(c_{p,s} + m_{str}c_{p,str}\big) + \sum_i n_i c_{p,g,i}
 + m_s\big(q_c c_{p,a,CO_2} + q_w c_{p,a,H_2O}\big)$$

The structure term uses the same mass-ratio form as the Simplified DAC
energy model (`mstr*Cpstr`): m_str is kg of contactor structure per kg
of sorbent, with its own c_p. Default m_str = 0 (bare sorbent); a wall
thermal mass m_W c_W can be added separately and is charged to the
sensible duty during heating.

LDF kinetics (eqn 3), wall (eqn 4, prescribed first-order approach):

$$\frac{dq_i}{dt} = k_i\,(q_i^* - q_i), \qquad
\frac{dT_W}{dt} = \frac{T_{set} - T_W}{\tau} \tag{3,4}$$

where T_set = T_des with τ_heat during heating and T_amb with τ_cool
otherwise. The wall has no thermal resistance of its own; it couples to
the bed through the volumetric coefficient U_a (W m⁻³ K⁻¹). An optional
wall thermal mass adds m_W c_W ΔT_W to the sensible duty.

The LDF coefficients can carry an Arrhenius temperature dependence
(eqn 5):

$$k_i(T) = k_{i,20}\,
  \exp\!\left[-\frac{E_{a,i}}{R}\left(\frac{1}{T} -
  \frac{1}{293.15\ \mathrm{K}}\right)\right] \tag{5}$$

The k sliders set the 20 °C values; E_a is entered per component in
kJ/mol (0–80, default 0 = the constant-k model of the presets). With
E_a > 0 desorption at T_des is faster than capture at ambient — e.g.
E_a = 30 kJ/mol gives k(120 °C) ≈ 23× k(20 °C) — which typically
raises recovery and lowers thermal energy because more CO₂ leaves the
sorbent within the same heating window. E_a(CO₂) is also available as
a sweep axis.

Per-step constraint: constant-P steps slave N through
dN/dt = −(N/T)(dT/dt) and solve the resulting linear equation for dT/dt;
the evacuation step advances N directly with
ṅ_O = (S_pump/3600)·N/V_C and includes the V_C dP/dt expansion term.

### Equilibria

CO₂ — Toth with temperature-dependent b and ω, modified for water
following the Stampi-Bombelli protocol (paper eqns 5–9):

$$q^*_{CO2} = \frac{q_s b P_{CO2}}{[1+(bP_{CO2})^\omega]^{1/\omega}},\quad
q_s = \frac{q_{s,0}e^{\chi(1-T/T_0)}}{1-\psi q_{H2O}},\quad
b = b_0 e^{\frac{\Delta H_0}{RT_0}(\frac{T_0}{T}-1)}(1+\beta q_{H2O}) \tag{6}$$

H₂O — GAB on the relative-humidity fraction x (no T dependence beyond
p_sat; paper eqn 10):

$$q^*_{H2O} = \frac{C_m C_G K x}{(1-Kx)(1+(C_G-1)Kx)} \tag{7}$$

Water activity uses Antoine's equation with the paper's constants
(log₁₀P[bar] = 4.6543 − 1435.264/(T − 64.848)). The coupling loading is
clamped per the source paper (NbOFFIVE: 75 % RH, SIFSIX: 90 % RH) to
avoid the negative-loading region of the fits.

**Custom sorbent.** The preset menu ends with *Custom*, which exposes
every parameter of eqns (6)–(7) as editable fields, seeded from the
last-selected preset (Reset restores the seed; a Reseed button pulls a
different preset's numbers). Water handling is a toggle — GAB (eqn 7)
or inert — and the CO₂–H₂O coupling offers three modes: *none*,
the *modulated* form of eqn (6) (ψ, β — positive values co-operative,
negative competitive, with the same x_max RH clamp mechanism the
presets use), or *site competition*, an extended Toth in which the
water partial pressure joins the shared denominator (eqn 8):

$$q^*_{CO2} = \frac{q_s\, b\, p_{CO2}}
  {\left[1 + \left(b\, p_{CO2} + b_w\, p_{H2O}\right)^{\omega}\right]^{1/\omega}},
  \qquad b_w = b_{w0}\,
  e^{\frac{\Delta H_{H2O}}{R T_0}\left(\frac{T_0}{T} - 1\right)} \tag{8}$$

Because p_H₂O is a gas-phase quantity, site competition also works with
water set to inert. −ΔH(CO₂) in the property boxes remains the Toth ΔH₀
(unless unlinked), and −ΔH(H₂O) doubles as the b_w Arrhenius heat in
competition mode. Nonphysical inputs (GAB K ≥ 1, Toth exponent ≤ 0 at
T_des, ψ·q_w near 1) raise red warning badges rather than being
silently clamped. Custom parameters are yours, not a published fit —
check the isotherm panel against your data before trusting the KPIs.

### Energy (Simplified-DAC forms)

- **Thermal** = jacket duty ∫U_a(T_W − T)⁺dt over the heating step
  (+ wall mass sensible). In the instant-kinetics limit this collapses
  to E_sensible + Δq·ΔH terms of `dac.energyModel` (see validation).
- **Fan (electrical)** — isentropic blower over ΔP (eqn 9):

$$\dot W_{fan} = \frac{\gamma}{\gamma-1}\frac{P_{amb}}{\eta_{fan}}
\left[\Big(\tfrac{P_{amb}+\Delta P}{P_{amb}}\Big)^{(\gamma-1)/\gamma}-1\right]\dot V_I \tag{9}$$

- **Vacuum pump (electrical)** — isothermal (default, Simplified-model
  convention, eqn 10) or adiabatic (Glaser's form, eqn 11), applied to
  the actual outflow at the instantaneous column pressure:

$$\dot W_{vac}^{iso} = \frac{\dot n_O R T_{amb}}{\eta_{vac}}\ln\frac{P_{amb}}{P},\qquad
\dot W_{vac}^{adia} = \frac{\dot n_O R T}{\eta_{vac}}\frac{\gamma}{\gamma-1}
\left[\Big(\tfrac{P_{amb}}{P}\Big)^{(\gamma-1)/\gamma}-1\right] \tag{10,11}$$

### KPIs (at CSS)

Purity (dry basis, product step only), recovery (captured / CO₂ fed over
the whole cycle), specific energies (MJ/kg CO₂ **and** kWh/tonne CO₂
— ×1000/3.6 — thermal | fan | vacuum), productivity per m³ of
**contactor** and per m³ of **sorbent** (contactor × 1/(1−ε_b)), working
capacity Δq_CO₂ (loading at the end of adsorption minus at the end of
heating) with Δq_H₂O and the bed temperature swing, water co-desorbed
per CO₂, cycles to CSS.
CSS = start-of-cycle state repeats to 10⁻⁵ (relative).

All five plots carry a live crosshair: hover to read every curve at that
abscissa (the profile plots also name the active step).

### Exporting figures, data and reports

Every figure carries **⬇ PNG** and **⬇ CSV** buttons. The CSV of a plotted
figure is derived from the *same declarative spec that was drawn*, so an
export cannot disagree with the picture; markers (the A and B cycle end
states) are written into the file header as comments. The sweep map
exports its grid in long format, one row per grid point, with every KPI
as a column and a blank row where a point failed to converge. Buttons are
disabled until the corresponding result exists.

Exports **force the light palette on an opaque white background**
regardless of the screen theme, so a figure saved from a dark-mode
machine is still publication-ready. (Canvases are cleared transparent and
the plotting code reads the CSS variables, so a naive `toDataURL` in dark
mode yields pale ink on a transparent field — invisible on white paper.)

**Generate report** (top right, in the same place in every app of the
suite that has it) produces a self-contained HTML document: introduction,
objectives, the full equation set typeset in MathML with a nomenclature
table, the numerical method, every input parameter with its units, the
KPI tables, and each selected figure with a numbered caption. Sections and
figures are individually selectable and an optional case-notes field
appears as an abstract. Equations and figures are cross-referenced by name
and numbered at build time, so the numbering cannot collide or drift. The
file embeds its figures as data URIs and makes no network requests; print
it to PDF for a paper-style document. Figures that have not been computed
in the session are named in place rather than silently omitted.

### Simplified-model comparison

The three energy cards carry a small second line: the same quantity from
the group's Simplified DAC energy model (`dac.energyModel`), evaluated
with **this cycle's own** Δq, ΔT and capture fraction η_cap:

$$E_{th}^{simp}=\frac{\big(c_{p,s}+m_{str}c_{p,str}+\Delta q_{CO_2}c_{p,a,CO_2}
+\Delta q_{H_2O}c_{p,a,H_2O}\big)\Delta T+\Delta q_{CO_2}\Delta H_{CO_2}
+\Delta q_{H_2O}\Delta H_{H_2O}}{\Delta q_{CO_2}M_{CO_2}} \tag{12}$$

$$E_{fan}^{simp}=\frac{\dot w_{air}}{C_{CO_2,air}\,\eta_{cap}M_{CO_2}},\qquad
E_{vac}^{simp}=\frac{RT_{amb}}{\eta_{vac}M_{CO_2}}\ln\frac{P_{amb}}{P_{regen}}
\left(1+\frac{\Delta q_{H_2O}}{\Delta q_{CO_2}}\right) \tag{13,14}$$

Feeding it the dynamic Δq/ΔT/η_cap makes the comparison a test of the
*energy forms alone*, so the gap is attributable to what the simplified
model omits: the void-gas inventory, sensible heat swept out with the
desorbing gas, and CO₂ released during evacuation instead of in the
product step. Across the ten anchor cases the thermal duty agrees to
0.2–8 % and the electrical to better than 0.5 % — except when the
adiabatic pump is selected, where the simplified model (isothermal by
construction) is ~8 % low, as it must be. A wall thermal mass, which has
no analogue in the simplified model, is added to its sensible term so a
non-zero wall does not by itself open a gap.

## Sorbents

| | q_s0 | b₀ [Pa⁻¹] | −ΔH₀ [kJ/mol] | ω₀ | α | T₀ [K] | GAB (Cm/C_G/K) | ψ, β [g/mmol] |
|---|---|---|---|---|---|---|---|---|
| TRI-PE-MCM-41 | 2.7 | 0.38776 | 90 | 0.3123 | 0.3472 | 298 | — (H₂O inert) | — |
| SI-AEATPMS | 1.0 | 2.0123 | 90 | 0.2696 | 0.08069 | 298 | — (H₂O inert) | — |
| NbOFFIVE-1-Ni | 2.22 | 0.17567 | 50 | 1.166 | −0.4937 | 273 | 6.159/36.81/0.4063 | 0, −0.11707 |
| SIFSIX-18-Ni-β | 3.0 | 0.04572 | 52 | 0.5455 | 1.592 | 273 | 0.4067/93.13/0.7859 | 0, −0.6772 |
| APDES-NFC-FD-S | 2.2 | 0.373 | 60 | 0.4247 | −0.4921 | 296 | 36.48/0.1489/0.5751 | 0.00958, 3.448 |
| Lewatit VP OC 1065 | 3.9387 | 0.20112 | 70.999 | 0.26906 | 0.045015 | 303.15 | 3.864/1.795/0.7683 | 0.0047, 0.569 |

### Editable properties

ε_b, ε_p, ρ_p, c_p,s, −ΔH(CO₂) and −ΔH(H₂O) are number boxes seeded from
the preset; switching sorbent (or pressing **Reset to preset**) restores
the published values, and any edited field is outlined and badged so a
KPI is never silently attributed to the published material. The
isotherm parameters of the six *presets* stay fixed; to change those,
switch to the **Custom** sorbent (see Equilibria above), which starts
from the current preset's numbers.

−ΔH(CO₂) is not a free parameter: it is also the Toth temperature
coefficient ΔH₀ in eqn (6), because the two are the same physical
quantity (Clausius–Clapeyron). Editing it therefore moves the isotherm
*and* the energy balance together. The **unlink** checkbox breaks that
tie for sensitivity testing — useful, but thermodynamically inconsistent,
so the app badges it in red while it is on. For the base case, raising
−ΔH(CO₂) from 71 to 80 kJ/mol gives:

| | recovery | E thermal |
|---|---|---|
| preset | 34.5 % | 11.54 MJ/kg |
| linked (isotherm moves too) | 44.0 % | 9.61 MJ/kg |
| unlinked (energy balance only) | 34.5 % | 11.74 MJ/kg |

−ΔH(H₂O) has no equilibrium counterpart (the GAB fit is on relative
humidity), so it varies freely; it is disabled for the two sorbents whose
water is treated as inert.

First five from the group's five-sorbent DAC screening paper (SI Tables
S1–S5); Lewatit from Balasubramaniam, Forbes, Picard, Sawada &
Rajendran, *Adsorption* 32(5):23, 2026 (same parameter set as the
FullModelLDF project; b₀ converted from bar⁻¹). Water on TMCM-41 and
TPMS is treated as inert (no pure-H₂O data; the paper's base case).
Lewatit ΔH_H2O = 43.8 kJ/mol and c_p,s = 1580 J/kg/K are assumed
(documented in HANDOFF.md). ρ_p, ε_b, c_p,s per SI Table S5
(APDES ε_b = 0.092; Lewatit ε_b = 0.37, ε_p = 0.35).

## Numerics

Fully implicit TR-BDF2 (Newton + finite-difference Jacobian, adaptive
steps, error-controlled), event location for the evacuation end
pressure, trapezoidal quadrature of flows/powers on the accepted steps.
Simulation runs in a Web Worker; a typical cycle takes ~50–100 ms and
CSS a few cycles, so the app responds near-interactively. The
parametric sweep runs every grid point to CSS.

## Validation (see `validation/README.md`)

- 373 cross-checks against an independent Python/scipy (BDF, rtol 1e-10)
  reference over 14 anchor cases spanning all six sorbents, both pump
  models, RH clamps, wall mass, structure mass, edited sorbent
  properties (linked and unlinked ΔH), Arrhenius LDF temperature
  dependence (wet and dry) and the custom site-competition isotherm
  (GAB and inert water): **worst deviation 6.4·10⁻⁶**. Custom seeded
  with a preset's exact numbers is checked to reproduce that preset
  bit-for-bit. The simplified-model comparison values (Δq, ΔT, η_cap and
  the three energies) are checked in both implementations too.
- Mole balances close to ~10⁻¹³ (Python, exact quadrature) and
  ≤ 2·10⁻⁵ (JS, trapezoid).
- Instant-kinetics limit collapses onto MATLAB `dac.energyModel`
  term by term: blower exact (7·10⁻¹⁶), all Δq-based terms ≤ 0.43 %
  with the residual fully attributed to the void-gas inventory that the
  simplified model deliberately neglects.

## Known modelling limits

- A 0-D CSTR under-predicts adsorption-step capture (no plug-flow
  front), so recovery is conservative vs 1-D models; purity is likewise
  diluted by the mixed residual N₂ at the start of desorption.
- Repressurisation is instantaneous and isothermal; its ~seconds-scale
  dynamics and compression heating are neglected.
- The pump is assumed able to hold P_regen during heating regardless of
  the desorption rate; steam-assisted regeneration is not modelled.
