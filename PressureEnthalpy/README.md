# Pressure–Enthalpy (P–h) Diagram

Single-file interactive P–h diagram for **CO₂, water, and ethane**:
`PressureEnthalpy.html` (~790 kB). Works offline; share by sending the file.

## Features

- **Axes**: specific enthalpy (kJ/kg) on x; pressure (atm) on a log-scale
  y-axis. Axis ranges are editable ("Reset axes" restores the fluid's
  full data range). Fluid chosen from a pull-down menu.
- **Two-phase dome** always drawn (bold), with the critical point marked.
- **Property-line families** (toggleable, colour-coded):
  isotherms (°C, orange), isentropes (kJ/kg·K, violet),
  isochores (kg/m³, green), and vapor-quality lines x = 0.1–0.9
  (cyan) inside the dome.
- **Crosshair readout**: P, T, h, s, ρ, v, quality x, and phase
  (liquid / vapor / two-phase / supercritical) at the cursor.
- **State points A–D**: click to place in turn, drag to move, or type
  P + T (single phase) / P + x (two-phase). Property table with one
  column per point, and a chained process path with Δh, Δs, ΔT, and
  pressure ratio per segment — handy for compression, throttling,
  and refrigeration-cycle work.
- **Cycle builder** — computes a full cycle and places its states as
  A–D with a closed process loop on the chart:
  - *Vapor-compression*: T_evap, T_cond, superheat, subcooling, and
    compressor isentropic efficiency η (h₂ = h₁ + (h₂ₛ − h₁)/η).
    Reports q_evap, w_comp, q_rej, **COP (cooling and heating)**,
    discharge T, pressure ratio, quality after the valve, and
    volumetric capacity (kJ/m³).
  - *Transcritical (gas cooler)*: for CO₂ heat rejection above the
    critical point — inputs are gas-cooler pressure and exit T.
  - *Rankine (power)*: boiler/condenser P, turbine-inlet T, turbine
    and pump efficiencies → turbine/pump work, thermal efficiency,
    turbine-exit quality, back-work ratio.
  Manually moving or editing a point switches the cycle off so the
  points become free again.
- **Export**: copy the point table (tab-separated) or download the
  chart as PNG.

## Data and accuracy

Property data are **precomputed with CoolProp** and embedded in the
HTML (no internet needed at use time):

- CO₂ — Span & Wagner (1996) reference EOS
- Water — IAPWS-95 (Wagner & Pruß)
- Ethane — Bücker & Wagner (2006)

Enthalpy/entropy reference states follow each EOS's published
convention (the NIST WebBook defaults), so values can be cross-checked
against webbook.nist.gov directly. The crosshair uses a 120×180
quantized grid (bilinear interpolation; ~0.01 °C in T, ~0.02 % in ρ),
and iso-curves are exact EOS polylines. Spot-checked against direct
CoolProp flashes for all three fluids in single-phase, two-phase, and
supercritical regions.

States below the CO₂ triple pressure (5.18 atm) where only solid +
vapor exist are blanked (the EOS does not cover solids).

Ranges: CO₂ 1–300 atm, −50…200 °C; water 0.01–300 atm, 1…700 °C;
ethane 0.1–100 atm, −180…200 °C.

## Regenerating / extending the data

`generate_data.py` rebuilds the embedded data (needs Python +
`pip install CoolProp`) and injects it into the HTML in place:

    python3 generate_data.py

To add a fluid, append an entry to `FLUIDS` in that script (any fluid
CoolProp supports — e.g. N₂, CH₄, NH₃, propane) and re-run; the app
picks it up automatically in the pull-down.
