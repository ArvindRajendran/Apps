# Psychrometric Chart

Single-file interactive psychrometric (moist-air) chart:
`PsychrometricChart.html`. Works offline; share by sending the file.

## Features

- **Axes**: dry-bulb temperature (°C) on x; humidity on the right y-axis,
  switchable between humidity ratio (g water / kg dry air) and water mole
  fraction (mol/mol). Default window 0–50 °C, 0–35 g/kg; temperature
  range, humidity limit, and barometric pressure (default 101.325 kPa)
  are all editable.
- **Property-line families** (each toggleable, colour-coded):
  relative humidity (green), thermodynamic wet-bulb temperature (blue),
  specific volume (purple), specific enthalpy (orange). The saturation
  curve (100% RH) is always drawn; the supersaturated "fog" region above
  it is shaded.
- **Crosshair readout**: hover anywhere for T dry-bulb, W, mole fraction,
  RH, wet-bulb, dew point, enthalpy, specific volume, and vapour pressure.
- **Anchor points A, B, C, D**: click to place in turn, drag to move,
  or type T + RH. Shows a property table with one column per placed
  point, plus a process path chaining the points in order with ΔT, ΔW,
  Δh and sensible heat ratio for every segment (A→B, B→C, C→D).
- **Export**: copy point data (tab-separated, pastes into Excel),
  download the chart as PNG, or download a fresh standalone copy of the
  app itself.
- Settings and anchors persist in the browser (localStorage). Light and
  dark themes follow the system.

## Formulations

ASHRAE Fundamentals-style ideal-gas moist-air relations, per kg dry air:

1. Saturation pressure: Hyland–Wexler correlations (separate
   coefficients over water, t ≥ 0 °C, and over ice, t < 0 °C).
2. Humidity ratio: W = 0.621945·p_w/(P − p_w)
3. Enthalpy: h = 1.006·t + W·(2501 + 1.86·t) kJ/kg d.a.
4. Specific volume: v = 0.287042·(t + 273.15)·(1 + 1.607858·W)/P
5. Wet bulb: implicit ASHRAE relation solved by bisection (ice-phase
   variant below 0 °C); dew point by inverting the saturation curve.
