# Gas Constant (R)

Single-file reference and converter for the universal gas constant:
`GasConstant.html`. Works offline; share by sending the file.

## Features

- **Standard values list** — ~15 common forms (J, kJ, kPa·L, bar·L,
  atm·L, Torr·L, bar·cm³, atm·cm³, cal, J/kmol, psi·ft³/lb-mol·°R,
  atm·ft³/lb-mol·°R, BTU/lb-mol·°R, ft·lbf/lb-mol·°R, …), each with a
  copy button that copies **the number only**.
- **Build your own** — two modes:
  - **P·V/(n·T)**: pick pressure (Pa, kPa, MPa, bar, atm, Torr, psi),
    volume (m³, L, cm³, ft³, in³), amount (mol, kmol, lb-mol), and
    temperature (K, °R).
  - **Energy/(n·T)**: pick energy (J, kJ, cal, kcal, BTU, ft·lbf, erg),
    amount, and temperature.
- **Significant figures** selector (4 / 6 / 9 / 15) applies to the list,
  the builder, and what gets copied.

## Provenance

Nothing is hardcoded: every value is computed at runtime from the exact
SI definition R = N_A·k_B = 8.31446261815324 J mol⁻¹ K⁻¹ (2019 SI) and
exact conversion factors (atm = 101325 Pa; Torr = 101325/760 Pa;
psi = 4.4482216152605 N / (0.0254 m)²; lb-mol = 453.59237 mol;
°R = 5/9 K; thermochemical cal = 4.184 J; IT BTU = 1055.05585262 J;
ft·lbf = 0.3048 m × 4.4482216152605 N; ft³ = 0.3048³ m³).

Temperature offers only K and °R because a gas constant requires an
absolute temperature scale (°C and °F intervals match K and °R, but the
constant is only meaningful with absolute units).
