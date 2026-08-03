# Unit Converter

`UnitConverter.html` — plain-text unit conversion with live output. Type
`2.7 atm in Pa` and the answer appears as you type; no buttons, no
"convert" keyword needed.

## How it works

Every unit is stored as a scale factor to SI plus a vector of base-SI
exponents (m, kg, s, K, mol, A, cd). A conversion parses both sides into
(factor, vector) pairs, checks the vectors match, and divides the
factors — so compound units (`W/m2K`, `kJ/mol`, `Pa.s`,
`BTU/h.ft2.F`), SI prefixes on any metric unit, and exponents
(`m2`, `m^2`, `s-1`, `m²`) all work with no special cases. The physical
quantity shown in the badge is recognized from the exponent vector
itself (e.g. [−1, 1, −2] → pressure).

### Input grammar

- `<number> <units> in|to|→|=  <units>[, <units>…]` — several targets
  separated by commas.
- The number can be an arithmetic expression: `10*10*10 Pa in atm`,
  `(1+2)/4 bar in kPa`, `2^3 atm`, `2.7e5`, `101,325`. Evaluated with
  a built-in parser (+ − * / ^ and parentheses, standard precedence) —
  no eval(), fully offline.
- Juxtaposed shorthand is understood: `W/m2K` = W m⁻² K⁻¹ (longest-
  valid-prefix splitting), successive `/` denominate (`J/mol/K` =
  `J/molK` = J mol⁻¹ K⁻¹).
- No target yet → the app names the quantity and suggests common
  targets (click a chip or press ↓).
- **Temperatures** (°C/C, °F/F, K, °R) convert absolutely when alone
  (`77 F in K` → 298.15) and as differences inside compound units
  (`J/molK`, `W/m2K`), which is the standard convention.
- **Gauge pressure**: `psig` carries the +1 atm offset when converted
  standalone (`5 psig in kPa` → 135.799).
- **mol ↔ mass crossing**: append `(MW 44)` (g/mol) to convert e.g.
  `mmol/g → mg/g` or `kWh/t → J/mol`; the result is labelled with the
  MW used. Without it, the app asks for one.
- Dimension mismatches are reported with both dimensions spelled out.

### Conventions

- mmHg and torr are both 101 325/760 Pa (they differ by 2×10⁻⁷
  relative; the torr definition is used for both).
- cal is thermochemical (4.184 J exactly); BTU is the international
  table value (1055.055 852 62 J); gal is US (3.785 411 784 L);
  yr = 365.25 d; M = mol/L (prefixable: mM, µM).
- ppm/ppb/% are dimensionless factors.

## Validation

70-case automated test suite (run headlessly against the same engine
script the app ships with) — all pass, including:

| Input | Result | Check |
|---|---|---|
| `2.7 atm in Pa` | 273 577.5 | exact (atm := 101 325 Pa) |
| `77 F in K` | 298.15 | exact offset chain |
| `-40 C in F` | −40 | crossover point |
| `1 W/m2K in BTU/h.ft2.F` | 0.176 110 | vs. published factor 0.176 110 2 |
| `5 psig in kPa` | 135.799 | gauge offset |
| `3.2 mmol/g in mg/g (MW 44)` | 140.8 | MW bridge |
| `10 kWh/t in J/mol (MW 44)` | 1584 | MW bridge, compound |
| `8.314 J/molK in cal/molK` | 1.987 | juxtaposed shorthand |
| `1 atm in mmHg` | 760 | exact by definition |
| `1 m2/s in cSt` | 10⁶ | prefix on non-SI unit |
| `10*10*10 Pa in atm` | 9.869 × 10⁻³ | arithmetic input |
| `(1+2)/4 bar in kPa` | 75 | parentheses + precedence |
| `1/2 in in cm` | 1.27 | fraction vs. "in" separator |

Settings (significant figures) and the recent-conversion history
persist in localStorage. Works offline — single file, share freely.
