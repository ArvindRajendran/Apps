# Unit → SI Notation Converter

Single-file, offline web app that converts slash-style units to
negative-exponent SI notation, e.g. `mmol/g` → `mmol g⁻¹`.

## Use

Open `UnitSIConverter.html` in any browser (double-click). Type or paste a
unit expression; the converted form appears instantly with two outputs:

- **LaTeX** — e.g. `J mol$^{-1}$ K$^{-1}$` — with a copy button
  (`µ` → `$\mu$`, `°` → `$^{\circ}$`, `%` → `\%` handled automatically).
- **Rich text (Word/RTF)** — copies formatted text with real superscripts
  that paste correctly into Word, Google Docs, PowerPoint, email, etc.

## Understood input

- Slashes, including chained and parenthesized: `W/m^2/K`, `J/(mol K)`,
  `((J/mol)/K)`
- Exponents as `^2`, `^{-1}`, bare digits (`s2`, `g-1`), or unicode (`m³`)
- Multiplication by space, `·`, `*`, or `.`
- Repeated units are combined: `m/s/s` → `m s⁻²`
- Symbols may include `µ`, `°`, `Ω`, `%`, `Å`

The file is fully self-contained (~15 kB) — share it freely.
