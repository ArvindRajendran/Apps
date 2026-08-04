# Research Apps

A collection of single-file, browser-based tools for chemical
engineering research and teaching, built by **Arvind Rajendran**
(University of Alberta) with the assistance of Claude Code (Anthropic).

Every app is one standalone HTML file: no installation, no internet
connection needed, no data uploaded or collected — all calculations run
entirely in your browser. Share an app by sending the one file, or
share the whole collection by sharing this folder.

**Open `index.html`** (or visit the deployed site) for the launcher —
hover an icon for a description, click to open. (`Apps.html` still
works — it's a redirect stub kept for old bookmarks/links.)

## The apps

### General purpose

| App | Description |
|---|---|
| 📚 **Journal Abbreviator** | Full journal name → ISO 4 abbreviation. 40,000+ embedded titles, rule-based fallback for unlisted journals, custom add/edit. |
| ⚖️ **Unit → SI Notation** | `mmol/g` → mmol g⁻¹ — negative-exponent SI style with LaTeX and rich-text (Word) copy buttons. |
| 🔢 **Text Counter** | Live word/character/sentence/paragraph counts (Microsoft Word conventions), including for a selected portion of the text. |
| 🔄 **Unit Converter** | Type a conversion in plain text — `2.7 atm in Pa` — and the answer appears as you type. Dimensional analysis recognizes the physical quantity; compound units, SI prefixes, multiple targets, mol↔mass crossing via molar mass. |

### Chemical engineering

| App | Description |
|---|---|
| 🌡️ **Psychrometric Chart** | Interactive moist-air chart: crosshair readout, RH / wet-bulb / specific-volume / enthalpy lines, adjustable pressure and ranges, A→D state points with process analysis, PNG export. |
| ♨️ **P–h Diagram** | Real-fluid pressure–enthalpy charts for CO₂, water, and ethane (CoolProp-precomputed reference EOS data): isotherms, isentropes, isochores, quality lines, state points, and a refrigeration-cycle builder with COP. |
| 🎈 **Gas Constant (R)** | Standard values of the universal gas constant with one-click copy, plus a builder for any P·V/(n·T) or energy/(n·T) unit combination. |
| ⚗️ **McCabe–Thiele Designer** | Binary distillation design at constant relative volatility: live sliders for α, compositions, feed quality, and reflux; stage-by-stage animation; Rmin and Nmin limits; Murphree efficiency. |

### Adsorption & chromatography

| App | Description |
|---|---|
| 🌊 **Equilibrium Chromatography** | Single-component equilibrium theory (Rhee–Aris–Amundson): envelope construction on any isotherm → waves and shocks, colour-mapped physical plane, column profiles and breakthrough curves. |
| 🔭 **Equilibrium Chromatography — Pulse** | Finite pulse injection by wavefront tracking: interacting fans, decaying shocks, collision events re-resolved by envelope construction, peak metrics and mass-balance check. |
| ⚛️ **Equilibrium Chromatography — 2×2 Systems** | Riemann problems for coupled conservation laws: isothermal binary chromatography with hodograph-plane analysis, exact constant-selectivity Langmuir construction plus a numerical reference solution for arbitrary isotherms. |
| 🌬️ **Equilibrium Chromatography — Sorption Effect** | Single adsorbable component whose adsorption/desorption changes the fluid velocity (isobaric, fixed inlet velocity): analytical Riemann solution via envelopes on νq vs y/(1−y), twin colour-mapped physical planes for composition and velocity with a linked crosshair. |
| 🔬 **SMB Triangle Designer** | Simulated moving bed design by triangle theory: complete-separation regions for linear and generalized-Langmuir isotherms, reduced-purity pentagon (linear), clickable (m₂,m₃)/(m₁,m₄) planes with live feed-concentration morphing, and open-loop flow-rate design from a pressure-drop limit. |
| 🏭 **PSA Simulator** | 1D non-isothermal pressure/vacuum swing adsorption cycle simulator (LDF kinetics, Ergun momentum, WENO finite volumes, stiff integration in a Web Worker): cycle builder, run to cyclic steady state, purity/recovery/energy KPIs, JSON/CSV round-trip. |

More apps are added over time.

## Design principles

- **One file per app** — inline CSS and JS, no external requests, works
  from a double-click on any machine.
- **Validation first** — every solver is checked against analytical
  results or an independent reference implementation before release;
  the validation cases are documented in each app folder's README.
- **Private by construction** — nothing leaves your browser; settings
  persist only in your browser's local storage.

## Structure

- `index.html` — the launcher (compact icon grid); this is what GitHub
  Pages (or any static host) serves at the site root automatically.
  `Apps.html` is a redirect stub to `index.html`, kept for old links.
- One folder per app, containing the app's single HTML file and a
  README documenting the model, data sources, and validation.

To deploy your own copy: push this folder to a GitHub repository and
enable GitHub Pages (Settings → Pages → deploy from branch). All links
are relative, so the site works at any URL.

## License & disclaimer

© 2026 Arvind Rajendran, licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (see
`LICENSE`) — share and adapt freely with attribution.

These tools are provided **"as is"**, for educational and research
purposes, without warranty of any kind. Verify results independently
before use in any critical application — **use at your own risk**. If
these tools are useful in your work, an acknowledgement is appreciated.
