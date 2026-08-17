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
| 💧 **Vapour Pressure Plotter** | Saturation curves for 653 components from published Antoine, Wagner, DIPPR-101 and PPDS coefficients plus corresponding-states estimates: overlay up to ten, click two temperatures to read every curve at both, switch to Cox or reduced coordinates, compare every correlation for one component, get α(T) for a pair and ΔHvap from the Clapeyron slope. |

### Adsorption & chromatography

#### Equilibria — fit and predict

| App | Description |
|---|---|
| 📐 **Isotherm Fitter** | Fit pure-component adsorption data to 16 isotherm models + custom expressions: paste (T, p, q) with unit conversion, Levenberg–Marquardt multi-start, AICc ranking with parameter-count penalty, ±95% CI and over-parametrization warnings, per-temperature or global van 't Hoff fits (ΔH, model-based isosteric heat), overlay/residual/parity plots, one-click parameter hand-off to the other apps; validated against scipy.optimize.least_squares. |
| 🥣 **IAST — Binary Mixtures** | Ideal Adsorbed Solution Theory (Myers & Prausnitz): binary mixture loadings from two pure isotherms — 11 models with analytic spreading pressures, guaranteed-unique monotone solver, consistency badges (Henry limit, monotonicity), live extrapolation warnings for the hypothetical pressures p_i⁰, ψ-construction and x–y diagrams, selectivity heatmaps + slices, selectivity-reversal detection, extended-Langmuir comparison, 3D surface, CSV export; cross-validated against pyIAST to machine precision. |

#### Kinetics — uptake and transport

| App | Description |
|---|---|
| ⏳ **Uptake Curves** | Gravimetric uptake predictor: step change at the surface of a bidisperse spherical adsorbent — micropore, macropore, macro+micro (Ruckenstein-type) or LDF transport, Darken correction Γ = dlnp/dlnq for any isotherm in the menu, Bosanquet Dₘ+D_K+τₚ macropore assembly, intracrystalline/intraparticle profiles, mₜ/m∞ with t ↔ √t axis toggle. |

#### Chromatography — equilibrium theory and SMB

| App | Description |
|---|---|
| 🌊 **Equilibrium Chromatography** | Single-component equilibrium theory (Rhee–Aris–Amundson): envelope construction on any isotherm → waves and shocks, colour-mapped physical plane, column profiles and breakthrough curves. |
| 🔭 **Equilibrium Chromatography — Pulse** | Finite pulse injection by wavefront tracking: interacting fans, decaying shocks, collision events re-resolved by envelope construction, peak metrics and mass-balance check. |
| ⚛️ **Equilibrium Chromatography — 2×2 Systems** | Riemann problems for coupled conservation laws: isothermal binary chromatography with hodograph-plane analysis, exact constant-selectivity Langmuir construction plus a numerical reference solution for arbitrary isotherms. |
| 🌬️ **Equilibrium Chromatography — Sorption Effect** | Single adsorbable component whose adsorption/desorption changes the fluid velocity (isobaric, fixed inlet velocity): analytical Riemann solution via envelopes on νq vs y/(1−y), twin colour-mapped physical planes for composition and velocity with a linked crosshair. |
| 🔬 **SMB Triangle Designer** | Simulated moving bed design by triangle theory: complete-separation regions for linear and generalized-Langmuir isotherms, reduced-purity pentagon (linear), clickable (m₂,m₃)/(m₁,m₄) planes with live feed-concentration morphing, and open-loop flow-rate design from a pressure-drop limit. |

#### Cycles — screening to full simulation

| App | Description |
|---|---|
| 🏷️ **BAAM — Cycle Screening** | Batch Adsorber Analogue Model (Balashankar 2019 + e-BAAM Liske 2026): 4-step VSA/PVSA CO₂/N₂ cycles with LPP or feed pressurization solved from one composition-path ODE + algebraic balances — cycle transitions drawn on the competitive isotherms with a live crosshair, three adjustable pressure levels (P_H up to 5 bar, with feed-compression work), design-space maps over (P_LOW, P_INT), purity–recovery Paretos with US-DOE classification (r_cut = 124.5) and full-model energy calibration (×1.58), extended dual-site Langmuir or IAST mixture equilibria, adsorbent presets from both papers, Isotherm Fitter hand-off; PNG/CSV on every figure and a **Generate report** button writing a self-contained, printable HTML document (model equations in MathML, all inputs, cycle results, design-space summary and classification verdict, selected figures with numbered captions); validated against the authors' original MATLAB (paper Table 2 reproduced exactly, converged agreement ≤10⁻⁵). |
| 🏭 **PSA Simulator** | 1D non-isothermal pressure/vacuum swing adsorption cycle simulator (LDF kinetics, Ergun momentum, WENO finite volumes, stiff integration in a Web Worker): cycle builder, run to cyclic steady state, purity/recovery/energy KPIs, JSON/CSV round-trip. |
| 🌁 **DAC — 0-D Cycle Simulator** | Dynamic 0-D (CSTR) simulation of a 4-step temperature–vacuum swing direct-air-capture cycle run to cyclic steady state: six sorbent presets (the group's five-sorbent DAC screening paper + Lewatit VP OC 1065) with Toth/GAB/modified-Toth binary CO₂–H₂O equilibria, LDF kinetics, first-order wall heating and volumetric bed–wall heat transfer, evacuation by pump speed (t_evac an outcome), Simplified-DAC energy forms (isentropic fan, isothermal/adiabatic vacuum pump), optional structure and wall thermal mass, user-editable sorbent properties (ε_b, ε_p, ρ_p, c_p,s, ΔH — with ΔH(CO₂) driving the Toth temperature coefficient unless explicitly unlinked), a fully custom sorbent (every Toth/GAB parameter editable, seeded from any preset, with modulated ψ/β or shared-denominator site-competition water coupling and GAB-vs-inert water toggle), optional Arrhenius temperature dependence of the LDF coefficients (per-component E_a, k referenced to 20 °C), purity/recovery/energy KPIs (MJ/kg and kWh/tonne) with productivity on both sorbent and contactor bases, live CSS cycle profiles with hover crosshairs, cycle-state markers on the isotherms, working capacity Δq, a side-by-side comparison with the Simplified DAC energy model on every energy card, and a full-CSS parametric sweep; PNG and CSV download on every figure (exported in the light palette on opaque white whatever the screen theme) and a **Generate report** button that writes a self-contained HTML document — introduction, objectives, the equation set in MathML with nomenclature, numerical method, every input parameter, KPI tables and the selected figures with numbered captions, printable to PDF; validated against an independent Python/scipy reference (373 checks over 14 anchor cases, ≤6·10⁻⁶) and term-by-term against MATLAB `dac.energyModel` in the instant-kinetics limit. |
| 🍁 **MAPLE — ANN PVSA Emulator** | The MAPLE framework (Pai, Prasad & Rajendran, IECR 2020 + ACS SCE 2021) as an offline app: the trained 12-input neural networks of github.com/ArvindRajendran/MAPLE ported to JavaScript, emulating a 4-step PVSA cycle (LPP or FP pressurization) for CO₂/N₂ capture at cyclic steady state — purity, recovery, energy (incl. the 100 %-pump-efficiency variant) and productivity in microseconds for any single-site-Langmuir adsorbent. Editable sorbent properties with data-driven trained-range guards (extrapolation badges), presets and a 36-material screening library from the papers' SI, Isotherm Fitter hand-off, pure-component isotherm comparison at 30 °C against known sorbents, instant 1-D sweeps and full 2-D KPI maps with DOE-feasible cells marked, and PNG/CSV on every figure and a **Generate report** button writing a self-contained, printable HTML document (surrogate and optimizer formulation, inputs with trained-range flags, predicted indicators, optimization outcome, selected figures), and in-browser NSGA-II Pareto optimization (purity–recovery, energy–recovery, productivity–recovery, and energy–productivity under editable purity/recovery constraints; narrowable decision-variable ranges, tunable population/generations/seed, live per-generation front animation with a hypervolume convergence plot, per-point tooltips and decision-variable trend plots; optional VSA mode) with multi-sorbent comparison tables. Validated: port exact to 10⁻¹² (2,120 anchor checks), R² 0.994–0.9997 against the ~96k detailed-model labels, and the 216 published optimum points of the Limits paper reproduced at the surrogate's own accuracy (purity 0.0 %, energy 1.8 % median). |

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
- `_shared/` — source of the few blocks that are genuinely shared.
  `report_core.js` is the figure-export and report module; because the
  apps are single-file by design it is *pasted* into each app between
  sentinel comments rather than linked. `sync_report_core.sh` re-pastes
  it, `check_report_core.sh` verifies every copy is byte-identical, and
  `_shared/README.md` documents the API, the integration checklist and
  the per-app gotchas. Carried by DAC, BAAM and MAPLE so far.

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
