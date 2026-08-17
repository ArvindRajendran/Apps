# `_shared/` — the figure-export and report module

`report_core.js` gives an app ⬇ PNG / ⬇ CSV buttons on every figure and a
**Generate report** button that writes a self-contained HTML document.

The apps are single-file by design, so the core is **pasted** into each
one between sentinel comments rather than linked:

```
/* == report-core v1 ==   …   /* == /report-core == */
```

| Script | Purpose |
|---|---|
| `sync_report_core.sh <file>…` | re-paste the canonical core into the named files |
| `check_report_core.sh` | diff every app's embedded copy against the canonical file |

Never edit an embedded copy. Edit `report_core.js`, run `sync`, run
`check`. For an app assembled by concatenation (DAC), sync the **source
part**, not the built file — the next build would overwrite it.

Carrying it: **DAC**, **BAAM**, **MAPLE**.

## Two behaviours that are the point of the module

**Exports force the light palette on opaque white.** Canvases across the
suite are cleared transparent and the plotting code reads CSS variables,
so a naive `toDataURL` on a dark-mode machine yields pale ink on
transparency — invisible on white paper. `withLight()` temporarily sets
every `prefers-color-scheme: dark` media rule to `not all`, calls the
host's `redraw`, captures onto an opaque white canvas, then restores and
redraws. This works because the apps define the light palette on a bare
`:root` and override it inside the media query — **keep that structure**.
A re-entrancy guard makes a whole report one flip, not one per figure.

**Equation and figure numbers are assigned at build time.** Prose writes
`{{eq:toth}}` and `{{fig:isoC}}`; the core numbers every equation in
document order and substitutes. Do not hand-number equations — that is
exactly the collision that had to be unpicked by hand in DAC's README.

## API

```js
Report.init({ app, version, filePrefix, redraw, meta, button })
Report.figure({ id, file, title, caption, ready, csv })
Report.spec(canvasId, spec)      // or (id, null) to clear
Report.mount()                   // after the DOM and figures exist
Report.refresh()                 // whenever results appear or go stale
Report.content({ … })            // the report prose, equations and tables
```

- `redraw` **must repaint every canvas** — it is what produces the
  light-palette figures. If it repaints only some, the rest are captured
  in the screen palette.
- `meta()` returns comment lines prepended to every CSV.
- `button: false` suppresses the header button (nothing uses this yet).
- `ready()` must be false whenever the canvas has never been drawn, not
  merely when the data is missing. A canvas inside a closed `<details>`
  has never been drawn; see the MAPLE note below.
- `csv()` may return `null` to **defer to the plotted spec** — that is how
  a canvas that is a line plot in one mode and a heatmap in another
  serves the right table in each.

### Where the buttons go

The **report button is mounted by the core**, not by the app's markup: it
finds `.wrap > h2`, wraps it in a flex `.apphd` row and appends the
button. Every app therefore gets it in the same place — top right, level
with the title — with zero per-app HTML. This is deliberate; do not
hand-place it.

**Per-figure buttons** depend on how many registered figures a card
holds. One figure → the buttons go in that card's `.phead`, beside the
title. Several → each canvas gets its own right-aligned `.figrow`
underneath. MAPLE stacks eight canvases in one card, which is why this
exists; do not assume one header per card.

## Adding the core to an app

1. Paste the sentinel block into its own `<script>` tag, **before** the
   app's main script and **outside** any block that is re-read as text.
   (DAC and BAAM build a Web Worker from a script tag's `textContent`;
   the core touches the DOM and must not land in there.)
2. If the app has a declarative plotter, end it with
   `Report.spec(cv.id, spec)` and CSV comes free for every figure drawn
   through it. Otherwise write a `csv()` per figure returning
   `{head, rows, notes}`.
3. **Clear the spec** on any canvas drawn by something other than the
   plotter — otherwise its CSV exports whatever was last plotted there
   under the new figure's name. This is a real bug that bit twice: DAC's
   water isotherm when the sorbent's water is inert, and MAPLE's sweep
   canvas when it switches to a heatmap.
4. `Report.init(...)`, then `Report.figure(...)` per figure, then
   `Report.content(...)`, then `Report.mount()` once the DOM exists.
5. Call `Report.refresh()` wherever results appear or become stale —
   after the main draw, and in the completion callback of anything
   asynchronous.
6. Write the report prose from the app's own README; those equations and
   validation statements are already vetted, and lifting them keeps the
   two from drifting.

## Verifying an integration

Beyond "it looks right", these are the checks that have actually caught
things:

- Every KPI in the report equals the on-screen card, digit for digit.
- The report's own tables close whatever balance the app reports.
- Generated report contains **no** `http` references (data URIs only).
- Equation/figure cross-references all resolve — no `(?)` or `Figure ?`.
- Uncomputed figures are named in place, never embedded blank.
- With a dark-mode screen, an exported PNG has a white corner pixel and
  dark ink; the app's own theme is restored afterwards.
- The report body does not scroll horizontally at a 420 px viewport.
- The app's numerical harness is unchanged — this work touches no engine.
- `check_report_core.sh` reports every copy identical.

A browser-pane trap while doing this: console history persists across
navigations in a tab, so an error thrown by your own console snippet
keeps reappearing and reads as an app bug. Confirm in a fresh tab.

## Per-app notes

**DAC** — spec-driven; `isoPlots` clears the water-isotherm spec in the
inert branch. Sync `dac_ui.html`, then rebuild `DAC.html`.

**BAAM** — `E` is already the engine handle, so its MathML builders are
named `ML`. Its four CSV builders were kept (they export more than the
plotted series) but converted to return tables instead of writing files.

**MAPLE** — the five decision-variable trends render lazily on a
`<details>` toggle, so their `ready()` checks `$("mp-dvwrap").open`;
without it the report embeds five blank frames. The report clamps purity
and recovery to 100 % exactly as the KPI cards do, since the surrogate
can return ~101.8 %.

**Not yet carrying it**: PSA Simulator (descoped), SMB Triangle and the
four Equilibrium Chromatography apps. All of those draw bespoke, so each
figure needs its own `csv()`; the four chromatography apps share most of
their theory prose.
