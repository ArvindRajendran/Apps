# UptakeCurves validation assets

Independent cross-validation pair for the engine embedded in
`../UptakeCurves.html`. Two implementations of the same semi-discrete
system (conservative FV + TR-BDF2), written in different languages with
different linear algebra, checked against exact solutions and against each
other.

| File | What it is |
|---|---|
| `uptake_ref.py` | Independent Python reference (numpy only). Checks A–F: exact sphere series, effective-D series, LDF, bidisperse limits, Darken/Γ table vs exact Langmuir, ads/des anchors, bidisperse cross-check anchors. |
| `uptake_engine.js` | Snapshot of the engine exactly as shipped inside the app's `<script id="uc-engine">` block (kept in sync manually — if the app engine changes, re-copy it here). |
| `uptake_test.js` | jsc test harness mirroring the Python checks, with the Python anchor numbers baked in. |

## Run

Python reference (~3 min, prints PASS/FAIL per check):

```bash
python3 uptake_ref.py
```

JS engine (macOS JavaScriptCore, no node needed; ~17 s):

```bash
/System/Library/Frameworks/JavaScriptCore.framework/Versions/Current/Helpers/jsc uptake_engine.js uptake_test.js
```

Both must end with `ALL PASS`.

## Recorded anchor values (2026-08-05)

- Sphere series U(τ = 0.01/0.05/0.1/0.2/0.4) = 0.308514 / 0.606940 /
  0.770479 / 0.915496 / 0.988269; both engines reproduce to ≤ 3.6·10⁻⁵.
- Darken Langmuir micro-only (q_s 5, b 10, y 0↔0.5, D_μ0/r_c² = 10⁻³ s⁻¹):
  ads t₅₀/t₉₀ = 11.24/59.69 s, des = 19.97/140.85 s (Python); JS within
  2.3·10⁻⁴ relative.
- Bidisperse linear anchors U(5, 20, 60, 150, 400 s) = 0.599232 /
  0.912529 / 0.998233 / 1 / 1 (JS within 6·10⁻⁴ — residual
  time-integration difference, both dt-converged).
- Bidisperse Langmuir + Darken anchors U(5, 20, 60, 150 s) = 0.841206 /
  0.999919 / 1 / 1 (JS within 1.3·10⁻⁵).
