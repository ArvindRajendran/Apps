# IAST validation assets

Cross-validation pair for the engine embedded in `../IAST.html`, using
deliberately different algorithms (Python: bisection on p₁⁰ with adaptive
Simpson ψ; JS: bisection+Newton on x₁ with analytic/Newton–Cotes-table ψ).

| File | What it is |
|---|---|
| `iast_ref.py` | Independent Python reference. Checks A–H incl. a **pyIAST cross-validation** (`pip install --user pyiast scipy pandas`; skips gracefully if absent). Run: `MPLCONFIGDIR=/tmp/mpl python3 iast_ref.py` |
| `iast_engine.js` | Snapshot of the engine exactly as shipped in the app's `<script id="ia-engine">` block (keep in sync). |
| `iast_test.js` | jsc harness with the Python anchors baked in. Run: `jsc iast_engine.js iast_test.js` (macOS JavaScriptCore; no node needed). |

Both suites must end `ALL PASS`.

## Recorded anchors (2026-08-11)

- Equal-q_s Langmuir pair ≡ extended Langmuir: 1.1e-13 (py) / 1.3e-15 (JS).
- pyIAST (Langmuir qs 3, b 20 vs qs 5, b 0.7; fitted params recovered
  exactly): max relative difference 1.654e-13.
- Henry-limit selectivity → H₁/H₂ = 15.4 exactly.
- DSL(3,50,4,1) vs Sips(5,2,0.85) anchors, (P, y₁) → q₁, q₂, x₁:
  (0.1,0.3) → 1.69484730, 0.41146444, 0.80465169;
  (1,0.5) → 3.86052802, 0.47496954, 0.89044636;
  (10,0.15) → 3.36197662, 2.02236558, 0.62439877;
  (100,0.85) → 6.89883293, 0.03895783, 0.99438469. JS within 1.6e-7 rel.
- Selectivity reversal for Langmuir (2,20)/(6,0.8): p* = 4.28054 bar;
  S(p*, y₁∈{0.2,0.5,0.8}) = 1 to 6.7e-16; no x=y crossing at P ≠ p*
  (binary IAST admits no composition azeotropes — S = p₂⁰/p₁⁰ depends on
  ψ alone).
- Degenerate pair (identical quad–quad): `findSReversal` returns
  `{ps: [], degenerate: true}` (not 400 bogus crossings), x₁ = y₁ to
  machine precision; genuine quad–quad crossing (b: 5,10 vs 8,3) found at
  p* = 3/7 bar exactly.
