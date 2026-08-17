# Isotherm Fitter validation assets

Cross-validation pair for the engine embedded in `../IsothermFitter.html`.
The Python reference uses scipy.optimize.least_squares (trf, native
bounds); the JS engine uses its own Levenberg–Marquardt on logit-bounded
parameters with multi-start. The seeded noisy datasets are GENERATED AND
PRINTED by the Python run and pasted verbatim into the jsc harness, so
both engines fit byte-identical data.

| File | Run |
|---|---|
| `fitter_ref.py` | `python3 fitter_ref.py` (needs scipy) |
| `fitter_engine.js` | snapshot of the embedded engine — keep in sync |
| `fitter_test.js` | `jsc fitter_engine.js fitter_test.js` (macOS JavaScriptCore) |

Both must end `ALL PASS`.

## Recorded anchors (2026-08-11)

- Noisy Langmuir (seed 7, 3% rel noise): scipy qs=3.97843437, b=8.10611476,
  SSE_w=1.0186424078e-2, AICc=−138.34; JS within 2.2e-8 / 2.5e-10.
- AICc ordering: Langmuir −138.3 beats DSL −135.2 on Langmuir data in both
  engines (parameter penalty working).
- DSL overfit valley: JS (30 starts) found SSE_w=8.394e-3, deeper than
  scipy's 8.712e-3.
- Toth (seed 11, 5%): rel-weight optimum (4.1614, 8.2519, 0.42282) vs
  abs-weight (4.2116, 10.1951, 0.40533) — both matched to ≤6e-7.
- Global van 't Hoff DSL (3 T, seed 3): truth (qs 3/4, U 30/12 kJ/mol);
  both engines recover within 4% with identical SSE_w to 1.4e-9;
  q_st(q) runs 29.6 → 12.9 kJ/mol across loading.
- Custom expression `qs*b*p/(1+b*p)` recovers (4, 8) to 2.9e-11.
