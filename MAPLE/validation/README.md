# MAPLE app — validation

The app is a port of trained MATLAB networks, so validation means:
(A) the port is bit-faithful, (B) it reproduces the published accuracy
claim against the detailed-model labels, (D) it reproduces published
numbers from the papers themselves. Section C covers app-side sanity.

## Files

- `export_maple_weights.m` — run in a clone of
  github.com/ArvindRajendran/MAPLE:
  `matlab -batch "run('export_maple_weights.m')"` → `export/maple_weights.json`
  + a MATLAB-side self-check (hand-rolled forward pass vs the network
  objects, ≤ 9·10⁻¹⁵).
- `make_maple_anchors.m` — same place; writes `maple_anchors.json`
  (200 LHS rows + 12 preset×condition anchors per cycle, MATLAB truth),
  the full label exports, and prints the reference R² values.
- `net_role_grid.m` — the diagnostic that established the net roles
  (R² of every net vs every candidate label column).
- `maple_weights.js / maple_anchors.js / maple_labels_*.js /
  maple_paper_pts.js` — the same data as JS literals for jsc.
- `maple_engine.js` — snapshot of the `maple-engine` block in
  `MAPLE.html`. **Keep in sync with the app.**
- `maple_test.js` — harness, 2,189 checks:
  `jsc maple_weights.js maple_engine.js maple_anchors.js maple_labels_LPP.js maple_labels_FP.js maple_paper_pts.js maple_test.js`

## A. Port fidelity (JS vs MATLAB nets)

200 labelled LHS rows + 12 named preset anchors per cycle, all five
outputs each: **worst relative deviation ≤ 10⁻¹²** (weights serialized
at %.15g; the forward pass is the same arithmetic).

## B. Published accuracy claim (vs detailed-model labels)

All ~48k labelled rows per cycle, JS predictions vs the detailed-model
KPIs, R² identical to MATLAB's to 10⁻⁶ and above the papers' ≥ 0.99:

| KPI | LPP | FP |
|---|---|---|
| purity | 0.99932 | 0.99974 |
| recovery | 0.99841 | 0.99890 |
| energy | 0.99609 | 0.99849 |
| productivity | 0.99710 | 0.99410 |
| energy (η=100%) | 0.99605 | 0.99521 |

## C. App-side sanity

Guard box (in/out cases, ratio + ordering constraints), SSL Henry-limit
and binary-reduction identities, b(298 K) identity, NSGA-II smoke test
(front non-dominated, reaches the high-purity/high-recovery corner).
Since v1.1: the stepwise GA API (`nsga2init/step/front`, used by the UI
for live per-generation drawing) is checked to reproduce the monolithic
`nsga2` **exactly** (same seed → identical front, deviation 0); the
energy–productivity mode yields a non-empty, constraint-feasible front
with sane energies; and `optProblem` without an explicit constraint
object reproduces the default 95 %/90 % constraint values bit-for-bit.
Since v1.2 the constraint penalty is linear + quadratic,
P = Σ 20 v + 5000 v² — a regression check verifies that a VSA
energy–productivity run delivers *strictly feasible* front points
(with the previous pure-quadratic penalty the GA parked the whole
front at Pu ≈ 94.9 %, just outside the constraint, and the results
table rejected every sorbent as "infeasible").

## D. External anchor — published optimum points

`MAPLE2Data(PuReAdded).xlsx` in the repo carries the Fig-3 curve data of
the Limits paper: 216 optimum points (curves E5/P5 = IISERP-MOF2,
E6/P6 = UTSA-16, E7/P7 = Zeolite 13X, 4-step LPP) whose KPIs come from
the **detailed dynamic model** re-run at the MAPLE-Opt conditions.
Evaluating the JS engine at those exact inputs:

| KPI | median | p90 | max |
|---|---|---|---|
| purity | 0.00 % | 0.11 % | 0.40 % |
| recovery | 0.02 % | 0.96 % | 3.26 % |
| energy | 1.80 % | 6.08 % | 18.9 % |
| productivity | 0.93 % | 4.33 % | 7.8 % |

This is the surrogate-vs-detailed-model gap at published conditions —
i.e. the app reproduces the papers' own accuracy, from the papers' own
numbers. (Purity is near-exact because the ≥ 95 % constraint pins it.)

## Findings about the public repo (documented, not "fixed")

1. **The LHS-samples xlsx has swapped Pu/Re column headers.** Its
   "Pu [%]" column is recovery and "Re [%]" is purity. Evidence: at the
   SI optimum points the purity constraint (≥ 95 %) pins the true
   purity, and it is net2 (the wrapper's `Pu`) that reproduces the
   pinned column, while matching the *other* LHS column. The `MAPLE.m`
   wrapper's role assignment (net1 = recovery, net2 = purity) is
   therefore correct; only the spreadsheet headers are wrong. All
   mappings in this folder use the corrected semantics.
2. **net6 matches nothing shipped** (R² < 0.3 against every label
   column, linear or log): an unused leftover in both .mat files. The
   app ports nets 1–5 only (net5 = log₁₀ energy at 100 % pump
   efficiency, used for the KPI card's second line).

## The absolute-pressure support lesson (v1.3)

The LHS training data sampled **absolute** pressures — P_I ∈ [0.07, 4.0]
and P_L ∈ [0.01, 1.0] bar, identical in both cycles — while the network
inputs are the **ratios** log₁₀(P_I/P_H), log₁₀(P_L/P_H). The global
ratio ranges (e.g. log₁₀(P_L/P_H) down to −2.695) are only the envelope
over P_H ∈ [1, 5]: the −2.695 corner exists solely at P_H = 5. Bounding
the optimizer by ratios alone therefore admitted joint data voids —
P_L ≈ 2 mbar at P_H = 1 has **zero** of the 48k training rows (the
sampled VSA floor is ~8–10 mbar) — and the GA loved those voids:

- the VSA Pu–Re front for UTSA-16 came out "flat" (Pu 99.4–100 across
  all recoveries), triple-pinned at the P_I, P_L and v_F bounds;
- the FP cycle appeared to beat LPP on energy (13X "148 kWh/t" at
  P_L = 2.8 mbar), inverting the physically expected cycle ranking.

Both were pure unsupported extrapolation *inside* the marginal guard
box. Since v1.3 the engine enforces the absolute windows (`PABS`) in
`guards()` and in the optimizer's validity check; harness checks
`guard.*-abs*` and `nsga.support-*` guard the regression. With the
floors in place the UTSA-16 VSA Pu–Re front is a genuine curve
(Pu 91→100 trade-off) and FP correctly costs *more* energy than LPP
for every preset. Moral: for surrogates trained on transformed inputs,
range-check the **sampled** variables, not just the transformed ones —
marginal boxes do not imply joint support.

## Fragmented Pu–Re fronts (v1.4)

Reported symptom: Pu–Re fronts looked smooth for Zeolite 13X and
Mg-MOF-74 but broken/gappy for UTSA-16 and IISERP-MOF2. Two causes,
both in this app rather than in the published networks:

1. **Rank-only tournament.** The selection operator compared rank only
   (an acknowledged shortcut). Once the whole population reaches rank 0
   — which happens early on a near-flat front — selection degenerates
   to a random walk with no pressure toward under-populated regions.
   v1.4 uses the standard crowded-comparison operator (rank, then
   crowding distance).
2. **Ratio-space DV sampling.** P_I and P_L were searched as
   log₁₀ ratios of P_H; the LHS sampled them as *absolute* pressures.
   v1.4 searches absolute log₁₀ P_I and log₁₀ P_L (with the ratio
   bounds applied as rejection constraints), matching the design.

Effect on the largest recovery gap in the front, LPP/VSA, y_F = 0.15:

| sorbent | before | after |
|---|---|---|
| IISERP-MOF2 | 16.6 pp | 2.3 pp |
| Zeolite 13X | 3.1 pp | 2.5 pp |
| UTSA-16 | 42.7 pp | 44.7 pp (isolated branch, see below) |

**Correction to an earlier diagnosis.** The fragmentation was at first
explained as the high-selectivity fronts being "below surrogate
resolution" — purity varying only ~0.01 pp across recovery 80–95 %.
That number was wrong: it came from a diagnostic that printed just the
first four points of the band. Measured over the full
band the purity span is 2.10 pp (UTSA-16), 2.29 pp (IISERP-MOF2) and
5.31 pp (13X) — smaller for the high-selectivity pair, but well above
the surrogate's ±0.15 pp purity accuracy. The fronts are resolvable;
the gaps were mostly algorithmic.

UTSA-16 retains one genuine discontinuity: an isolated high-purity
branch near Re ≈ 64 separated from the main branch. A 40,000-point
random DV sweep finds only 4 non-dominated points in that interval, so
the attainable set really is sparse there — not a solver artifact.

The `resnote` badge in the app still flags a front whose purity span
falls below ±0.15 pp; after v1.4 it rarely triggers for the shipped
presets, which is the intended outcome.

## Notes for future edits

- The input convention is: [q_sat·ρ, log₁₀ b_CO2(298 K),
  log₁₀ b_N2(298 K), ΔU_CO2, ΔU_N2, ρ, y_F, t_ADS, P_H,
  log₁₀(P_I/P_H), log₁₀(P_L/P_H), v_F], z-scored with the stored
  mue/sig. b(298) uses T = 298.0 exactly (their code), not 298.15.
- Optimizer results are surrogate-optima; the papers re-ran them
  through the detailed model. Expect a few per cent of optimism near
  sparse corners (e.g. v_F at its bound) — that is a property of the
  published networks, not of this port.
- If the upstream repo re-trains the networks, re-run both .m scripts
  and rebuild the .js wrappers (`var NAME = <json>;`).
