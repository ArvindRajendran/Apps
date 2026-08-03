# PSA / VSA Cycle Simulator

Single-file 1D non-isothermal PSA/VSA cycle simulator:
`PSASimulator.html`. JavaScript port of the validated MATLAB simulator in
`Claude/Projects/PSA modelling/` (which remains the reference). Works
offline; the numerical engine runs in an inlined Web Worker so the UI
stays live during multi-minute CSS runs.

## Model

1D axial-dispersed plug flow, ideal gas, 3 components (all may adsorb),
LDF mass transfer, Ergun momentum, non-isothermal, non-isobaric:

1. ∂(c·y_i)/∂t = ∂/∂z(D_L·c·∂y_i/∂z) − ∂(v·c·y_i)/∂z − φ·ρ_p·∂q_i/∂t
2. ∂q_i/∂t = k_i(q_i* − q_i)
3. (ε·C_vg/R)·∂P/∂t + (1−ε)·ρ_p·(C_ps + C_pa·Σq)·∂T/∂t
   = K_z·∂²T/∂z² − ε·C_pg·∂(v·c·T)/∂z + (1−ε)·ρ_p·Σ(ΔH_i·∂q_i/∂t) − wall loss
4. Ergun: −∂P/∂z = 150μ(1−ε)²/(ε³d_p²)·u + 1.75(1−ε)ρ_g/(ε³d_p)·u|u|

Thermal options: isothermal (dT/dt = 0), lumped wall loss (2U/r_in), or
the full wall PDE (paper Eq. 10, T_w(z) state). Per cell, energy + the
differentiated EOS give a 2×2 solve for (dT/dt, dP/dt).

Isotherm: competitive dual-site Langmuir, concentration basis, with the
Wilkins & Rajendran (2019) Table 2 CO₂/N₂-on-13X parameters built in
(EES and UES N₂ sets) and all values editable.

## Numerics

- Finite volumes (default N = 15, editable), 2-candidate upwind-biased
  WENO face reconstruction of (y₁, y₂, T), Ergun face velocities from
  cell-centre ΔP, central dispersion/conduction.
- **TR-BDF2** (γ = 2−√2, L-stable) with modified Newton; colored
  finite-difference Jacobian (bandwidth-aware, ~3n_f RHS evaluations),
  dense LU with caching (Jacobian refreshed on Newton failure or every
  50 steps; factorization reused within a dt deadband). Error control
  via the divided-difference third-derivative estimate with the
  Hosea–Shampine stiff scaling; deadband step controller.
- Quadrature states (per-end species mole integrals + pump work)
  integrate with the same TR-BDF2 stage weights but stay outside the
  Newton system — stream accounting is discretization-consistent.
- Events (pressure triggers, pressure-floor guards) by sign check +
  bisection on the accepted step.
- Steps are BC pairs on a generic bed integrator: closed / inlet
  (Danckwerts) / outletP with exponential vacuum ramp
  P_bc = P_set + (P_start−P_set)e^(−λt), λ = 6/t_step.
- CSS: max|Δ(dimensionless states)| < tol for 3 consecutive cycles;
  per-component cyclic mass-balance closure reported (normalized per
  component, floored at 1e-5 of throughput).

## Validation against the MATLAB reference

The MATLAB energy balance had a bug (found during this port, fixed in
`core/bedRHS.m`): the convective enthalpy term applied C_pg twice
(ε·C_pg²·∂(vcT)/∂z). Both codes now use the corrected form; see
`PSA modelling/HTML_SIMULATOR_SPEC.md` §5.

Cross-checks (Wilkins operating point, wall model, cold start, N = 30):

- Isotherm anchors: pure CO₂ 1 bar 25 °C → 5.3594 mol/kg; pure N₂
  0.85 bar → 0.3174; competitive 15/85 → 3.8544 — identical to MATLAB
  to 4 decimals (paper measured 3.83–3.85).
- Single cold cycle, JS vs MATLAB (ode15s): every stream quadrature
  matches to 4–5 significant figures — e.g. feed CO₂ 4.8968/4.8968 mol,
  BLO waste N₂ 15.4274/15.4271, EVAC extract CO₂ 1.4010/1.4012,
  vacuum work 108.9/109.2 kJ; LPP pressure-event duration
  89.24 s vs 89.241 s; end-of-cycle bed T range [298.09, 332.27] K in
  both. Mass-balance closure ~3e-7 per component.
- ADS-step exotherm 298.2→354.6 K vs MATLAB 298.15→355.09 K, and the
  swing scales correctly with ΔH (the fixed physics).
- CSS targets for the default case (N = 15, MATLAB):
  wall model (300 cycles): purity 89.79 %, recovery 88.33 %,
  152.0 kWh/t, t_cycle 428.5 s; isothermal (100 cycles): 84.98 %,
  96.33 %, 188.9 kWh/t. A browser run-to-CSS takes a few hours at
  these settings — compare against these numbers. (MATLAB N = 30
  trends the same: ~91 % / 88 % / 152 after 200 cycles; the paper
  reports 96.8 / 90.8 / 213 with loading-dependent ΔH_CO2 and the
  exact recorded-raffinate LPP source.)

## Performance & tips

- N = 15 runs a cycle in roughly 30–60 s of wall time (browser,
  worker); N = 30 is ~4× slower. Run-to-CSS on the default case takes
  tens of minutes — start it and leave the tab open (progress and
  convergence are shown live; cancellable).
- The bed state is kept warm between runs: after a CSS run, further
  runs continue from the converged state (changing thermal model or
  isotherm set resets it).
- Export/import the full setup as JSON; results (KPIs + profiles) as CSV.
