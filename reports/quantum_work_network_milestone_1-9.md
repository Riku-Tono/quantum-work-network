# Quantum Work Network

## Project overview

Quantum Work Network is a Rust reference implementation for numerically computing how much of the
energy flowing through a small quantum network arrives at a receiver (the load) as **extractable
work (ergotropy)**.

The physical system consists of several two-level sites (qubits) connected in a chain, plus one
receiver (a 3-level load) attached at one end. One end site is driven by an external pulse to inject
energy, and the way that energy propagates along the chain and accumulates in the load is tracked as
the time evolution of a density matrix. The purpose is to compare conditions with and without phase
noise (dephasing) and to **directly measure**, within a finite-size, finite-time scope, how the
locally extractable work changes when noise is present in the same situation.

What this repository computes is strictly "the values directly observed for this physical model,
these initial conditions, this finite time, and this Rust implementation." It does not prove any
quantum advantage or any universal scaling law. Explanations of cause (why something happens) and
observed results (what happened) are kept distinct.

Development proceeded in stages from Milestone 1 through 9c (and its validation). Each stage was
built by adding new implementations or comparisons on top of the previous stage without changing the
earlier physical model, CSVs, or reports.

---

## Model and numerical conventions

The following conventions were frozen in Milestone 1 and maintained in all subsequent stages.

**Physical model**

- Site configuration: N two-level sites (`|0>` = empty, `|1>` = excited) plus one 3-level load.
- Default parameters (ℏ=1): inter-site coupling `J=1`, site–load coupling `g=0.25`, angular
  frequency of each site, the load, and the drive = 1.
- Chain onsite Hamiltonian: `omega * sum_i |1><1|_i`.
- Drive site = 0 (one end), load coupling site = N-1 (the opposite end).
- Coherent drive: `H_drive(t) = Omega f(t){ exp(-i omega t) sigma_1^+ + exp(+i omega t) sigma_1^- }`,
  envelope `f(t)=sin^2(pi t / tau)` (for `0<=t<=tau`, else 0), defaults `tau=3.2`, `Omega=0.2`.
- Phase noise (dephasing): `L_phi,j = sqrt(gamma_phi/2) sigma_z,j` on each chain site. No direct
  noise is placed on the load.
- Time evolution: the Lindblad master equation `d rho/dt = -i[H(t), rho] + sum_k D[L_k] rho`.

**Basis order and vectorization**

- For general N the tensor order is `|q1, q2, ..., qN, load>`. The rightmost load index varies
  fastest.
- For N=3 this is `|q1, q2, q3, load>`.
- Density matrices use **column-major vectorization**:
  `vec(rho) = [rho(0,0), rho(1,0), ..., rho(0,1), ...]^T`. Therefore
  `vec(A rho B) = (B^T ⊗ A) vec(rho)`.
- For N=3, a `24 x 24` density matrix → a length-576 vector, and the Liouvillian is `576 x 576`.

**Liouvillian convention**

```
L = -i (I ⊗ H - H^T ⊗ I)
    + sum_k [ L_k* ⊗ L_k
            - 1/2 I ⊗ (L_k^dagger L_k)
            - 1/2 (L_k^dagger L_k)^T ⊗ I ]
```

Collapse operators are passed **with their coefficients already included** (e.g.,
`sqrt(gamma) sigma_minus`, `sqrt(gamma_phi/2) sigma_z`).

**Local ergotropy (extractable work)**

```
W(rho_L) = Tr(rho_L H_L) - min_U Tr(U rho_L U^dagger H_L)
```

This is the maximum work locally extractable from the load, computed from the load's reduced state.

---

## Milestone list

For each stage, the "purpose," "what was added/changed," "what was run or verified," and "what was
not yet confirmed at that stage" are given briefly. For numerical details, see the individual
MILESTONE reports and CSVs.

### Milestone 1
- **Purpose**: correctly implement the static modules (operators / partial_trace / ergotropy)
  needed before any time evolution can be trusted.
- **Added/changed**: froze the conventions for basis order, load operators, ergotropy computation,
  etc.
- **Run/verified**: unit tests for each module (the environment used to generate this project had no
  Rust toolchain, so the contemporaneous materials explicitly note that the source and tests were
  produced but not compiled there).
- **Not confirmed**: the correctness of the time evolution (Liouvillian) is out of scope at this
  stage.

### Milestone 2
- **Purpose**: add the two time-evolution modules.
- **Added/changed**: `liouvillian` (column-major vectorization and construction of the
  Hamiltonian/Lindblad superoperators) and `propagator` (accuracy-first propagation via a dense
  matrix exponential at each requested time).
- **Run/verified**: 26 regular tests (vectorization, `t=0` identity, trace/Hermiticity/positivity
  preservation, agreement with closed-system unitary evolution, analytic amplitude-damping and pure
  dephasing, etc.).
- **Not confirmed**: diagnostics, protocol matching, parameter sweeps, and plotting are not yet
  introduced.

### Milestone 2.1
- **Purpose**: confirm that the complete 24-dimensional model composes correctly from a downstream
  user's perspective.
- **Added/changed**: `tests/full_24d_short_time.rs` (an opt-in integration smoke test that builds
  the `576 x 576` Liouvillian and propagates from `t=0` to `t=0.001`).
- **Run/verified**: with the injection collapse operator `sqrt(0.1) sigma_1_plus` on
  `ModelParams::default()`, verified dimensions, trace, Hermiticity, positivity, finiteness, and a
  nonzero change from the vacuum state.
- **Not confirmed**: no diagnostics or matching yet.

### Milestone 3
- **Purpose**: add state diagnostics and signed-power accounting.
- **Added/changed**: `diagnostics` (reduced load energy and ergotropy, energy decomposition, load
  current, source/dephasing power, physicality metrics, signed trapezoidal power integration).
- **Run/verified**: existing tests maintained with the diagnostics present. Protocol matching,
  efficiency claims, parameter sweeps, and plotting are deliberately not introduced.
- **Not confirmed**: efficiency claims, matching, and sweeps are not done at this stage.

### Milestone 3.1
- **Purpose**: make the treatment of the signed-power integral explicit and correct (not a new
  physical experiment separate from Milestone 3, but housekeeping that preserves the diagnostics).
- **Added/changed**: in `integrate_signed_power`, made explicit the handling that, when the power
  changes sign across the endpoints of a trapezoidal interval, splits the interval at the linearly
  interpolated zero crossing. Example: `[(0.0, 1.0), (1.0, -1.0)]` yields `energy_net=0`,
  `energy_in=0.25`, `energy_out=0.25`.
- **Run/verified**: confirmed that `energy_in`/`energy_out` are correctly separated by the
  zero-crossing split.
- **Not confirmed**: since the physical model itself is unchanged, there is no new physical claim.

### Milestone 4
- **Purpose**: using single-shot transport from an initial coherent state, compare the ergotropy of
  A (no noise) and B (with phase noise) when matched to the same time and the same load energy.
- **Added/changed**: the coherent-input experiment and `MILESTONE_4_RESULT.md`, which freezes the
  confirmed results.
- **Run/verified**: energy match across all 16 conditions of evaluation times `3.0/5.0/7.9/10.0` and
  B dephasing strengths `0.1/0.2/0.5/1.0`. In all conditions A ergotropy > B ergotropy (matched
  ratio `1.229`–`49.318`). Diagonal ergotropy was zero in all conditions, so the difference is
  coherence-derived. In the equal-input (`p_B=p_A`) comparison, A > B in all 16 conditions as well.
  Zero physical-check failures.
- **Not confirmed**: this is single-shot transport, not creation from vacuum or continuous supply.
  Some conditions are not equal-input-cost comparisons. Classical comparison and quantum advantage
  are untested.

### Milestone 5a
- **Purpose**: implement and verify the RK4 propagator for the time-dependent Lindblad equation.
- **Added/changed**: `src/time_dependent.rs` (fixed-max-step RK4, time-dependent `H(t)` and collapse
  operators, save schedule), public module exposure, and a single-qubit sanity-check binary.
- **Run/verified**: on a time-invariant reference problem, monotonic convergence toward the dense
  exponential solution under step halving (error `7.80e-7 → 4.99e-8 → 3.19e-9`, error ratio ≈15.6,
  close to the RK4 expectation of 16). `cargo test --release`: 47 passed / 0 failed / 1 ignored.
- **Not confirmed**: time-dependent driving on the production network, A/B comparison, energy
  matching, continuous supply, and whether RK4 preserves positivity for arbitrarily large steps and
  strong drives.

### Milestone 5b
- **Purpose**: confirm (as a sanity check) that a finite pulse drive from the vacuum initial state
  generates load coherence and nonzero ergotropy.
- **Added/changed**: a single `sin^2` pulse applied to the existing 24-dimensional `H0`. A:
  `gamma_phi=0`; B: `gamma_phi=0.5` (all three sites).
- **Run/verified**: maximum load ergotropy A `5.5424e-2` (t=9.48), B `3.0302e-3` (t=5.63). The
  7-item same-time comparison success conditions PASS. Physical checks, energy ledger, and dt-halving
  convergence confirmed.
- **Not confirmed**: fair comparison at the same time and same load energy, energy matching, search,
  continuous driving, work extraction, classical comparison, quantum advantage. Each maximum may
  occur at a different time, so this is for sanity-check use.

### Milestone 5c (the central comparison result)
- **Purpose**: at t=10, conditionally match A/B load energy and compare their ergotropy.
- **Added/changed**: sweep `Omega_B` over `0.2–1.0` (81 points) and refine the sign-change interval
  by bisection. Monotonicity is not assumed.
- **Run/verified**: from the unique sign-change interval, `Omega_B=0.431953125`. At t=10, load energy
  matched to relative difference `4.001e-5` (`<1e-4`). Ergotropy A `5.2798e-2` / B `8.2846e-3`,
  **A/B ratio = 6.373**. All 10 success conditions PASS. Direction unchanged under dt halving.
- **Not confirmed**: only the comparison time, final load energy, model, pulse shape, and drive
  frequency are matched. **The drive strength Omega and the total input energy are not matched**, so
  this is neither an equal-input-cost comparison nor a causal comparison varying noise alone.
  Continuous operation, work extraction, classical comparison, and quantum advantage are unconfirmed.

### Milestone 6a
- **Purpose**: from the 5c final state, actually recover the work using the ideal local unitary
  corresponding to the ergotropy definition — an implementation cross-check (not a new physical
  discovery).
- **Added/changed**: deterministically reconstruct the 5c setup at `dt=0.0025`, cross-check
  `H_drive(10)=0`, and after a sudden switch-off apply an ideal, instantaneous load-local unitary.
- **Run/verified**: recovered gross extracted work = load ergotropy (A `5.2798e-2`, B `8.2846e-3`),
  with post-extraction ergotropy zero in both. Switch work is numerically zero. 18/18 items pass for
  both A and B. Gross-work A-B ratio `6.373`.
- **Not confirmed**: control cost of the extraction unitary, finite-time switch-off/extraction,
  repeated discharge, continuous operation, global extraction from correlations, classical
  comparison, quantum advantage. The switch-work definition is idealized and does not fully account
  for device cost.

### Milestone 7a
- **Purpose**: in the fixed N=3 model, compare, by position, the effect on the load of placing phase
  noise on only one site.
- **Added/changed**: an entry point to specify the set of sites that receive noise (the existing
  all-three-site API is retained as a wrapper).
- **Run/verified**: at `gamma_phi=0.5`, compared the three placements (entrance/middle/exit) against
  noise-free. At t=10, the minimum W was at `site1` (entrance) and the minimum usable fraction at
  `site3` (exit). Middle noise causes less loss than the other two placements. All physical checks
  PASS; the minimum-W position is unchanged under dt halving.
- **Not confirmed**: other `gamma_phi`, Omega, pulse, longer chains, longer times, continuous
  operation, extraction, classical comparison. No universal law of noise position or causal
  assertion is made.

### Milestone 7b
- **Purpose**: from the fixed 7a time series, describe "when" the difference between the noisy
  conditions and noise-free appears persistently (no new time evolution).
- **Added/changed**: a CSV-reading, analysis-only binary (it does not call the Hamiltonian/RK4/etc.).
  Computes persistent onset, time of maximum loss, site population at onset, time-window comparisons,
  and rank switches.
- **Run/verified**: all 40 checks PASS. E and W begin at a diagnostic-level t=2.25 for all three
  conditions and three thresholds. Usable fraction is threshold-dependent (at medium, exit 1.27 <
  entrance 1.57 < middle 1.62). It presents a "clue" that the reason middle's damage is lighter
  cannot be fully explained by maximum population or time-area alone.
- **Not confirmed**: causal decomposition into specific physical processes, explanation by population
  alone, coherence/inter-site current, recovery via protection, generalization.

### Milestone 7c
- **Purpose**: in a counterfactual condition where, from all-three-site noise, only the phase-noise
  operator of a specified site is ideally removed, compare the upper bound of load-quantity recovery.
- **Added/changed**: five conditions `all_noisy=[0,1,2]`, `protect_entrance=[1,2]`,
  `protect_exit=[0,1]`, `protect_both_ends=[1]`, `noise_free=[]`. This is ideal protection (complete
  removal of the collapse operator), not real devices, control, cost, or error correction.
- **Run/verified**: all 75 checks PASS. At t=10, `protect_both_ends` (leaving only middle noise) has
  the largest recovery in W, usable fraction, and E. Entrance and exit protection have similar
  recovery amounts, and their ranking swaps over time. The non-additivity (synergy) of protection is
  classified as positive_nonadditivity (synergy is the non-additivity of the observable response,
  not an interaction energy).
- **Not confirmed**: realistic implementation, cost, imperfect protection, other parameters, long
  times, causal mechanism. No claim that middle noise is harmless.

### Milestone 7d
- **Purpose**: fixing the middle-site gamma at 0.5, study the recovery curve as only the end-site
  gammas are lowered from `0.5→0`.
- **Added/changed**: a per-site gamma API (rejects negative/non-finite values, excludes gamma=0
  collapse). Sweep points `0.50/0.40/0.30/0.20/0.15/0.10/0.05/0.00`.
- **Run/verified**: all 152 checks PASS. Endpoints reproduce 7c within absolute error `1e-9`. Both
  W_max and usable fraction are discretely monotonic non-decreasing. The maximum-sensitivity interval
  is `0.05→0.00` in every case. Curvature and sensitivity are computed from discrete points only and
  are not called critical exponents or phase transitions.
- **Not confirmed**: the implementable protection strength required, physical susceptibility /
  critical exponents / phase transitions, other parameters, long times / continuous operation,
  causal mechanism.

### Milestone 8a
- **Purpose**: change only the chain length from N=3 to N=5 and compare load energy, ergotropy, and
  usable fraction.
- **Added/changed**: operator-construction and drive-execution APIs that take a `chain_length`
  argument (the existing N=3 API is preserved). `dim=2^N*3`.
- **Run/verified**: reproduces the N=3 regression within absolute error `2e-9`. Confirmed load
  ergotropy generation for N=5 noise-free. W_max ratio N5/N3 is free `0.3965`, noisy `0.3620`. At
  t=10 and at the individual peaks, N=5 W < N=3 W and noisy < free (same ordering conclusion, but
  different ratios). The dt-halving consistency PASS for the three specified conditions (the N=3
  noisy halving was not run).
- **Not confirmed**: N>5, continuous N sweep, position-specific weaknesses, general scaling, total-
  noise matching, protection cost. Changing chain length simultaneously changes distance, dimension,
  bond count, and (in all-site) total noise (not an equal-total-dissipation comparison).

### Milestone 8b
- **Purpose**: assess, with only a short-time probe, whether N=7 dense Lindblad RK4 is feasible in
  the current compute environment.
- **Added/changed**: a feasibility-probe binary doing construction-only, one step, and short-time
  probes.
- **Run/verified**: Hilbert dimension 384, density matrix `384x384`. Numerical quality of
  construction, one step, and short-time probe PASS. From the `t=0.1` noisy step time `21.32 s/step`,
  the t=10 estimate is ≈23.7 hours, establishing **infeasible_with_current_dense_method** for the
  current dense approach.
- **Not confirmed**: the N=7 t=10 final values, halved-step measurement, long-time stability,
  post-optimization performance. This is feasibility "in the current environment and current
  implementation," not physical possibility.

### Milestone 8c
- **Purpose**: exactly accelerate only the Lindblad term of local sigma_z phase noise, without
  changing the physical model.
- **Added/changed**: `DiagonalDephasingKernel` (adds `-Gamma[a,b] rho[a,b]` to each element — an
  exact component-wise representation of the same Lindblad dissipator, not a physical approximation).
  The old dense path is retained.
- **Run/verified**: checks CSV **140 PASS / 0 FAIL**. Agrees with the dense path and existing values
  within tolerance for N=3/N=5. N=7 noisy median `0.7116 s/step`, **29.96x** vs the old value. Re-
  estimated t=10 ≈0.791 hours, classification **feasible_candidate**.
- **Not confirmed**: the N=7 t=10 final physical quantities, long-time measured performance,
  N scaling, application to other Lindblad operators.

### Milestone 9a
- **Purpose**: run N=7 noise-free to t=10 at `dt=0.0025` and confirm arrival and the difference vs
  N=3/N=5.
- **Added/changed**: `n7_noise_free_full.rs` and the associated artifacts.
- **Run/verified**: all checks PASS (max trace `2.109e-15`, min eigenvalue `-6.578e-12`, etc.). At
  t=10: E `1.3961e-2`, W `1.3085e-2`, usable `0.9373`. W_max `2.2436e-2` (t=7.71). W_max ratio
  N7/N5=`1.0210`, N7/N3=`0.4048`. Classification `peak_resolved`.
- **Not confirmed**: N=7 noisy, fine steps, t>10, the final arrival upper bound, N>7, scaling,
  real-device performance.

### Milestone 9b
- **Purpose**: run N=7 all-site noisy with `gamma_phi=0.5` on all 7 sites to t=10.
- **Added/changed**: `n7_all_site_noisy_full.rs` (using the 8c kernel) and the associated artifacts.
- **Run/verified**: all checks PASS. At t=10: E `3.3041e-3`, W `3.1157e-4`. W_max `4.0437e-4`
  (t=7.65), classification `peak_resolved`. In fixed-per-site noisy, W_max is N3 `3.0302e-3` >
  N5 `1.0968e-3` > N7 `4.0437e-4` (N7/N5=`0.3687`). **The "N7 W_max > N5 W_max" feature seen in
  noise-free did not survive in this fixed-per-site noisy condition.**
- **Not confirmed**: dt halving, t>10, position-specific noise, protection, gamma/Omega sweep, N>7,
  real-device performance. Judgment **completed_comparison**. Note that at this stage, distance and
  the number of noisy sites (hence total noise) increase together with N (not an equal-total-noise
  comparison).

### Milestone 9c
- **Purpose**: fix the simple sum of all site gammas at `TOTAL_GAMMA=1.5` to partially disentangle
  the confound of "increasing N" and "increasing total noise" (fixed-total-noise comparison).
- **Added/changed**: newly compute N=5 (gamma_site=0.3) and N=7 (gamma_site=1.5/7). N=3 references
  the existing gamma_site=0.5 result.
- **Run/verified (interim status at this stage)**: N=5 checks=true. For N=7, one minimum-eigenvalue
  diagnostic point at `t=0.02` became `NaN`, and the `positivity` and `finite_values` checks FAILed.
  **The report for this stage, `MILESTONE_9C_REPORT.md`, recorded `numerical_issue_stop` as an
  interim judgment** (as explained below, this is not the final conclusion). The minimum eigenvalue
  at the other 1000 saved points is finite, with a minimum of `-5.278e-18`, and trace/Hermiticity/
  ledger were normal.
- **Not confirmed (at this stage)**: isolating the cause of the `NaN`, and whether the comparison
  results may be formally adopted, were incomplete at this stage (resolved by the diagnostic and
  validation).

### Milestone 9c diagnostic
- **Purpose**: isolate the N=7 `t=0.02` minimum-eigenvalue `NaN` using only a short-time
  recomputation (the `t=10` full run is not re-executed).
- **Added/changed**: a diagnostic binary that recomputes `t=0` to `t=0.03` over 12 RK4 steps, three
  times.
- **Run/verified**: at `t=0.02` the density matrix rho is **finite in all elements**, trace `≈1`,
  Hermiticity error 0, with reproducibility across the three runs. The cause is on the
  **eigenvalue-solver side**: raw / Hermitianized SymmetricEigen returned non-finite eigenvalues, and
  the 9c CSV formatter recorded these uniformly as `NaN`. An independent solver (Complex Schur)
  returned all-finite eigenvalues with minimum `-3.237e-24`. Maximum absolute difference of 0 across
  24 scalar comparisons.
- **Not confirmed (at this stage)**: reconstruction beyond `t=0.03`, other dt/gamma, t=10
  recomputation. **This diagnostic is merely a supplementary record and does not overwrite the
  original report. The diagnostic label is `state_level_numerical_issue` (a label from the stopping
  rule); it does not mean rho itself had a NaN, a trace anomaly, a Hermiticity anomaly, or
  nondeterminism.**

### Milestone 9c validation (the final source of truth for 9c)
- **Purpose**: re-verify N=7 fixed-total-noise to t=10 with a robust positivity diagnostic that
  judges state finiteness and solver finiteness independently.
- **Added/changed**: a solver policy with Hermitianized SymmetricEigen as primary, falling back to
  Hermitianized Complex Schur only on failure. The physical model and time evolution are unchanged
  (the fallback is in the diagnostic layer only).
- **Run/verified**: conditions N=7, total gamma=1.5, gamma_site=1.5/7, `dt=0.0025`, 4000 RK4 steps,
  1001 saved points. Validation values appear in "Current main results" below. **Final judgment
  `completed_comparison_with_fallback_diagnostic`. The physical comparison results may be formally
  adopted.**
- **Not confirmed**: other dt, other gamma, t>10, N>7, other external solver crates. The "10% band"
  is not a statistically significant difference and not a proof of distance-only causation.

---

## Milestone 1–3: foundational implementation

Milestones 1–3.1 lay the groundwork before the time evolution can be trusted. First the static
modules (operators, partial trace, ergotropy) were fixed (M1); then the Liouvillian based on
column-major vectorization and the dense matrix-exponential propagator were added (M2); and an
integration test confirmed that the complete 24-dimensional model composes (M2.1). After that, state
diagnostics and signed-power accounting were introduced (M3), and the treatment of signed-power
integration — splitting sign-reversal intervals at the zero crossing — was made explicit (M3.1).

The propagator at this stage (`DenseExponentialPropagator`) computes `exp(L t) vec(rho(0))`
independently at each time; it is an accuracy-first implementation. The `576 x 576` matrix
exponential is heavy, but it serves as a correctness baseline. Efficient computation is deliberately
deferred to later stages.

---

## Milestone 4–6: comparison experiments and work extraction

This is where the core comparison begins. In Milestone 4, using single-shot transport from an
initial coherent state, noise-free A exceeded phase-noisy B in all 16 conditions matched to the same
time and the same load energy. Milestone 5a verified the time-dependent RK4 propagator, and 5b
confirmed that a finite pulse from vacuum generates nonzero load ergotropy.

The central result is Milestone 5c. When A/B load energy at t=10 was matched to relative error
`4.001e-5`, the ergotropy ratio was **A/B=6.373**. However, this match is a conditional match on the
"final load energy," and the drive strength Omega and total input energy differ between A and B. It
is therefore neither an equal-input-cost comparison nor a causal comparison varying noise alone.

Milestone 6a is not a new physical discovery but an **implementation cross-check** showing that,
by constructing the ideal local unitary corresponding to the ergotropy definition, the predicted
amount of work can be recovered.

(A summary of this 4–6a segment is also collected in `Milestone_4-6a_研究結果ノート.pdf`.)

---

## Milestone 7: noise position and partial protection

This segment, in the fixed N=3 model, studies "where to place" and "where to remove" noise. 7a is a
harmful-placement comparison putting noise on only one site; at t=10 the minimum W was at the
entrance (site1) and the minimum usable fraction at the exit (site3). 7b describes, from that time
series, when the difference appears persistently (no new time evolution).

7c does the reverse: comparing the recovery upper bound when noise on a specific site is ideally
removed from the all-site noise, with both-ends protection (leaving only middle noise) giving the
largest recovery in W, usable fraction, and E. 7d fixes the middle gamma and lowers the end gammas
from `0.5→0`, confirming that recovery connects in a discretely monotonic non-decreasing way. All of
these are ideal noise removal, not real protection devices, cost, or imperfect protection.

---

## Milestone 8: chain-length generalization and N=7 feasibility

This segment generalizes the chain length. 8a changes only the chain length from N=3 to N=5,
confirming load-ergotropy arrival and the shrinking of W_max. 8b probes the feasibility of N=7 with
the dense method and judges it infeasible (≈23.7 hours to t=10 with the current dense
implementation). 8c, without changing the physical model, replaces only the sigma_z dephasing term
with an exact component-wise representation (`DiagonalDephasingKernel`); after confirming it agrees
with the dense path for N=3/N=5, it obtains an ≈30x speedup and updates N=7 t=10 to
feasible_candidate.

Note that changing chain length simultaneously changes propagation distance, Hilbert dimension, bond
count, and (in all-site) total noise. No distance-only causation or scaling law is claimed from a
2–3 point finite-length comparison.

---

## Milestone 9: N=7 full runs and fixed-total-noise comparison

Using the 8c speedup, this segment proceeds to the N=7 full runs.

9a runs N=7 noise-free to t=10, confirming load-energy and ergotropy arrival. 9b runs N=7 all-site
noisy (gamma=0.5 per site, fixed-per-site). In this fixed-per-site condition the total noise grows
with chain length (N3=1.5, N5=2.5, N7=3.5), so the "N7 W_max > N5 W_max" feature seen in noise-free
did not survive (N7/N5=0.369).

9c performs a comparison that fixes the simple sum of all site gammas at `TOTAL_GAMMA=1.5`
(fixed-total-noise), to partially disentangle the confound of "increasing N" and "increasing total
noise." For 9c, the history is organized in the following order:

1. **An initial diagnostic found a failure of the eigenvalue solver.** In the 9c full run for N=7,
   one minimum-eigenvalue diagnostic point at `t=0.02` became `NaN`, and the `positivity` and
   `finite_values` checks FAILed. For this reason `MILESTONE_9C_REPORT.md` recorded
   `numerical_issue_stop` as an interim status.

2. **An exact diagnostic was performed.** In `MILESTONE_9C_DIAGNOSTIC.md`, a short-time
   recomputation of `t=0.02` found the density matrix rho finite in all elements, trace `≈1`,
   Hermiticity error 0, and agreement across three runs. The `NaN` originates from non-finite output
   **on the eigenvalue-solver side, not the state**, which the CSV formatter recorded uniformly as
   `NaN`. An independent Complex Schur gave all-finite eigenvalues.

3. **At the two points where the symmetric eigenvalue computation failed, a Schur-decomposition
   fallback was used.** In `MILESTONE_9C_VALIDATION.md`, a robust diagnostic was introduced with
   Hermitianized SymmetricEigen as primary, falling back to Hermitianized Complex Schur only at the
   two failing times (`t=0.01`, `t=0.02`). The fallback is in the diagnostic layer only and does not
   change the physical time evolution.

4. **Positivity was determined at all 1001 times.** Primary success 999, failure 2, fallback success
   2/2, solver_failure 0. Worst selected minimum eigenvalue `-5.278e-18`.

5. **The comparison against the existing 9c trajectory had a maximum difference of 0.** The specified
   physical quantities were compared over 1001 times at tolerance `1e-12`, and all maximum absolute
   differences were 0 (all checks=true).

6. **The final state is `completed_comparison_with_fallback_diagnostic`.**

7. **Therefore, the N=3, N=5, N=7 comparison under fixed total phase noise may be adopted as a formal
   comparison result.**

From the above, the **final source of truth for 9c is `MILESTONE_9C_VALIDATION.md`**. The
`numerical_issue_stop` in `MILESTONE_9C_REPORT.md` is an interim status, not the final conclusion.
`MILESTONE_9C_DIAGNOSTIC.md` is the record of the interim diagnostic.

---

## Current main results

All of these are "values directly confirmed for this model, these initial conditions, this finite
time, and this Rust implementation." They are not general laws.

**Milestone 5c (N=3, conditional load-energy match, t=10)**
- Ergotropy A/B ratio: `6.373` (A `5.2798e-2` / B `8.2846e-3`), load-energy relative difference
  `4.001e-5`.
- Only the final load energy is matched. Omega and total input energy are not matched.

**Milestone 9a/9b (N=7, fixed-per-site, for reference)**
- Noise-free W_max: `2.2436e-2` (t=7.71).
- All-site noisy (gamma=0.5) W_max: `4.0437e-4` (t=7.65).
- Fixed-per-site W_max: N3 `3.0302e-3` > N5 `1.0968e-3` > N7 `4.0437e-4`. The noise-free "N7 > N5"
  does not survive in this condition.

**Milestone 9c validation (N=7, fixed-total-noise, total gamma=1.5, final source of truth)**
- Positivity diagnostic completed: **1001 / 1001 times**, solver_failure=0.
- SymmetricEigen (primary) successes: **999**, failures: **2** (`t=0.01`, `t=0.02`).
- Complex Schur fallback successes: **2 / 2**.
- Maximum difference vs the existing 9c trajectory: **0** (1001 times, tolerance `1e-12`).
- N=7 W_max: **0.0038080717406769921** (t=7.70).
- Fixed-total W_max: N3 `3.0302e-3`, N5 `3.4876e-3`, N7 `3.8081e-3`.
- N=7 / N=5 W_max ratio: **1.0918825139**. This ratio is **within the descriptive "10% band."**

Note that this "10% band" is not a theoretical law but merely a **descriptive guideline for this
finite-size comparison**. It is neither a statistically significant difference nor a proof of
distance-only causation.

---

## Build and test

Run in an environment with a Rust toolchain (`Cargo.toml`: edition 2021, dependencies `nalgebra`,
`num-complex`, `thiserror`, dev-dependency `approx`).

```bash
# Format check
cargo fmt --all -- --check

# Regular tests (release, offline)
cargo test --release --offline

# Build
cargo build --release --offline
```

Run the deliberately-ignored `576 x 576` 24-dimensional smoke test explicitly:

```bash
cargo test --release full_24d_short_time_smoke_test -- --ignored --nocapture
```

---

## Main run commands

Each milestone runs via a dedicated binary (representative examples; for detailed arguments and
outputs, see the "generated files" list at the end of each report).

```bash
# 5a: time-dependent RK4 sanity check
cargo run --release --offline --bin time_dependent_sanity

# 8c: exact dephasing-kernel benchmark and equivalence
cargo run --release --offline --bin dephasing_kernel_benchmark

# 9a: N=7 noise-free full run
cargo run --release --offline --bin n7_noise_free_full

# 9b: N=7 all-site noisy full run
cargo run --release --offline --bin n7_all_site_noisy_full

# 9c: fixed-total-noise comparison
cargo run --release --offline --bin fixed_total_noise_comparison

# 9c diagnostic: isolating the t=0.02 eigenvalue NaN
cargo run --release --offline --bin n7_t002_eigen_diagnostic

# 9c validation: final re-verification via robust positivity diagnostic
cargo run --release --offline --bin n7_fixed_total_validation
```

The N=7 full runs take on the order of tens of minutes (9a total ≈ 2953s, 9b total ≈ 2899s, 9c
validation total ≈ 2703s).

---

## Output files

Each stage generates a report (`MILESTONE_*.md`) and several CSVs. This README is an entry point; for
detailed tables, see the reports and CSVs. Representative ones:

- **5c**: `coherent_drive_match_{grid,roots,comparison,timeseries,convergence}.csv`
- **6a**: `explicit_load_extraction_{results,checks,mapping}.csv`
- **7a/7b**: `local_noise_placement_*.csv` / `local_noise_damage_*.csv`
- **7c/7d**: `ideal_partial_protection_*.csv` / `partial_end_protection_*.csv`
- **8a/8b/8c**: `chain_length_reachability_*.csv` / `n7_feasibility_*.csv` / `dephasing_kernel_*.csv`
- **9a/9b**: N=7 noise-free / all-site noisy timeseries, summary, checks, etc.
- **9c validation (source of truth)**:
  - `n7_fixed_total_validation_summary.csv` (t=10 values, W_max, solver accounting, final judgment)
  - `n7_fixed_total_validation_checks.csv` (per-check PASS record)
  - `n7_fixed_total_validation_trajectory_comparison.csv` (difference vs the existing 9c trajectory,
    all 0)
  - `n7_fixed_total_validation_eigen_diagnostics.csv`, `..._timeseries.csv`, `..._performance.csv`
  - `fixed_total_noise_final_comparison.csv` (N3/N5/N7 W_max comparison)
  - `robust_eigen_diagnostic_unit_checks.csv`

**Key points for reading the CSVs**: `W_time_area` and `E_time_area` are time-areas of state
quantities, not cumulative inflow energy or cumulative extracted work. `W/Ein` is not an overall
efficiency including control cost. `usable_fraction` is the ratio of ergotropy to load energy.

---

## Numerical notes

- The phase-noise time evolution is computed with dense density-matrix RK4. The
  `DiagonalDephasingKernel` introduced in 8c is an **exact component-wise representation** of the
  sigma_z dephasing term, not a physical approximation (agreement with the dense path is confirmed
  within tolerance).
- Because RK4 does not strictly guarantee a completely positive map, the minimum eigenvalue is
  checked explicitly. Observed negative values are within rounding error (e.g., the 9c validation
  worst `-5.278e-18`), and no eigenvalue correction is applied.
- Even when the eigenvalue **diagnostic** produces a non-finite value, that does not mean the density
  matrix itself is anomalous. In 9c, state finiteness and solver finiteness are separated, and using
  a Schur fallback only on primary failure, positivity was determined at all 1001 times (the fallback
  is in the diagnostic layer only).
- Because the ledger-residual denominator is nearly zero, the absolute difference (rather than the
  relative difference) is used as the primary criterion in some places.
- Changing chain length N simultaneously changes distance, dimension, bond count, and (in all-site)
  total noise. Distance-only causation or scaling laws cannot be derived from a comparison of
  finitely many N.

---

## What has not been confirmed

The following are things this project **has not yet shown**. They are also things that must not be
asserted in this README.

- Quantum advantage / quantum supremacy, or being higher-performing than a classical method.
- Universal scaling laws, exponential/power-law decay, the thermodynamic limit, or N→∞ behavior.
- A "general mechanism by which noise improves performance," or universal laws of noise position or
  protection strength.
- Real-world power-transmission / storage / device efficiency, or effectiveness on real hardware.
- Control cost of the extraction unitary, finite-time switch-off/extraction, repeated cycles,
  continuous supply, steady-state output, long-time stability.
- Other dt, other gamma, t>10, N>7, or other external solver crates (out of scope for 9c
  validation).
- A fair baseline comparison against classical wave or stochastic models.
- Novelty / priority relative to the literature (no literature review has been done).
- Causal explanations not written in the reports (e.g., asserting that noise on a specific site
  breaks "only injection/transport/handoff").

In the 9c validation, all 47 runtime-check items in `n7_fixed_total_validation_checks.csv` PASSed.
The 10 required items of the robust eigenvalue diagnostic recorded in
`robust_eigen_diagnostic_unit_checks.csv` also all PASSed. A total of "90/90" that is referenced
elsewhere for the whole Cargo test suite cannot be confirmed from the materials used to write this
README, so it is not adopted as a verification count here.

---

## Repository layout

The main modules exposed by `Cargo.toml` (`quantum_work_network`, edition 2021) and `src/lib.rs`:

```
src/
  operators.rs            # operators (M1)
  partial_trace.rs        # partial trace (M1)
  ergotropy.rs            # local ergotropy (M1)
  matrix.rs               # ComplexMatrix / C64
  error.rs                # PhysicsError
  liouvillian.rs          # column-major vectorization / superoperators (M2)
  propagator.rs           # dense matrix-exponential propagation (M2)
  diagnostics.rs          # state diagnostics / signed-power accounting (M3 / M3.1)
  time_dependent.rs       # time-dependent RK4 propagator (M5a)
  coherent_drive.rs       # coherent drive (M5b/5c line)
  coherent_drive_matching.rs
  matching.rs / protocol.rs / experiment.rs
  load_extraction.rs      # ideal local-unitary extraction (M6a)
  dephasing_kernel.rs     # exact dephasing kernel (M8c)
  bin/
    time_dependent_sanity.rs        # M5a
    local_noise_placement.rs        # M7a
    local_noise_damage_analysis.rs  # M7b
    ideal_partial_protection.rs     # M7c
    partial_end_protection.rs       # M7d
    chain_length_reachability.rs    # M8a
    n7_feasibility_probe.rs         # M8b
    dephasing_kernel_benchmark.rs   # M8c
    n7_noise_free_full.rs           # M9a
    n7_all_site_noisy_full.rs       # M9b
    fixed_total_noise_comparison.rs # M9c
    n7_t002_eigen_diagnostic.rs   # M9c diagnostic
    n7_fixed_total_validation.rs  # M9c validation (final source of truth)
tests/
  full_24d_short_time.rs            # M2.1 (ignored smoke test)
MILESTONE_*.md                      # per-stage reports
*.csv                               # per-stage outputs
```

The `src/lib.rs` module declarations above are based on the actual files. The mapping inside `bin/`
is organized from the file names listed in each report's "generated files" section. For a complete
list of individual CSVs, see the end of the corresponding MILESTONE report.

---

## Milestone 10: fixed-total comparison and the XGamma diagnostic

In Milestone 10, after organizing the existing results, fixed total phase noise `TOTAL_GAMMA=1.5`
and `3.0` were compared for N=3, 5, 7 under the same diagnostic system.

- **Milestone 10a**: without running new physical computations, organized and compared existing
  results and made explicit the missing values in the fixed-total conditions.
- **Milestone 10b**: computed `TOTAL_GAMMA=3.0` for N=3, 5, 7 and first introduced XGamma.
- **Milestone 10c**: recomputed `TOTAL_GAMMA=1.5` with XGamma for N=3, 5, 7, filling in the gaps
  from 10a. The N=7 physical trajectory agreed with the 9c source of truth across 1001 times.

For both fixed-total 1.5 and 3.0, `W_max`, `W(t=10)`, and usable fraction were `N=7 > N=5 > N=3`,
`W_time_area` was `N=3 > N=5 > N=7`, and ergotropy arrival was fastest for N=3. Thus, rather than
deciding a single overall ranking of chain lengths, this is a finite-condition result in which
**the ranking reverses depending on the metric**. `W_time_area` is the time-integral of the
ergotropy state quantity, not cumulative extracted work.

XGamma is the following dephasing-kernel-weighted coherence exposure:

```text
x_gamma(t) = sum_ab Gamma[a,b] |rho[a,b](t)|^2
XGamma(T) = integral_0^T x_gamma(t) dt
```

XGamma is a diagnostic quantity — not lost work, dissipated energy, dephasing power, heat, entropy
production, efficiency, or damage. This finite comparison does not show a general advantage of longer
chains, a scaling law, a causal mechanism of XGamma, a universal multiplier for doubling gamma, or
quantum advantage.
