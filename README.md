# Quantum Work Network

A Rust reference implementation for numerically computing how much of the energy driven through a
small quantum network arrives at a receiver (the *load*) as **extractable work (ergotropy)**, and how
that changes under phase noise. The project is developed as a single, incremental research and
implementation history from Milestone 1 through Milestone 11M, followed by finite-chain (N=3..7)
event-structure diagnostics.

> Intended to be committed as the repository-root `README.md`.

---

## Project overview

The physical system is a chain of two-level sites (qubits) with one 3-level load attached at one end.
One end site is driven by an external pulse to inject energy; the propagation of that energy along the
chain and its accumulation in the load are tracked as the time evolution of a density matrix under a
Lindblad master equation. The recurring question across all milestones is narrow and concrete: for the
*same situation*, how does the locally extractable work at the load change when phase noise
(dephasing) is present — measured directly, within a finite-size, finite-time scope.

Each milestone builds on the previous one without modifying earlier code, CSVs, or reports. Later
stages generalize the chain length (N=3 -> 5 -> 7), separate two noise conventions (fixed-per-site vs
fixed-total), and finally compare N=3 and N=7 at **equal drive input energy** and analyze the shape of
the resulting work curves.

---

## Scope and limitations

What this repository provides are "values directly observed for this physical model, these initial
conditions, this finite time, and this Rust implementation." In particular, this project **does not**:

- prove any quantum advantage or that quantum beats a classical method;
- establish universal scaling laws, exponential/power-law decay, thermodynamic limits, or N->infinity
  behavior;
- claim a general advantage of longer chains, or universal laws of noise position or protection
  strength;
- treat curve-transform fits or decay fits as proofs of a physical mechanism — these are
  **descriptive shape diagnostics** only;
- treat `W_time_area` as cumulative extracted work, or `usable_fraction` as an independent performance
  gain, or XGamma as lost work / dissipated energy / efficiency / damage;
- assert causal explanations not written in the reports.

Throughout, **observed results** (what happened) and **candidate interpretations** (why it might
happen) are kept distinct, and the fixed-total-dephasing comparison (Milestone 10) is not conflated
with the equal-input comparison (Milestone 11).

---

## Physical model and numerical conventions

Frozen in Milestone 1 and maintained thereafter.

**Model.** N two-level sites (`|0>` empty, `|1>` excited) plus one 3-level load. Defaults (hbar=1):
inter-site coupling `J=1`, site-load coupling `g=0.25`, all angular frequencies = 1. Chain onsite
Hamiltonian `omega * sum_i |1><1|_i`. Drive site = 0 (one end), load coupling site = N-1 (other end).

**Coherent drive.**
`H_drive(t) = Omega f(t){ exp(-i omega t) sigma_1^+ + exp(+i omega t) sigma_1^- }`,
envelope `f(t)=sin^2(pi t / tau)` for `0<=t<=tau` (else 0), defaults `tau=3.2`, `Omega=0.2`.

**Phase noise.** `L_phi,j = sqrt(gamma_phi/2) sigma_z,j` on each chain site; no direct noise on the
load. Time evolution: `d rho/dt = -i[H(t), rho] + sum_k D[L_k] rho`.

**Basis and vectorization.** Tensor order `|q1, q2, ..., qN, load>` (rightmost load index varies
fastest; for N=3, `|q1, q2, q3, load>`). Column-major vectorization
`vec(rho) = [rho(0,0), rho(1,0), ..., rho(0,1), ...]^T`, so `vec(A rho B) = (B^T kron A) vec(rho)`. For
N=3 the `24 x 24` density matrix maps to a length-576 vector and a `576 x 576` Liouvillian.

**Liouvillian.**

```
L = -i (I kron H - H^T kron I)
    + sum_k [ L_k* kron L_k
            - 1/2 I kron (L_k^dagger L_k)
            - 1/2 (L_k^dagger L_k)^T kron I ]
```

Collapse operators are passed with coefficients already included (e.g. `sqrt(gamma_phi/2) sigma_z`).

**Local ergotropy (extractable work).**

```
W(rho_L) = Tr(rho_L H_L) - min_U Tr(U rho_L U^dagger H_L)
```

computed from the load's reduced state.

---

## Milestone roadmap

Numerical details for every stage live in the individual `MILESTONE_*.md` reports and CSVs; this
section gives only the arc and the representative outcomes.

### Milestones 1-3: foundations

Static modules first (operators, partial trace, ergotropy) with conventions frozen (**M1**); then the
Liouvillian and an accuracy-first dense matrix-exponential propagator (**M2**, 26 regular tests); an
opt-in `576 x 576` integration smoke test confirming the full 24-dimensional model composes (**M2.1**);
state diagnostics and signed-power accounting (**M3**); and an explicit fix for signed-power
integration that splits sign-reversal intervals at the linearly interpolated zero crossing (**M3.1**,
e.g. `[(0,1),(1,-1)] -> energy_net=0, energy_in=0.25, energy_out=0.25`). The `DenseExponentialPropagator`
is a correctness baseline; efficiency is deferred.

### Milestones 4-6: comparison and extraction

**M4** — single-shot transport from an initial coherent state: across all 16 conditions matched to the
same time and same load energy, noise-free A > phase-noisy B in ergotropy (matched ratio
`1.229`-`49.318`); differences are coherence-derived; zero physical-check failures. **M5a** implements
and verifies the time-dependent RK4 propagator (step-halving convergence `7.80e-7 -> 4.99e-8 -> 3.19e-9`,
ratio ~15.6; `47 passed / 0 failed / 1 ignored`). **M5b** confirms a finite pulse from vacuum generates
nonzero load ergotropy (A `5.5424e-2`, B `3.0302e-3`; sanity check).

**M5c (central comparison).** With A/B load energy at t=10 matched to relative difference `4.001e-5`
(`Omega_B=0.431953125`), ergotropy A `5.2798e-2` / B `8.2846e-3`, **A/B = 6.373**; all 10 success
conditions PASS. Only the final load energy is matched — **Omega and total input energy are not** — so
this is neither an equal-input-cost comparison nor a noise-only causal comparison.

**M6a** is an implementation cross-check (not a new physical discovery): reconstruct the 5c final state
and recover gross extracted work = load ergotropy via the ideal local unitary (18/18 checks pass;
post-extraction ergotropy zero; gross-work A-B ratio `6.373`).

### Milestone 7: noise location and protection (fixed N=3)

**M7a** places noise on a single site: at t=10 the minimum W is at the entrance (`site1`) and the
minimum usable fraction at the exit (`site3`); middle noise costs less than either end. **M7b** (no new
time evolution) describes *when* the noise-vs-noise-free difference appears persistently (E and W begin
at a diagnostic-level t=2.25; usable fraction is threshold-dependent). **M7c** ideally removes noise
from selected sites; both-ends protection (leaving only middle noise) gives the largest recovery in W,
usable fraction, and E, with protection non-additivity classified as `positive_nonadditivity`. **M7d**
sweeps end-site gamma `0.5->0` (center fixed): recovery is discretely monotonic non-decreasing, maximum
sensitivity interval `0.05->0.00`. All of 7c/7d use **ideal** noise removal — not real devices, cost, or
imperfect protection; curvatures are not called critical exponents or phase transitions.

### Milestones 8-9: chain-length extension and validation

**M8a** changes only chain length N=3 -> N=5 (regression within `2e-9`); W_max ratio N5/N3 free `0.3965`,
noisy `0.3620`. **M8b** probes N=7 dense feasibility and finds it infeasible with the current dense
method (~23.7 h to t=10). **M8c** replaces only the sigma_z dephasing term with an exact component-wise
`DiagonalDephasingKernel` (not an approximation; `140 PASS / 0 FAIL`), giving ~29.96x speedup and
updating N=7 t=10 to `feasible_candidate`.

**M9a** runs N=7 noise-free to t=10 (W_max `2.2436e-2` at t=7.71; `peak_resolved`). **M9b** runs N=7
all-site noisy (gamma=0.5 per site, **fixed-per-site**, so total noise grows with N: N3=1.5, N5=2.5,
N7=3.5): the noise-free "N7 W_max > N5 W_max" feature does **not** survive here (N7/N5=0.369).

**M9c (fixed-total-noise)** fixes the summed site gammas at `TOTAL_GAMMA=1.5` to partially disentangle
"increasing N" from "increasing total noise." Its history, in order:

1. An initial diagnostic found an **eigenvalue-solver** failure: one minimum-eigenvalue point at
   `t=0.02` became `NaN`, failing `positivity` and `finite_values`, so `MILESTONE_9C_REPORT.md`
   recorded `numerical_issue_stop` **as an interim status**.
2. An exact diagnostic (`MILESTONE_9C_DIAGNOSTIC.md`) showed rho itself finite at `t=0.02` (trace `~1`,
   Hermiticity error 0, reproducible); the `NaN` came from non-finite **solver** output recorded
   uniformly by the CSV formatter. Independent Complex Schur gave all-finite eigenvalues.
3. The robust re-verification (`MILESTONE_9C_VALIDATION.md`) uses Hermitianized SymmetricEigen as
   primary and falls back to Hermitianized Complex Schur only at the two failing times (`t=0.01`,
   `t=0.02`); the fallback is in the diagnostic layer only and does not change the time evolution.
4. Positivity was determined at all **1001/1001** times (primary 999 success / 2 failure, fallback
   2/2, solver_failure 0; worst selected minimum eigenvalue `-5.278e-18`).
5. The comparison against the existing 9c trajectory had **maximum difference 0** (1001 times,
   tolerance `1e-12`).
6. Final state `completed_comparison_with_fallback_diagnostic`.
7. The N=3, N=5, N=7 comparison under fixed total phase noise may therefore be adopted as a formal
   result.

**The final source of truth for 9c is `MILESTONE_9C_VALIDATION.md`.** The `numerical_issue_stop` in
`MILESTONE_9C_REPORT.md` is interim; `MILESTONE_9C_DIAGNOSTIC.md` is the interim diagnostic record.

### Milestone 10: fixed-total-dephasing comparison and XGamma

Milestone 10 organizes the fixed-total comparison the later work builds on, distinguishing two noise
conventions: **fixed-per-site** (per-site rate constant, summed total grows with N) and **fixed-total**
(summed rate held at `TOTAL_GAMMA`, per-site rate `TOTAL_GAMMA/N`). For N=5, 7 these are not the same
total-noise condition.

- **M10a** compares existing results (no new time evolution), separates the two conventions, and leaves
  missing fixed-total N=3/N=5 entries as `not_available` rather than guessing
  (`completed_with_explicit_missing_values`).
- **M10b** newly computes `TOTAL_GAMMA=3.0` for N=3, 5, 7 and introduces **XGamma**
  (`completed_fixed_total_gamma_3_comparison`).
- **M10c** recomputes `TOTAL_GAMMA=1.5` under the same XGamma diagnostics; its N=7 trajectory matches
  the 9c source of truth to maximum difference 0 across 7 quantities x 1001 times. **This 10c N=3
  series is the reference reused as the N=3 side of the Milestone 11 comparison**
  (`completed_with_fallback_diagnostic`).
- **M10 Final** reads only the 10a/10b/10c artifacts. Under **both** `TOTAL_GAMMA = 1.5` and `3.0`, the
  same finite-condition rankings held:

| Metric | Ranking (both TOTAL_GAMMA = 1.5 and 3.0) |
|---|---|
| `W_max` | N=7 > N=5 > N=3 |
| `W(t=10)` | N=7 > N=5 > N=3 |
| `usable_fraction` | N=7 > N=5 > N=3 |
| `W_time_area` | N=3 > N=5 > N=7 |
| ergotropy arrival | N=3 fastest, N=5 middle, N=7 slowest |
| XGamma | N=7 > N=5 > N=3 |

The reading is deliberately **not** "longer chains are better": under fixed total dephasing the metric
rankings differ (longer chains have higher peak/final ergotropy and usable fraction; shorter chains
have larger `W` time-area and earlier arrival). This is a finite-condition descriptive result, not a
general law. `W_time_area`/`E_time_area` are time-integrals of state quantities, not cumulative work or
input; the two `TOTAL_GAMMA` points are not used to infer any functional form.

XGamma is a dephasing-kernel-weighted coherence-exposure **diagnostic**,
`x_gamma(t) = sum_ab Gamma[a,b] |rho[a,b](t)|^2`, time-integrated — **not** lost work, dissipated
energy, dephasing power, heat, entropy production, efficiency, or damage.

**Why Milestone 11.** In Milestone 10 the total dephasing was matched but the **integrated drive input
energy** of N=3 vs N=7 was not, so ranking differences alone cannot isolate the effect of chain length.
Milestone 11 fixes that. Note the Milestone 11 matched N=7 trajectory (`Omega ~ 0.18748`) is a
*different* trajectory from the Milestone 10c N=7 run (`Omega = 0.2`); the two must not be conflated.

### Milestone 11: equal-input N=3 vs N=7 analysis

**Equal-input matching (11c-11f).** The N=7 drive strength `Omega` is chosen so its integrated drive
input energy equals the N=3 reference (`target_E_drive_in = 5.9618618770136536e-2`), by an auditable,
one-trajectory-at-a-time procedure: a weak-drive initial guess (11c, no new evolution), an upper probe
(11d, `lower_probe_required`), a lower probe confirming a bracket `[0.18700000, 0.18770762191709489]`
(11e, `local_input_bracket_confirmed`), and evaluation of the interpolated point (11f). At
**`Omega = 0.18748395731510084`** the measured input matches the target to **relative mismatch
`2.8157011099880636e-6`**, with passing numerical-quality checks
(`matched_input_found_with_fallback_diagnostic`). This matches **drive input energy only**; chain
length, dimension, bond count, and geometry still differ, and no unique root, global uniqueness,
monotonicity, or `dt` convergence of the matching is claimed.

**Shape at equal input (11g).** Comparing the two matched trajectories (N=3 from the 10c series, N=7
from 11f): N=3 `W_max = 3.0302e-3` at t=5.63; N=7 `W_max = 3.3853e-3` at t=7.70 (ratio ~ **1.11717**).
N=7 has the **higher, later, narrower** peak but a **smaller** 0->10 `W` time-area (the `W` curves cross
twice). `W_time_area` is a state-quantity time-area, not cumulative extracted work; the higher N=7
usable fraction is a descriptive ratio, not an independent performance gain.

**Descriptive shape fits (11h, 11j, 11k).** These describe how much of the N=7 curve can be *described*
by transforming the N=3 curve; they are not mechanism proofs, and `AIC-like`/`BIC-like` are descriptive
complexity comparisons only.

- **11h** single amplitude/shift/scale transform: normalized RMSE ~ **0.05485**
  (`A=1.02343645, delta=2.49, s=0.91`), leaving structured residuals.
- **11j** asymmetric rise/fall transform: normalized RMSE ~ **0.03899**
  (`A=1.05674517, delta=2.32, s_rise=0.965, s_fall=0.515`); `asymmetric_time_scaling_partially_supported`.
- **11k** post-peak decay: best low-parameter model is a two-stage exponential, normalized RMSE ~
  **0.00810** (`switch=8.50, lambda_1=0.15, lambda_2=0.20`); the remaining 11j residual is
  concentrated late (**~59.6%** of its absolute area in t=9->10). Final judgment
  **`late_tail_structure_remains`**; auxiliary (internal to the decay fit)
  `two_stage_tail_decay_supported`.

**Design and audit (11i).** Milestone 11i is a **design and audit**, not a measurement. It confirmed
that up to that point only aggregate load quantities were stored (no site-resolved populations/currents,
full density matrices, mutual information, negativity, or mode occupation), and specified a minimal
follow-up: re-run the single matched N=7 trajectory once, saving site-resolved diagnostics
(`completed_design_with_targeted_recomputation_required`). 11i itself introduced no new time evolution.
The minimal site-resolved diagnostic it specified was subsequently executed in Milestone 11M (below);
11i remains the **design**, 11M is its **single-trajectory execution**. The candidate mechanisms
themselves (group velocity, boundary reflection, mode beating, entanglement/correlation fronts) remain
undetermined — 11M records time structure descriptively and does not establish any of them.

**dt-halving convergence of the matched N=7 condition (11L).** For the equal-input matched N=7 condition
(`Omega = 0.18748395731510084`, fixed, **not** re-matched), the internal time step was halved from
`dt=0.0025` to `dt=0.00125` (8000 RK4 steps, t=0..10, 1001 saved points; one new full trajectory). The
11f formal coarse condition was not recomputed and its SHA-256 matched before and after. Two separate
verdicts are reported and kept distinct:

- **Matching preservation: PASS.** With `Omega` held fixed, fine `E_drive_in = 5.9618450901723317e-2`,
  relative mismatch `2.8157045010698592e-6`; fine-vs-coarse input change `-2.0217e-13` (relative
  `3.3911e-12`). No re-matching was performed and **no new matching root is claimed**.
- **Physical convergence: PASS.** Scalar metrics **23/23 PASS** (e.g. W_max coarse/fine
  `3.3852501213e-3`/`3.3852501213e-3`, W(t=10) `2.3264679067e-3` both, XGamma `5.3283047839e-2` vs
  `5.3283047838e-2`); same-time trajectory quantities **10/10 PASS**; W-peak shape PASS (peak time 7.70,
  half-max width 3.85, both grids). Numerical quality: finite PASS, max trace error `3.775e-15`,
  Hermiticity error 0, `solver_failure = 0`, with one robust-positivity fallback point
  (primary 1000 success / 1 failure, fallback 1/1).

Final judgment **`matched_condition_dt_halving_passed_with_fallback_diagnostic`**. This is convergence of
**this one matched N=7 trajectory** under a single dt-halving; it does not claim convergence at arbitrary
dt, for other matched conditions, for a fine-dt-specific matching root, or for any other N/Omega/noise.

**Site-resolved transport and load-local change (11M).** Executing the 11i design on a single new
trajectory (N=7, `TOTAL_GAMMA=1.5`, `Omega = 0.18748395731510084`, `dt=0.0025`, t=0..10 — no fine-dt
re-run, since 11L already confirmed dt convergence), 11M saves per-site populations
`n_j = <sigma_j^+ sigma_j^->`, nearest-neighbor bond currents (`I_{j->j+1} = -2 J Im(z_j)`, positive =
drive->load), a site7->load **energy** current `I_7_to_load = <i[H_interaction, H_load]>`, and the load's
E / W / passive energy, all on the same 1001 times. Continuity and unit checks PASS (max continuity
residual `1.373e-5`); the 8 aggregate 11f quantities reproduce (8/8 PASS, max abs difference 0);
formal-input SHA-256 matched before and after. Observed time structure around the W peak:

| Event | Time |
|---|---|
| W peak | t = 7.70 (W_peak `3.3852501213e-3`) |
| site7->load current, decline onset | t = 7.5 |
| site7->load current, zero crossing after peak | t ~ 8.15 |
| site7->load current, persistent negative onset | t = 8.16 |
| load passive-energy increase onset | t = 7.70 |

From the W peak to t=10: load `Delta E = -1.86345075e-4`, `Delta W = -1.05878221e-3`,
`Delta passive = +8.72437140e-4`. Final classification
**`mixed_transport_and_load_local_change_with_fallback_diagnostic`**. Both the transport-side signature
(current decline / zero crossing / persistent negative direction) and the load-local signature (W
decrease with concurrent passive increase) are logged with `status = observed` only. **No causal claim is
made**: the current zero crossing is not called a proof of boundary reflection, the current change is not
said to cause the W decrease, and the passive-energy increase is not called thermalization or
decoherence. Not examined: causality, boundary reflection, group velocity, modes, mutual information,
negativity, entanglement, t>10, N>7, other Omega.

### Finite-chain (N=3..7) event-structure diagnostics

Two follow-up analyses read only stored / formally adopted trajectories and add no new mechanism claims.

**N=3 W-peak branch check.** Testing whether the noise-dependent N=7 W-peak candidate generalizes to
N=3, only the three stored N=3, `Omega=0.2` trajectories (eta = 0, 1.5, 3.0) were read — no new time
evolution, smoothing, or fitting. Under fixed extraction rules, the resolved peak sits in different time
windows across the branches: eta=0 at **t=9.48** (late window), eta=1.5 at **t=5.63** and eta=3.0 at
**t=5.39** (early window). They cannot be aligned as a single time-window peak branch. Verdict
**`peak_branch_correspondence_incomplete`** (secondary `none`): neither an early- nor a late-window
correspondence is supported on the stored grid, and the incompleteness comes from absent local peaks in
the predefined windows, not from grid resolution. The earlier N=7 preregistered W-peak prediction
(post-11k experiment 3) is **retained as-is**; it is **not** generalized into an N-common law, and no
claim is made about eigenmodes, phase transitions, exponential-law validity, or causal noise selection.

**Equal-input N=3..7 event comparison.** N=3, 4, 5, 6, 7 are compared at equal drive input energy, each
using its own matched `Omega` (N=7 reuses the 11M trajectory; N=8 is **not** included). Per-chain
categorical signatures:

| N | parity | t_W_peak_end | current negative onset | signature |
|---:|---|---:|---:|---|
| 3 | odd | 5.63 | not observed | `negative_not_observed` |
| 4 | even | 6.04 | 8.69 | `W->passive->negative` |
| 5 | odd | 6.62 | not observed | `negative_not_observed` |
| 6 | even | 7.16 | 7.72 | `W->passive->negative` |
| 7 | odd | 7.70 | 8.16 | `W->passive->negative` |

The formal artifact classification is **`odd_even_alternation_candidate`**, driven by the *continuous*
side of the precommitment: three quantities (`backflow_amount_post_peak`, `post_peak_W_loss_fraction`,
`Delta_passive_peak_to_t10`) show non-zero alternating local differences `A_4, A_5, A_6`. However the
**categorical** odd-even condition is **not** met: N=7 (odd) shares the even-group `W->passive->negative`
signature rather than the `negative_not_observed` signature of N=3 and N=5, so N=7 breaks the odd-group
commonality. `t_W_peak_end` also varies smoothly with chain length rather than alternating. The careful
statement is therefore: a simple odd/even rule was **not** confirmed; only some continuous quantities
show local alternating differences. This is **not** a parity mechanism, a universal parity order
parameter, or a statistically significant difference. Numerical quality across N=3..7: 1001 points each,
state/solver nonfinite 0, max continuity residual `1.373e-5` (< 5e-4), max trace error `2.220e-15`, max
energy-ledger residual `5.197e-7`; N=7 carries one robust-positivity fallback point (consistent with the
reused 11M trajectory).

### Post-11k mini-experiments: input dependence on stored trajectories

Building on the stored time series and the tail structure obtained up to Milestone 11k — and
keeping mechanism hypotheses separate — three small questions were examined on the *existing*
stored trajectories before any larger computation. They propose no new theory or mechanism. Only
experiment 3 ran a new trajectory (a single eta=0.75 run); experiments 1 and 2 reuse stored
trajectories with no new time evolution. All results are observations for this finite model and
these finite conditions.

- **Event order (Case B).** Across total-noise inputs eta = 0, 1.5, 3.0, the per-metric event
  order did not agree (the main difference is the position of the usable event). A simple
  input-independent common event order is not supported. Stop here for this question.
- **Coherence vs W (Case D, not comparable).** The three trajectories share no common absolute W
  range (common lower bound `1.3085e-2` exceeds common upper bound `6.9917e-4`), so a raw C(W)
  comparison at equal W cannot be performed. This is neither support nor non-support of a common
  C(W) relation; "not comparable" is not a negative result.
- **W-peak vs total phase noise (Case A).** Using only the existing eta = 0/1.5/3.0 points, an
  exponential-decay candidate `W_peak(eta) = A exp(-k eta)` and its prediction were fixed in
  advance (A, k, thresholds, and eta = 0.75 committed in `PRECOMMIT.md`). A single new eta = 0.75
  trajectory was then run (1001 saved points, numerical checks passed,
  `completed_with_fallback_diagnostic`): predicted `9.3038620128203936e-3`, observed
  `9.1509856038800939e-3`, relative error `1.643150%`, verdict
  `phase_noise_W_peak_exponential_candidate_retained`. The result is agreement with one added point
  evaluated after precommitment, within the pre-specified threshold; the exponential-decay
  candidate is provisionally retained. This is not establishment of an exponential law. The 3-point
  R-squared was not used for the verdict; A, k, and thresholds were not changed after seeing
  eta = 0.75; no second eta and no alternative model were run.

**What can be said:** for this finite model and these finite conditions, event order is
input-dependent (Case B); the stored trajectories are not comparable at equal W (Case D); and one
preregistered added point is consistent, within the pre-specified threshold, with a fixed
exponential-decay candidate that is provisionally retained (Case A).

**Not claimed:** no universal law, no new quantum channel, no discovered causal mechanism, no
claim that coherence determines W, no proof of an exponential law, no general scaling law.
`W_time_area` remains a state-quantity time-area, not cumulative extracted work; usable fraction
and passive energy are not treated as independent performance quantities.

---

## Current main results

Representative values only; all are directly-observed results for this model, these conditions, this
finite time, and this implementation — not general laws. Full tables live in the reports and CSVs.

- **M5c (N=3, conditional load-energy match, t=10):** ergotropy A/B = `6.373` (A `5.2798e-2` /
  B `8.2846e-3`); load-energy relative difference `4.001e-5`. Omega and total input energy not matched.
- **M9a/9b (N=7, fixed-per-site, reference):** noise-free W_max `2.2436e-2` (t=7.71); all-site-noisy
  (gamma=0.5) W_max `4.0437e-4` (t=7.65); fixed-per-site W_max N3 `3.0302e-3` > N5 `1.0968e-3` >
  N7 `4.0437e-4` (noise-free "N7 > N5" does not survive here).
- **M9c validation (N=7, fixed-total, total gamma=1.5 — source of truth):** positivity 1001/1001,
  solver_failure 0; SymmetricEigen 999 success / 2 failure; Schur fallback 2/2; max difference vs
  existing 9c trajectory `0`; N=7 W_max **`0.0038080717406769921`** (t=7.70); fixed-total W_max
  N3 `3.0302e-3`, N5 `3.4876e-3`, N7 `3.8081e-3`; N7/N5 W_max ratio **`1.0918825139`**, within the
  descriptive "10% band" (not a statistically significant difference, not distance-only causation).
- **M11g (N=3 vs N=7, equal drive input energy — distinct from the fixed-total M9c run above):**
  matched at `Omega=0.18748395731510084` (relative input mismatch `2.8157e-6`); N=7 `W_max 3.3853e-3`
  at t=7.70 vs N=3 `3.0302e-3` at t=5.63 (ratio ~ `1.11717`); N=7 higher/later/narrower peak but
  smaller `W` time-area.
- **M11L (matched N=7, dt-halving of the same trajectory):** `Omega` fixed at `0.18748395731510084`;
  `dt` halved `0.0025 -> 0.00125`. Matching preserved (relative mismatch `2.8157045e-6`, input change
  `~3.39e-12` relative); scalar 23/23 PASS, trajectory 10/10 PASS, W-peak shape PASS (peak t=7.70,
  W_max `3.3853e-3` both grids); `solver_failure = 0`, one positivity fallback point.
  `matched_condition_dt_halving_passed_with_fallback_diagnostic`.
- **M11M (matched N=7, site-resolved single-trajectory diagnostic):** W peak t=7.70; site7->load energy
  current zero crossing after peak ~ t=8.15, persistent negative onset t=8.16; load passive-increase
  onset t=7.70. Peak->t10: load `Delta E = -1.86345075e-4`, `Delta W = -1.05878221e-3`,
  `Delta passive = +8.72437140e-4`. `mixed_transport_and_load_local_change_with_fallback_diagnostic`
  (transport-side and load-local signatures both `observed`, no causal claim).
- **N=3 W-peak branch check (stored N=3, Omega=0.2):** resolved peaks at eta=0 t=9.48, eta=1.5 t=5.63,
  eta=3.0 t=5.39 — not alignable as one time-window branch; `peak_branch_correspondence_incomplete`.
- **Equal-input N=3..7 event comparison:** signatures N3 `negative_not_observed`, N4
  `W->passive->negative`, N5 `negative_not_observed`, N6 `W->passive->negative`, N7
  `W->passive->negative`. Formal `odd_even_alternation_candidate`, but the categorical odd/even condition
  is **not** met (N=7 breaks the odd-group commonality); only some continuous quantities alternate. Not a
  parity mechanism, universal rule, or statistically significant difference. N=8 excluded.

---

## Directly validated findings

- The equal-input match of N=7 to N=3 at `Omega=0.18748395731510084` (relative mismatch `2.8157e-6`),
  with passing numerical-quality checks (M11f).
- Peak values/timing N=3 vs N=7 and the two `W` crossings; N=7 higher/later/narrower peak, smaller
  0->10 `W` time-area (M11g).
- The 9c fixed-total positivity accounting (1001/1001, fallback 2/2, solver_failure 0) and the max-0
  agreement with the existing 9c trajectory (M9c validation).
- Milestone 10 fixed-total metric rankings under both `TOTAL_GAMMA = 1.5` and `3.0` (M10 Final).
- Numerical-quality checks (per-run) and, for 11j/11k, SHA-256 input-integrity checks; the 11i checks
  are design/audit checks.
- dt-halving (`0.0025 -> 0.00125`) of the *same* matched N=7 trajectory with `Omega` fixed: matching
  preserved and scalar 23/23 / trajectory 10/10 / W-peak shape all PASS, `solver_failure = 0` (M11L).
- Site-resolved single-trajectory diagnostic for the matched N=7 run: per-site populations, bond
  currents, site7->load energy current, nearest-neighbor coherences, and load E/W/passive on 1001 times,
  with the W-peak-relative time ordering above; continuity/unit checks PASS, 11f aggregates reproduced
  8/8 (M11M).
- N=3 stored-trajectory W-peak branch windows (eta=0 t=9.48; eta=1.5 t=5.63; eta=3.0 t=5.39) and their
  non-correspondence as a single time-window branch (`peak_branch_correspondence_incomplete`).
- Equal-input N=3..7 per-chain categorical signatures and numerical-quality checks; the categorical
  odd/even condition is not met (N=7 breaks the odd-group commonality).

## Descriptive model support

- Single amplitude/shift/scale transform, normalized RMSE ~ `0.05485` (M11h).
- Asymmetric rise/fall transform, normalized RMSE ~ `0.03899`,
  `asymmetric_time_scaling_partially_supported` (M11j).
- Post-peak two-stage exponential, normalized RMSE ~ `0.00810`, `two_stage_tail_decay_supported`;
  late-tail concentration ~ `59.6%` in t=9->10, `late_tail_structure_remains` (M11k).

These are descriptive shape diagnostics, not evidence of any physical mechanism.

---

## What has not been confirmed

- Quantum advantage/supremacy or superiority over a classical method.
- Universal scaling laws, exponential/power-law decay, thermodynamic limits, N->infinity behavior.
- Any general advantage of longer chains, or universal laws of noise position/protection strength.
- Any physical mechanism behind the observed shapes — group velocity, boundary reflection, mode
  beating, entanglement/correlation fronts, or any causal relationship — and that transform/decay fits
  correspond to specific physical processes.
- Control cost of extraction, finite-time switch-off/extraction, repeated cycles, continuous supply,
  steady output, long-time stability.
- `dt`-halving convergence *beyond the single matched N=7 trajectory checked in 11L* (arbitrary dt, other
  matched conditions, a fine-dt-specific matching root); matching-root uniqueness/monotonicity.
- Equal-input results for `TOTAL_GAMMA=3.0`, for other `Omega` roots, for N>7, or for t>10; N=8 is not
  included in the formal chain-length comparison; other external solver crates (out of 9c-validation
  scope). (Equal-input N=4/5/6 now exist only as finite-chain event-structure diagnostics, not as a
  scaling law.)
- Any parity mechanism, universal odd/even effect, parity order parameter, or statistically significant
  parity difference; the N=3..7 `odd_even_alternation_candidate` is a descriptive finite-size candidate
  only, and N=7 breaks the categorical odd-group commonality.
- Any causal reading of the 11M site-resolved diagnostic: that the current change causes the W decrease,
  that the current zero crossing / reversal proves boundary reflection, or that the passive-energy
  increase proves thermalization or decoherence.
- Fair baselines against classical wave/stochastic models; novelty/priority vs the literature (no
  literature review done).

---

## Reproduction and verification

Requires a Rust toolchain (`Cargo.toml`: edition 2021; deps `nalgebra`, `num-complex`, `thiserror`;
dev-dep `approx`).

```bash
cargo fmt --all -- --check
cargo test --release --offline
cargo build --release --offline
# opt-in 24-dimensional smoke test:
cargo test --release full_24d_short_time_smoke_test -- --ignored --nocapture
```

Representative milestone binaries:

```bash
cargo run --release --offline --bin time_dependent_sanity          # M5a
cargo run --release --offline --bin dephasing_kernel_benchmark     # M8c
cargo run --release --offline --bin n7_noise_free_full             # M9a
cargo run --release --offline --bin n7_all_site_noisy_full         # M9b
cargo run --release --offline --bin fixed_total_noise_comparison   # M9c
cargo run --release --offline --bin n7_t002_eigen_diagnostic       # M9c diagnostic
cargo run --release --offline --bin n7_fixed_total_validation      # M9c validation (source of truth)
```

The N=7 full runs take on the order of tens of minutes (e.g. 9a ~ 2953s, 9b ~ 2899s, 9c validation ~
2703s). Per-milestone reports record `cargo fmt` PASS, the analysis binary PASS, and release-test
counts that grow across the line (e.g. M2 26 tests; M5a `47 passed / 0 failed / 1 ignored`; M10b 101,
M10c 104, M10 Final 107; M11c 110, M11d 113, M11e 116, M11f and M11h-M11k 119, each `0 failed / 1
ignored`). The 9c-validation runtime checks are 47/47 PASS and the robust eigenvalue diagnostic's
required checks are 10/10 PASS. 11j and 11k additionally verify their inputs by SHA-256 before and
after analysis.

---

## Repository layout

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
    n7_t002_eigen_diagnostic.rs     # M9c diagnostic
    n7_fixed_total_validation.rs    # M9c validation (final source of truth)
    matched_n7_dt_halving_validation.rs    # M11L
    matched_n7_site_resolved_transport.rs  # M11M
tests/
  full_24d_short_time.rs            # M2.1 (ignored smoke test)
MILESTONE_*.md                      # per-stage reports
*.csv                               # per-stage outputs
```

The `src/lib.rs` module declarations are based on the actual files; the `bin/` mapping and the Milestone
10-11 binaries are organized from the "generated files" sections of the respective reports. For a
complete list of individual CSVs, see the end of each MILESTONE report.

---

## Detailed reports and data

Minimal pointers for the load-bearing results (each report lists its own full CSV set):

- **Foundations & central comparison:** `MILESTONE_4_RESULT.md`, `MILESTONE_5A/5B/5C`,
  `MILESTONE_6A_REPORT.md`; `Milestone_4-6a_研究結果ノート.pdf` summarizes 4-6a.
- **Noise location / protection:** `MILESTONE_7A/7B/7C/7D_REPORT.md`.
- **Chain length / feasibility:** `MILESTONE_8A/8B/8C_REPORT.md`.
- **N=7 full runs:** `MILESTONE_9A/9B_REPORT.md`.
- **9c (source of truth + records):** `MILESTONE_9C_VALIDATION.md` (final),
  `MILESTONE_9C_DIAGNOSTIC.md`, `MILESTONE_9C_REPORT.md`; key CSVs
  `n7_fixed_total_validation_summary.csv`, `..._checks.csv`, `..._trajectory_comparison.csv`,
  `fixed_total_noise_final_comparison.csv`, `robust_eigen_diagnostic_unit_checks.csv`.
- **Fixed-total / XGamma:** `MILESTONE_10A/10B/10C_REPORT.md`, `MILESTONE_10_FINAL_REPORT.md`.
- **Equal-input N=3 vs N=7:** `MILESTONE_11C_PRECHECK_REPORT.md`, `MILESTONE_11D/11E_REPORT.md`,
  `MILESTONE_11F/11G/11H_REPORT.md`, `MILESTONE_11I_MOUNTAIN_MECHANISM_DESIGN.md`,
  `MILESTONE_11J_ASYMMETRIC_TRANSFORM_REPORT.md`, `MILESTONE_11K_POST_PEAK_DECAY_REPORT.md`; key CSVs
  `input_matching_interpolated_trial_summary.csv`, `equal_input_timeseries_shape_summary.csv`,
  `equal_input_peak_widths.csv`, `equal_input_curve_transform_models.csv`,
  `equal_input_asymmetric_transform_models.csv`, `equal_input_post_peak_decay_models.csv`.
- **Matched N=7 dt-halving (11L):** `MILESTONE_11L_MATCHED_N7_DT_HALVING_REPORT.md` (scalar/trajectory
  comparison CSVs referenced therein).
- **Matched N=7 site-resolved transport (11M):** `MILESTONE_11M_SITE_RESOLVED_TRANSPORT_REPORT.md`; key
  CSV `matched_n7_site_resolved_events.csv`.
- **N=3 W-peak branch check:** `MILESTONE_N3_W_PEAK_BRANCH_REPORT.md`; key CSV
  `n3_W_peak_branch_comparison.csv`.
- **Equal-input N=3..7 event comparison:** `MILESTONE_N3_7_EVENT_COMPARISON_REPORT.md` (provided as
  `FINAL_REPORT.md`); key CSVs `event_summary.csv`, `odd_even_temporal_comparison.csv`,
  `numerical_checks.csv`.
- **Post-11k mini-experiments** (paths TODO — confirm folder layout before linking):
  - Event order: `qwn_event_order_mini_REPORT.md` (CSV `event_orders.csv`).
  - Coherence vs W: `qwn_coherence_W_response_REPORT.md` (CSV `common_W_range.csv`).
  - W-peak prediction: `qwn_phase_noise_W_peak_prediction_REPORT.md` and `PRECOMMIT.md`
    (CSV `prediction_test_result.csv`).

Reading conventions for the CSVs: `W_time_area`/`E_time_area` are time-areas of state quantities (not
cumulative work or input); `W/Ein` is not an overall efficiency; `usable_fraction` is ergotropy over
load energy; XGamma is a diagnostic, not a loss/efficiency.

---

## Items requiring confirmation

No numerical or judgment contradictions were found between the two source READMEs when merging;
Milestone 10 appeared in both and has been unified. The following points could not be resolved from the
provided materials and are left open rather than guessed:

- **Whole-suite Cargo test count "90/90".** The 9c-validation runtime checks are 47/47 PASS and the
  robust eigenvalue diagnostic's required checks are 10/10 PASS, but a total of "90/90" referenced
  elsewhere for the entire Cargo test suite cannot be confirmed from the materials and is not adopted
  here.
- **Milestone 11g release-test count.** The 11g report records formatting and analysis execution as
  PASS but does not state a numerical release-test count.
- **Two distinct N=7 W_max values are not a contradiction.** `3.8081e-3` (M9c fixed-total,
  `Omega=0.2`) and `3.3853e-3` (M11g equal-input, `Omega~0.18748`) belong to *different* trajectories
  and conditions; readers should not conflate them. Noted here to avoid misreading.

---

## Citation and license

License: MIT (see `Cargo.toml` / `LICENSE`). When citing, please reference this repository and the
specific `MILESTONE_*.md` report(s) underlying the result used, since each report defines the exact
conditions and reservations for its numbers.
