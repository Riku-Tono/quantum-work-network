# Quantum Work Network

A Rust reference implementation for numerically computing how much of the energy driven through a
small quantum network arrives at a receiver (the *load*) as **extractable work (ergotropy)**, and how
that changes under phase noise. The project is developed as a single, incremental research and
implementation history from Milestone 1 through Milestone 15A, via finite-chain (N=3..7)
event-structure diagnostics, a descriptive decomposition of the load's coherent ergotropy (M12), a
short-horizon counterfactual branch experiment (M13), coupling-control switching-work accounting (M14),
and a single end-point ideal local extraction (M15A). After M15A, a supplementary comparison extended
the C2 protocol to t=10 and compared it against a long-window C0 baseline on the common saved grid; it
is not assigned a new Milestone number.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21435208.svg)](https://doi.org/10.5281/zenodo.21435208)
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

### Milestone 12: what the post-peak W decrease looks like (fixed matched N=7)

Milestone 12 keeps the M11M matched N=7 condition fixed (`TOTAL_GAMMA=1.5`, `Omega =
0.18748395731510084`, `dt=0.0025`, t=0..10, 1001 saved points, vacuum initial state) and asks three
*descriptive* questions about the post-peak decrease of the load's extractable work: which component
changes (12A), in what temporal order (12B), and how load-chain correlation behaves around it (12C).
No matching, re-optimization, root search, fine-dt run, or other N/gamma/Omega was performed.

**12A — local ergotropy decomposition.** One new trajectory under the same fixed condition. Using the
existing energy-basis convention (`H_L = diag(0,1,2)`), the load ergotropy is split as
`W = W_pop + W_coh`, with `W_pop = W(Delta_H(rho_L))` the diagonal (population) part. At **every** saved
time `W_pop = 0` and therefore `W = W_coh`; over `t = 7.70 -> 10.00`,
`Delta W = Delta W_coh = -1.0587822145716654e-3`, `Delta W_pop = 0`,
`Delta E_load = -1.8634507496082930e-4`, `Delta E_passive = 8.7243713961083613e-4`. The observed signed
W decrease is thus *described* entirely as a `W_coh` decrease. The pre-registered ratio rule, however,
required **both** components to decrease positively, and `L_pop = 0` made the ratio undefined, so the
formal classification is **`component_change_inconclusive`**. That verdict reflects the pre-committed
rule being inapplicable — it does **not** mean the decomposition failed; the decomposition itself is
exact (`max |W - (W_pop + W_coh)| = 0`). Numerical quality **20/20 PASS**, M11M aggregates reproduced
**8/8 PASS** (max abs difference 0), formal input SHA-256 7/7 unchanged. Population and coherent parts
are not independent conserved quantities or separate physical energies.

**12B — temporal relation (no new time evolution).** Reading only the 12A saved 1001-point series under
a pre-registered five-interval persistence rule (`M12B_PRECOMMIT.md`; no smoothing, fitting, FFT, or
extremum interpolation), the saved-grid order is:

| Event | Time |
|---|---|
| site7->load current decline onset | t = 7.50 |
| `W_coh` peak / persistent decline | t = 7.70 / 7.70 |
| `Cl1` (coherence l1-norm) peak / persistent decline | t = 7.73 / 7.73 |
| current zero crossing | t = 8.15 |
| current persistent negative onset | t = 8.16 |
| current persistent positive return | t = 9.68 |

Formal classification **`coherence_and_ergotropy_decline_before_negative_current`** (3/5/7-interval
sensitivity all agree). The persistent negative current begins **0.46 later** than the `W_coh` decline
onset, so the simple account "the turn to negative current initiated the `W_coh` decrease" is **not
supported by the time ordering**. This is a statement about ordering only; **no causal relation between
current and `W_coh` is established** in either direction. Checks 15/15 PASS, `W = W_coh` max difference
0, input SHA-256 3/3 unchanged. The current integrals are time-integrals of a state-derived quantity,
not work, efficiency, loss, dissipation, or heat.

**12C — load-chain correlation.** One new full trajectory under the same fixed condition, evaluating
load-chain mutual information `I(C:L) = S_chain + S_load - S_total` (natural log, nats), the individual
entropies, and purities at all 1001 saved times. The MI global peak and its five-interval persistent
decline both fall at **t = 7.51**, with peak value **`4.5985207125390781e-3` nats** — **0.19 earlier**
than the `W_coh` peak/decline at t=7.70. Formal classification
**`mutual_information_peaks_before_Wcoh_decline`**, unchanged across the 3/5/7-interval sensitivity
grid. Numerical quality **29/29 PASS**; existing formal artifacts SHA-256 **30/30 unchanged**; M12A
common quantities reproduced to max abs difference 0. Mutual information is a **total-correlation**
observable including classical and quantum correlation — it is **not** work, ergotropy, or
entanglement. No claim is made that the MI decrease caused the `W_coh` decrease, or that `W` was
converted into MI.

### Milestone 13: short-horizon counterfactual branches and their temporal resolution limit

**13A — formal counterfactual branch experiment.** Starting from the *same* stored total-system state
as M12C and changing only the **generator after the branch time**, short branches were propagated on a
0.01 saved grid to t=8.50 from two branch starts (t=7.50 and t=7.70; 728 formal saved rows total):

- **B0** baseline continuation (reproduces M12C exactly, max abs difference 0);
- **B1** chain dephasing set to 0 after the branch;
- **B3** chain-load coupling halved after the branch (`g: 0.25 -> 0.125`);
- **B4** coupling set to 0, a boundary control (max invariant drift `1.9984014443252818e-15`).

**B3 gives the same formal classification at both branch starts:
`weaker_coupling_weakens_Wcoh_decline`.**

| Branch start | endpoint `W_coh` relative difference | `W_coh` time-area relative difference |
|---|---:|---:|
| 7.50 | +6.1424321151409798 % | +1.4670925175957812 % |
| 7.70 | +6.1471886834895455 % | +2.1027954259932112 % |

B1 differs between starts and was **not** averaged: `dephasing_counterfactual_mixed` at 7.50,
`future_dephasing_change_has_small_effect` at 7.70. B4 is
`zero_coupling_boundary_control_pass` at both starts. All numerical-quality rows PASS (max trace error
`1.776e-15`, Hermiticity error 0, no fallback used, drive after branch 0).

The admissible summary is exactly: **from the same state, halving the chain-load coupling after the
branch produced a consistently weaker `W_coh` decrease within the finite window, at both branch starts,
relative to the baseline that keeps the normal coupling.** Explicitly **not** claimed: that coupling is
the sole cause of the `W_coh` decrease; that weaker coupling is always better; that an optimal coupling
was found; that halving the coupling was shown to protect `W_coh` by some mechanism. The B1 and B3
percentages are **not** compared head-to-head to rank which knob "matters more". `W_coh` time-area is a
time-integral of a state quantity, not cumulative extracted work.

**13B — temporal response decomposition (no new computation).** Using only the formal 13A saved CSVs
(new trajectories 0, new branches 0, new parameter conditions 0), 13B asks which saved observable first
develops a persistent B3-minus-B0 difference. Under the frozen five-interval rule:

- **start 7.50:** current at t=7.50, then all five other quantities (`Cl1`, MI, `W_coh`, load energy,
  passive energy) together at t=7.51;
- **start 7.70:** current at t=7.70, then MI / load energy / passive energy / `W_coh` at t=7.71, then
  `Cl1` at t=7.75.

However the three-interval rule places the main state quantities simultaneous with `W_coh`, while the
five- and seven-interval rules produce a start-time dependence; across the full 3/5/7 x epsilon
(x0.5/x1/x2) grid the classification is not stable (`Wcoh_and_other_differences_emerge_together` vs
`branch_response_order_depends_on_start_time`). The formal classification is therefore
**`temporal_response_inconclusive`**. **This is a resolution limit, not a computational failure, and it
is not a negation of the 13A B3 result**: the numerical audit is **35/35 PASS** with **11/11 SHA-256
unchanged** for the formal M13A inputs, and the M13A B3 classifications and relative differences
reproduce within 1e-12. The immediate current difference at the branch time is **structural** — the
current expression itself contains the coupling coefficient that was changed (audited as
`current_B3 = 0.5 * current_B0` at branch start) — and it is **not** read as the current causing the
later `W_coh` difference. Also not shown: that passive energy converted into coherent ergotropy, or
that correlation suppression protected work.

**Where this leaves the project.** M12 described the *component* (`W = W_coh`, 12A), the *temporal
order* (12B), and the *load-chain correlation* (12C) of the post-peak `W_coh` decrease. M13A examined
the finite-time response to a generator change from an identical state and found, at both branch
starts, that halving the coupling weakens the `W_coh` decrease. M13B could not uniquely resolve the
internal response ordering of that branch at the available saved resolution. **Coupling sensitivity is
therefore a confirmed observation, while the specific physical mechanism and any optimal condition
remain undetermined.**

### Milestone 14: coupling-control design and switching-work accounting

Milestone 14 asks what it costs, *within the model's Hamiltonian ledger*, to make the M13A coupling
change — and keeps that ledger strictly separate from the state ledger. Fixed protocol names used from
here on:

| Protocol | Definition | Corresponds to |
|---|---|---|
| **C0** | `g = 0.25` throughout | M13A B0 |
| **C1** | `g: 0.25 -> 0.125` at t=7.50, held to 8.50 | M13A B3, start 7.50 |
| **C2** | `g: 0.25 -> 0.125` at t=7.70, held to 8.50 | M13A B3, start 7.70 |

No `g=0`, additional coupling value, switch-back, ramp, feedback, or multiple switch is admitted.

**14A — design only (no formal control trajectory).** Milestone 14A is a **design and audit
Milestone**: it ran **no formal control trajectory, no coupling sweep, no optimization, and no ramp**.
It audited the implemented interaction term (`H_CL(g) = g V_CL`, `V_CL = sigma_plus_N b +
sigma_minus_N b_dagger`, exactly linear in `g`), confirmed the drive envelope is exactly zero for
t>tau=3.2 so that post-switch protocol differences are **not** drive-input differences, fixed the
ideal-quench convention `W_switch_on_system = Tr[rho(t_s)(H_after - H_before)]` (positive = supplied by
the external controller to the modeled system; state unchanged across the quench), and fixed **three
separate ledgers** — state quantities, switching work, and post-switch energy — which are **not**
summed or subtracted into any net quantity. Its no-propagation unit checks are **11/11 PASS**. The
reuse audit found that the M13A artifacts do not persist full density matrices or `Tr[rho V_CL]`, so
the verdict was `single_prefix_recomputation_required`: re-run one baseline prefix only. Final design
verdict **`coupling_control_experiment_ready`**, recommended plan **Plan B
(`minimal_state_recomputation_for_switching_accounting`)**. The state-level trade-off preview quoted
there is explicitly *not* a formal control verdict, because the switching work was still missing.

**14 Plan B — formal switching-work accounting (executed).** Exactly one baseline prefix was propagated
from t=0 to t=7.70, full density matrices were saved at t=7.50 and t=7.70, and the ideal same-state
switch was evaluated at both. New post-switch branches, coupling values, sweeps, ramps, extraction,
optimization, fine-dt runs, and t>8.50 propagation: **0**. The existing M13A C0/C1/C2 post-switch state
trajectories were reused unchanged.

| Protocol | switch time | `W_switch_on_system` | `W_switch_out` |
|---|---:|---:|---:|
| C1 | 7.50 | `1.6980822985791084e-15` | `-1.6980822985791084e-15` |
| C2 | 7.70 | `8.8624029002963674e-16` | `-8.8624029002963674e-16` |

Both ideal-quench energy jumps fall below the fixed `1e-12` audit tolerance, so the correct description
is **numerically zero at this tolerance** (status `numerically_zero_within_1e-12_audit_tolerance`).
**A numerically zero Hamiltonian-jump work does not mean a real controller is free**: the model has no
actuator Hamiltonian, finite ramp, bandwidth, or device-cost ledger. Reusing the M13A state values, the
formal composite classification for **both** C1 and C2 is
**`Wcoh_improves_with_load_energy_tradeoff + switching_accounting_available`** — endpoint `W_coh` and
`W_coh` area each improve by at least 1%, while endpoint load-energy retention stays below the fixed
99% threshold (C1: `W_coh` endpoint +6.142%, area +1.467%, retention 97.159%; **C2: `W_coh` endpoint
+6.147%, area +2.103%, retention 98.961%**). Numerical checks **22/22 PASS**; M13A B0 quantities
reproduce at both switch times to max abs difference `3.839e-13`; source/formal input hashes 15/15
unchanged. The reported switching-work ratios are scale diagnostics, **not efficiencies**, and the
switching ledger is not converted into a profit, a "net protected `W_coh`", or a control-efficiency
claim.

### Milestone 15A: end-point extraction accounting after coupling control

Milestone 15A propagates only C0 and C2 from 7.70 to 8.50 and performs **one** ideal local extraction
at the single fixed evaluation time **t = 8.50**. The formal comparison is **C2 - C0 only**. Verdict
**`state_improvement_preserved_after_extraction`**, mandatory numerical failures **0**.

| Quantity at t=8.50 | C0 | C2 | C2 - C0 |
|---|---:|---:|---:|
| `W_coh` before extraction | `3.0626774684866464e-3` | `3.2509460312412416e-3` | `1.8826856275459521e-4` |
| load energy before extraction | `6.601877630500e-3` | `6.533284638378e-3` | `-6.8592992121339208e-5` |
| `W_load_gross` | `3.0626774684866503e-3` | `3.2509460312412546e-3` | `1.8826856275460432e-4` |
| `W_operation_out_total` | `3.062677468536e-3` | `3.250946031265e-3` | `1.8826856272917414e-4` |
| interaction jump | `-4.957529045397e-14` | `-2.416320556056e-14` | `2.5412084893412970e-14` |
| load ergotropy after extraction | `0` | `0` | `0` |

The gross-work difference reproduces the pre-extraction `W_coh` difference with residual
`9.1072982488782372e-18` (ratio `1.000000000000`). Because the interaction-jump difference is only
`2.5412084893412970e-14`, the bare-load and total-Hamiltonian accountings give the **same sign** for the
comparison. Numerical and provenance checks **40/40 PASS**, including trace, Hermiticity, positivity,
post-extraction zero ergotropy, M13A endpoint reuse, M14 state SHA-256, and switching-accounting reuse.

The admissible conclusion is exactly: **at t=8.50, a single ideal local extraction on C2 yielded a
larger gross extracted work than the baseline C0; the C2-C0 gross-work difference agreed with the
pre-extraction `W_coh` difference to numerical precision; and including the interaction term in a
total-Hamiltonian accounting did not change the sign of the comparison.** In short, the state-level
improvement was **preserved** through this one ideal local extraction.

**Ledgers stay separate.** `W_coh`, load energy, switching work, `W_load_gross`, and
`W_operation_out_total` are **not** summed into any single "benefit". The M14 Plan B C2 switching work
(`8.8624029002963674e-16`) is carried over for reference only and is **not** combined with either
extraction-work column. Not performed in 15A: new coupling conditions, sweeps, ramps, optimization,
feedback, measurement, multiple extractions, other N/gamma/Omega, or any propagation beyond t=8.50.
Nothing here is optimal control, a free increase of work, an efficiency improvement, a confirmed net
benefit, a zero controller cost, an experimental-hardware advantage, or a quantum advantage.

### C2 extension to t=10: long-window C0 comparison

This supplementary comparison asks whether the C2-minus-C0 bare-load ergotropy difference observed at
t=8.50 is an artifact of that window's endpoint. It is not a new Milestone. The formally saved M14 Plan
B full-system state at t=7.70 was propagated **once and continuously** with C2 (`g=0.125`) to t=10.00 —
no extraction at t=8.50, no switch-back, no rematching, and no split formal trajectory. Fixed
conditions: N=7, `TOTAL_GAMMA=1.5`, per-site gamma `1.5/7`, `Omega=0.18748395731510084`, restart time
7.70, final time 10.00, internal `dt=0.0025`, saved interval 0.01, 920 internal steps, 231 saved
points, drive exactly 0. The C0 long-window baseline is the M12A CSV, fixed by an explicit artifact
audit. Primary classification **`C2_advantage_positive_at_t10`**.

| t | `W_C0` | `W_C2` | `Delta W` | relative difference |
|---:|---:|---:|---:|---:|
| 8.50 | `3.062677468487e-3` | `3.250946031241e-3` | `+1.882685627546e-4` | +6.147189% |
| 9.68 | `2.417441223901e-3` | `2.999355919229e-3` | `+5.819146953273e-4` | +24.071514% |
| 10.00 | `2.326467906728e-3` | `2.972492411029e-3` | `+6.460245043013e-4` | +27.768468% |

After t=8.50 the difference increased in **all 150 saved intervals** (decrease count 0, near-flat count
0), never reached near-zero, and never turned negative; its global minimum is at t=8.50 and its global
maximum at t=10.00, with no interior local extrema. **This is a statement about the C2-C0 difference,
not about C2 itself: C2 formed no new W peak after t=8.50**, its in-window maximum W remains at t=7.70
(`3.3852501212993710e-3`, the same value as C0's maximum), and the difference grows because C0 declines
faster — t=10 retention is 87.807172% for C2 versus 68.723664% for C0. Over 7.70..10.00 the ergotropy
time-areas are `A_C0 = 6.5655097929267292e-3`, `A_C2 = 7.2916111984482872e-3`,
`A_Delta = 7.2610140552155800e-4`, relative difference +11.059330%; these are time-areas of the
instantaneous bare-load ergotropy and are **not** cumulative extracted work. At t=10 the accompanying
energy differences are `Delta E_load = +2.106309466093e-5` and `Delta E_passive = -6.249614096404e-4`,
with energy-identity residual 0. The site7-to-load currents are recorded as observations under the
existing definition only and are **not** used to explain the ergotropy difference.

Audit: mandatory numerical checks **28/28 PASS**, existing-artifact reproduction checks **15/15 PASS**,
SHA audit PASS. The C0 short/long observable residual is 0, the C2 t=7.70..8.50 reproduction residual is
0, and the M15A t=8.50 values reproduce to a maximum residual of `9.107e-18`. Eigensolver fallback count
0; formal runtime 682.479 s.

Scope: this is a finite-condition, saved-grid comparison from one propagation of one stored state
against one fixed baseline. **No extraction was performed in this extension** — the only formal
extraction in the project remains the single ideal local extraction at t=8.50 in M15A. It does not
establish general weak-coupling or C2 superiority, an optimal coupling, switch time, or extraction time,
any causal mechanism, any conversion of passive energy into ergotropy, any net benefit combining
switching and extraction work, repeated or cyclic operation, steady-state output, behavior beyond t=10,
other N/gamma/Omega/g, efficiency improvement, or quantum advantage. Where the word *advantage* is used
here, it means only the sign of the bare-load ergotropy difference under these fixed finite conditions.

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
- **M12A (matched N=7, load ergotropy decomposition):** `W_pop = 0` at every saved time, so `W = W_coh`;
  over t=7.70->10.00 `Delta W = Delta W_coh = -1.0587822146e-3`, `Delta W_pop = 0`,
  `Delta E_load = -1.8634507496e-4`, `Delta E_passive = 8.7243713961e-4`. Formal
  `component_change_inconclusive` (the pre-committed both-components ratio rule was inapplicable, not a
  failed decomposition); 20/20 PASS, M11M aggregates 8/8 PASS.
- **M12B (same series, no new evolution):** current decline onset t=7.50; `W_coh` peak/decline t=7.70;
  `Cl1` peak/decline t=7.73; current zero crossing t=8.15; persistent negative onset t=8.16; persistent
  positive return t=9.68. `coherence_and_ergotropy_decline_before_negative_current`; the persistent
  negative current is 0.46 *later* than the `W_coh` decline onset.
- **M12C (same condition, one new trajectory):** load-chain mutual information peaks and begins its
  persistent decline at t=7.51, peak `4.5985207125e-3` nats — 0.19 earlier than the `W_coh` peak.
  `mutual_information_peaks_before_Wcoh_decline`; 29/29 PASS, 30/30 SHA-256 unchanged. MI is a
  total-correlation observable, not work or ergotropy.
- **M13A (short-horizon counterfactual branches from the M12C state):** halving the chain-load coupling
  after the branch (B3) gives `weaker_coupling_weakens_Wcoh_decline` at **both** starts — endpoint
  `W_coh` relative difference `+6.1424321151%` (start 7.50) and `+6.1471886835%` (start 7.70);
  `W_coh` time-area relative difference `+1.4670925176%` and `+2.1027954260%`. B1 start-dependent
  (`dephasing_counterfactual_mixed` / `future_dephasing_change_has_small_effect`); B4
  `zero_coupling_boundary_control_pass`, max invariant drift `1.998e-15`.
- **M13B (13A CSVs only, no new computation):** the B3-B0 persistent-difference ordering is not stable
  across the 3/5/7 x epsilon sensitivity grid — `temporal_response_inconclusive`; audit 35/35 PASS,
  11/11 SHA-256 unchanged, M13A values reproduced within 1e-12. A resolution limit, not a failure, and
  not a negation of the M13A B3 result.
- **M14A (design only):** no formal control trajectory, sweep, ramp, or optimization was run. Fixed the
  protocol set C0/C1/C2, the ideal-quench convention
  `W_switch_on_system = Tr[rho(t_s)(H_after-H_before)]`, and three separate ledgers; unit checks 11/11
  PASS. Verdict `coupling_control_experiment_ready`, recommended Plan B.
- **M14 Plan B (executed):** one baseline prefix to t=7.70; ideal same-state switch evaluated at both
  switch times. `W_switch_on_system` = `1.6980822985791084e-15` (C1) and `8.8624029002963674e-16` (C2)
  — **numerically zero within the 1e-12 audit tolerance**, which does **not** imply a real controller is
  free. Both protocols: `Wcoh_improves_with_load_energy_tradeoff + switching_accounting_available`
  (C2 endpoint `W_coh` relative difference +6.147%, endpoint load-energy retention 98.961%). Checks
  22/22 PASS; hashes 15/15 unchanged.
- **M15A (executed, C2-C0 at t=8.50 only):** one ideal local extraction.
  `W_coh` before `3.0626774684866464e-3` (C0) vs `3.2509460312412416e-3` (C2);
  `Delta W_coh = 1.8826856275459521e-4`; `W_load_gross` `3.0626774684866503e-3` vs
  `3.2509460312412546e-3`, `Delta = 1.8826856275460432e-4`;
  `Delta W_operation_out_total = 1.8826856272917414e-4`;
  `Delta` interaction jump `2.5412084893412970e-14`; `Delta` load energy before
  `-6.8592992121339208e-5`; gross-work-vs-`W_coh` residual `9.1072982488782372e-18`; post-extraction
  load ergotropy 0 for both. Checks 40/40 PASS. Verdict
  `state_improvement_preserved_after_extraction`. The ledgers are not summed into a net benefit.
- **C2 extension to t=10 (supplementary, no new Milestone number, no extraction):** propagating the
  saved t=7.70 state once with C2 to t=10.00 and comparing against the fixed M12A C0 long-window
  baseline, the C2-C0 bare-load ergotropy difference is `+1.882685627546e-4` (+6.147189%) at t=8.50,
  `+5.819146953273e-4` (+24.071514%) at t=9.68, and `+6.460245043013e-4` (+27.768468%) at t=10.00,
  increasing in all 150 saved intervals after t=8.50 with no near-zero or negative crossing. C2 formed
  **no new W peak** after t=8.50 — the difference grew because C0 declined faster (t=10 retention
  87.807172% vs 68.723664%). Ergotropy time-area over 7.70..10.00: `A_Delta = 7.2610140552155800e-4`
  (+11.059330%), a time-area of instantaneous ergotropy, not cumulative extracted work. Classification
  `C2_advantage_positive_at_t10`; 28/28 and 15/15 checks PASS.

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
- For the matched N=7 condition, that the load's diagonal (population) ergotropy is exactly 0 at every
  saved time, so `W = W_coh`, with an exact decomposition residual of 0 and the signed peak->t10 changes
  above (M12A).
- The saved-grid event ordering `current decline (7.50) -> W_coh decline (7.70) -> Cl1 decline (7.73) ->
  current zero crossing (8.15) -> persistent negative current (8.16) -> persistent positive return
  (9.68)`, stable across 3/5/7-interval rules (M12B).
- The load-chain mutual information peak/persistent-decline time t=7.51 and peak value
  `4.5985207125e-3` nats, 0.19 earlier than the `W_coh` peak, stable across 3/5/7-interval rules (M12C).
- That, from an identical stored state, halving the chain-load coupling after the branch weakens the
  `W_coh` decrease within the finite window at **both** branch starts, with the endpoint and time-area
  relative differences above; plus the B0 reproduction (max abs difference 0) and the B4 zero-coupling
  boundary control (max invariant drift `1.998e-15`) (M13A).
- The numerical/SHA-256 audits of the re-analysis stage: 35/35 checks PASS, 11/11 formal M13A inputs
  unchanged, M13A B3 values reproduced within 1e-12 (M13B).
- That the implemented interaction term is exactly linear in `g` and that the drive is exactly zero at
  t=7.50, 7.70, and 8.50, so post-switch protocol differences are not drive-input differences; unit
  checks 11/11 PASS with zero residuals (M14A, no time evolution).
- The ideal-quench Hamiltonian-jump work for both coupling-down protocols,
  `1.6980822985791084e-15` (C1) and `8.8624029002963674e-16` (C2), each below the 1e-12 audit
  tolerance, with the state exactly unchanged across the evaluated switch and M13A B0 quantities
  reproduced to `3.839e-13`; checks 22/22 PASS (M14 Plan B).
- That at t=8.50 a single ideal local extraction gives `W_load_gross` larger for C2 than for C0 by
  `1.8826856275460432e-4`, agreeing with the pre-extraction `W_coh` difference to a residual of
  `9.1072982488782372e-18`, with the total-Hamiltonian accounting giving the same sign and
  post-extraction load ergotropy exactly 0 in both protocols; checks 40/40 PASS (M15A).
- That, extending C2 from the formally saved t=7.70 state to t=10 and comparing it with the fixed C0
  baseline on the common saved grid, the C2-C0 bare-load ergotropy difference is still positive at t=10
  and increased across every saved interval after t=8.50, under these fixed finite conditions. This
  concerns the difference between the two protocols, not growth of C2's own ergotropy, which formed no
  new peak after t=8.50; checks 28/28 and 15/15 PASS (C2 extension, no extraction performed).

## Descriptive model support

- Single amplitude/shift/scale transform, normalized RMSE ~ `0.05485` (M11h).
- Asymmetric rise/fall transform, normalized RMSE ~ `0.03899`,
  `asymmetric_time_scaling_partially_supported` (M11j).
- Post-peak two-stage exponential, normalized RMSE ~ `0.00810`, `two_stage_tail_decay_supported`;
  late-tail concentration ~ `59.6%` in t=9->10, `late_tail_structure_remains` (M11k).

These are descriptive shape diagnostics, not evidence of any physical mechanism.

**Counterfactual branch response (M13A).** Not a curve fit but the same category of evidence: a
descriptive, finite-window comparison between branches that share an identical state and differ only in
the generator applied after the branch time. The B3 (coupling halved) branch shows a weaker `W_coh`
decrease than the B0 baseline at both starts (endpoint `W_coh` relative difference `+6.142%` / `+6.147%`;
time-area relative difference `+1.467%` / `+2.103%`), with the B4 zero-coupling boundary control passing.
This establishes a **sensitivity of the observed `W_coh` decline to the chain-load coupling within this
window** — it is not a mechanism, not a protection scheme, not an optimum, and the B1 (dephasing) and B3
(coupling) percentages are not ranked against each other. M13B, which re-analyzed only these saved
values, could not uniquely order the internal response
(`temporal_response_inconclusive`); that is a limit of the saved temporal resolution and does not
withdraw the B3 result above.

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
- The internal temporal ordering of the weaker-coupling branch response: under the 3/5/7-interval x
  epsilon sensitivity grid the B3-B0 persistent-difference order was **not stable**, so the order in
  which current, load energy, passive energy, `Cl1`, MI, and `W_coh` respond could not be uniquely
  determined at the current saved interval and persistence rules (M13B,
  `temporal_response_inconclusive`). The immediate current difference at the branch time is structural
  (the current contains the changed coupling coefficient) and is not evidence of a response order.
- The physical mechanism behind any of the M12/M13 observations. Coupling sensitivity is observed
  (M13A), but no mechanism is identified: it is not shown that coupling is the sole cause of the `W_coh`
  decrease, that weaker coupling is generally better, that an optimal coupling exists or was found, that
  passive energy converted into coherent ergotropy, that correlation suppression protected work, or that
  the mutual-information decrease caused (or was converted into) the `W_coh` decrease.
- Generalization of any M12/M13 result beyond this single fixed matched N=7 condition, these two branch
  start times, the finite window (branch t<=8.50, trajectory t<=10), and this saved resolution; within
  M12/M13 themselves no extraction, recovery, control, or repeated-cycle operation was performed (the
  later C2 extension to t=10 is a separate supplementary check and likewise performed no extraction).
- Anything beyond the ideal, instantaneous idealizations used in M14/M15A. Specifically not confirmed:
  a **finite-time switch or ramp** (the modeled quench has mathematically zero duration, no controller
  Hamiltonian and no bandwidth limit); the **implementation cost of any actuator or controller** (a
  numerically zero Hamiltonian-jump work is not a zero device cost); the **implementation cost of the
  extraction unitary**; any **net benefit combining switching work and extraction work** (the ledgers
  are deliberately kept separate and are never summed); **repeated extraction or cycle operation**;
  **general long-time stability, steady state, or behavior beyond t=10** — note that the formal windows
  and the single extraction time of M13A and M15A themselves remain t=8.50, while a separate
  supplementary check has propagated C2 to t=10, so saved-grid comparisons at t=9.68 and t=10 have been
  carried out (without extraction); **parameter generalization** to other coupling values, sweeps, N,
  gamma, Omega, branch times, or conditions beyond those fixed here; **optimal control** of any kind (no
  optimum was searched for or found, and weaker coupling is not shown to be generally better); any
  **experimental advantage** on real hardware; and any **quantum advantage**.
- Fair baselines against classical wave/stochastic models; novelty/priority vs the literature (no
  literature review done).

---

## Reproduction and verification

Complete reproduction packages are stored as Milestone-specific ZIP archives in the repository's
`zip/` directory. Source code, reports, CSV outputs, audit files, and any saved states required by a
given stage are contained in the corresponding archive; they are not all exposed as individual files
at the GitHub repository root.

Download and extract the relevant archive, then run the commands below from the extracted directory
that contains its `Cargo.toml`. A Rust toolchain is required (`Cargo.toml`: edition 2021; deps
`nalgebra`, `num-complex`, `thiserror`; dev-dep `approx`). Archive contents vary by Milestone, so a
binary or data file listed below may be present only in the package for that stage.

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

Later stages report their own runtimes and audits: the M14A design checks needed no time evolution
(3.455 s after build, 11/11 PASS); M14 Plan B propagated a single 3080-step baseline prefix to t=7.70 in
2177.4 s (36.29 min; construction 5.2 s, prefix 2169.8 s, two-state diagnostics/write 2.4 s) with
22/22 checks PASS and 15/15 hashes unchanged; M15A reports 40/40 numerical and provenance checks PASS.
The two saved 384x384 states from Plan B are stored as row-major little-endian complex-f64 binaries
with `QWNRHO1` metadata (2,359,320 bytes each), individually hashed in `m14_plan_b_state_files.csv`.
The supplementary C2 extension to t=10 reports a formal runtime of 682.479 s (11.37 min) with 28/28
mandatory numerical checks PASS, 15/15 existing-artifact reproduction checks PASS, SHA audit PASS, and
eigensolver fallback count 0; it saves two unextracted C2 states (`c2_pre_extraction_state_t850.bin`,
`c2_pre_extraction_state_t1000.bin`) in the same `QWNRHO1` format, each with a recorded SHA-256.

---

## Repository and archive layout

The public GitHub repository is organized as a small entry point plus Milestone-specific reproduction
archives:

```
README.md                 # integrated project overview
LICENSE                   # MIT license
reports/                  # optional browser-readable report copies, where provided
zip/                      # Milestone-specific complete reproduction packages
  *.zip
```

The `reports/` directory is a convenience for reading selected reports in the browser. The
reproduction packages in `zip/` are the self-contained records used for downloading and reproducing a
stage. A report may therefore appear both as a browser-readable copy and inside its corresponding ZIP.

After extracting the relevant ZIP, a typical Rust reproduction package has the following internal
layout. Not every archive contains every later-stage file or binary:

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

This tree describes archive contents, not the GitHub repository root. The `src/lib.rs` module
declarations are based on the actual files; the `bin/` mapping and the Milestone 10-11 binaries are
organized from the "generated files" sections of the respective reports. For a complete list of
individual CSVs, see the end of each MILESTONE report in the relevant ZIP.

---

## Detailed reports and data

The filenames below refer to files inside the relevant Milestone reproduction ZIP. Some reports may
also be mirrored under `reports/` for browser reading. Each report lists its own full CSV set.

Minimal pointers for the load-bearing results:

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
- **Post-peak ergotropy decomposition (12A):**
  `MILESTONE_12A_LOCAL_ERGOTROPY_DECOMPOSITION_REPORT.md` (fixed-condition trajectory and
  decomposition/quality CSVs listed therein).
- **Temporal relation (12B):** `MILESTONE_12B_TEMPORAL_RELATION_REPORT.md` and `M12B_PRECOMMIT.md`
  (event, interval-change, and normalized-comparison CSVs listed therein).
- **Load-chain correlation (12C):** `MILESTONE_12C_LOAD_CHAIN_CORRELATION_REPORT.md`; key CSV
  `m12c_interval_changes.csv`.
- **Counterfactual branch experiment (13A):**
  `MILESTONE_13A_FORMAL_SHORT_HORIZON_COUNTERFACTUAL_BRANCH_REPORT.md`; key CSV
  `m13a_formal_branch_timeseries.csv` (728 formal rows).
- **Branch temporal response (13B):** `MILESTONE_13B_WEAKER_COUPLING_TEMPORAL_RESPONSE_REPORT.md` and
  `M13B_PRECOMMIT.md`; key CSVs `m13b_input_audit.csv`, `m13b_maximum_differences.csv`,
  `m13b_interval_summary.csv`, `m13b_current_event_comparison.csv`, `m13b_energy_decomposition.csv`,
  `m13b_normalized_delta_timeseries.csv`, `m13b_start_time_comparison.csv`.
- **Coupling-control design (14A, design only):**
  `MILESTONE_14A_COUPLING_CONTROL_AND_SWITCHING_WORK_DESIGN.md` (no formal control trajectory; the
  11-row no-propagation unit-check CSV is described therein).
- **Switching-work accounting (14 Plan B):** `MILESTONE_14_PLAN_B_SWITCHING_ACCOUNTING_REPORT.md`; key
  CSVs `m14_plan_b_switching_work.csv`, `m14_plan_b_classification.csv`,
  `m14_plan_b_scale_diagnostics.csv`, `m14_plan_b_state_files.csv`.
- **End-point extraction accounting (15A):** `MILESTONE_15A_ENDPOINT_EXTRACTION_ACCOUNTING.md`; key
  CSVs `endpoint_extraction_comparison.csv`, `gross_work.csv`, `energy_ledgers.csv`,
  `numerical_checks.csv`.
- **C2 extension to t=10 (supplementary, no new Milestone number):**
  `MILESTONE_C2_EXTENSION_TO_T10_COMPARISON.md` and `C2_EXTENSION_INPUT_AUDIT.md`; key CSVs
  `c0_c2_extended_summary.csv`, `c0_c2_extended_comparison_timeseries.csv`,
  `c0_c2_extended_time_area_summary.csv`, `c0_c2_extended_difference_events.csv`,
  `c2_extension_numerical_checks.csv`, `c2_extension_input_manifest.csv`.


Reading conventions for the CSVs: `W_time_area`/`E_time_area` are time-areas of state quantities (not
cumulative work or input); `W/Ein` is not an overall efficiency; `usable_fraction` is ergotropy over
load energy; XGamma is a diagnostic, not a loss/efficiency. `W_coh`/`W_pop` are a descriptive
energy-basis split of the load ergotropy, not independent conserved quantities; mutual information is a
total-correlation observable in nats, not work, ergotropy, or entanglement; `W_coh` time-area is a
time-integral of a state quantity, not cumulative extracted work. Switching work
(`W_switch_on_system`/`W_switch_out`) is an ideal instantaneous Hamiltonian-jump ledger, not an actuator
or controller cost; `W_load_gross` and `W_operation_out_total` are single-shot ideal-extraction ledgers.
The state, switching, and extraction ledgers are reported separately and are never added or subtracted
into a net benefit.

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
- **M14/M15A binary names and repository paths.** The 14A, 14 Plan B, and 15A reports identify their
  source and output paths on the author's working machine but do not state repository-relative binary
  names, so no entries were added to the `bin/` map in *Repository layout* rather than guessing.
- **The M14A unit-check CSV and the M15A `numerical_checks.csv` are distinct artifacts.** 14A records
  11/11 no-propagation unit-check rows; M15A's `numerical_checks.csv` holds the 40/40 extraction-stage
  checks. An earlier-milestone file of the same generic name exists for the N=3..7 event comparison; the
  three should not be conflated when collecting CSVs.
- **A numerically zero switching work is not a cost statement.** `W_switch_on_system` values of order
  `1e-15`/`1e-16` are below the 1e-12 audit tolerance for the *modeled Hamiltonian jump only*. No
  actuator, ramp, bandwidth, or device-cost model exists in the project, so no controller-cost
  conclusion can be drawn from them.

---

## Citation and license

License: MIT (see `Cargo.toml` / `LICENSE`). When citing, please reference this repository and the
specific `MILESTONE_*.md` report(s) underlying the result used, since each report defines the exact
conditions and reservations for its numbers.
