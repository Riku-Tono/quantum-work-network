# From fixed-total comparison to equal-input N=3 vs N=7 analysis

This README summarizes two connected lines of work in the Quantum Work Network project:
Milestone 10 (a fixed-total-dephasing comparison of N=3, 5, 7) and Milestones 11c–11k (an
equal-input comparison of N=3 vs N=7). Milestone 10 compared chains at equal *total* dephasing but
**not** equal drive input energy; that gap is exactly what motivated the equal-input matching in
Milestone 11.

Everything below is drawn from the attached official reports and CSVs. No new time
evolution, parameter scan, or numerical estimate was produced for this README. Where a
statement is a directly measured value it is written as such; where it is a candidate
interpretation it is flagged as descriptive only. Curve-transform fits and decay fits are
**descriptive shape diagnostics**, not proofs of any physical mechanism.

---

## 1. Purpose and scope

The goal is to let a reader follow, in order:

1. What Milestone 10 fixed and compared, and why its result alone could not isolate the effect of
   chain length.
2. What was then held fixed for the equal-input comparison, and how that matched point was reached.
3. What was directly confirmed about the N=3 vs N=7 shape difference.
4. How much of that difference can be *described* by transforming the N=3 curve.
5. What remains unconfirmed, and what is worth examining next.

The intended reading order is Milestone 10 (fixed-total context) → 11c–11f (fix the equal-input
point) → 11g (shape) → 11h (single-scale transform) → 11j (asymmetric transform) →
11k (post-peak decay) → 11i (a diagnostic *design* that has not been executed).

---

## 2. Milestone 10: fixed-total comparison and the motivation for equal-input matching

Milestone 10 organized the fixed-total-dephasing comparison that the equal-input work builds on. It
distinguished two noise conventions that are easy to conflate: **fixed-per-site** (each site keeps
the same rate, so the summed total grows with N) and **fixed-total** (the summed rate is held at a
constant `TOTAL_GAMMA`, so the per-site rate is `TOTAL_GAMMA/N`). For N=5 and N=7 these are not the
same total-noise condition.

- **Milestone 10a** compared existing results across N=3, 5, 7 without running any new time
  evolution. It separated the fixed-per-site and fixed-total conditions and, where the official
  sources lacked fixed-total N=3/N=5 values beyond `W_max`, left them as `not_available` rather than
  filling them with guesses or interpolation. Judgment: `completed_with_explicit_missing_values`.
- **Milestone 10b** newly computed `TOTAL_GAMMA = 3.0` for N=3, 5, 7 (one trajectory each) and
  introduced **XGamma**, a dephasing-kernel-weighted coherence-exposure diagnostic
  `x_gamma(t) = sum_ab Gamma[a,b] |rho[a,b](t)|^2` (time-integrated). XGamma is a diagnostic
  quantity only — **not** lost work, dissipated energy, dephasing power, heat, entropy production,
  damage, or efficiency. Judgment: `completed_fixed_total_gamma_3_comparison`.
- **Milestone 10c** re-computed `TOTAL_GAMMA = 1.5` for N=3, 5, 7 under the same XGamma diagnostics,
  filling the fixed-total gaps from 10a with official values. Its N=7 trajectory matched the 9c
  primary record to a maximum difference of 0 across 7 quantities × 1001 times. **This 10c N=3
  series is the reference trajectory later reused as the N=3 side of the Milestone 11 comparison.**
  Judgment: `completed_with_fallback_diagnostic`.
- **Milestone 10 Final** read only the 10a/10b/10c official artifacts (no new time evolution). Under
  both `TOTAL_GAMMA = 1.5` and `3.0`, the same finite-condition rankings held:

| Metric | Ranking (both TOTAL_GAMMA = 1.5 and 3.0) |
|---|---|
| `W_max` | N=7 > N=5 > N=3 |
| `W(t=10)` | N=7 > N=5 > N=3 |
| `usable_fraction` | N=7 > N=5 > N=3 |
| `W_time_area` | N=3 > N=5 > N=7 |
| ergotropy arrival | N=3 fastest, N=5 middle, N=7 slowest |
| XGamma | N=7 > N=5 > N=3 |

The reading is deliberately not "longer chains are better." Under fixed total dephasing, the
**metric rankings differ**: longer chains have the higher peak/final ergotropy and usable fraction,
while shorter chains have the larger `W` time-area and the earlier arrival. This is a
finite-condition descriptive result for this model, not a general law. `W_time_area` and
`E_time_area` are time-integrals of state quantities, not cumulative extracted work or cumulative
input energy; the higher N=7 usable fraction is a descriptive ratio, not an independent performance
gain; and the two `TOTAL_GAMMA` points are not used to infer any functional form or universal
multiplier. Judgment: `completed_with_explicit_missing_values`.

**Why this led to Milestone 11.** In Milestone 10 the total dephasing was matched, but the
**integrated drive input energy** of N=3 and N=7 was not. So the ranking differences alone could not
isolate the effect of chain length from the effect of unequal input. Milestones 11c–11f therefore
adjusted the N=7 `Omega` until its drive input energy matched N=3's, and 11g onward compared the
time-series shape at that equal-input condition. Note that the Milestone 11 matched N=7 trajectory
(`Omega ≈ 0.18748`) is a different trajectory from the Milestone 10c N=7 run (`Omega = 0.2`); the
two must not be conflated.

---

## 3. Physical conditions and equal-input matching

The physical model is unchanged from the earlier milestones. The same model construction, drive
protocol, RK4 scheme, dt, load, and vacuum preparation were used. Chain length, Hilbert dimension,
bond count, and per-site dephasing rate differ between N=3 and N=7. What this line of work adds is
a **matched drive input**: the N=7 drive strength `Omega` is chosen so that the integrated
drive input energy `E_drive_in` equals the N=3 reference value.

| Fixed / matched quantity | Value |
|---|---|
| Chains compared | N=3 vs N=7 |
| Total dephasing `TOTAL_GAMMA` | 1.5 (distributed across chain sites) |
| Time step `dt` | 0.0025 |
| Horizon `T` | 10 (4000 RK4 steps, 1001 saved points) |
| N=3 reference `E_drive_in` (target) | 5.9618618770136536e-2 |
| N=7 matched `Omega` | **0.18748395731510084** |
| N=7 measured `E_drive_in` | 5.9618450901925489e-2 |
| Absolute input mismatch | 1.6786821104702865e-7 |
| **Relative input mismatch** | **2.8157011099880636e-6** |

The match is a match on **drive input energy only**. Chain length, Hilbert dimension, bond
count, and geometry all still differ between N=3 and N=7; this is not an equal-cost,
equal-efficiency, or single-variable causal comparison. `W/E_drive_in` and similar ratios are
input-normalized descriptors, **not** device or overall efficiency.

---

## 4. Milestones 11c–11f: fixing the equal-input point

The equal-input point was reached by an auditable, one-trajectory-at-a-time procedure rather
than an opaque solve. Milestones 11c–11e narrowed and bracketed the N=7 drive strength `Omega`,
and Milestone 11f evaluated the final interpolated point. The full per-step search numbers live
in the individual reports; the essential chain is:

- **11c (precheck).** Using only the two stored official trajectories (no new time evolution),
  the drive input energy was re-integrated with the existing signed-power integrator. The N=3
  reference was `target_E_drive_in = 5.9618618770136536e-2`; the N=7 trajectory at the earlier
  `Omega = 0.2` integrated *above* target, so a lower `Omega` was indicated. A weak-drive
  (`E_drive_in ∝ Omega^2`) relation was used **only** to place an initial guess
  `Omega ≈ 0.18770762` (a `quadratic_response_initial_guess`, explicitly not a root). Judgment:
  `completed_input_matching_precheck`.
- **11d.** One trajectory was run at that guess `Omega = 0.18770762191709489`. Its measured input
  was slightly above target (matching residual `+1.397e-4`, tolerance FAIL), so a lower probe was
  required to form a bracket. Judgment: `lower_probe_required`.
- **11e.** One trajectory was run at `Omega = 0.1870`. Its measured input was below target
  (residual `-3.023e-4`, tolerance FAIL). The sign change between 0.1870 and 0.18770762 confirmed
  a local bracket `[0.18700000, 0.18770762191709489]`, and the linear interpolation inside it gave
  the candidate `Omega_interpolated = 0.18748395731510084`. Judgment:
  `local_input_bracket_confirmed`.
- **11f.** One trajectory was run at that interpolated point,
  `Omega = 0.18748395731510084`. The measured input matched the target to relative mismatch
  `2.8157011099880636e-6`, within tolerance, with passing numerical-quality checks (worst
  eigenvalue `-5.061e-18`, primary solver success/failure `1000/1`, one fallback attempt
  succeeding, solver failure `0`, ledger max `5.197e-7`). Final classification:
  `matched_input_found_with_fallback_diagnostic`.

Each of 11d, 11e, and 11f ran exactly one new N=7 trajectory; 11c ran none. The reports are
consistent that this procedure establishes a **matched input at one bracketed point**, and does
**not** establish a unique root, global uniqueness, monotonicity of the matching function in
`Omega`, or `dt` convergence of the matching condition. All physical quantities saved by 11d/11e
carry the reservation that they are not to be used for a fair performance conclusion until the
match is complete — which is why the shape analysis in Section 5 onward uses the 11f matched
trajectory.

---

## 5. Milestone 11g: time-series shape and peak width

Milestone 11g analyzed only the two stored, input-matched trajectories (N=3 from the Milestone 10c
official series, N=7 from the 11f matched series). No new time evolution was run.

| Quantity | N=3 | N=7 | N7 / N3 |
|---|---|---|---|
| `W_max` | 3.0302005931286463e-3 | 3.3852501212993710e-3 | **1.1171703** |
| `t` at `W_max` | 5.63 | 7.70 | 1.3676732 |
| Time above half-max (longest) | 6.1210465 | 3.8540023 | 0.6296313 |
| Time above 90% of max | 3.4755391 | 1.4792893 | 0.4256287 |
| `W` time-area (0→10) | 1.7361323e-2 | 1.2014985e-2 | 0.6920547 |
| `usable_fraction` at t=10 | 1.8776464e-1 | 3.6900619e-1 | 1.9652592 |

The `W` curves cross twice; N=7 is above N=3 over one interval (largest positive-difference
interval ≈ 6.8717–9.3786, area ≈ 8.937e-4). The negative-difference area (6.240e-3) is larger,
so the net `W` time-area difference is negative: **N=7 has the higher peak, but a narrower
peak and a smaller 0→10 `W` time-area than N=3.**

Two cautions carried by the report itself: `W_time_area` is the time-area of a state quantity,
**not** cumulative extracted work; and the higher N=7 `usable_fraction` is a descriptive ratio,
**not** an independent performance gain. Reported `W`–coherence Pearson correlations are
descriptive, not causal.

---

## 6. Milestone 11h: single time-scale transform

Milestone 11h fit the N=7 curve as a simple transform of the N=3 curve,
`W_model(t) = A · W3((t − delta)/s)`, using only stored trajectories. Models 0–3 were searched
on a deterministic coarse→fine grid; extrapolation was excluded.

The best post-arrival single-scale transform (model 3):

| Parameter | Value |
|---|---|
| `A` (amplitude) | 1.02343645 |
| `delta` (shift) | 2.49 |
| `s` (time scale) | 0.91 |
| Normalized RMSE | **0.05485** |
| Max absolute residual | 4.673381e-4 (at t=10.00) |
| Residual sign changes | 3 |

Classification: `partial_transform_with_structured_residual`. A single amplitude/shift/scale
transform captures most of the curve but leaves structured residuals. The report records two
distinct decompositions of the usable-fraction difference: (1) an exact algebraic split into a
W-difference term and an energy-denominator term for a chosen ordering, and (2) a separate
first-order finite-difference approximation with a nonlinear residual. Neither decomposition is
causal. `AIC-like`/`BIC-like` values used here are descriptive complexity comparisons only.

---

## 7. Milestone 11j: asymmetric rise/fall time transform

Milestone 11j asked whether letting the rising side and falling side of the curve have
**different** time scales improves the description. The boundary between rise and fall is tied
to the N=3 peak (`t3 = 5.63`) and is continuous at that boundary, so it adds no free parameter;
the free parameters are `A`, `delta`, `s_rise`, `s_fall`. The comparison interval is the same
post-arrival window (t=3.83→10, 618 points) as 11h.

| Model | `A` | `delta` | `s_rise` | `s_fall` | Normalized RMSE | BIC-like |
|---|---|---|---|---|---|---|
| Baseline (single scale) | 1.02343645 | 2.49 | 0.91 | 0.91 | 0.05485 | −10599.80 |
| **Asymmetric** | 1.05674517 | 2.32 | 0.965 | 0.515 | **0.03899** | −11015.26 |

The scale difference `|s_rise − s_fall| = 0.45` is larger than the fine search step (0.005).
Residuals: 3 sign changes, max positive `1.847604e-4` (t=7.65), max negative `-2.826439e-4`
(t=9.50), with residual MAE growing from early to late window.

Rule-based classification: **`asymmetric_time_scaling_partially_supported`**.
Complexity-penalized (descriptive) classification: `improvement_supported_after_complexity_penalty`.
The report is explicit that `AIC-like`/`BIC-like` are descriptive comparisons, not proof of a
generative model or a physical mechanism, and that reflection, group velocity, mode beating,
and entanglement were **not** assessed. All 7 input files matched the pre-analysis 11h artifacts
by SHA-256.

---

## 8. Milestone 11k: post-peak decay and late-tail residual

Milestone 11k characterized the N=7 curve **after** its peak (main interval t=7.70→10.00) and
examined where the 11j residual concentrates. Low-parameter decay models were compared.

| Decay model | Free params | Normalized RMSE | BIC-like |
|---|---|---|---|
| A — single exponential | 2 | 0.01014 | −4738.20 |
| **B — two-stage exponential** | 3 | **0.00810** | **−4836.73** |
| C — exponential × linear | 3 | 0.00913 | −4791.18 |

Best post-peak model (BIC-like): **two-stage exponential**, with
`switch = 8.50`, `lambda_1 = 0.15`, `lambda_2 = 0.20`, `|lambda_2 − lambda_1| = 0.05`.
Two-stage decay is **descriptively** supported (BIC-like improves by ≈98.5 over single
exponential) but this does **not** imply two physical processes.

Tail of the 11j residual: absolute area ≈ 4.075e-4, with **≈ 59.6%** of that absolute area
concentrated in **t = 9→10**; 1 sign change; max positive `1.839769e-4` (t=7.70), max negative
`-2.826439e-4` (t=9.50). Classification of the tail: `tail_residual_concentrated_late`.

Component normalization at t=10 (relative to t=7.70): `E = 0.9713`, `W = 0.6872`,
`passive = 1.2809`, `coherence = 0.8289`. Over this window `W` falls faster than total `E`,
and the passive part does not fall but rises. The N7−N3 `usable_fraction` difference stays
positive but shrinks (0.2444 at 7.70 → 0.1812 at t=10); the report does **not** read this as an
independent usable-fraction gain, and does **not** claim coherence drives the `W` decay.

Final judgment: **`late_tail_structure_remains`**.
Auxiliary judgment (internal to the decay fit): **`two_stage_tail_decay_supported`**.

---

## 9. Milestone 11i: unexecuted physical-mechanism diagnostic design

Milestone 11i is a **design and audit**, not a measurement. It audited the stored 11h artifacts
and confirmed that only aggregate load quantities (energy, ergotropy, coherence), power/quality
diagnostics, and comparison summaries are saved — **site-resolved populations, bond and
site–load currents, site–site coherences, full density matrices, mutual information, negativity,
and mode occupation are not stored.** It therefore concluded that the ≈5.5% residual cannot be
attributed to any mechanism from the existing load time series alone.

The design lays out candidate diagnostics (group velocity / finite-chain modes, boundary
reflection, correlation/entanglement fronts via mutual information and negativity, mode beating),
fixes the pre-selected inspection times, and specifies the minimal follow-up: **re-run the single
N=7 matched trajectory once** under identical physical conditions and additionally save
site-resolved currents/populations plus a few reduced states at selected times.

Status: `completed_design_with_targeted_recomputation_required`. **This site-resolved
re-computation has not been executed.** Every candidate mechanism (group velocity, dispersion,
boundary reflection, correlation/entanglement front, mode beating) is explicitly listed as
undetermined.

---

## 10. What has been confirmed so far

Directly measured, within this fixed model, this finite horizon, and this implementation:

- The N=7 trajectory at `Omega = 0.18748395731510084` matches the N=3 drive input to relative
  mismatch `2.8157e-6`, with passing numerical-quality checks (11f).
- Peak values and timing: N=3 `W_max = 3.0302e-3` at t=5.63; N=7 `W_max = 3.3853e-3` at t=7.70;
  ratio ≈ 1.11717 (11g).
- N=7 has the higher peak but a narrower peak and a **smaller** 0→10 `W` time-area than N=3;
  the `W` curves cross twice (11g).
- A single amplitude/shift/scale transform describes the post-arrival curve with normalized
  RMSE ≈ 0.05485, leaving structured residuals (11h).
- An asymmetric rise/fall transform lowers normalized RMSE to ≈ 0.03899
  (`asymmetric_time_scaling_partially_supported`) (11j).
- After the peak, a two-stage exponential is the best descriptive decay model
  (normalized RMSE ≈ 0.00810; switch 8.50, `lambda_1 = 0.15`, `lambda_2 = 0.20`) (11k).
- The remaining 11j residual is concentrated late: ≈ 59.6% of its absolute area falls in
  t = 9→10; final judgment `late_tail_structure_remains` (11k).

**Overall (at the intended strength), in two stages.**

*Milestone 10.* Under fixed total dephasing, metric rankings differ: longer chains have higher
peak and final ergotropy and higher usable fraction, while shorter chains have larger `W` time-area
and earlier arrival. This is a finite-condition descriptive result, not a general law.

*Milestone 11.* After matching drive input energy between N=3 and N=7, the N=7 curve still has a
higher, later, narrower peak, but a smaller `W` time-area. Transform and decay fits describe much of
the shape difference, while a late-tail residual remains and its physical mechanism is unresolved
from the stored load time series alone.

---

## 11. What has not yet been confirmed

- Any physical mechanism behind the observed shape, including group velocity, boundary
  reflection, mode beating, entanglement/correlation fronts, and any causal relationship.
- That the transform fits or the decay fits correspond to specific physical processes (they are
  descriptive shape fits).
- `dt`-halving convergence of the matched N=7 condition.
- Uniqueness of the matching root, global uniqueness, or monotonicity of the matching function.
- Equal-input results for N=5, for `TOTAL_GAMMA = 3.0`, for other `Omega` roots, or for N > 7.
- Behavior for t > 10, extraction cycles, continuous operation, or real-device performance.
- Any quantum advantage or general advantage of longer chains.
- That the higher N=7 `usable_fraction` is an independent performance gain (it is a descriptive
  ratio).

---

## 12. Candidate next steps

These are candidates recorded by the reports; none is executed automatically.

- Decide whether a `dt`-halving verification of the matched N=7 condition is warranted (raised as
  a question by 11f, 11g, and 11h).
- If the late-tail structure is to be pursued, run the **single** N=7 matched trajectory once more
  under identical physical conditions and additionally save the site-resolved currents/populations
  and selected reduced states specified in the 11i design — then compare those against the late
  residual. This is the 11i plan, still unexecuted.

---

## 13. Reproduction and verification information

All milestones in these lines report `cargo fmt --all -- --check`: PASS and the per-milestone
analysis binary: PASS. For the Milestone 10 line, the reports list 101 passed for 10b, 104 for 10c,
and 107 for 10 Final, each with 0 failed / 1 ignored; 10a records 16/16 required checks PASS but
does not state a numerical release-test count. For the Milestone 11 line, the reports explicitly
list 110 passed for 11c, 113 for 11d, 116 for 11e, and 119 for 11f and 11h–11k, each with
0 failed / 1 ignored. The 11g report records formatting and analysis execution as PASS but does not
state the numerical release-test count.

Numerical-quality and integrity notes recorded by the reports:

- Single-trajectory runs 11d/11e/11f each passed their checks (worst selected minimum eigenvalue
  on the order of `-5e-18`, solver failure `0`, ledger max on the order of `5.2e-7`). 11f: primary
  solver success/failure `1000/1`, fallback success/attempt `1/1`.
- 11i audit: 15/15 checks PASS. 11j: 23/23 checks PASS. 11k: 13/13 checks PASS.
- 11j and 11k verified their official inputs by SHA-256 both before and after analysis; the
  hashes are listed in the respective reports.

Only Milestones 11d, 11e, and 11f performed new N=7 time evolution (one trajectory each);
11c, 11g, 11h, 11i, 11j, and 11k operate on stored trajectories only.

---

## 14. Referenced artifacts

Reports (primary sources used here):

- `MILESTONE_10A_REPORT.md`, `MILESTONE_10B_REPORT.md`, `MILESTONE_10C_REPORT.md`,
  `MILESTONE_10_FINAL_REPORT.md` (fixed-total context and the N=3 reference trajectory).
- `MILESTONE_11C_PRECHECK_REPORT.md`, `MILESTONE_11D_REPORT.md`, `MILESTONE_11E_REPORT.md`
  (fixing the equal-input point).
- `MILESTONE_11F_REPORT.md`, `MILESTONE_11G_REPORT.md`, `MILESTONE_11H_REPORT.md`,
  `MILESTONE_11I_MOUNTAIN_MECHANISM_DESIGN.md`, `MILESTONE_11J_ASYMMETRIC_TRANSFORM_REPORT.md`,
  `MILESTONE_11K_POST_PEAK_DECAY_REPORT.md`.

CSVs:

- `fixed_total_gamma_1_5_trajectory_comparison.csv` (10c N=7 vs 9c primary record),
  `fixed_total_gamma_three_point_comparison.csv` (10b, TOTAL_GAMMA = 0/1.5/3.0, available values only)
- `input_matching_precheck_integrals.csv`, `input_matching_precheck_guess.csv` (11c re-integration
  and initial guess)
- `input_matching_interpolated_trial_summary.csv` (11f matched point and quality)
- `equal_input_timeseries_shape_summary.csv`, `equal_input_peak_widths.csv` (11g)
- `equal_input_curve_transform_models.csv` (11h single-scale transforms)
- `equal_input_asymmetric_transform_models.csv`,
  `equal_input_asymmetric_transform_residual_summary.csv` (11j)
- `equal_input_post_peak_decay_models.csv`, `equal_input_post_peak_decay_summary.csv`,
  `equal_input_post_peak_component_comparison.csv` (11k)

Detailed tables live in these reports and CSVs; this README is an entry point, not a reprint.

---

## Summary lists

**Directly validated**
- N=7 input match to N=3 at `Omega = 0.18748395731510084`, relative mismatch `2.8157e-6` (11f).
- Peaks: N=3 `W_max = 3.0302e-3` at t=5.63; N=7 `W_max = 3.3853e-3` at t=7.70; ratio ≈ 1.11717 (11g).
- N=7 higher peak, narrower peak, smaller 0→10 `W` time-area; two `W` crossings (11g).
- Numerical-quality checks for 11f, design/audit checks for 11i, and SHA-256 integrity checks for
  11j and 11k.

**Descriptive model support**
- Single amplitude/shift/scale transform: normalized RMSE ≈ 0.05485
  (`A = 1.02343645`, `delta = 2.49`, `s = 0.91`) (11h).
- Asymmetric rise/fall transform: normalized RMSE ≈ 0.03899
  (`A = 1.05674517`, `delta = 2.32`, `s_rise = 0.965`, `s_fall = 0.515`);
  `asymmetric_time_scaling_partially_supported` (11j).
- Post-peak two-stage exponential: normalized RMSE ≈ 0.00810
  (switch 8.50, `lambda_1 = 0.15`, `lambda_2 = 0.20`); `two_stage_tail_decay_supported` (11k).
- Late-tail concentration: ≈ 59.6% of the 11j residual absolute area in t = 9→10;
  `late_tail_structure_remains` (11k).

**Not yet tested**
- `dt`-halving convergence of the matched N=7 condition; root uniqueness/monotonicity.
- Equal-input N=5, `TOTAL_GAMMA = 3.0`, other `Omega` roots, N > 7, t > 10.
- Any physical mechanism: group velocity, boundary reflection, mode beating, entanglement,
  causal relationships.
- Quantum advantage or a general long-chain advantage.

**Proposed but not executed**
- The Milestone 11i site-resolved diagnostic: one additional N=7 matched trajectory under
  identical conditions, saving site-resolved currents/populations and selected reduced states.
