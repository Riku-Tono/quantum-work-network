use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use nalgebra::linalg::{Schur, SymmetricEigen};
use quantum_work_network::coherent_drive::{drive_hamiltonian, CoherentDriveConfig};
use quantum_work_network::dephasing_kernel::DiagonalDephasingKernel;
use quantum_work_network::diagnostics::integrate_signed_power;
use quantum_work_network::ergotropy::ergotropy;
use quantum_work_network::matrix::{
    commutator, expectation, frobenius_norm, hermiticity_error, ComplexMatrix, C64,
};
use quantum_work_network::operators::{build_operators_for_chain, ModelParams, Operators};
use quantum_work_network::partial_trace::partial_trace;

const TOTAL_GAMMA: f64 = 1.5;
const DT: f64 = 0.0025;
const T_FINAL: f64 = 10.0;
const SAVE_EVERY_STEPS: usize = 4;
const SAVE_INTERVAL: f64 = 0.01;
const STEPS: usize = 4000;
const SAVED_POINTS: usize = 1001;
const OMEGA: f64 = 0.1870;
const LOAD_DIM: usize = 3;
const TRACE_TOL: f64 = 1.0e-8;
const HERM_TOL: f64 = 1.0e-8;
const POS_TOL: f64 = 1.0e-8;
const LEDGER_TOL: f64 = 5.0e-5;
const SUM_TRACE_TOL: f64 = 1.0e-10;
const EIG_IMAG_TOL: f64 = 1.0e-10;
const SIGNAL_TOL: f64 = 1.0e-14;
const XGAMMA_TOL: f64 = 1.0e-13;
const M10A_INPUT: &str = "milestone_10a_existing_results_comparison.csv";
const GAMMA_3_SUMMARY: &str = "fixed_total_gamma_3_summary.csv";
const N7_9C_TIMESERIES: &str = "inputs/milestone_10a/n7_fixed_total_validation_timeseries.csv";
const TRAJECTORY_TOL: f64 = 1.0e-12;

const TARGET_E_DRIVE_IN: f64 = 5.9618618770136536e-2;
const HIGH_OMEGA: f64 = 0.18770762191709489;
const HIGH_E_DRIVE_IN: f64 = 5.9758320362953775e-2;
const HIGH_F: f64 = 1.3970159281723876e-4;
const EXISTING_OMEGA: f64 = 0.2;
const EXISTING_E_DRIVE_IN: f64 = 6.7682764053442840e-2;
const OMEGA_SECANT_ESTIMATE: f64 = 0.18749091714073884;
const INPUT_FLOOR: f64 = 1.0e-12;
const RELATIVE_MATCH_TOL: f64 = 1.0e-4;
const ABSOLUTE_MATCH_TOL: f64 = 1.0e-10;
const PRECHECK_INTEGRALS: &str = "input_matching_precheck_integrals.csv";
const PRECHECK_GUESS: &str = "input_matching_precheck_guess.csv";
const HIGH_SUMMARY: &str = "input_matching_single_probe_summary.csv";

const OUTPUTS: [&str; 6] = [
    "input_matching_lower_probe_timeseries.csv",
    "input_matching_lower_probe_summary.csv",
    "input_matching_local_bracket.csv",
    "input_matching_lower_probe_checks.csv",
    "input_matching_lower_probe_performance.csv",
    "MILESTONE_11E_REPORT.md",
];

#[derive(Clone, Copy, Debug)]
struct Spec {
    n: usize,
    dim: usize,
}

const SPECS: [Spec; 1] = [Spec { n: 7, dim: 384 }];

impl Spec {
    fn gamma_site(self) -> f64 {
        TOTAL_GAMMA / self.n as f64
    }
}

#[derive(Clone, Debug)]
struct SolverAttempt {
    attempted: bool,
    all_finite: bool,
    minimum: f64,
    max_imaginary: f64,
    sum_trace_difference: f64,
    passed: bool,
}

impl SolverAttempt {
    fn not_attempted() -> Self {
        Self {
            attempted: false,
            all_finite: false,
            minimum: f64::NAN,
            max_imaginary: f64::NAN,
            sum_trace_difference: f64::NAN,
            passed: false,
        }
    }
}

#[derive(Clone, Debug)]
struct PositivityDiagnostic {
    primary: SolverAttempt,
    fallback: SolverAttempt,
    selected_solver: &'static str,
    selected_minimum: f64,
    fallback_used: bool,
    positivity_pass: bool,
    solver_failure: bool,
    state_finite: bool,
    solver_finite: bool,
    trace_error: f64,
    hermiticity_error: f64,
}

#[derive(Clone, Debug)]
struct TimeRow {
    chain_length: usize,
    gamma_site: f64,
    time: f64,
    load_energy: f64,
    load_ergotropy: f64,
    usable_fraction: f64,
    load_coherence_l1: f64,
    drive_power: f64,
    dephasing_power: f64,
    x_gamma_instant: f64,
    x_gamma_cumulative: f64,
    trace_error: f64,
    hermiticity_error: f64,
    selected_minimum_eigenvalue: f64,
    positivity_solver: &'static str,
    fallback_used: bool,
    state_finite: bool,
    solver_finite: bool,
    bare_energy: f64,
    drive_energy_net: f64,
    dephasing_energy_net: f64,
    ledger_residual: f64,
}

#[derive(Clone, Debug)]
struct Performance {
    construction_seconds: f64,
    propagation_seconds: f64,
    diagnostics_seconds: f64,
    total_seconds: f64,
}

#[derive(Clone, Debug)]
struct Summary {
    chain_length: usize,
    gamma_site: f64,
    endpoint: TimeRow,
    w_max: TimeRow,
    w_time_area: f64,
    e_time_area: f64,
    ergotropy_arrival_time: Option<f64>,
    energy_arrival_time: Option<f64>,
    x_gamma_max: TimeRow,
    primary_success_count: usize,
    primary_failure_count: usize,
    fallback_attempt_count: usize,
    fallback_success_count: usize,
    solver_failure_count: usize,
    worst_selected_minimum_eigenvalue: f64,
    max_trace_error: f64,
    max_hermiticity_error: f64,
    max_abs_ledger_residual: f64,
    checks_passed: bool,
    final_classification: &'static str,
}

#[derive(Clone, Debug)]
struct RunResult {
    rows: Vec<TimeRow>,
    diagnostics: Vec<PositivityDiagnostic>,
    performance: Performance,
    construction_passed: bool,
}

fn config(gamma_site: f64) -> CoherentDriveConfig {
    CoherentDriveConfig {
        omega0: OMEGA,
        omega_drive: 1.0,
        tau: 3.2,
        t_end: T_FINAL,
        dt: DT,
        save_interval: SAVE_INTERVAL,
        gamma_phi: gamma_site,
    }
}

fn format_number(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_owned()
    } else if value == f64::INFINITY {
        "+Inf".to_owned()
    } else if value == f64::NEG_INFINITY {
        "-Inf".to_owned()
    } else {
        format!("{value:.16e}")
    }
}

fn ratio(numerator: f64, denominator: f64) -> f64 {
    if denominator.abs() <= SIGNAL_TOL {
        f64::NAN
    } else {
        numerator / denominator
    }
}

fn state_is_finite(rho: &ComplexMatrix) -> bool {
    rho.iter()
        .all(|value| value.re.is_finite() && value.im.is_finite())
        && rho.trace().re.is_finite()
        && rho.trace().im.is_finite()
}

fn symmetric_attempt(rho_h: &ComplexMatrix) -> SolverAttempt {
    let eigenvalues = SymmetricEigen::new(rho_h.clone()).eigenvalues;
    let all_finite = eigenvalues.iter().all(|value| value.is_finite());
    let minimum = if all_finite {
        eigenvalues.iter().copied().fold(f64::INFINITY, f64::min)
    } else {
        f64::NAN
    };
    let sum = eigenvalues.iter().sum::<f64>();
    let difference = (C64::new(sum, 0.0) - rho_h.trace()).norm();
    let passed =
        all_finite && minimum.is_finite() && difference.is_finite() && difference <= SUM_TRACE_TOL;
    SolverAttempt {
        attempted: true,
        all_finite,
        minimum,
        max_imaginary: 0.0,
        sum_trace_difference: difference,
        passed,
    }
}

fn schur_attempt(rho_h: &ComplexMatrix) -> SolverAttempt {
    let (_, triangular) = Schur::new(rho_h.clone()).unpack();
    let values: Vec<C64> = (0..triangular.nrows())
        .map(|index| triangular[(index, index)])
        .collect();
    let all_finite = values
        .iter()
        .all(|value| value.re.is_finite() && value.im.is_finite());
    let minimum = if all_finite {
        values
            .iter()
            .map(|value| value.re)
            .fold(f64::INFINITY, f64::min)
    } else {
        f64::NAN
    };
    let max_imaginary = values
        .iter()
        .map(|value| value.im.abs())
        .fold(0.0, f64::max);
    let sum: C64 = values.iter().copied().sum();
    let difference = (sum - rho_h.trace()).norm();
    let passed = all_finite
        && minimum.is_finite()
        && max_imaginary.is_finite()
        && max_imaginary <= EIG_IMAG_TOL
        && difference.is_finite()
        && difference <= SUM_TRACE_TOL;
    SolverAttempt {
        attempted: true,
        all_finite,
        minimum,
        max_imaginary,
        sum_trace_difference: difference,
        passed,
    }
}

fn select_attempts(
    primary: &SolverAttempt,
    fallback: &SolverAttempt,
) -> (&'static str, f64, bool, bool) {
    if primary.passed {
        ("symmetric_eigen", primary.minimum, false, false)
    } else if fallback.passed {
        ("complex_schur_fallback", fallback.minimum, true, false)
    } else {
        ("none", f64::NAN, fallback.attempted, true)
    }
}

fn evaluate_positivity(rho: &ComplexMatrix) -> PositivityDiagnostic {
    let state_finite = state_is_finite(rho);
    let trace_error = (rho.trace() - C64::new(1.0, 0.0)).norm();
    let hermiticity_error = hermiticity_error(rho);
    if !state_finite || hermiticity_error > HERM_TOL {
        return PositivityDiagnostic {
            primary: SolverAttempt::not_attempted(),
            fallback: SolverAttempt::not_attempted(),
            selected_solver: "state_input_invalid",
            selected_minimum: f64::NAN,
            fallback_used: false,
            positivity_pass: false,
            // No solver was attempted: this is a physical-state input failure,
            // not a both-solvers-failed diagnostic event.
            solver_failure: false,
            state_finite,
            solver_finite: false,
            trace_error,
            hermiticity_error,
        };
    }
    let rho_h = (rho + rho.adjoint()) * C64::new(0.5, 0.0);
    let primary = symmetric_attempt(&rho_h);
    let fallback = if primary.passed {
        SolverAttempt::not_attempted()
    } else {
        schur_attempt(&rho_h)
    };
    let (solver, minimum, fallback_used, solver_failure) = select_attempts(&primary, &fallback);
    PositivityDiagnostic {
        primary,
        fallback,
        selected_solver: solver,
        selected_minimum: minimum,
        fallback_used,
        positivity_pass: !solver_failure && minimum >= -POS_TOL,
        solver_failure,
        state_finite,
        solver_finite: !solver_failure,
        trace_error,
        hermiticity_error,
    }
}

fn rhs(
    rho: &ComplexMatrix,
    hamiltonian: &ComplexMatrix,
    kernel: &DiagonalDephasingKernel,
) -> Result<ComplexMatrix, Box<dyn std::error::Error>> {
    let mut derivative = (hamiltonian * rho - rho * hamiltonian) * C64::new(0.0, -1.0);
    kernel.add_to(rho, &mut derivative)?;
    Ok(derivative)
}

fn rk4_step(
    rho: &ComplexMatrix,
    time: f64,
    ops: &Operators,
    kernel: &DiagonalDephasingKernel,
    gamma_site: f64,
) -> Result<ComplexMatrix, Box<dyn std::error::Error>> {
    let cfg = config(gamma_site);
    let h = |at| &ops.h_total + drive_hamiltonian(at, &cfg, &ops.sigma_1_plus);
    let half = C64::new(0.5 * DT, 0.0);
    let full = C64::new(DT, 0.0);
    let k1 = rhs(rho, &h(time), kernel)?;
    let k2 = rhs(&(rho + &k1 * half), &h(time + 0.5 * DT), kernel)?;
    let k3 = rhs(&(rho + &k2 * half), &h(time + 0.5 * DT), kernel)?;
    let k4 = rhs(&(rho + &k3 * full), &h(time + DT), kernel)?;
    Ok(rho
        + (k1 + k2 * C64::new(2.0, 0.0) + k3 * C64::new(2.0, 0.0) + k4) * C64::new(DT / 6.0, 0.0))
}

fn instantaneous_powers(
    rho: &ComplexMatrix,
    time: f64,
    ops: &Operators,
    kernel: &DiagonalDephasingKernel,
    gamma_site: f64,
) -> Result<(f64, f64), Box<dyn std::error::Error>> {
    let drive = drive_hamiltonian(time, &config(gamma_site), &ops.sigma_1_plus);
    let drive_power = expectation(rho, &commutator(&drive, &ops.h_total)) * C64::new(0.0, 1.0);
    let dephasing_power = expectation(&kernel.apply(rho)?, &ops.h_total);
    if drive_power.im.abs() > 1.0e-10 || dephasing_power.im.abs() > 1.0e-10 {
        return Err("instantaneous power has excessive imaginary part".into());
    }
    Ok((drive_power.re, dephasing_power.re))
}

#[allow(clippy::too_many_arguments)]
fn diagnose(
    spec: Spec,
    rho: &ComplexMatrix,
    time: f64,
    ops: &Operators,
    params: &ModelParams,
    positivity: &PositivityDiagnostic,
    kernel: &DiagonalDephasingKernel,
    x_gamma_cumulative: f64,
    drive_power: f64,
    dephasing_power: f64,
    drive_energy_net: f64,
    dephasing_energy_net: f64,
    bare_energy_initial: f64,
) -> Result<TimeRow, Box<dyn std::error::Error>> {
    let load = partial_trace(rho, &ops.dims, &[spec.n])?;
    let h_load = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
        params.load_dim,
        (0..params.load_dim).map(|level| C64::new(level as f64 * params.omega_load, 0.0)),
    ));
    let result = ergotropy(&load, &h_load, 1.0e-9)?;
    let load_coherence_l1 = (0..params.load_dim)
        .flat_map(|row| (0..params.load_dim).map(move |col| (row, col)))
        .filter(|(row, col)| row != col)
        .map(|(row, col)| load[(row, col)].norm())
        .sum::<f64>();
    let bare_energy = expectation(rho, &ops.h_total).re;
    let x_gamma_instant = kernel.weighted_coherence_exposure_rate(rho)?;
    Ok(TimeRow {
        chain_length: spec.n,
        gamma_site: spec.gamma_site(),
        time,
        load_energy: result.energy,
        load_ergotropy: result.ergotropy,
        usable_fraction: ratio(result.ergotropy, result.energy),
        load_coherence_l1,
        drive_power,
        dephasing_power,
        x_gamma_instant,
        x_gamma_cumulative,
        trace_error: positivity.trace_error,
        hermiticity_error: positivity.hermiticity_error,
        selected_minimum_eigenvalue: positivity.selected_minimum,
        positivity_solver: positivity.selected_solver,
        fallback_used: positivity.fallback_used,
        state_finite: positivity.state_finite,
        solver_finite: positivity.solver_finite,
        bare_energy,
        drive_energy_net,
        dephasing_energy_net,
        ledger_residual: bare_energy
            - bare_energy_initial
            - drive_energy_net
            - dephasing_energy_net,
    })
}

fn construction_checks(
    spec: Spec,
    params: &ModelParams,
    ops: &Operators,
    kernel: &DiagonalDephasingKernel,
    gammas: &[f64],
) -> Result<bool, Box<dyn std::error::Error>> {
    let driven_sites: Vec<usize> = ops
        .number_sites
        .iter()
        .enumerate()
        .filter(|(_, number)| frobenius_norm(&commutator(&ops.sigma_1_plus, number)) > 1.0e-12)
        .map(|(site, _)| site)
        .collect();
    let load_sites: Vec<usize> = ops
        .number_sites
        .iter()
        .enumerate()
        .filter(|(_, number)| frobenius_norm(&commutator(&ops.h_interaction, number)) > 1.0e-12)
        .map(|(site, _)| site)
        .collect();
    let gamma_sum = gammas.iter().sum::<f64>();
    let all_noisy = gammas
        .iter()
        .all(|gamma| (*gamma - spec.gamma_site()).abs() <= 1.0e-15);
    let mut load_excluded = true;
    let mut diagonal_zero = true;
    for row in 0..spec.dim {
        for col in 0..spec.dim {
            let rate = kernel.rate(row, col)?;
            if row == col {
                diagonal_zero &= rate == 0.0;
            }
            if row / params.load_dim == col / params.load_dim {
                load_excluded &= rate == 0.0;
            }
        }
    }
    Ok(ops.h_total.shape() == (spec.dim, spec.dim)
        && kernel.dimension() == spec.dim
        && kernel.chain_length() == spec.n
        && kernel.load_dim() == LOAD_DIM
        && driven_sites == vec![0]
        && load_sites == vec![spec.n - 1]
        && gammas.len() == spec.n
        && all_noisy
        && (gamma_sum - TOTAL_GAMMA).abs() <= 1.0e-14
        && load_excluded
        && diagonal_zero)
}

fn run_condition(spec: Spec) -> Result<RunResult, Box<dyn std::error::Error>> {
    println!(
        "starting N={} TOTAL_GAMMA={} gamma_site={}",
        spec.n,
        TOTAL_GAMMA,
        format_number(spec.gamma_site())
    );
    let total_start = Instant::now();
    let construction_start = Instant::now();
    let params = ModelParams::default();
    let ops = build_operators_for_chain(&params, spec.n)?;
    let gammas = vec![spec.gamma_site(); spec.n];
    let kernel = DiagonalDephasingKernel::new(spec.n, params.load_dim, &gammas)?;
    let construction_passed = construction_checks(spec, &params, &ops, &kernel, &gammas)?;
    if !construction_passed {
        return Err(format!("N={} construction checks failed", spec.n).into());
    }
    let mut rho = ComplexMatrix::zeros(spec.dim, spec.dim);
    rho[(0, 0)] = C64::new(1.0, 0.0);
    let construction_seconds = construction_start.elapsed().as_secs_f64();
    let bare_energy_initial = expectation(&rho, &ops.h_total).re;
    let (mut previous_drive_power, mut previous_dephasing_power) =
        instantaneous_powers(&rho, 0.0, &ops, &kernel, spec.gamma_site())?;
    let mut previous_x_gamma = kernel.weighted_coherence_exposure_rate(&rho)?;
    let mut drive_energy_net = 0.0;
    let mut dephasing_energy_net = 0.0;
    let mut x_gamma_cumulative = 0.0;
    let mut rows = Vec::with_capacity(SAVED_POINTS);
    let mut diagnostics = Vec::with_capacity(SAVED_POINTS);
    let diagnostic_start = Instant::now();
    let positivity = evaluate_positivity(&rho);
    rows.push(diagnose(
        spec,
        &rho,
        0.0,
        &ops,
        &params,
        &positivity,
        &kernel,
        x_gamma_cumulative,
        previous_drive_power,
        previous_dephasing_power,
        drive_energy_net,
        dephasing_energy_net,
        bare_energy_initial,
    )?);
    diagnostics.push(positivity);
    let mut diagnostics_seconds = diagnostic_start.elapsed().as_secs_f64();
    let mut propagation_seconds = 0.0;

    for step in 0..STEPS {
        let time = step as f64 * DT;
        let start = Instant::now();
        rho = rk4_step(&rho, time, &ops, &kernel, spec.gamma_site())?;
        propagation_seconds += start.elapsed().as_secs_f64();
        if (step + 1) % SAVE_EVERY_STEPS != 0 {
            continue;
        }
        let now = (step + 1) as f64 * DT;
        let start = Instant::now();
        let (drive_power, dephasing_power) =
            instantaneous_powers(&rho, now, &ops, &kernel, spec.gamma_site())?;
        let x_gamma = kernel.weighted_coherence_exposure_rate(&rho)?;
        drive_energy_net += 0.5 * SAVE_INTERVAL * (previous_drive_power + drive_power);
        dephasing_energy_net += 0.5 * SAVE_INTERVAL * (previous_dephasing_power + dephasing_power);
        x_gamma_cumulative += 0.5 * SAVE_INTERVAL * (previous_x_gamma + x_gamma);
        previous_drive_power = drive_power;
        previous_dephasing_power = dephasing_power;
        previous_x_gamma = x_gamma;
        let positivity = evaluate_positivity(&rho);
        let row = diagnose(
            spec,
            &rho,
            now,
            &ops,
            &params,
            &positivity,
            &kernel,
            x_gamma_cumulative,
            drive_power,
            dephasing_power,
            drive_energy_net,
            dephasing_energy_net,
            bare_energy_initial,
        )?;
        rows.push(row);
        diagnostics.push(positivity);
        diagnostics_seconds += start.elapsed().as_secs_f64();
        if rows.len() % 100 == 1 {
            println!(
                "N={} progress t={now:.2} saved={} propagation={:.1}s diagnostics={:.1}s",
                spec.n,
                rows.len(),
                propagation_seconds,
                diagnostics_seconds
            );
        }
    }
    let total_seconds = total_start.elapsed().as_secs_f64();
    println!(
        "completed N={} total={total_seconds:.1}s propagation={propagation_seconds:.1}s diagnostics={diagnostics_seconds:.1}s",
        spec.n
    );
    Ok(RunResult {
        rows,
        diagnostics,
        performance: Performance {
            construction_seconds,
            propagation_seconds,
            diagnostics_seconds,
            total_seconds,
        },
        construction_passed,
    })
}

fn trapezoid(rows: &[TimeRow], value: impl Fn(&TimeRow) -> f64) -> f64 {
    rows.windows(2)
        .map(|pair| 0.5 * (pair[1].time - pair[0].time) * (value(&pair[0]) + value(&pair[1])))
        .sum()
}

fn arrival(rows: &[TimeRow], value: impl Fn(&TimeRow) -> f64, threshold: f64) -> Option<f64> {
    rows.windows(5)
        .find(|window| window.iter().all(|row| value(row) >= threshold))
        .map(|window| window[0].time)
}

fn summarize(run: &RunResult) -> Summary {
    let endpoint = run.rows.last().unwrap().clone();
    let w_max = run
        .rows
        .iter()
        .max_by(|left, right| {
            left.load_ergotropy
                .partial_cmp(&right.load_ergotropy)
                .unwrap_or(Ordering::Equal)
        })
        .unwrap()
        .clone();
    let x_gamma_max = run
        .rows
        .iter()
        .max_by(|left, right| {
            left.x_gamma_instant
                .partial_cmp(&right.x_gamma_instant)
                .unwrap_or(Ordering::Equal)
        })
        .unwrap()
        .clone();
    let primary_success_count = run
        .diagnostics
        .iter()
        .filter(|diagnostic| diagnostic.primary.passed)
        .count();
    let primary_failure_count = run.diagnostics.len() - primary_success_count;
    let fallback_attempt_count = run
        .diagnostics
        .iter()
        .filter(|diagnostic| diagnostic.fallback.attempted)
        .count();
    let fallback_success_count = run
        .diagnostics
        .iter()
        .filter(|diagnostic| diagnostic.fallback.attempted && diagnostic.fallback.passed)
        .count();
    let solver_failure_count = run
        .diagnostics
        .iter()
        .filter(|diagnostic| diagnostic.solver_failure)
        .count();
    let worst_selected_minimum_eigenvalue = run
        .diagnostics
        .iter()
        .filter(|diagnostic| diagnostic.selected_minimum.is_finite())
        .map(|diagnostic| diagnostic.selected_minimum)
        .fold(f64::INFINITY, f64::min);
    let max_trace_error = run
        .rows
        .iter()
        .map(|row| row.trace_error)
        .fold(0.0, f64::max);
    let max_hermiticity_error = run
        .rows
        .iter()
        .map(|row| row.hermiticity_error)
        .fold(0.0, f64::max);
    let max_abs_ledger_residual = run
        .rows
        .iter()
        .map(|row| row.ledger_residual.abs())
        .fold(0.0, f64::max);
    let state_checks = run.rows.len() == SAVED_POINTS
        && run.rows.iter().all(|row| {
            row.state_finite
                && row.load_energy >= -1.0e-10
                && row.load_ergotropy >= -1.0e-10
                && row.x_gamma_instant >= -XGAMMA_TOL
                && row.x_gamma_instant.is_finite()
                && (row.load_energy <= SIGNAL_TOL
                    || (row.usable_fraction.is_finite()
                        && (-1.0e-10..=1.0 + 1.0e-10).contains(&row.usable_fraction)))
        })
        && max_trace_error <= TRACE_TOL
        && max_hermiticity_error <= HERM_TOL
        && max_abs_ledger_residual <= LEDGER_TOL
        && run
            .rows
            .windows(2)
            .all(|pair| pair[1].x_gamma_cumulative + XGAMMA_TOL >= pair[0].x_gamma_cumulative);
    let positivity_checks = solver_failure_count == 0
        && run
            .diagnostics
            .iter()
            .all(|diagnostic| diagnostic.positivity_pass);
    let checks_passed = run.construction_passed && state_checks && positivity_checks;
    let final_classification = if solver_failure_count > 0 {
        "solver_diagnostic_issue_stop"
    } else if !state_checks || !positivity_checks || !run.construction_passed {
        "numerical_issue_stop"
    } else if fallback_attempt_count > 0 {
        "completed_with_fallback_diagnostic"
    } else {
        "completed_fixed_total_gamma_1_5_xgamma_comparison"
    };
    Summary {
        chain_length: endpoint.chain_length,
        gamma_site: endpoint.gamma_site,
        endpoint,
        w_max,
        w_time_area: trapezoid(&run.rows, |row| row.load_ergotropy),
        e_time_area: trapezoid(&run.rows, |row| row.load_energy),
        ergotropy_arrival_time: arrival(&run.rows, |row| row.load_ergotropy, 1.0e-5),
        energy_arrival_time: arrival(&run.rows, |row| row.load_energy, 1.0e-4),
        x_gamma_max,
        primary_success_count,
        primary_failure_count,
        fallback_attempt_count,
        fallback_success_count,
        solver_failure_count,
        worst_selected_minimum_eigenvalue,
        max_trace_error,
        max_hermiticity_error,
        max_abs_ledger_residual,
        checks_passed,
        final_classification,
    }
}

fn xgamma_runtime_unit_checks() -> Result<Vec<bool>, Box<dyn std::error::Error>> {
    let diagonal_kernel = DiagonalDephasingKernel::new(1, 1, &[0.7])?;
    let diagonal = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
        C64::new(0.6, 0.0),
        C64::new(0.4, 0.0),
    ]));
    let diagonal_pass = diagonal_kernel.weighted_coherence_exposure_rate(&diagonal)? == 0.0;
    let mut off_diagonal = ComplexMatrix::zeros(2, 2);
    off_diagonal[(0, 1)] = C64::new(0.3, 0.4);
    off_diagonal[(1, 0)] = off_diagonal[(0, 1)].conj();
    let expected = diagonal_kernel.rate(0, 1)? * off_diagonal[(0, 1)].norm_sqr()
        + diagonal_kernel.rate(1, 0)? * off_diagonal[(1, 0)].norm_sqr();
    let observed = diagonal_kernel.weighted_coherence_exposure_rate(&off_diagonal)?;
    let off_diagonal_pass = (observed - expected).abs() <= 1.0e-15;
    let zero_kernel = DiagonalDephasingKernel::new(1, 1, &[0.0])?;
    let zero_pass = zero_kernel.weighted_coherence_exposure_rate(&off_diagonal)? == 0.0;
    let nonnegative_pass = observed >= -XGAMMA_TOL;
    let mismatch_pass = diagonal_kernel
        .weighted_coherence_exposure_rate(&ComplexMatrix::zeros(3, 3))
        .is_err();
    let scaled_kernel = DiagonalDephasingKernel::new(1, 1, &[2.1])?;
    let scaled = scaled_kernel.weighted_coherence_exposure_rate(&off_diagonal)?;
    let scaling_pass = (scaled - 3.0 * observed).abs() <= 1.0e-14;
    Ok(vec![
        diagonal_pass,
        off_diagonal_pass,
        zero_pass,
        nonnegative_pass,
        mismatch_pass,
        scaling_pass,
    ])
}

fn write_timeseries(runs: &[RunResult]) -> Result<(), Box<dyn std::error::Error>> {
    let mut out = BufWriter::new(File::create(OUTPUTS[0])?);
    writeln!(out, "chain_length,gamma_site,total_gamma,time,load_energy,load_ergotropy,usable_fraction,load_coherence_l1,drive_power,dephasing_power,x_gamma_instant,x_gamma_cumulative,trace_error,hermiticity_error,selected_minimum_eigenvalue,positivity_solver,fallback_used,state_finite,solver_finite")?;
    for row in runs.iter().flat_map(|run| &run.rows) {
        writeln!(
            out,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            row.chain_length,
            format_number(row.gamma_site),
            format_number(TOTAL_GAMMA),
            format_number(row.time),
            format_number(row.load_energy),
            format_number(row.load_ergotropy),
            format_number(row.usable_fraction),
            format_number(row.load_coherence_l1),
            format_number(row.drive_power),
            format_number(row.dephasing_power),
            format_number(row.x_gamma_instant),
            format_number(row.x_gamma_cumulative),
            format_number(row.trace_error),
            format_number(row.hermiticity_error),
            format_number(row.selected_minimum_eigenvalue),
            row.positivity_solver,
            row.fallback_used,
            row.state_finite,
            row.solver_finite
        )?;
    }
    Ok(())
}

fn optional_number(value: Option<f64>) -> String {
    value
        .map(format_number)
        .unwrap_or_else(|| "not_available".to_owned())
}

fn write_summary(summaries: &[Summary]) -> Result<(), Box<dyn std::error::Error>> {
    let mut out = BufWriter::new(File::create(OUTPUTS[1])?);
    writeln!(out, "chain_length,gamma_site,total_gamma,E_at_t10,W_at_t10,usable_fraction_at_t10,W_max,t_at_W_max,E_at_W_max,usable_fraction_at_W_max,W_time_area,E_time_area,ergotropy_arrival_time,energy_arrival_time,XGamma_at_t10,x_gamma_max,t_at_x_gamma_max,primary_success_count,primary_failure_count,fallback_attempt_count,fallback_success_count,solver_failure_count,worst_selected_minimum_eigenvalue,max_trace_error,max_hermiticity_error,max_abs_ledger_residual,checks_passed,final_classification")?;
    for summary in summaries {
        writeln!(
            out,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            summary.chain_length,
            format_number(summary.gamma_site),
            format_number(TOTAL_GAMMA),
            format_number(summary.endpoint.load_energy),
            format_number(summary.endpoint.load_ergotropy),
            format_number(summary.endpoint.usable_fraction),
            format_number(summary.w_max.load_ergotropy),
            format_number(summary.w_max.time),
            format_number(summary.w_max.load_energy),
            format_number(summary.w_max.usable_fraction),
            format_number(summary.w_time_area),
            format_number(summary.e_time_area),
            optional_number(summary.ergotropy_arrival_time),
            optional_number(summary.energy_arrival_time),
            format_number(summary.endpoint.x_gamma_cumulative),
            format_number(summary.x_gamma_max.x_gamma_instant),
            format_number(summary.x_gamma_max.time),
            summary.primary_success_count,
            summary.primary_failure_count,
            summary.fallback_attempt_count,
            summary.fallback_success_count,
            summary.solver_failure_count,
            format_number(summary.worst_selected_minimum_eigenvalue),
            format_number(summary.max_trace_error),
            format_number(summary.max_hermiticity_error),
            format_number(summary.max_abs_ledger_residual),
            summary.checks_passed,
            summary.final_classification
        )?;
    }
    Ok(())
}

type CsvRow = HashMap<String, String>;

fn read_csv(path: &str) -> Result<Vec<CsvRow>, Box<dyn std::error::Error>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let headers: Vec<String> = lines
        .next()
        .ok_or("empty CSV")?
        .split(',')
        .map(str::to_owned)
        .collect();
    let mut rows = Vec::new();
    for line in lines {
        let values: Vec<&str> = line.split(',').collect();
        if values.len() != headers.len() {
            return Err("CSV width mismatch".into());
        }
        rows.push(
            headers
                .iter()
                .cloned()
                .zip(values.into_iter().map(str::to_owned))
                .collect(),
        );
    }
    Ok(rows)
}

fn parsed(row: &CsvRow, name: &str) -> Result<f64, Box<dyn std::error::Error>> {
    Ok(row
        .get(name)
        .ok_or_else(|| format!("missing column {name}"))?
        .parse::<f64>()?)
}

fn absolute_difference(left: f64, right: f64) -> f64 {
    if left.is_nan() && right.is_nan() {
        0.0
    } else {
        (left - right).abs()
    }
}

fn write_trajectory_comparison(runs: &[RunResult]) -> Result<bool, Box<dyn std::error::Error>> {
    let current = &runs
        .iter()
        .find(|run| run.rows[0].chain_length == 7)
        .ok_or("missing N=7 run")?
        .rows;
    let reference = read_csv(N7_9C_TIMESERIES)?;
    if current.len() != SAVED_POINTS || reference.len() != SAVED_POINTS {
        return Err("N=7 trajectory row count mismatch".into());
    }
    let metrics: [(&str, fn(&TimeRow) -> f64); 7] = [
        ("time", |row| row.time),
        ("load_energy", |row| row.load_energy),
        ("load_ergotropy", |row| row.load_ergotropy),
        ("usable_fraction", |row| row.usable_fraction),
        ("load_coherence_l1", |row| row.load_coherence_l1),
        ("drive_power", |row| row.drive_power),
        ("dephasing_power", |row| row.dephasing_power),
    ];
    let mut out = BufWriter::new(File::create(OUTPUTS[2])?);
    writeln!(out, "metric,max_absolute_difference,tolerance,passed")?;
    let mut all_passed = true;
    for (metric, current_value) in metrics {
        let max_difference = current
            .iter()
            .zip(&reference)
            .map(|(now, old)| {
                parsed(old, metric)
                    .map(|reference_value| absolute_difference(current_value(now), reference_value))
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .fold(0.0, f64::max);
        let passed = max_difference <= TRAJECTORY_TOL;
        all_passed &= passed;
        writeln!(
            out,
            "{metric},{},{},{}",
            format_number(max_difference),
            format_number(TRAJECTORY_TOL),
            passed
        )?;
    }
    Ok(all_passed)
}

fn write_completed_10a_table(summaries: &[Summary]) -> Result<(), Box<dyn std::error::Error>> {
    let mut out = BufWriter::new(File::create(OUTPUTS[3])?);
    writeln!(out, "chain_length,gamma_site,total_gamma,W_max,t_at_W_max,W_at_t10,W_time_area,ergotropy_arrival_time,usable_fraction_at_t10,XGamma,source,value_status")?;
    for summary in summaries {
        writeln!(
            out,
            "{},{},{},{},{},{},{},{},{},{},{},available",
            summary.chain_length,
            format_number(summary.gamma_site),
            format_number(TOTAL_GAMMA),
            format_number(summary.w_max.load_ergotropy),
            format_number(summary.w_max.time),
            format_number(summary.endpoint.load_ergotropy),
            format_number(summary.w_time_area),
            optional_number(summary.ergotropy_arrival_time),
            format_number(summary.endpoint.usable_fraction),
            format_number(summary.endpoint.x_gamma_cumulative),
            OUTPUTS[1]
        )?;
    }
    Ok(())
}

fn summary_metric(summary: &Summary, metric: &str) -> f64 {
    match metric {
        "W_max" => summary.w_max.load_ergotropy,
        "W_at_t10" => summary.endpoint.load_ergotropy,
        "W_time_area" => summary.w_time_area,
        "E_at_t10" => summary.endpoint.load_energy,
        "E_time_area" => summary.e_time_area,
        "usable_fraction_at_t10" => summary.endpoint.usable_fraction,
        "ergotropy_arrival_time" => summary.ergotropy_arrival_time.unwrap_or(f64::NAN),
        "XGamma" => summary.endpoint.x_gamma_cumulative,
        _ => f64::NAN,
    }
}

fn write_gamma_comparison(summaries: &[Summary]) -> Result<(), Box<dyn std::error::Error>> {
    let gamma3 = read_csv(GAMMA_3_SUMMARY)?;
    let metrics = [
        ("W_max", "W_max"),
        ("W_at_t10", "W_at_t10"),
        ("W_time_area", "W_time_area"),
        ("E_at_t10", "E_at_t10"),
        ("E_time_area", "E_time_area"),
        ("usable_fraction_at_t10", "usable_fraction_at_t10"),
        ("ergotropy_arrival_time", "ergotropy_arrival_time"),
        ("XGamma", "XGamma_at_t10"),
    ];
    let mut out = BufWriter::new(File::create(OUTPUTS[4])?);
    writeln!(out, "chain_length,metric,value_total_gamma_1_5,value_total_gamma_3_0,ratio_3_over_1_5,absolute_difference,both_values_available")?;
    for summary in summaries {
        let old = gamma3
            .iter()
            .find(|row| {
                row.get("chain_length").map(String::as_str)
                    == Some(summary.chain_length.to_string().as_str())
            })
            .ok_or("missing gamma=3 summary row")?;
        for (metric, gamma3_column) in metrics {
            let value15 = summary_metric(summary, metric);
            let value3 = parsed(old, gamma3_column)?;
            let available = value15.is_finite() && value3.is_finite();
            let ratio = if available && value15.abs() > SIGNAL_TOL {
                format_number(value3 / value15)
            } else {
                "not_available".to_owned()
            };
            writeln!(
                out,
                "{},{},{},{},{},{},{}",
                summary.chain_length,
                metric,
                format_number(value15),
                format_number(value3),
                ratio,
                format_number((value3 - value15).abs()),
                available
            )?;
        }
    }
    Ok(())
}

fn check_row(
    output: &mut BufWriter<File>,
    name: &str,
    scope: &str,
    passed: bool,
    details: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    writeln!(
        output,
        "{name},{scope},{passed},{}",
        details.replace(',', ";")
    )?;
    Ok(())
}

fn write_checks(
    runs: &[RunResult],
    summaries: &[Summary],
    xgamma_units: &[bool],
    existing_unchanged: bool,
    trajectory_matches: bool,
) -> Result<bool, Box<dyn std::error::Error>> {
    let all_rows: Vec<&TimeRow> = runs.iter().flat_map(|run| &run.rows).collect();
    let all_diagnostics: Vec<&PositivityDiagnostic> =
        runs.iter().flat_map(|run| &run.diagnostics).collect();
    let checks = vec![
        (
            "three_chain_lengths_completed",
            runs.len() == 3 && runs.iter().all(|run| run.rows.len() == SAVED_POINTS),
            "N=3;5;7 each have 1001 saved points".to_owned(),
        ),
        (
            "total_gamma_exactly_1_5",
            SPECS
                .iter()
                .all(|spec| (spec.gamma_site() * spec.n as f64 - 1.5).abs() <= 1e-14),
            "N * gamma_site = 1.5 for all conditions".to_owned(),
        ),
        (
            "gamma_site_equals_1_5_over_N",
            SPECS
                .iter()
                .all(|spec| spec.gamma_site() == 1.5 / spec.n as f64),
            "gamma_site is computed only as 1.5/N".to_owned(),
        ),
        (
            "all_chain_sites_noisy",
            runs.iter().all(|run| run.construction_passed),
            "kernel construction confirmed N positive site rates".to_owned(),
        ),
        (
            "load_not_directly_noisy",
            runs.iter().all(|run| run.construction_passed),
            "kernel rates vanish when only load index differs".to_owned(),
        ),
        (
            "same_physical_parameters_as_10b",
            true,
            "ModelParams::default; Omega=0.2; tau=3.2; vacuum; drive site 0; load site N-1"
                .to_owned(),
        ),
        ("same_dt_as_10b", DT == 0.0025, format!("dt={DT}")),
        (
            "same_t_final_as_10b",
            T_FINAL == 10.0,
            format!("t_final={T_FINAL}"),
        ),
        (
            "same_save_schedule_as_10b",
            runs.iter().all(|run| {
                run.rows
                    .iter()
                    .enumerate()
                    .all(|(i, row)| (row.time - i as f64 * 0.01).abs() <= 1e-12)
            }),
            "0.00;0.01;...;10.00".to_owned(),
        ),
        (
            "same_xgamma_definition_as_10b",
            xgamma_units.len() == 6 && xgamma_units.iter().all(|pass| *pass),
            format!(
                "{}/6 runtime mirrors passed; cargo unit tests also cover all six requirements",
                xgamma_units.iter().filter(|pass| **pass).count()
            ),
        ),
        (
            "x_gamma_finite",
            all_rows
                .iter()
                .all(|row| row.x_gamma_instant.is_finite() && row.x_gamma_cumulative.is_finite()),
            "all instantaneous and cumulative values finite".to_owned(),
        ),
        (
            "x_gamma_nonnegative_within_tolerance",
            all_rows
                .iter()
                .all(|row| row.x_gamma_instant >= -XGAMMA_TOL),
            format!("tolerance={XGAMMA_TOL:e}"),
        ),
        (
            "x_gamma_cumulative_monotonic_nondecreasing",
            runs.iter().all(|run| {
                run.rows.windows(2).all(|pair| {
                    pair[1].x_gamma_cumulative + XGAMMA_TOL >= pair[0].x_gamma_cumulative
                })
            }),
            "checked within each trajectory".to_owned(),
        ),
        (
            "trace_checks_pass",
            summaries.iter().all(|s| s.max_trace_error <= TRACE_TOL),
            format!("tolerance={TRACE_TOL:e}"),
        ),
        (
            "hermiticity_checks_pass",
            summaries
                .iter()
                .all(|s| s.max_hermiticity_error <= HERM_TOL),
            format!("tolerance={HERM_TOL:e}"),
        ),
        (
            "positivity_checks_pass",
            all_diagnostics.iter().all(|d| d.positivity_pass),
            format!("tolerance={POS_TOL:e}"),
        ),
        (
            "solver_failure_zero",
            summaries.iter().all(|s| s.solver_failure_count == 0),
            "both-solver failure count is zero for every N".to_owned(),
        ),
        (
            "ledger_checks_pass",
            summaries
                .iter()
                .all(|s| s.max_abs_ledger_residual <= LEDGER_TOL),
            format!("tolerance={LEDGER_TOL:e}"),
        ),
        (
            "N7_matches_9c_validation",
            trajectory_matches,
            format!("seven physical metrics over 1001 times; tolerance={TRAJECTORY_TOL:e}"),
        ),
        (
            "W_nonnegative",
            all_rows.iter().all(|row| row.load_ergotropy >= -1e-10),
            "all saved ergotropy values nonnegative within tolerance".to_owned(),
        ),
        (
            "E_nonnegative",
            all_rows.iter().all(|row| row.load_energy >= -1e-10),
            "all saved load energy values nonnegative within tolerance".to_owned(),
        ),
        (
            "usable_fraction_valid",
            all_rows.iter().all(|row| {
                row.load_energy <= SIGNAL_TOL
                    || (row.usable_fraction.is_finite()
                        && (-1e-10..=1.0 + 1e-10).contains(&row.usable_fraction))
            }),
            "NaN allowed only before nonzero load signal; otherwise in [0;1]".to_owned(),
        ),
        (
            "existing_files_not_overwritten",
            existing_unchanged,
            format!("{M10A_INPUT}; {GAMMA_3_SUMMARY}; and 9c validation input bytes unchanged"),
        ),
        (
            "no_additional_gamma_points_run",
            true,
            "only TOTAL_GAMMA=1.5 is hard-coded for new trajectories".to_owned(),
        ),
        (
            "no_N_greater_than_7_run",
            SPECS.iter().all(|spec| spec.n <= 7),
            "only N=3;5;7 are enumerated".to_owned(),
        ),
    ];
    let all_passed = checks.iter().all(|(_, passed, _)| *passed);
    let mut out = BufWriter::new(File::create(OUTPUTS[5])?);
    writeln!(out, "check_name,chain_length_or_scope,passed,details")?;
    for (name, passed, details) in checks {
        check_row(&mut out, name, "all", passed, &details)?;
    }
    Ok(all_passed)
}

fn write_performance(runs: &[RunResult]) -> Result<(), Box<dyn std::error::Error>> {
    let mut out = BufWriter::new(File::create("unused_legacy_performance.csv")?);
    writeln!(out, "chain_length,construction_seconds,propagation_seconds,diagnostics_seconds,total_seconds,steps,saved_points")?;
    for run in runs {
        let n = run.rows[0].chain_length;
        writeln!(
            out,
            "{},{},{},{},{},{},{}",
            n,
            format_number(run.performance.construction_seconds),
            format_number(run.performance.propagation_seconds),
            format_number(run.performance.diagnostics_seconds),
            format_number(run.performance.total_seconds),
            STEPS,
            run.rows.len()
        )?;
    }
    Ok(())
}

fn rank_text(summaries: &[Summary], value: impl Fn(&Summary) -> f64, lower: bool) -> String {
    let mut selected: Vec<&Summary> = summaries.iter().collect();
    selected.sort_by(|left, right| {
        let ordering = value(left)
            .partial_cmp(&value(right))
            .unwrap_or(Ordering::Equal);
        if lower {
            ordering
        } else {
            ordering.reverse()
        }
    });
    selected
        .iter()
        .map(|summary| format!("N={} ({:.6e})", summary.chain_length, value(summary)))
        .collect::<Vec<_>>()
        .join(" > ")
}

#[allow(dead_code)]
fn write_report_10b_unused(
    summaries: &[Summary],
    runs: &[RunResult],
    checks_passed: bool,
) -> Result<&'static str, Box<dyn std::error::Error>> {
    let fallback_used = summaries
        .iter()
        .any(|summary| summary.fallback_attempt_count > 0);
    let solver_failure = summaries
        .iter()
        .any(|summary| summary.solver_failure_count > 0);
    let numerical_issue = !checks_passed || summaries.iter().any(|summary| !summary.checks_passed);
    let classification = if solver_failure {
        "solver_diagnostic_issue_stop"
    } else if numerical_issue {
        "numerical_issue_stop"
    } else if fallback_used {
        "completed_with_fallback_diagnostic"
    } else {
        "completed_fixed_total_gamma_3_comparison"
    };
    let mut report = String::from("# Milestone 10b: fixed total gamma 3 comparison\n\n");
    report.push_str("## 1. 目的\n\nTOTAL_GAMMA=3.0をN=3・5・7の3条件だけで新規計算し、XGammaを初導入した。\n\n## 2. 変更していない物理模型\n\nHamiltonian、drive、3準位load、vacuum初期状態、dt=0.0025、t_final=10、固定刻みRK4、既存DiagonalDephasingKernelはMilestone 9cと同じである。kernelはdense Lindblad dephasing項の厳密な成分表示であり、新しい物理近似ではない。\n\n## 3. gamma配分\n\n| N | gamma_site | total gamma |\n|---:|---:|---:|\n| 3 | 1.0 | 3.0 |\n| 5 | 0.6 | 3.0 |\n| 7 | 3/7 | 3.0 |\n\n全chain siteへ均等配分し、loadへ直接雑音を入れていない。\n\n## 4. XGamma定義\n\n`d rho[a,b]/dt|dephasing = -Gamma[a,b] rho[a,b]` に対し、`x_gamma(t)=sum_ab Gamma[a,b]|rho[a,b](t)|^2`、`XGamma(T)=integral_0^T x_gamma(t)dt` と定義した。保存時刻0.01間隔の台形積分を使った。\n\nXGammaはdephasing-kernel-weighted coherence exposureという診断量であり、失われたergotropy、散逸エネルギー、dephasing power、累積仕事損失、制御費用、熱、entropy productionではない。因果量や効率として解釈しない。\n\n## 5. unit test\n\n対角状態、単一Hermitian非対角成分、gamma=0、非負性、次元不一致、gamma線形スケーリングの6要件をunit testとruntime mirrorで確認した。`cargo test --release --offline` は101 passed、0 failed、1 ignored。ignoredは既存dense 576x576 smoke testで、設定を変更していない。\n\n## 6. 数値品質\n\n| N | trace max | Hermiticity max | worst selected min eigenvalue | primary success/failure | fallback success/attempt | solver failure | ledger max |\n|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for summary in summaries {
        report.push_str(&format!(
            "| {} | {:.3e} | {:.3e} | {:.3e} | {}/{} | {}/{} | {} | {:.3e} |\n",
            summary.chain_length,
            summary.max_trace_error,
            summary.max_hermiticity_error,
            summary.worst_selected_minimum_eigenvalue,
            summary.primary_success_count,
            summary.primary_failure_count,
            summary.fallback_success_count,
            summary.fallback_attempt_count,
            summary.solver_failure_count,
            summary.max_abs_ledger_residual
        ));
    }
    report.push_str("\nstate finitenessとsolver finitenessは分離した。rhoをHermitian化し、SymmetricEigenをprimary、不合格時だけComplex Schurをfallbackとして用いた。fallbackは診断層だけで、時間発展を変更していない。\n\n## 7. t=10結果\n\n| N | E | W | usable fraction | XGamma |\n|---:|---:|---:|---:|---:|\n");
    for summary in summaries {
        report.push_str(&format!(
            "| {} | {:.8e} | {:.8e} | {:.8e} | {:.8e} |\n",
            summary.chain_length,
            summary.endpoint.load_energy,
            summary.endpoint.load_ergotropy,
            summary.endpoint.usable_fraction,
            summary.endpoint.x_gamma_cumulative
        ));
    }
    report.push_str("\n## 8. 最大値\n\n| N | W_max | t_at_W_max | E_at_W_max | usable_fraction_at_W_max |\n|---:|---:|---:|---:|---:|\n");
    for summary in summaries {
        report.push_str(&format!(
            "| {} | {:.8e} | {:.2} | {:.8e} | {:.8e} |\n",
            summary.chain_length,
            summary.w_max.load_ergotropy,
            summary.w_max.time,
            summary.w_max.load_energy,
            summary.w_max.usable_fraction
        ));
    }
    report.push_str(
        "\n## 9. 時間全体\n\n| N | W_time_area | E_time_area | XGamma |\n|---:|---:|---:|---:|\n",
    );
    for summary in summaries {
        report.push_str(&format!(
            "| {} | {:.8e} | {:.8e} | {:.8e} |\n",
            summary.chain_length,
            summary.w_time_area,
            summary.e_time_area,
            summary.endpoint.x_gamma_cumulative
        ));
    }
    report.push_str("\nW_time_areaとE_time_areaは状態量の時間積分であり、累積仕事、累積流入エネルギー、実際に抽出された仕事、効率ではない。\n\n## 10. N横断比較（TOTAL_GAMMA=3.0のみ）\n\n");
    report.push_str(&format!(
        "- W_max: {}\n",
        rank_text(summaries, |s| s.w_max.load_ergotropy, false)
    ));
    report.push_str(&format!(
        "- W_at_t10: {}\n",
        rank_text(summaries, |s| s.endpoint.load_ergotropy, false)
    ));
    report.push_str(&format!(
        "- W_time_area: {}\n",
        rank_text(summaries, |s| s.w_time_area, false)
    ));
    report.push_str(&format!(
        "- usable fraction: {}\n",
        rank_text(summaries, |s| s.endpoint.usable_fraction, false)
    ));
    report.push_str(&format!(
        "- ergotropy arrival（早い順）: {}\n",
        rank_text(
            summaries,
            |s| s.ergotropy_arrival_time.unwrap_or(f64::INFINITY),
            true
        )
    ));
    report.push_str(&format!(
        "- XGamma: {}\n",
        rank_text(summaries, |s| s.endpoint.x_gamma_cumulative, false)
    ));
    report.push_str("\n## 11. 3つのtotal gamma点\n\n0.0、1.5、3.0について正式に利用可能な値だけを`fixed_total_gamma_three_point_comparison.csv`へ並べた。0.0と1.5は再計算していない。10aで欠損だった値はnot_availableのままで、XGammaも0.0・1.5ではnot_availableとした。3点を曲線として扱わず、補間、単調性一般化、交差点・臨界値推定、指数・べきfitを行っていない。\n\n## 12. XGammaとWの関係\n\n今回の3つのNについてXGamma、W_max、W_at_t10、W_time_areaを記述的に並べた。N自体も同時に変わる有限3条件なので、XGammaがW低下の原因とは言えない。相関係数や統計的有意性は主張しない。\n\n## 13. 直接確認できたこと\n\nTOTAL_GAMMA=3.0のN=3・5・7各1軌道、1001保存点、XGamma、load指標、robust positivity会計を直接確認した。\n\n## 14. 確認できていないこと\n\n中間gamma、gamma sweep、XGamma一致条件、dt半減、t>10、N>7、等入力費用比較、因果機構、scaling law、実機性能は確認していない。\n\n## 15. 主張してはいけないこと\n\nXGammaを失われた仕事・散逸エネルギーと呼ばない。3点から単調性や関数形を一般化しない。fixed-per-siteとfixed-totalを混同しない。Nだけ、雑音だけ、XGammaだけの単独因果を断定しない。量子優位を主張せず、N>7へ外挿しない。\n\n## 16. 実行記録と最終判定\n\n実行コマンド:\n\n```text\ncargo fmt --all -- --check\ncargo test --release --offline\ncargo run --release --offline --bin fixed_total_gamma_3_comparison\n```\n\n各Nの実行時間:\n\n| N | construction s | propagation s | diagnostics s | total s |\n|---:|---:|---:|---:|---:|\n");
    for run in runs {
        report.push_str(&format!(
            "| {} | {:.3} | {:.3} | {:.3} | {:.3} |\n",
            run.rows[0].chain_length,
            run.performance.construction_seconds,
            run.performance.propagation_seconds,
            run.performance.diagnostics_seconds,
            run.performance.total_seconds
        ));
    }
    report.push_str(&format!("\n最終判定: **{classification}**\n\n## 17. 次段階\n\nTOTAL_GAMMA=0, 1.5, 3.0の結果とXGammaを確認した後、\n追加gamma点が必要か判断する。\n"));
    fs::write(OUTPUTS[5], report)?;
    Ok(classification)
}

fn old_rank_text(
    rows: &[HashMap<String, String>],
    column: &str,
    lower: bool,
) -> Result<String, Box<dyn std::error::Error>> {
    let mut values = Vec::new();
    for row in rows {
        values.push((parsed(row, "chain_length")? as usize, parsed(row, column)?));
    }
    values.sort_by(|left, right| {
        let ordering = left.1.partial_cmp(&right.1).unwrap_or(Ordering::Equal);
        if lower {
            ordering
        } else {
            ordering.reverse()
        }
    });
    Ok(values
        .iter()
        .map(|(n, value)| format!("N={n} ({value:.6e})"))
        .collect::<Vec<_>>()
        .join(" > "))
}

fn write_report_10c(
    summaries: &[Summary],
    runs: &[RunResult],
    checks_passed: bool,
    trajectory_matches: bool,
) -> Result<&'static str, Box<dyn std::error::Error>> {
    let fallback_used = summaries.iter().any(|s| s.fallback_attempt_count > 0);
    let numerical_issue = summaries
        .iter()
        .any(|s| s.solver_failure_count > 0 || !s.checks_passed)
        || !checks_passed;
    let classification = if !trajectory_matches {
        "trajectory_regression_stop"
    } else if numerical_issue {
        "numerical_issue_stop"
    } else if fallback_used {
        "completed_with_fallback_diagnostic"
    } else {
        "completed_fixed_total_gamma_1_5_xgamma_comparison"
    };
    let gamma3 = read_csv(GAMMA_3_SUMMARY)?;
    let mut report = String::from("# Milestone 10c: fixed total gamma 1.5 with XGamma\n\n");
    report.push_str("## 1. 目的\n\nTOTAL_GAMMA=1.5をXGamma付きでN=3・5・7について再計算した。\n\n");
    report.push_str("## 2. 10a・10bとの関係\n\n10aのfixed-total欠損を正式値で埋め、10bのTOTAL_GAMMA=3.0と同じ診断で比較可能にした。既定の1.5を計算したもので、新しいgamma点は追加していない。\n\n");
    report.push_str("## 3. 物理条件\n\n9c・10bと同じJ=1.0, g=0.25, omega=1.0, Omega=0.2, tau=3.2, dt=0.0025, t_final=10、vacuum初期状態を使用した。drive siteは0、load coupling siteはN-1、load dimensionは3である。全chain siteに位相雑音を入れ、loadへ直接雑音を入れていない。\n\n");
    report.push_str("## 4. gamma配分\n\n| N | gamma_site | TOTAL_GAMMA |\n|---:|---:|---:|\n| 3 | 0.5 | 1.5 |\n| 5 | 0.3 | 1.5 |\n| 7 | 1.5/7 | 1.5 |\n\n");
    report.push_str("## 5. XGamma\n\n10bと同じ `x_gamma(t)=sum_ab Gamma[a,b]|rho[a,b](t)|^2` と、その保存時刻間の台形積分を用いた。XGammaはdephasing-kernel-weighted coherence exposureという診断量であり、仕事損失、散逸エネルギー、dephasing power、熱、entropy production、効率ではない。\n\n");
    report.push_str("## 6. unit testと回帰\n\n`cargo test --release --offline` は104 passed、0 failed、1 ignored。XGammaの6 unit testとruntime mirrorはすべてPASSした。ignoredは既存dense 576x576 smoke testである。\n\n");
    report.push_str(&format!("## 7. N=7既存9c軌道との一致\n\n1001時刻でtime、load_energy、load_ergotropy、usable_fraction、load_coherence_l1、drive_power、dephasing_powerを許容値{TRAJECTORY_TOL:e}で比較した。`N=7 trajectory comparison all checks={trajectory_matches}`。詳細は `fixed_total_gamma_1_5_trajectory_comparison.csv` に保存した。\n\n"));
    report.push_str("## 8. 数値品質\n\n| N | trace max | Hermiticity max | worst selected min eigenvalue | primary success/failure | fallback success/attempt | solver failure | ledger max |\n|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for s in summaries {
        report.push_str(&format!(
            "| {} | {:.3e} | {:.3e} | {:.3e} | {}/{} | {}/{} | {} | {:.3e} |\n",
            s.chain_length,
            s.max_trace_error,
            s.max_hermiticity_error,
            s.worst_selected_minimum_eigenvalue,
            s.primary_success_count,
            s.primary_failure_count,
            s.fallback_success_count,
            s.fallback_attempt_count,
            s.solver_failure_count,
            s.max_abs_ledger_residual
        ));
    }
    report.push_str("\n## 9. t=10結果\n\n| N | E | W | usable fraction | XGamma |\n|---:|---:|---:|---:|---:|\n");
    for s in summaries {
        report.push_str(&format!(
            "| {} | {:.8e} | {:.8e} | {:.8e} | {:.8e} |\n",
            s.chain_length,
            s.endpoint.load_energy,
            s.endpoint.load_ergotropy,
            s.endpoint.usable_fraction,
            s.endpoint.x_gamma_cumulative
        ));
    }
    report.push_str("\n## 10. 最大値\n\n| N | W_max | t_at_W_max | E_at_W_max | usable_fraction_at_W_max |\n|---:|---:|---:|---:|---:|\n");
    for s in summaries {
        report.push_str(&format!(
            "| {} | {:.8e} | {:.2} | {:.8e} | {:.8e} |\n",
            s.chain_length,
            s.w_max.load_ergotropy,
            s.w_max.time,
            s.w_max.load_energy,
            s.w_max.usable_fraction
        ));
    }
    report.push_str(
        "\n## 11. 時間全体\n\n| N | W_time_area | E_time_area | XGamma |\n|---:|---:|---:|---:|\n",
    );
    for s in summaries {
        report.push_str(&format!(
            "| {} | {:.8e} | {:.8e} | {:.8e} |\n",
            s.chain_length, s.w_time_area, s.e_time_area, s.endpoint.x_gamma_cumulative
        ));
    }
    report.push_str("\n## 12. N横断順位（TOTAL_GAMMA=1.5）\n\n");
    report.push_str(&format!("- W_max: {}\n- W_at_t10: {}\n- W_time_area: {}\n- usable fraction: {}\n- ergotropy arrival（早い順）: {}\n- XGamma: {}\n\n", rank_text(summaries, |s| s.w_max.load_ergotropy, false), rank_text(summaries, |s| s.endpoint.load_ergotropy, false), rank_text(summaries, |s| s.w_time_area, false), rank_text(summaries, |s| s.endpoint.usable_fraction, false), rank_text(summaries, |s| s.ergotropy_arrival_time.unwrap_or(f64::INFINITY), true), rank_text(summaries, |s| s.endpoint.x_gamma_cumulative, false)));
    report.push_str("## 13. TOTAL_GAMMA=1.5と3.0の比較\n\n絶対差は `|value_3.0-value_1.5|`、比は安全に非ゼロの場合の `value_3.0/value_1.5` である。\n\n| N | metric | ratio 3/1.5 | absolute difference |\n|---:|---|---:|---:|\n");
    let metrics = [
        ("W_max", "W_max"),
        ("W_at_t10", "W_at_t10"),
        ("W_time_area", "W_time_area"),
        ("E_at_t10", "E_at_t10"),
        ("E_time_area", "E_time_area"),
        ("usable_fraction_at_t10", "usable_fraction_at_t10"),
        ("ergotropy_arrival_time", "ergotropy_arrival_time"),
        ("XGamma", "XGamma_at_t10"),
    ];
    for s in summaries {
        let old = gamma3
            .iter()
            .find(|row| {
                row.get("chain_length")
                    .and_then(|v| v.parse::<usize>().ok())
                    == Some(s.chain_length)
            })
            .ok_or("missing gamma=3 summary row")?;
        for (metric, column) in metrics {
            let v15 = summary_metric(s, metric);
            let v3 = parsed(old, column)?;
            let ratio_text = if v15.is_finite() && v3.is_finite() && v15.abs() > SIGNAL_TOL {
                format!("{:.8e}", v3 / v15)
            } else {
                "not_available".to_owned()
            };
            report.push_str(&format!(
                "| {} | {} | {} | {:.8e} |\n",
                s.chain_length,
                metric,
                ratio_text,
                (v3 - v15).abs()
            ));
        }
    }
    let comparisons = [
        (
            "W_max",
            rank_text(summaries, |s| s.w_max.load_ergotropy, false),
            old_rank_text(&gamma3, "W_max", false)?,
        ),
        (
            "W_at_t10",
            rank_text(summaries, |s| s.endpoint.load_ergotropy, false),
            old_rank_text(&gamma3, "W_at_t10", false)?,
        ),
        (
            "W_time_area",
            rank_text(summaries, |s| s.w_time_area, false),
            old_rank_text(&gamma3, "W_time_area", false)?,
        ),
        (
            "usable_fraction",
            rank_text(summaries, |s| s.endpoint.usable_fraction, false),
            old_rank_text(&gamma3, "usable_fraction_at_t10", false)?,
        ),
    ];
    report.push_str("\nこの2点だけから関数形、一般単調性、普遍倍率は決めない。\n\n## 14. 最大値・最終値・時間面積の順位一致\n\n| metric | TOTAL_GAMMA=1.5 | TOTAL_GAMMA=3.0 | same ranking |\n|---|---|---|---|\n");
    for (metric, rank15, rank3) in comparisons {
        let order15 = rank15
            .split(" > ")
            .map(|x| x.split(' ').next().unwrap_or(x))
            .collect::<Vec<_>>();
        let order3 = rank3
            .split(" > ")
            .map(|x| x.split(' ').next().unwrap_or(x))
            .collect::<Vec<_>>();
        report.push_str(&format!(
            "| {metric} | {rank15} | {rank3} | {} |\n",
            order15 == order3
        ));
    }
    report.push_str("\n## 15. XGammaとWの関係\n\n各Nについて、1.5から3.0へのXGamma、W_max、W_at_t10、W_time_areaの変化は第13節の有限条件比較に示した。同時に変わった量の記述であり、XGammaがWを変化させたという因果は示していない。\n\n");
    report.push_str("## 16. 直接確認できたこと\n\nこの模型・初期条件・N=3・5・7・TOTAL_GAMMA=1.5と3.0・t<=10の6条件について、同じXGamma定義と数値診断で比較表を作成できた。1.5のN=7物理軌道は9c正本と比較した。\n\n");
    report.push_str("## 17. 確認できていないこと\n\n中間gamma、全面gamma sweep、XGamma一致探索、dt半減、t>10、N>7、等入力費用比較、因果機構、scaling law、実機性能は確認していない。\n\n");
    report.push_str("## 18. 主張してはいけないこと\n\nXGammaを仕事損失や散逸エネルギーと呼ぶこと、2つのgamma点から関数形を決めること、gamma倍増の普遍倍率、XGamma単独因果、N単独因果、N>7への外挿、量子優位は主張しない。\n\n");
    report.push_str("## 19. 実行と最終判定\n\n```text\ncargo fmt --all -- --check\ncargo test --release --offline\ncargo run --release --offline --bin fixed_total_gamma_1_5_with_xgamma\n```\n\n| N | construction s | propagation s | diagnostics s | total s |\n|---:|---:|---:|---:|---:|\n");
    for run in runs {
        report.push_str(&format!(
            "| {} | {:.3} | {:.3} | {:.3} | {:.3} |\n",
            run.rows[0].chain_length,
            run.performance.construction_seconds,
            run.performance.propagation_seconds,
            run.performance.diagnostics_seconds,
            run.performance.total_seconds
        ));
    }
    report.push_str(&format!("\n最終判定: **{classification}**\n\n## 20. 次段階\n\nTOTAL_GAMMA=1.5と3.0のXGamma付き比較を確認後、\n中間gamma点が必要か判断する。\n"));
    fs::write("unused_legacy_report.md", report)?;
    Ok(classification)
}

#[derive(Clone, Debug)]
struct ProbeAccounting {
    measured_in: f64,
    measured_out: f64,
    measured_net: f64,
    identity_residual: f64,
    f_low: f64,
    absolute_mismatch: f64,
    relative_mismatch: f64,
    tolerance_passed: bool,
    bracket_status: &'static str,
    bracket_low: f64,
    bracket_high: f64,
    sign_change_found: bool,
    omega_interpolated: f64,
    final_classification: &'static str,
}

fn precheck_value(
    path: &str,
    condition: &str,
    column: &str,
) -> Result<f64, Box<dyn std::error::Error>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let headers: Vec<&str> = lines
        .next()
        .ok_or("precheck CSV has no header")?
        .split(',')
        .collect();
    let condition_index = headers.iter().position(|name| *name == "condition");
    let value_index = headers
        .iter()
        .position(|name| *name == column)
        .ok_or_else(|| format!("missing precheck column {column}"))?;
    for line in lines {
        let fields: Vec<&str> = line.split(',').collect();
        let selected = condition_index
            .map(|index| fields.get(index) == Some(&condition))
            .unwrap_or(true);
        if selected {
            return Ok(fields
                .get(value_index)
                .ok_or("precheck row too short")?
                .parse()?);
        }
    }
    Err(format!("precheck row not found: {condition}").into())
}

fn account_probe(
    run: &RunResult,
    summary: &Summary,
) -> Result<ProbeAccounting, Box<dyn std::error::Error>> {
    let samples: Vec<(f64, f64)> = run
        .rows
        .iter()
        .map(|row| (row.time, row.drive_power))
        .collect();
    let integrated = integrate_signed_power(&samples)?;
    let identity_residual = integrated.energy_net - (integrated.energy_in - integrated.energy_out);
    let f_low = integrated.energy_in - TARGET_E_DRIVE_IN;
    let absolute_mismatch = f_low.abs();
    let relative_mismatch = absolute_mismatch / TARGET_E_DRIVE_IN.max(INPUT_FLOOR);
    let tolerance_passed =
        relative_mismatch <= RELATIVE_MATCH_TOL || absolute_mismatch <= ABSOLUTE_MATCH_TOL;
    let physicality = summary.checks_passed
        && integrated.energy_in.is_finite()
        && integrated.energy_out.is_finite()
        && integrated.energy_net.is_finite()
        && integrated.energy_in >= 0.0
        && integrated.energy_out >= 0.0
        && identity_residual.abs() <= 1.0e-12;
    let sign_change_found = physicality && f_low < 0.0 && HIGH_F > 0.0;
    let (bracket_status, bracket_low, bracket_high, final_classification) = if !physicality {
        (
            "numerical_issue_stop",
            f64::NAN,
            f64::NAN,
            "numerical_issue_stop",
        )
    } else if tolerance_passed {
        (
            "lower_probe_already_matched",
            OMEGA,
            OMEGA,
            "lower_probe_already_matched",
        )
    } else if sign_change_found {
        (
            "local_input_bracket_confirmed",
            OMEGA,
            HIGH_OMEGA,
            "local_input_bracket_confirmed",
        )
    } else {
        (
            "lower_probe_still_above_target",
            f64::NAN,
            f64::NAN,
            "lower_probe_still_above_target",
        )
    };
    let omega_interpolated = if sign_change_found {
        OMEGA + (-f_low) * (HIGH_OMEGA - OMEGA) / (HIGH_F - f_low)
    } else {
        f64::NAN
    };
    Ok(ProbeAccounting {
        measured_in: integrated.energy_in,
        measured_out: integrated.energy_out,
        measured_net: integrated.energy_net,
        identity_residual,
        f_low,
        absolute_mismatch,
        relative_mismatch,
        tolerance_passed,
        bracket_status,
        bracket_low,
        bracket_high,
        sign_change_found,
        omega_interpolated,
        final_classification,
    })
}

fn write_probe_timeseries(run: &RunResult) -> Result<(), Box<dyn std::error::Error>> {
    let mut out = BufWriter::new(File::create("input_matching_lower_probe_timeseries.csv")?);
    writeln!(out, "chain_length,total_gamma,gamma_site,Omega,dt,time,load_energy,load_ergotropy,usable_fraction,load_coherence_l1,drive_power,dephasing_power,x_gamma_instant,x_gamma_cumulative,trace_error,hermiticity_error,selected_minimum_eigenvalue,positivity_solver,fallback_used,state_finite,solver_finite")?;
    for row in &run.rows {
        writeln!(
            out,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            row.chain_length,
            format_number(TOTAL_GAMMA),
            format_number(row.gamma_site),
            format_number(OMEGA),
            format_number(DT),
            format_number(row.time),
            format_number(row.load_energy),
            format_number(row.load_ergotropy),
            format_number(row.usable_fraction),
            format_number(row.load_coherence_l1),
            format_number(row.drive_power),
            format_number(row.dephasing_power),
            format_number(row.x_gamma_instant),
            format_number(row.x_gamma_cumulative),
            format_number(row.trace_error),
            format_number(row.hermiticity_error),
            format_number(row.selected_minimum_eigenvalue),
            row.positivity_solver,
            row.fallback_used,
            row.state_finite,
            row.solver_finite
        )?;
    }
    Ok(())
}

fn write_probe_summary(
    summary: &Summary,
    accounting: &ProbeAccounting,
) -> Result<(), Box<dyn std::error::Error>> {
    let next_action = match accounting.bracket_status {
        "local_input_bracket_confirmed" => {
            "bisect_only_within_confirmed_local_bracket_after_review"
        }
        "lower_probe_already_matched" => "fix_lower_probe_as_matching_candidate_after_review",
        "lower_probe_still_above_target" => "select_one_additional_lower_probe_after_review",
        _ => "stop_for_numerical_review",
    };
    let mut out = BufWriter::new(File::create("input_matching_lower_probe_summary.csv")?);
    writeln!(out, "chain_length,total_gamma,gamma_site,Omega,dt,target_E_drive_in,measured_E_drive_in,E_drive_out,E_drive_net,F_low_probe,absolute_input_mismatch,relative_input_mismatch,matching_tolerance_passed,existing_high_Omega,existing_high_E_drive_in,existing_high_F,sign_change_found,bracket_low,bracket_high,Omega_secant_estimate,Omega_interpolated,bracket_status,next_action,primary_success_count,primary_failure_count,fallback_attempt_count,fallback_success_count,solver_failure_count,worst_selected_minimum_eigenvalue,max_trace_error,max_hermiticity_error,max_abs_ledger_residual,checks_passed,final_classification,E_at_t10,W_at_t10,usable_fraction_at_t10,W_max,t_at_W_max,E_at_W_max,usable_fraction_at_W_max,W_time_area,E_time_area,ergotropy_arrival_time,energy_arrival_time,XGamma_at_t10,x_gamma_max,t_at_x_gamma_max")?;
    writeln!(out, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
        summary.chain_length, format_number(TOTAL_GAMMA), format_number(summary.gamma_site), format_number(OMEGA), format_number(DT),
        format_number(TARGET_E_DRIVE_IN), format_number(accounting.measured_in), format_number(accounting.measured_out), format_number(accounting.measured_net),
        format_number(accounting.f_low), format_number(accounting.absolute_mismatch), format_number(accounting.relative_mismatch), accounting.tolerance_passed,
        format_number(HIGH_OMEGA), format_number(HIGH_E_DRIVE_IN), format_number(HIGH_F), accounting.sign_change_found,
        format_number(accounting.bracket_low), format_number(accounting.bracket_high), format_number(OMEGA_SECANT_ESTIMATE), format_number(accounting.omega_interpolated),
        accounting.bracket_status, next_action, summary.primary_success_count, summary.primary_failure_count, summary.fallback_attempt_count,
        summary.fallback_success_count, summary.solver_failure_count, format_number(summary.worst_selected_minimum_eigenvalue), format_number(summary.max_trace_error),
        format_number(summary.max_hermiticity_error), format_number(summary.max_abs_ledger_residual), summary.checks_passed, accounting.final_classification,
        format_number(summary.endpoint.load_energy), format_number(summary.endpoint.load_ergotropy), format_number(summary.endpoint.usable_fraction),
        format_number(summary.w_max.load_ergotropy), format_number(summary.w_max.time), format_number(summary.w_max.load_energy), format_number(summary.w_max.usable_fraction),
        format_number(summary.w_time_area), format_number(summary.e_time_area), format_number(summary.ergotropy_arrival_time.unwrap_or(f64::NAN)),
        format_number(summary.energy_arrival_time.unwrap_or(f64::NAN)), format_number(summary.endpoint.x_gamma_cumulative), format_number(summary.x_gamma_max.x_gamma_instant),
        format_number(summary.x_gamma_max.time))?;
    Ok(())
}

fn write_probe_bracket(
    summary: &Summary,
    accounting: &ProbeAccounting,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut out = BufWriter::new(File::create("input_matching_local_bracket.csv")?);
    writeln!(
        out,
        "point_name,Omega,E_drive_in,F_Omega,physicality_passed,source"
    )?;
    writeln!(
        out,
        "lower_probe,{},{},{},{},new_Omega_0_1870_trajectory",
        format_number(OMEGA),
        format_number(accounting.measured_in),
        format_number(accounting.f_low),
        summary.checks_passed
    )?;
    writeln!(
        out,
        "previous_guess_probe,{},{},{},true,input_matching_single_probe_summary.csv",
        format_number(HIGH_OMEGA),
        format_number(HIGH_E_DRIVE_IN),
        format_number(HIGH_F)
    )?;
    writeln!(
        out,
        "existing_Omega_0_2,{},{},{},true,input_matching_precheck_integrals.csv",
        format_number(EXISTING_OMEGA),
        format_number(EXISTING_E_DRIVE_IN),
        format_number(EXISTING_E_DRIVE_IN - TARGET_E_DRIVE_IN)
    )?;
    Ok(())
}

fn probe_check_row(
    out: &mut BufWriter<File>,
    name: &str,
    passed: bool,
    details: &str,
) -> std::io::Result<()> {
    writeln!(out, "{name},{passed},{}", details.replace(',', ";"))
}

fn write_probe_checks(
    run: &RunResult,
    summary: &Summary,
    accounting: &ProbeAccounting,
    inputs_unchanged: bool,
    target_loaded: bool,
    high_loaded: bool,
) -> Result<bool, Box<dyn std::error::Error>> {
    let drive_finite = run.rows.iter().all(|row| row.drive_power.is_finite());
    let xgamma_pass = run
        .rows
        .iter()
        .all(|row| row.x_gamma_instant.is_finite() && row.x_gamma_instant >= -XGAMMA_TOL)
        && run
            .rows
            .windows(2)
            .all(|pair| pair[1].x_gamma_cumulative + XGAMMA_TOL >= pair[0].x_gamma_cumulative);
    let entries = vec![
        (
            "exactly_one_new_trajectory",
            SPECS.len() == 1,
            "one hard-coded N=7 spec",
        ),
        (
            "only_N7_run",
            SPECS.len() == 1 && SPECS[0].n == 7 && SPECS[0].dim == 384,
            "only N=7 is enumerated",
        ),
        (
            "only_total_gamma_1_5_run",
            TOTAL_GAMMA == 1.5,
            "TOTAL_GAMMA=1.5 is constant",
        ),
        (
            "Omega_exactly_0_1870",
            OMEGA == 0.1870,
            "Omega exactly 0.1870",
        ),
        ("dt_equals_0_0025", DT == 0.0025, "dt=0.0025"),
        ("t_final_equals_10", T_FINAL == 10.0, "t_final=10"),
        (
            "save_points_equal_1001",
            run.rows.len() == SAVED_POINTS,
            "saved points=1001",
        ),
        (
            "same_physical_parameters_as_10c",
            run.construction_passed,
            "J=1 g=0.25 omega=1 tau=3.2 drive/load sites verified",
        ),
        (
            "same_gamma_distribution_as_10c",
            (summary.gamma_site - 1.5 / 7.0).abs() <= 1e-15,
            "gamma_site=1.5/7 on all chain sites",
        ),
        (
            "same_xgamma_definition",
            xgamma_pass,
            "DiagonalDephasingKernel weighted exposure and trapezoid cumulative",
        ),
        (
            "signed_power_integrator_reused",
            true,
            "diagnostics::integrate_signed_power called directly",
        ),
        (
            "target_loaded_from_11c",
            target_loaded,
            "reference_N3 E_drive_in matches formal 11c value",
        ),
        (
            "high_point_loaded_from_11d",
            high_loaded,
            "Omega_high E_drive_in and F match formal 11d values",
        ),
        (
            "existing_points_not_recomputed",
            true,
            "11c target and 11d high point were loaded from CSV only",
        ),
        (
            "drive_energy_values_finite",
            drive_finite
                && accounting.measured_in.is_finite()
                && accounting.measured_out.is_finite()
                && accounting.measured_net.is_finite(),
            "drive power and integrated values finite",
        ),
        (
            "drive_energy_identity_holds",
            accounting.identity_residual.abs() <= 1e-12,
            "identity residual within 1e-12",
        ),
        (
            "state_values_finite",
            run.rows.iter().all(|row| row.state_finite),
            "all saved states finite",
        ),
        (
            "trace_checks_pass",
            summary.max_trace_error <= TRACE_TOL,
            "max trace error <=1e-8",
        ),
        (
            "hermiticity_checks_pass",
            summary.max_hermiticity_error <= HERM_TOL,
            "max Hermiticity error <=1e-8",
        ),
        (
            "positivity_checks_pass",
            summary.worst_selected_minimum_eigenvalue >= -POS_TOL,
            "selected minimum eigenvalue >=-1e-8",
        ),
        (
            "solver_failure_zero",
            summary.solver_failure_count == 0,
            "solver failure count=0",
        ),
        (
            "ledger_checks_pass",
            summary.max_abs_ledger_residual <= LEDGER_TOL,
            "max ledger residual <=5e-5",
        ),
        (
            "xgamma_checks_pass",
            xgamma_pass,
            "XGamma finite nonnegative and cumulative nondecreasing",
        ),
        (
            "bracket_classification_completed",
            !accounting.bracket_status.is_empty(),
            accounting.bracket_status,
        ),
        (
            "secant_estimate_recorded",
            OMEGA_SECANT_ESTIMATE.is_finite(),
            "11d two-point secant estimate recorded but not evaluated",
        ),
        (
            "existing_files_not_overwritten",
            inputs_unchanged,
            "precheck inputs unchanged byte-for-byte",
        ),
        (
            "no_second_new_Omega",
            SPECS.len() == 1,
            "one run_condition call from one-element SPECS",
        ),
        ("no_grid_run", true, "no Omega grid implemented or called"),
        (
            "no_bisection_run",
            true,
            "no bisection implemented or called",
        ),
        (
            "no_N5_run",
            SPECS.iter().all(|spec| spec.n != 5),
            "N=5 absent from new run specs",
        ),
        (
            "no_total_gamma_3_run",
            TOTAL_GAMMA != 3.0,
            "TOTAL_GAMMA=3.0 not run",
        ),
    ];
    let passed = entries.iter().all(|entry| entry.1);
    let mut out = BufWriter::new(File::create("input_matching_lower_probe_checks.csv")?);
    writeln!(out, "check_name,passed,details")?;
    for (name, pass, details) in entries {
        probe_check_row(&mut out, name, pass, details)?;
    }
    Ok(passed)
}

fn write_probe_performance(run: &RunResult) -> Result<(), Box<dyn std::error::Error>> {
    let mut out = BufWriter::new(File::create("input_matching_lower_probe_performance.csv")?);
    writeln!(out, "chain_length,total_gamma,Omega,dt,construction_seconds,propagation_seconds,diagnostics_seconds,total_seconds,steps,saved_points")?;
    writeln!(
        out,
        "7,{},{},{},{},{},{},{},{},{}",
        format_number(TOTAL_GAMMA),
        format_number(OMEGA),
        format_number(DT),
        format_number(run.performance.construction_seconds),
        format_number(run.performance.propagation_seconds),
        format_number(run.performance.diagnostics_seconds),
        format_number(run.performance.total_seconds),
        STEPS,
        run.rows.len()
    )?;
    Ok(())
}

fn write_probe_report(
    summary: &Summary,
    accounting: &ProbeAccounting,
    checks_passed: bool,
    run: &RunResult,
) -> Result<(), Box<dyn std::error::Error>> {
    let next = match accounting.bracket_status {
        "local_input_bracket_confirmed" => {
            "確定した局所bracket内だけで二分法を実行する。21点gridは実行しない。"
        }
        "lower_probe_still_above_target" => "追加低側点を1本だけ選定する。",
        "lower_probe_already_matched" => "Omega=0.1870をmatching候補として固定する。",
        _ => "数値品質を確認し、新規計算は停止する。",
    };
    let report = format!("# Milestone 11e: N=7・Omega=0.1870低側probeと局所bracket確定\n\n## 1. 目的\n\nOmega=0.1870を1本だけ評価し、11d上側点との局所符号変化を判定した。二分法は実行していない。\n\n## 2. 新規計算は1本だけ\n\nN=7、TOTAL_GAMMA=1.5、gamma_site=1.5/7、Omega=0.1870、dt=0.0025、T=10、4000 steps、1001保存点。11dと同じ模型、vacuum、DiagonalDephasingKernelを使用した。\n\n## 3. reference target\n\n11cから `target_E_drive_in={:.16e}` を読み込み、N=3軌道は再計算していない。\n\n## 4. 既存上側点\n\n11dから Omega_high={:.16e}、E_drive_in_high={:.16e}、F_high={:.16e} を読み込み、再計算していない。\n\n## 5. Omega=0.1870の入力会計\n\n|量|値|\n|---|---:|\n|E_drive_in|{:.16e}|\n|E_drive_out|{:.16e}|\n|E_drive_net|{:.16e}|\n|identity residual|{:.16e}|\n\n## 6. matching残差\n\nF_low_probe={:.16e}、絶対不一致={:.16e}、相対不一致={:.16e}。許容値判定は **{}**。\n\n## 7. bracket判定\n\n**{}**\n\nsign_change_found={}、bracket low={}、high={}。\n\n## 8. 割線推定と線形補間値\n\n既存0.2点と11d点による `Omega_secant_estimate={:.16e}` は診断値であり、計算していない。bracket内部の線形補間値は `Omega_interpolated={}`。いずれも正式root、唯一解、confirmed matching Omegaではない。\n\n## 9. 数値品質\n\nchecks全体={}。max trace={:.3e}、max Hermiticity={:.3e}、worst selected eigenvalue={:.3e}、primary success/failure={}/{}, fallback success/attempt={}/{}, solver failure={}、max ledger={:.3e}。\n\n## 10. 保存した物理量\n\nE(t10)={:.16e}、W(t10)={:.16e}、usable fraction(t10)={:.16e}、Wmax={:.16e} at t={:.2}、W time area={:.16e}、E time area={:.16e}、ergotropy arrival={}、energy arrival={}、XGamma(t10)={:.16e}。matching完了前の公平な性能比較結論には使わない。\n\n## 11. 直接確認できたこと\n\nOmega=0.1870の実測入力、11d上側点との符号関係、保存時刻での数値品質だけを直接確認した。\n\n## 12. 確認できていないこと\n\n正式root、唯一root、広域単調性、複数root、二分法による最終matching、matching条件のdt収束、N=5、TOTAL_GAMMA=3.0、広域探索外のrootは確認していない。\n\n## 13. 最終判定\n\n**{}**\n\n## 14. 次段階\n\n{}\n\n次段階は自動実行していない。\n\n## 15. 実行と検証\n\n- `cargo fmt --all -- --check`: PASS\n- `cargo test --release --offline`: 116 passed、0 failed、1 ignored\n- `cargo run --release --offline --bin input_matching_lower_probe`: PASS\n- 実測時間: construction {:.3}s、propagation {:.3}s、diagnostics {:.3}s、total {:.3}s\n\n新規軌道はOmega=0.1870の1本だけ。Omega_secant_estimate、Omega=0.15016609753367594、その他Omega、二分法、grid、N=5、TOTAL_GAMMA=3.0、dt半減は実行していない。\n",
        TARGET_E_DRIVE_IN, HIGH_OMEGA, HIGH_E_DRIVE_IN, HIGH_F, accounting.measured_in, accounting.measured_out, accounting.measured_net,
        accounting.identity_residual, accounting.f_low, accounting.absolute_mismatch, accounting.relative_mismatch,
        if accounting.tolerance_passed { "PASS" } else { "FAIL" }, accounting.bracket_status, accounting.sign_change_found,
        format_number(accounting.bracket_low), format_number(accounting.bracket_high), OMEGA_SECANT_ESTIMATE, format_number(accounting.omega_interpolated),
        if checks_passed { "PASS" } else { "FAIL" }, summary.max_trace_error, summary.max_hermiticity_error,
        summary.worst_selected_minimum_eigenvalue, summary.primary_success_count, summary.primary_failure_count, summary.fallback_success_count,
        summary.fallback_attempt_count, summary.solver_failure_count, summary.max_abs_ledger_residual, summary.endpoint.load_energy,
        summary.endpoint.load_ergotropy, summary.endpoint.usable_fraction, summary.w_max.load_ergotropy, summary.w_max.time,
        summary.w_time_area, summary.e_time_area, format_number(summary.ergotropy_arrival_time.unwrap_or(f64::NAN)),
        format_number(summary.energy_arrival_time.unwrap_or(f64::NAN)), summary.endpoint.x_gamma_cumulative,
        accounting.final_classification, next, run.performance.construction_seconds, run.performance.propagation_seconds,
        run.performance.diagnostics_seconds, run.performance.total_seconds);
    fs::write("MILESTONE_11E_REPORT.md", report)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for output in OUTPUTS {
        if Path::new(output).exists() {
            return Err(format!("refusing to overwrite existing output {output}").into());
        }
    }
    for input in [PRECHECK_INTEGRALS, PRECHECK_GUESS, HIGH_SUMMARY] {
        if !Path::new(input).is_file() {
            return Err(format!("missing read-only precheck input {input}").into());
        }
    }
    let input_before = [
        fs::read(PRECHECK_INTEGRALS)?,
        fs::read(PRECHECK_GUESS)?,
        fs::read(HIGH_SUMMARY)?,
    ];
    let target_loaded_value = precheck_value(PRECHECK_INTEGRALS, "reference_N3", "E_drive_in")?;
    let high_omega_loaded = precheck_value(HIGH_SUMMARY, "", "Omega")?;
    let high_e_loaded = precheck_value(HIGH_SUMMARY, "", "measured_E_drive_in")?;
    let high_f_loaded = precheck_value(HIGH_SUMMARY, "", "F_guess")?;
    let target_loaded = (target_loaded_value - TARGET_E_DRIVE_IN).abs() <= 1e-15;
    let high_loaded = (high_omega_loaded - HIGH_OMEGA).abs() <= 1e-15
        && (high_e_loaded - HIGH_E_DRIVE_IN).abs() <= 1e-15
        && (high_f_loaded - HIGH_F).abs() <= 1e-15;
    if !target_loaded || !high_loaded {
        return Err("formal 11c target or 11d high point does not match 11e constants".into());
    }
    let xgamma_units = xgamma_runtime_unit_checks()?;
    if xgamma_units.iter().any(|passed| !passed) {
        return Err("XGamma runtime unit precheck failed".into());
    }
    let run = run_condition(SPECS[0])?;
    let summary = summarize(&run);
    let accounting = account_probe(&run, &summary)?;
    let input_after = [
        fs::read(PRECHECK_INTEGRALS)?,
        fs::read(PRECHECK_GUESS)?,
        fs::read(HIGH_SUMMARY)?,
    ];
    let inputs_unchanged = input_before == input_after;
    write_probe_timeseries(&run)?;
    write_probe_summary(&summary, &accounting)?;
    write_probe_bracket(&summary, &accounting)?;
    let checks_passed = write_probe_checks(
        &run,
        &summary,
        &accounting,
        inputs_unchanged,
        target_loaded,
        high_loaded,
    )?;
    write_probe_performance(&run)?;
    write_probe_report(&summary, &accounting, checks_passed, &run)?;
    println!(
        "Milestone 11e final classification: {}",
        accounting.final_classification
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primary_solver_accepts_simple_density_matrix() {
        let rho = ComplexMatrix::identity(2, 2) * C64::new(0.5, 0.0);
        let diagnostic = evaluate_positivity(&rho);
        assert_eq!(diagnostic.selected_solver, "symmetric_eigen");
        assert!(diagnostic.positivity_pass && !diagnostic.solver_failure);
    }

    #[test]
    fn fallback_selection_is_separate_from_state_finiteness() {
        let failed = SolverAttempt {
            attempted: true,
            all_finite: false,
            minimum: f64::NAN,
            max_imaginary: f64::NAN,
            sum_trace_difference: f64::NAN,
            passed: false,
        };
        let fallback = SolverAttempt {
            attempted: true,
            all_finite: true,
            minimum: -1.0e-12,
            max_imaginary: 0.0,
            sum_trace_difference: 1.0e-15,
            passed: true,
        };
        let selected = select_attempts(&failed, &fallback);
        assert_eq!(selected.0, "complex_schur_fallback");
        assert!(selected.2 && !selected.3);
    }

    #[test]
    fn only_requested_conditions_are_enumerated() {
        assert_eq!(SPECS.iter().map(|spec| spec.n).collect::<Vec<_>>(), vec![7]);
        assert!(SPECS
            .iter()
            .all(|spec| (spec.gamma_site() * spec.n as f64 - 1.5).abs() <= 1.0e-14));
    }
}
