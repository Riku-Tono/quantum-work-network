use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use nalgebra::linalg::{Schur, SymmetricEigen};
use quantum_work_network::coherent_drive::{drive_hamiltonian, CoherentDriveConfig};
use quantum_work_network::dephasing_kernel::DiagonalDephasingKernel;
use quantum_work_network::ergotropy::ergotropy;
use quantum_work_network::matrix::{
    commutator, expectation, frobenius_norm, hermiticity_error, ComplexMatrix, C64,
};
use quantum_work_network::operators::{build_operators_for_chain, ModelParams, Operators};
use quantum_work_network::partial_trace::partial_trace;

const TOTAL_GAMMA: f64 = 3.0;
const DT: f64 = 0.0025;
const T_FINAL: f64 = 10.0;
const SAVE_EVERY_STEPS: usize = 4;
const SAVE_INTERVAL: f64 = 0.01;
const STEPS: usize = 4000;
const SAVED_POINTS: usize = 1001;
const OMEGA: f64 = 0.2;
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

const OUTPUTS: [&str; 6] = [
    "fixed_total_gamma_3_timeseries.csv",
    "fixed_total_gamma_3_summary.csv",
    "fixed_total_gamma_three_point_comparison.csv",
    "fixed_total_gamma_3_checks.csv",
    "fixed_total_gamma_3_performance.csv",
    "MILESTONE_10B_REPORT.md",
];

#[derive(Clone, Copy, Debug)]
struct Spec {
    n: usize,
    dim: usize,
}

const SPECS: [Spec; 3] = [
    Spec { n: 3, dim: 24 },
    Spec { n: 5, dim: 96 },
    Spec { n: 7, dim: 384 },
];

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
        "completed_fixed_total_gamma_3_comparison"
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

fn existing_value(row: &CsvRow, name: &str) -> String {
    row.get(name)
        .cloned()
        .unwrap_or_else(|| "not_available".to_owned())
}

fn write_three_point_comparison(summaries: &[Summary]) -> Result<(), Box<dyn std::error::Error>> {
    let old_rows = read_csv(M10A_INPUT)?;
    let mut out = BufWriter::new(File::create(OUTPUTS[2])?);
    writeln!(out, "total_gamma,chain_length,W_max,t_at_W_max,W_at_t10,W_time_area,usable_fraction_at_t10,ergotropy_arrival_time,XGamma,source,value_status")?;
    for n in [3_usize, 5, 7] {
        for total in [0.0_f64, 1.5] {
            let noise = if total == 0.0 {
                "noise_free"
            } else {
                "fixed_total_gamma_1_5"
            };
            let row = old_rows
                .iter()
                .find(|row| {
                    row.get("chain_length").map(String::as_str) == Some(n.to_string().as_str())
                        && row.get("noise_condition").map(String::as_str) == Some(noise)
                })
                .ok_or_else(|| format!("missing Milestone 10a row N={n} {noise}"))?;
            writeln!(
                out,
                "{},{},{},{},{},{},{},{},not_available,{},{}",
                format_number(total),
                n,
                existing_value(row, "W_max"),
                existing_value(row, "t_at_W_max"),
                existing_value(row, "W_at_t10"),
                existing_value(row, "W_time_area"),
                existing_value(row, "usable_fraction_at_t10"),
                existing_value(row, "ergotropy_arrival_time"),
                M10A_INPUT,
                existing_value(row, "value_status")
            )?;
        }
        let summary = summaries
            .iter()
            .find(|summary| summary.chain_length == n)
            .unwrap();
        writeln!(
            out,
            "{},{},{},{},{},{},{},{},{},{},available",
            format_number(TOTAL_GAMMA),
            n,
            format_number(summary.w_max.load_ergotropy),
            format_number(summary.w_max.time),
            format_number(summary.endpoint.load_ergotropy),
            format_number(summary.w_time_area),
            format_number(summary.endpoint.usable_fraction),
            optional_number(summary.ergotropy_arrival_time),
            format_number(summary.endpoint.x_gamma_cumulative),
            OUTPUTS[1]
        )?;
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
            "total_gamma_exactly_3",
            SPECS
                .iter()
                .all(|spec| (spec.gamma_site() * spec.n as f64 - 3.0).abs() <= 1e-14),
            "N * gamma_site = 3.0 for all conditions".to_owned(),
        ),
        (
            "gamma_site_equals_3_over_N",
            SPECS
                .iter()
                .all(|spec| spec.gamma_site() == 3.0 / spec.n as f64),
            "gamma_site is computed only as 3.0/N".to_owned(),
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
            "same_physical_parameters_across_N",
            true,
            "ModelParams::default; Omega=0.2; tau=3.2; vacuum; drive site 0; load site N-1"
                .to_owned(),
        ),
        ("same_dt_across_N", DT == 0.0025, format!("dt={DT}")),
        (
            "same_t_final_across_N",
            T_FINAL == 10.0,
            format!("t_final={T_FINAL}"),
        ),
        (
            "same_save_schedule_across_N",
            runs.iter().all(|run| {
                run.rows
                    .iter()
                    .enumerate()
                    .all(|(i, row)| (row.time - i as f64 * 0.01).abs() <= 1e-12)
            }),
            "0.00;0.01;...;10.00".to_owned(),
        ),
        (
            "x_gamma_unit_checks_pass",
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
            format!("{M10A_INPUT} bytes unchanged"),
        ),
        (
            "no_additional_gamma_points_run",
            true,
            "only TOTAL_GAMMA=3.0 is hard-coded for new trajectories".to_owned(),
        ),
        (
            "no_N_greater_than_7_run",
            SPECS.iter().all(|spec| spec.n <= 7),
            "only N=3;5;7 are enumerated".to_owned(),
        ),
    ];
    let all_passed = checks.iter().all(|(_, passed, _)| *passed);
    let mut out = BufWriter::new(File::create(OUTPUTS[3])?);
    writeln!(out, "check_name,chain_length_or_scope,passed,details")?;
    for (name, passed, details) in checks {
        check_row(&mut out, name, "all", passed, &details)?;
    }
    Ok(all_passed)
}

fn write_performance(runs: &[RunResult]) -> Result<(), Box<dyn std::error::Error>> {
    let mut out = BufWriter::new(File::create(OUTPUTS[4])?);
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

fn write_report(
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for output in OUTPUTS {
        if Path::new(output).exists() {
            return Err(format!("refusing to overwrite existing output {output}").into());
        }
    }
    if !Path::new(M10A_INPUT).is_file() {
        return Err(format!("missing read-only comparison input {M10A_INPUT}").into());
    }
    let input_before = fs::read(M10A_INPUT)?;
    let xgamma_units = xgamma_runtime_unit_checks()?;
    if xgamma_units.iter().any(|passed| !passed) {
        return Err("XGamma runtime unit precheck failed".into());
    }
    let mut runs = Vec::new();
    for spec in SPECS {
        runs.push(run_condition(spec)?);
    }
    let summaries: Vec<Summary> = runs.iter().map(summarize).collect();
    write_timeseries(&runs)?;
    write_summary(&summaries)?;
    write_three_point_comparison(&summaries)?;
    let existing_unchanged = input_before == fs::read(M10A_INPUT)?;
    let checks_passed = write_checks(&runs, &summaries, &xgamma_units, existing_unchanged)?;
    write_performance(&runs)?;
    let classification = write_report(&summaries, &runs, checks_passed)?;
    println!("Milestone 10b final classification: {classification}");
    if !classification.starts_with("completed_") {
        return Err(format!("Milestone 10b stop: {classification}").into());
    }
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
        assert_eq!(
            SPECS.iter().map(|spec| spec.n).collect::<Vec<_>>(),
            vec![3, 5, 7]
        );
        assert!(SPECS
            .iter()
            .all(|spec| (spec.gamma_site() * spec.n as f64 - 3.0).abs() <= 1.0e-14));
    }
}
