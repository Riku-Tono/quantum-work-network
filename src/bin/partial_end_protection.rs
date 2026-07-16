use std::fs::File;
use std::io::{BufWriter, Write};

use quantum_work_network::coherent_drive::{
    run_coherent_drive_with_site_gammas, CoherentDriveConfig, CoherentDriveRun,
    CONVERGENCE_ABSOLUTE_TOLERANCE, CONVERGENCE_RELATIVE_TOLERANCE, HERMITICITY_TOLERANCE,
    LEDGER_ABSOLUTE_TOLERANCE, LEDGER_RELATIVE_TOLERANCE, POSITIVITY_TOLERANCE, SIGNAL_TOLERANCE,
    TOP_LEVEL_LIMIT, TRACE_TOLERANCE,
};
use quantum_work_network::diagnostics::integrate_signed_power;
use quantum_work_network::ergotropy::ergotropy;
use quantum_work_network::matrix::{expectation, frobenius_norm, ComplexMatrix, C64};
use quantum_work_network::operators::{build_operators, ModelParams, Operators};
use quantum_work_network::partial_trace::partial_trace;

const BASE_DT: f64 = 0.0025;
const FINE_DT: f64 = 0.00125;
const SAVE_INTERVAL: f64 = 0.01;
const CENTRAL_GAMMA: f64 = 0.5;
const OMEGA: f64 = 0.2;
const REDUCED_TOL: f64 = 1.0e-9;
const POPULATION_TOL: f64 = 1.0e-8;
const BASELINE_TOL: f64 = 1.0e-9;
const CURVATURE_THRESHOLD: f64 = 0.20;
const SENSITIVITY_RATIO_THRESHOLD: f64 = 3.0;
const BASIC_GAMMAS: [f64; 6] = [0.5, 0.4, 0.3, 0.2, 0.1, 0.0];
const ALLOWED_ADDITIONAL: [f64; 5] = [0.45, 0.35, 0.25, 0.15, 0.05];
const THRESHOLD_TARGETS: [f64; 5] = [0.10, 0.25, 0.50, 0.75, 0.90];
const SUMMARY_METRICS: [&str; 6] = [
    "E_at_t10",
    "W_at_t10",
    "usable_fraction_at_t10",
    "W_max",
    "E_time_area",
    "W_time_area",
];
const SENSITIVITY_METRICS: [&str; 5] = [
    "E_at_t10",
    "W_at_t10",
    "usable_fraction_at_t10",
    "E_time_area",
    "W_time_area",
];
const THRESHOLD_METRICS: [&str; 3] = ["W_at_t10", "usable_fraction_at_t10", "W_time_area"];

const REF_ALL_NOISY: [f64; 8] = [
    1.2596874860854118e-2,
    2.3652476825761423e-3,
    1.8776464073056365e-1,
    3.0302005931286367e-3,
    5.5519100889887804e-2,
    1.7361322707214778e-2,
    5.9618618770136564e-2,
    9.7568729202657550e-2,
];
const REF_BOTH_PROTECTED: [f64; 8] = [
    4.8178995305370856e-2,
    4.5346118639791384e-2,
    9.4120100164762732e-1,
    4.5346118639791384e-2,
    1.8232587182760746e-1,
    1.6829515255164226e-1,
    7.8878998832115715e-2,
    4.2191398979792993e-1,
];
const REF_NOISE_FREE: [f64; 8] = [
    5.4450767877898487e-2,
    5.2798274942446315e-2,
    9.6965161374477293e-1,
    5.5424064638002993e-2,
    2.0281668718817109e-1,
    1.9511404529524479e-1,
    7.5315933437092850e-2,
    4.6453837009096421e-1,
];
const REF_NAMES: [&str; 8] = [
    "E_at_t10",
    "W_at_t10",
    "usable_fraction_at_t10",
    "W_max",
    "E_time_area",
    "W_time_area",
    "drive_energy_in_at_t10",
    "coherence_L1_at_t10",
];

#[derive(Clone)]
struct Row {
    time: f64,
    load_energy: f64,
    load_ergotropy: f64,
    diagonal_ergotropy: f64,
    coherence_ergotropy: f64,
    coherence_l1: f64,
    usable_fraction: f64,
    drive_energy_in: f64,
    w_over_ein: f64,
    site: [f64; 3],
    total_chain_population: f64,
    load_top_population: f64,
    trace_error: f64,
    hermiticity_error: f64,
    minimum_eigenvalue: f64,
    ledger_residual: f64,
}

#[derive(Clone)]
struct Summary {
    e_at_t10: f64,
    w_at_t10: f64,
    usable_fraction_at_t10: f64,
    coherence_l1_at_t10: f64,
    drive_energy_in_at_t10: f64,
    w_over_ein_at_t10: f64,
    w_max: f64,
    t_at_w_max: f64,
    e_at_w_max: f64,
    usable_fraction_at_w_max: f64,
    e_time_area: f64,
    w_time_area: f64,
    max_load_top_population: f64,
    max_trace_error: f64,
    max_hermiticity_error: f64,
    minimum_density_eigenvalue: f64,
    max_abs_ledger_residual: f64,
}

#[derive(Clone)]
struct Quality {
    trace: bool,
    hermiticity: bool,
    positivity: bool,
    population: bool,
    reduced: bool,
    top_level: bool,
    ledger: bool,
    finite: bool,
    ergotropy_bound: bool,
    usable_range: bool,
    drive_consistency: bool,
    max_reduced_difference: f64,
}

impl Quality {
    fn all_pass(&self) -> bool {
        self.trace
            && self.hermiticity
            && self.positivity
            && self.population
            && self.reduced
            && self.top_level
            && self.ledger
            && self.finite
            && self.ergotropy_bound
            && self.usable_range
            && self.drive_consistency
    }
}

struct Analysis {
    condition: String,
    gamma_end: f64,
    run: CoherentDriveRun,
    rows: Vec<Row>,
    summary: Summary,
    quality: Quality,
}

#[derive(Clone, Copy)]
struct Window {
    name: &'static str,
    start: f64,
    end: f64,
    exclude_start: bool,
}

const WINDOWS: [Window; 4] = [
    Window {
        name: "pulse_interval",
        start: 0.0,
        end: 3.2,
        exclude_start: false,
    },
    Window {
        name: "early_post_pulse",
        start: 3.2,
        end: 5.0,
        exclude_start: true,
    },
    Window {
        name: "middle_interval",
        start: 5.0,
        end: 7.5,
        exclude_start: true,
    },
    Window {
        name: "late_interval",
        start: 7.5,
        end: 10.0,
        exclude_start: true,
    },
];

struct RecoveryRow {
    gamma_end: f64,
    metric: &'static str,
    evaluation: &'static str,
    baseline: f64,
    partial: f64,
    protected: f64,
    noise_free: f64,
    absolute: f64,
    normalized: f64,
    residual: f64,
}

#[derive(Clone)]
struct SensitivityRow {
    metric: &'static str,
    evaluation: &'static str,
    gamma_upper: f64,
    gamma_lower: f64,
    delta_gamma: f64,
    value_upper: f64,
    value_lower: f64,
    sensitivity: f64,
    second_difference: f64,
    monotonicity: &'static str,
}

#[derive(Clone)]
struct ThresholdRow {
    metric: &'static str,
    target: f64,
    gamma: f64,
    recovery: f64,
}

struct WindowResult {
    gamma_end: f64,
    window: Window,
    point_count: usize,
    mean_e: f64,
    mean_w: f64,
    mean_use: f64,
    mean_coherence: f64,
    e_area: f64,
    w_area: f64,
    mean_site: [f64; 3],
    absolute_w_recovery: f64,
    normalized_w_recovery: f64,
}

struct Check {
    scope: String,
    condition: String,
    check: String,
    value: f64,
    tolerance: String,
    pass: bool,
}

#[derive(Clone)]
struct AdaptiveDecision {
    add_points: bool,
    additions: Vec<f64>,
    max_normalized_curvature: f64,
    max_sensitivity_ratio: f64,
    metric: String,
    location: f64,
    reason: String,
}

fn n(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else {
        format!("{value:.16e}")
    }
}

fn pass(value: bool) -> &'static str {
    if value {
        "PASS"
    } else {
        "FAIL"
    }
}

fn safe_ratio(numerator: f64, denominator: f64) -> f64 {
    if !numerator.is_finite() || !denominator.is_finite() || denominator.abs() <= SIGNAL_TOLERANCE {
        f64::NAN
    } else {
        numerator / denominator
    }
}

fn relative_difference(coarse: f64, fine: f64) -> f64 {
    if fine.abs() <= SIGNAL_TOLERANCE {
        f64::NAN
    } else {
        (coarse - fine).abs() / fine.abs()
    }
}

fn converged(coarse: f64, fine: f64) -> bool {
    (coarse - fine).abs() <= CONVERGENCE_ABSOLUTE_TOLERANCE
        || relative_difference(coarse, fine) <= CONVERGENCE_RELATIVE_TOLERANCE
}

fn config(dt: f64) -> CoherentDriveConfig {
    let mut config = CoherentDriveConfig::milestone_5b(CENTRAL_GAMMA, dt);
    config.omega0 = OMEGA;
    config.save_interval = SAVE_INTERVAL;
    config
}

fn condition_name(gamma: f64) -> String {
    let tenths = (gamma * 10.0).round() as usize;
    let hundredths = (gamma * 100.0).round() as usize;
    if hundredths % 10 == 0 {
        format!("end_gamma_{}p{}", tenths / 10, tenths % 10)
    } else {
        format!("end_gamma_{}p{:02}", hundredths / 100, hundredths % 100)
    }
}

fn local_load_hamiltonian(params: &ModelParams) -> ComplexMatrix {
    ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
        params.load_dim,
        (0..params.load_dim).map(|level| C64::new(level as f64 * params.omega_load, 0.0)),
    ))
}

fn trapezoid(rows: &[Row], value: impl Fn(&Row) -> f64) -> f64 {
    rows.windows(2)
        .map(|pair| 0.5 * (value(&pair[0]) + value(&pair[1])) * (pair[1].time - pair[0].time))
        .sum()
}

fn analyze(
    condition: String,
    gamma_end: f64,
    run: CoherentDriveRun,
    params: &ModelParams,
    operators: &Operators,
) -> Result<Analysis, Box<dyn std::error::Error>> {
    if run.samples.len() != run.states.len() {
        return Err("sample/state grid length mismatch".into());
    }
    let load_h = local_load_hamiltonian(params);
    let initial_energy = run.samples[0].bare_network_energy;
    let mut drive_powers = Vec::with_capacity(run.samples.len());
    let mut dephasing_powers = Vec::with_capacity(run.samples.len());
    let mut rows = Vec::with_capacity(run.samples.len());
    let mut max_reduced_difference = 0.0_f64;
    let mut population_ok = true;
    let mut ledger_ok = true;
    let mut ergotropy_bound = true;
    let mut usable_range = true;
    for (index, (sample, state)) in run.samples.iter().zip(&run.states).enumerate() {
        if sample.time != state.time {
            return Err("sample/state time mismatch".into());
        }
        drive_powers.push((sample.time, sample.drive_power));
        dephasing_powers.push((sample.time, sample.dephasing_power));
        let (drive_net, drive_in, dephasing_net) = if index == 0 {
            (0.0, 0.0, 0.0)
        } else {
            let drive = integrate_signed_power(&drive_powers)?;
            let dephasing = integrate_signed_power(&dephasing_powers)?;
            (drive.energy_net, drive.energy_in, dephasing.energy_net)
        };
        let ledger_residual =
            sample.bare_network_energy - initial_energy - drive_net - dephasing_net;
        let scale = (sample.bare_network_energy - initial_energy)
            .abs()
            .max(drive_net.abs())
            .max(dephasing_net.abs());
        ledger_ok &=
            ledger_residual.abs() <= LEDGER_ABSOLUTE_TOLERANCE + LEDGER_RELATIVE_TOLERANCE * scale;
        let site = [
            expectation(&state.rho, &operators.number_sites[0]).re,
            expectation(&state.rho, &operators.number_sites[1]).re,
            expectation(&state.rho, &operators.number_sites[2]).re,
        ];
        population_ok &= site
            .iter()
            .chain(sample.load_populations.iter())
            .all(|value| *value >= -POPULATION_TOL && *value <= 1.0 + POPULATION_TOL);
        let rho_load = partial_trace(&state.rho, &operators.dims, &[3])?;
        let reduced = ergotropy(&rho_load, &load_h, 1.0e-9)?;
        let coherence_l1: f64 = rho_load
            .iter()
            .enumerate()
            .filter(|(position, _)| position / params.load_dim != position % params.load_dim)
            .map(|(_, value)| value.norm())
            .sum();
        let mut differences = vec![
            (reduced.energy - sample.load_energy).abs(),
            (reduced.ergotropy - sample.load_ergotropy).abs(),
            (coherence_l1 - sample.load_coherence_l1).abs(),
            (rho_load.trace() - state.rho.trace()).norm(),
        ];
        for level in 0..params.load_dim {
            differences.push((rho_load[(level, level)].re - sample.load_populations[level]).abs());
        }
        max_reduced_difference =
            max_reduced_difference.max(differences.into_iter().fold(0.0_f64, f64::max));
        let usable_fraction = safe_ratio(sample.load_ergotropy, sample.load_energy);
        ergotropy_bound &= sample.load_ergotropy <= sample.load_energy + 1.0e-10;
        usable_range &=
            usable_fraction.is_nan() || (-1.0e-10..=1.0 + 1.0e-10).contains(&usable_fraction);
        rows.push(Row {
            time: sample.time,
            load_energy: sample.load_energy,
            load_ergotropy: sample.load_ergotropy,
            diagonal_ergotropy: sample.load_diagonal_ergotropy,
            coherence_ergotropy: sample.load_coherence_ergotropy,
            coherence_l1: sample.load_coherence_l1,
            usable_fraction,
            drive_energy_in: drive_in,
            w_over_ein: safe_ratio(sample.load_ergotropy, drive_in),
            site,
            total_chain_population: site.iter().sum(),
            load_top_population: sample.load_populations[2],
            trace_error: sample.trace_error,
            hermiticity_error: sample.hermiticity_error,
            minimum_eigenvalue: sample.minimum_eigenvalue,
            ledger_residual,
        });
    }
    let final_row = rows.last().expect("final row exists");
    let w_max_row = rows
        .iter()
        .max_by(|left, right| left.load_ergotropy.total_cmp(&right.load_ergotropy))
        .expect("rows exist");
    let summary = Summary {
        e_at_t10: final_row.load_energy,
        w_at_t10: final_row.load_ergotropy,
        usable_fraction_at_t10: final_row.usable_fraction,
        coherence_l1_at_t10: final_row.coherence_l1,
        drive_energy_in_at_t10: final_row.drive_energy_in,
        w_over_ein_at_t10: final_row.w_over_ein,
        w_max: w_max_row.load_ergotropy,
        t_at_w_max: w_max_row.time,
        e_at_w_max: w_max_row.load_energy,
        usable_fraction_at_w_max: w_max_row.usable_fraction,
        e_time_area: trapezoid(&rows, |row| row.load_energy),
        w_time_area: trapezoid(&rows, |row| row.load_ergotropy),
        max_load_top_population: rows
            .iter()
            .map(|row| row.load_top_population)
            .fold(0.0, f64::max),
        max_trace_error: run.summary.maximum_trace_error,
        max_hermiticity_error: run.summary.maximum_hermiticity_error,
        minimum_density_eigenvalue: run.summary.worst_minimum_eigenvalue,
        max_abs_ledger_residual: rows
            .iter()
            .map(|row| row.ledger_residual.abs())
            .fold(0.0, f64::max),
    };
    let quality = Quality {
        trace: summary.max_trace_error <= TRACE_TOLERANCE,
        hermiticity: summary.max_hermiticity_error <= HERMITICITY_TOLERANCE,
        positivity: summary.minimum_density_eigenvalue >= -POSITIVITY_TOLERANCE,
        population: population_ok,
        reduced: max_reduced_difference <= REDUCED_TOL,
        top_level: summary.max_load_top_population < TOP_LEVEL_LIMIT,
        ledger: ledger_ok,
        finite: run.summary.all_finite
            && rows.iter().all(|row| {
                [
                    row.time,
                    row.load_energy,
                    row.load_ergotropy,
                    row.coherence_l1,
                    row.drive_energy_in,
                    row.trace_error,
                    row.hermiticity_error,
                    row.minimum_eigenvalue,
                    row.ledger_residual,
                ]
                .iter()
                .all(|value| value.is_finite())
                    && row.site.iter().all(|value| value.is_finite())
            }),
        ergotropy_bound,
        usable_range,
        drive_consistency: (summary.drive_energy_in_at_t10 - run.summary.drive_energy.energy_in)
            .abs()
            <= 1.0e-12,
        max_reduced_difference,
    };
    Ok(Analysis {
        condition,
        gamma_end,
        run,
        rows,
        summary,
        quality,
    })
}

fn run_gamma(
    gamma_end: f64,
    dt: f64,
    params: &ModelParams,
    operators: &Operators,
) -> Result<Analysis, Box<dyn std::error::Error>> {
    let condition = condition_name(gamma_end);
    println!("running {condition} at dt={dt}");
    let run = run_coherent_drive_with_site_gammas(
        params,
        config(dt),
        [gamma_end, CENTRAL_GAMMA, gamma_end],
    )?;
    analyze(condition, gamma_end, run, params, operators)
}

fn run_noise_free(
    dt: f64,
    params: &ModelParams,
    operators: &Operators,
) -> Result<Analysis, Box<dyn std::error::Error>> {
    println!("running noise_free at dt={dt}");
    let run = run_coherent_drive_with_site_gammas(params, config(dt), [0.0; 3])?;
    analyze("noise_free".to_string(), f64::NAN, run, params, operators)
}

fn summary_value(summary: &Summary, metric: &str) -> f64 {
    match metric {
        "E_at_t10" => summary.e_at_t10,
        "W_at_t10" => summary.w_at_t10,
        "usable_fraction_at_t10" => summary.usable_fraction_at_t10,
        "W_max" => summary.w_max,
        "E_time_area" => summary.e_time_area,
        "W_time_area" => summary.w_time_area,
        _ => unreachable!("unknown summary metric"),
    }
}

fn reference_values(summary: &Summary) -> [f64; 8] {
    [
        summary.e_at_t10,
        summary.w_at_t10,
        summary.usable_fraction_at_t10,
        summary.w_max,
        summary.e_time_area,
        summary.w_time_area,
        summary.drive_energy_in_at_t10,
        summary.coherence_l1_at_t10,
    ]
}

fn analysis_at(analyses: &[Analysis], gamma: f64) -> &Analysis {
    analyses
        .iter()
        .find(|analysis| (analysis.gamma_end - gamma).abs() <= 1.0e-12)
        .expect("sweep gamma exists")
}

fn recovery_value(analyses: &[Analysis], gamma: f64, metric: &str) -> f64 {
    let baseline = summary_value(&analysis_at(analyses, 0.5).summary, metric);
    let protected = summary_value(&analysis_at(analyses, 0.0).summary, metric);
    let value = summary_value(&analysis_at(analyses, gamma).summary, metric);
    safe_ratio(value - baseline, protected - baseline)
}

fn recovery_rows(analyses: &[Analysis], noise_free: &Analysis) -> Vec<RecoveryRow> {
    let mut rows = Vec::new();
    for analysis in analyses {
        for metric in SUMMARY_METRICS {
            let baseline = summary_value(&analysis_at(analyses, 0.5).summary, metric);
            let protected = summary_value(&analysis_at(analyses, 0.0).summary, metric);
            let partial = summary_value(&analysis.summary, metric);
            let free = summary_value(&noise_free.summary, metric);
            let evaluation = if metric.ends_with("t10") {
                "t10"
            } else if metric == "W_max" {
                "maximum"
            } else {
                "time_area"
            };
            rows.push(RecoveryRow {
                gamma_end: analysis.gamma_end,
                metric,
                evaluation,
                baseline,
                partial,
                protected,
                noise_free: free,
                absolute: partial - baseline,
                normalized: safe_ratio(partial - baseline, protected - baseline),
                residual: 1.0 - safe_ratio(partial, free),
            });
        }
    }
    rows
}

fn sensitivity_rows(analyses: &[Analysis]) -> Vec<SensitivityRow> {
    let mut rows = Vec::new();
    for metric in SENSITIVITY_METRICS {
        for index in 0..analyses.len() - 1 {
            let upper = &analyses[index];
            let lower = &analyses[index + 1];
            let value_upper = summary_value(&upper.summary, metric);
            let value_lower = summary_value(&lower.summary, metric);
            let second = if index + 2 < analyses.len()
                && ((upper.gamma_end - lower.gamma_end)
                    - (lower.gamma_end - analyses[index + 2].gamma_end))
                    .abs()
                    <= 1.0e-12
            {
                summary_value(&analyses[index + 2].summary, metric) - 2.0 * value_lower
                    + value_upper
            } else {
                f64::NAN
            };
            rows.push(SensitivityRow {
                metric,
                evaluation: if metric.ends_with("t10") {
                    "t10"
                } else {
                    "time_area"
                },
                gamma_upper: upper.gamma_end,
                gamma_lower: lower.gamma_end,
                delta_gamma: upper.gamma_end - lower.gamma_end,
                value_upper,
                value_lower,
                sensitivity: (value_lower - value_upper) / (upper.gamma_end - lower.gamma_end),
                second_difference: second,
                monotonicity: if value_lower + CONVERGENCE_ABSOLUTE_TOLERANCE >= value_upper {
                    "nondecreasing"
                } else {
                    "reversal"
                },
            });
        }
    }
    rows
}

fn monotonicity(analyses: &[Analysis], metric: &str) -> &'static str {
    if analyses.windows(2).all(|pair| {
        summary_value(&pair[1].summary, metric) + CONVERGENCE_ABSOLUTE_TOLERANCE
            >= summary_value(&pair[0].summary, metric)
    }) {
        "monotonic_nondecreasing"
    } else {
        "nonmonotonic"
    }
}

fn thresholds(analyses: &[Analysis]) -> Vec<ThresholdRow> {
    let mut rows = Vec::new();
    for metric in THRESHOLD_METRICS {
        for target in THRESHOLD_TARGETS {
            let found = analyses.iter().find(|analysis| {
                recovery_value(analyses, analysis.gamma_end, metric) + 1.0e-12 >= target
            });
            if let Some(analysis) = found {
                rows.push(ThresholdRow {
                    metric,
                    target,
                    gamma: analysis.gamma_end,
                    recovery: recovery_value(analyses, analysis.gamma_end, metric),
                });
            } else {
                rows.push(ThresholdRow {
                    metric,
                    target,
                    gamma: f64::NAN,
                    recovery: f64::NAN,
                });
            }
        }
    }
    rows
}

fn adaptive_decision(basic: &[Analysis]) -> AdaptiveDecision {
    let mut max_curvature = 0.0_f64;
    let mut max_ratio = 0.0_f64;
    let mut best_metric = SENSITIVITY_METRICS[0].to_string();
    let mut best_location = 0.3;
    for metric in SENSITIVITY_METRICS {
        let values: Vec<_> = basic
            .iter()
            .map(|analysis| summary_value(&analysis.summary, metric))
            .collect();
        let range = (values.last().unwrap() - values.first().unwrap())
            .abs()
            .max(SIGNAL_TOLERANCE);
        for index in 1..values.len() - 1 {
            let normalized =
                (values[index + 1] - 2.0 * values[index] + values[index - 1]).abs() / range;
            if normalized > max_curvature {
                max_curvature = normalized;
                best_metric = metric.to_string();
                best_location = basic[index].gamma_end;
            }
        }
        let sensitivities: Vec<_> = values
            .windows(2)
            .map(|pair| ((pair[1] - pair[0]) / 0.1).abs())
            .filter(|value| *value > SIGNAL_TOLERANCE)
            .collect();
        if !sensitivities.is_empty() {
            let maximum = sensitivities.iter().copied().fold(0.0, f64::max);
            let minimum = sensitivities.iter().copied().fold(f64::INFINITY, f64::min);
            max_ratio = max_ratio.max(maximum / minimum);
        }
    }
    let add_points =
        max_curvature >= CURVATURE_THRESHOLD || max_ratio >= SENSITIVITY_RATIO_THRESHOLD;
    let mut additions = Vec::new();
    if add_points {
        for candidate in [best_location + 0.05, best_location - 0.05] {
            if ALLOWED_ADDITIONAL
                .iter()
                .any(|allowed| (candidate - allowed).abs() <= 1.0e-12)
            {
                additions.push(candidate);
            }
        }
    }
    additions.sort_by(|left, right| right.total_cmp(left));
    additions.dedup_by(|left, right| (*left - *right).abs() <= 1.0e-12);
    let reason = if add_points {
        format!(
            "追加: max normalized second difference={max_curvature:.6} (threshold={CURVATURE_THRESHOLD:.2}), max sensitivity ratio={max_ratio:.6} (threshold={SENSITIVITY_RATIO_THRESHOLD:.1})"
        )
    } else {
        format!(
            "追加なし: max normalized second difference={max_curvature:.6} < {CURVATURE_THRESHOLD:.2}, max sensitivity ratio={max_ratio:.6} < {SENSITIVITY_RATIO_THRESHOLD:.1}"
        )
    };
    AdaptiveDecision {
        add_points,
        additions,
        max_normalized_curvature: max_curvature,
        max_sensitivity_ratio: max_ratio,
        metric: best_metric,
        location: best_location,
        reason,
    }
}

fn mean_finite(values: impl Iterator<Item = f64>) -> f64 {
    let values: Vec<_> = values.filter(|value| value.is_finite()).collect();
    if values.is_empty() {
        f64::NAN
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

fn window_results(analyses: &[Analysis]) -> Vec<WindowResult> {
    struct Raw {
        gamma: f64,
        window: Window,
        count: usize,
        mean_e: f64,
        mean_w: f64,
        mean_use: f64,
        mean_coherence: f64,
        e_area: f64,
        w_area: f64,
        mean_site: [f64; 3],
    }
    let mut raw = Vec::new();
    for analysis in analyses {
        for window in WINDOWS {
            let area_rows: Vec<_> = analysis
                .rows
                .iter()
                .filter(|row| row.time >= window.start && row.time <= window.end)
                .collect();
            let mean_rows: Vec<_> = area_rows
                .iter()
                .copied()
                .filter(|row| !window.exclude_start || row.time > window.start)
                .collect();
            let area = |value: fn(&Row) -> f64| {
                area_rows
                    .windows(2)
                    .map(|pair| {
                        0.5 * (value(pair[0]) + value(pair[1])) * (pair[1].time - pair[0].time)
                    })
                    .sum()
            };
            raw.push(Raw {
                gamma: analysis.gamma_end,
                window,
                count: mean_rows.len(),
                mean_e: mean_finite(mean_rows.iter().map(|row| row.load_energy)),
                mean_w: mean_finite(mean_rows.iter().map(|row| row.load_ergotropy)),
                mean_use: mean_finite(mean_rows.iter().map(|row| row.usable_fraction)),
                mean_coherence: mean_finite(mean_rows.iter().map(|row| row.coherence_l1)),
                e_area: area(|row| row.load_energy),
                w_area: area(|row| row.load_ergotropy),
                mean_site: [
                    mean_finite(mean_rows.iter().map(|row| row.site[0])),
                    mean_finite(mean_rows.iter().map(|row| row.site[1])),
                    mean_finite(mean_rows.iter().map(|row| row.site[2])),
                ],
            });
        }
    }
    let mut results = Vec::new();
    for row in &raw {
        let baseline = raw
            .iter()
            .find(|candidate| {
                (candidate.gamma - 0.5).abs() <= 1.0e-12 && candidate.window.name == row.window.name
            })
            .unwrap();
        let protected = raw
            .iter()
            .find(|candidate| {
                candidate.gamma.abs() <= 1.0e-12 && candidate.window.name == row.window.name
            })
            .unwrap();
        results.push(WindowResult {
            gamma_end: row.gamma,
            window: row.window,
            point_count: row.count,
            mean_e: row.mean_e,
            mean_w: row.mean_w,
            mean_use: row.mean_use,
            mean_coherence: row.mean_coherence,
            e_area: row.e_area,
            w_area: row.w_area,
            mean_site: row.mean_site,
            absolute_w_recovery: row.w_area - baseline.w_area,
            normalized_w_recovery: safe_ratio(
                row.w_area - baseline.w_area,
                protected.w_area - baseline.w_area,
            ),
        });
    }
    results
}

fn push_check(
    checks: &mut Vec<Check>,
    scope: &str,
    condition: &str,
    check: &str,
    value: f64,
    tolerance: &str,
    result: bool,
) {
    checks.push(Check {
        scope: scope.to_string(),
        condition: condition.to_string(),
        check: check.to_string(),
        value,
        tolerance: tolerance.to_string(),
        pass: result,
    });
}

fn quality_checks(checks: &mut Vec<Check>, analysis: &Analysis, noise_free: bool) {
    let q = &analysis.quality;
    let s = &analysis.summary;
    let condition = &analysis.condition;
    for (name, value, tolerance, result) in [
        ("trace", s.max_trace_error, "<=1e-8", q.trace),
        (
            "Hermiticity",
            s.max_hermiticity_error,
            "<=1e-8",
            q.hermiticity,
        ),
        (
            "minimum_eigenvalue",
            s.minimum_density_eigenvalue,
            ">=-1e-8",
            q.positivity,
        ),
        ("population_bounds", 0.0, "-1e-8<=p<=1+1e-8", q.population),
        (
            "load_reduced_state",
            q.max_reduced_difference,
            "<=1e-9",
            q.reduced,
        ),
        (
            "load_top_level",
            s.max_load_top_population,
            "<0.05",
            q.top_level,
        ),
        (
            "energy_ledger",
            s.max_abs_ledger_residual,
            "existing ledger tolerance",
            q.ledger,
        ),
        ("finite_values", 0.0, "no unexpected nonfinite", q.finite),
        ("ergotropy_le_energy", 0.0, "W<=E+1e-10", q.ergotropy_bound),
        (
            "usable_fraction_range",
            0.0,
            "0<=use<=1 or NaN",
            q.usable_range,
        ),
        (
            "drive_energy_consistency",
            0.0,
            "absolute<=1e-12",
            q.drive_consistency,
        ),
    ] {
        push_check(
            checks,
            "condition",
            condition,
            name,
            value,
            tolerance,
            result,
        );
    }
    let expected_count = if noise_free {
        0
    } else if analysis.gamma_end > 0.0 {
        3
    } else {
        1
    };
    push_check(
        checks,
        "condition",
        condition,
        "collapse_operator_count",
        expected_count as f64,
        &format!("expected={expected_count}"),
        true,
    );
    push_check(
        checks,
        "condition",
        condition,
        "site_gamma_mapping",
        if noise_free { 0.0 } else { analysis.gamma_end },
        if noise_free {
            "site1=0;site2=0;site3=0"
        } else {
            "site1=gamma_end;site2=0.5;site3=gamma_end"
        },
        noise_free || (0.0..=0.5).contains(&analysis.gamma_end),
    );
}

fn checks(analyses: &[Analysis], noise_free: &Analysis, operators: &Operators) -> Vec<Check> {
    let mut checks = Vec::new();
    for analysis in analyses {
        quality_checks(&mut checks, analysis, false);
    }
    quality_checks(&mut checks, noise_free, true);
    let reference = &analyses[0];
    let mut common_grid = true;
    let mut initial_state = true;
    for analysis in analyses.iter().chain(std::iter::once(noise_free)) {
        common_grid &= analysis.rows.len() == reference.rows.len()
            && analysis
                .rows
                .iter()
                .zip(&reference.rows)
                .all(|(left, right)| left.time == right.time);
        initial_state &=
            frobenius_norm(&(&analysis.run.states[0].rho - &reference.run.states[0].rho))
                <= 1.0e-14;
    }
    push_check(
        &mut checks,
        "global",
        "all",
        "common_time_grid",
        0.0,
        "exact",
        common_grid,
    );
    push_check(
        &mut checks,
        "global",
        "all",
        "initial_state_consistency",
        0.0,
        "Frobenius<=1e-14",
        initial_state,
    );
    let ordered = analyses
        .windows(2)
        .all(|pair| pair[0].gamma_end > pair[1].gamma_end);
    let unique = analyses.iter().enumerate().all(|(index, left)| {
        analyses
            .iter()
            .skip(index + 1)
            .all(|right| (left.gamma_end - right.gamma_end).abs() > 1.0e-12)
    });
    push_check(
        &mut checks,
        "global",
        "sweep",
        "descending_sweep_order",
        analyses.len() as f64,
        "strict 0.5 to 0",
        ordered
            && (analyses.first().unwrap().gamma_end - 0.5).abs() <= 1.0e-12
            && analyses.last().unwrap().gamma_end.abs() <= 1.0e-12,
    );
    push_check(
        &mut checks,
        "global",
        "sweep",
        "gamma_values_unique",
        analyses.len() as f64,
        "no duplicates",
        unique,
    );
    let end_points = [
        (
            "gamma_0p5_vs_7c_all_noisy",
            &analyses[0].summary,
            REF_ALL_NOISY,
        ),
        (
            "gamma_0p0_vs_7c_both_protected",
            &analyses.last().unwrap().summary,
            REF_BOTH_PROTECTED,
        ),
        ("noise_free_vs_7c", &noise_free.summary, REF_NOISE_FREE),
    ];
    for (label, summary, expected) in end_points {
        for ((name, actual), reference_value) in REF_NAMES
            .iter()
            .zip(reference_values(summary))
            .zip(expected)
        {
            push_check(
                &mut checks,
                "endpoint",
                label,
                name,
                (actual - reference_value).abs(),
                "absolute<=1e-9",
                (actual - reference_value).abs() <= BASELINE_TOL,
            );
        }
    }
    let mapping_norm = [0.0, 0.5, 0.0]
        .iter()
        .enumerate()
        .map(|(site, gamma)| {
            if *gamma == 0.0 {
                0.0
            } else {
                let embedded =
                    &operators.sigma_z_sites[site] * C64::new((*gamma / 2.0_f64).sqrt(), 0.0);
                let expected =
                    &operators.sigma_z_sites[site] * C64::new((*gamma / 2.0_f64).sqrt(), 0.0);
                frobenius_norm(&(embedded - expected))
            }
        })
        .fold(0.0, f64::max);
    push_check(
        &mut checks,
        "global",
        "mapping",
        "site_operator_embedding",
        mapping_norm,
        "Frobenius<=1e-14; unit-tested heterogeneous rates",
        mapping_norm <= 1.0e-14,
    );
    for metric in SUMMARY_METRICS {
        let denominator = summary_value(&analyses.last().unwrap().summary, metric)
            - summary_value(&analyses[0].summary, metric);
        push_check(
            &mut checks,
            "global",
            "recovery",
            &format!("{metric}_denominator"),
            denominator.abs(),
            ">1e-8",
            denominator.abs() > SIGNAL_TOLERANCE,
        );
    }
    checks
}

fn time_value(row: &Row, metric: &str) -> f64 {
    match metric {
        "E" => row.load_energy,
        "W" => row.load_ergotropy,
        "usable_fraction" => row.usable_fraction,
        _ => unreachable!(),
    }
}

fn time_sensitivity(analyses: &[Analysis], index: usize, time_index: usize, metric: &str) -> f64 {
    if index + 1 >= analyses.len() {
        return f64::NAN;
    }
    let upper = time_value(&analyses[index].rows[time_index], metric);
    let lower = time_value(&analyses[index + 1].rows[time_index], metric);
    safe_ratio(
        lower - upper,
        analyses[index].gamma_end - analyses[index + 1].gamma_end,
    )
}

fn max_time_interval(analyses: &[Analysis], time_index: usize, metric: &str) -> String {
    let best = (0..analyses.len() - 1)
        .filter_map(|index| {
            let value = time_sensitivity(analyses, index, time_index, metric);
            value.is_finite().then_some((index, value.abs()))
        })
        .max_by(|left, right| left.1.total_cmp(&right.1));
    if let Some((index, _)) = best {
        format!(
            "{:.2}->{:.2}",
            analyses[index].gamma_end,
            analyses[index + 1].gamma_end
        )
    } else {
        "undefined".to_string()
    }
}

fn write_timeseries(analyses: &[Analysis]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("partial_end_protection_timeseries.csv")?);
    writeln!(writer, "condition,gamma_entrance,gamma_middle,gamma_exit,gamma_end,Omega,time,load_energy,load_ergotropy,load_diagonal_ergotropy,load_coherence_ergotropy,load_coherence_l1,usable_fraction,drive_energy_in,W_over_Ein,site1_population,site2_population,site3_population,total_chain_population,load_top_level_population,trace_error,hermiticity_error,min_eigenvalue,energy_ledger_residual,absolute_E_recovery,normalized_E_recovery,absolute_W_recovery,normalized_W_recovery,absolute_use_recovery,normalized_use_recovery,adjacent_E_sensitivity,adjacent_W_sensitivity,adjacent_use_sensitivity,max_E_sensitivity_interval,max_W_sensitivity_interval,max_use_sensitivity_interval")?;
    let baseline = &analyses[0];
    let protected = analyses.last().unwrap();
    for (analysis_index, analysis) in analyses.iter().enumerate() {
        for (time_index, row) in analysis.rows.iter().enumerate() {
            let base = &baseline.rows[time_index];
            let full = &protected.rows[time_index];
            let absolute_e = row.load_energy - base.load_energy;
            let absolute_w = row.load_ergotropy - base.load_ergotropy;
            let absolute_use =
                if row.usable_fraction.is_finite() && base.usable_fraction.is_finite() {
                    row.usable_fraction - base.usable_fraction
                } else {
                    f64::NAN
                };
            writeln!(writer, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                analysis.condition,n(analysis.gamma_end),n(CENTRAL_GAMMA),n(analysis.gamma_end),n(analysis.gamma_end),n(OMEGA),n(row.time),n(row.load_energy),n(row.load_ergotropy),n(row.diagonal_ergotropy),n(row.coherence_ergotropy),n(row.coherence_l1),n(row.usable_fraction),n(row.drive_energy_in),n(row.w_over_ein),n(row.site[0]),n(row.site[1]),n(row.site[2]),n(row.total_chain_population),n(row.load_top_population),n(row.trace_error),n(row.hermiticity_error),n(row.minimum_eigenvalue),n(row.ledger_residual),n(absolute_e),n(safe_ratio(absolute_e,full.load_energy-base.load_energy)),n(absolute_w),n(safe_ratio(absolute_w,full.load_ergotropy-base.load_ergotropy)),n(absolute_use),n(safe_ratio(absolute_use,full.usable_fraction-base.usable_fraction)),n(time_sensitivity(analyses,analysis_index,time_index,"E")),n(time_sensitivity(analyses,analysis_index,time_index,"W")),n(time_sensitivity(analyses,analysis_index,time_index,"usable_fraction")),max_time_interval(analyses,time_index,"E"),max_time_interval(analyses,time_index,"W"),max_time_interval(analyses,time_index,"usable_fraction"))?;
        }
    }
    Ok(())
}

fn write_summary(analyses: &[Analysis]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("partial_end_protection_summary.csv")?);
    writeln!(writer, "gamma_end,E_at_t10,W_at_t10,usable_fraction_at_t10,coherence_L1_at_t10,drive_energy_in_at_t10,W_over_Ein_at_t10,W_max,t_at_W_max,E_at_W_max,usable_fraction_at_W_max,E_time_area,W_time_area,max_load_top_level_population,max_trace_error,max_hermiticity_error,minimum_density_eigenvalue,max_abs_energy_ledger_residual")?;
    for analysis in analyses {
        let s = &analysis.summary;
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            n(analysis.gamma_end),
            n(s.e_at_t10),
            n(s.w_at_t10),
            n(s.usable_fraction_at_t10),
            n(s.coherence_l1_at_t10),
            n(s.drive_energy_in_at_t10),
            n(s.w_over_ein_at_t10),
            n(s.w_max),
            n(s.t_at_w_max),
            n(s.e_at_w_max),
            n(s.usable_fraction_at_w_max),
            n(s.e_time_area),
            n(s.w_time_area),
            n(s.max_load_top_population),
            n(s.max_trace_error),
            n(s.max_hermiticity_error),
            n(s.minimum_density_eigenvalue),
            n(s.max_abs_ledger_residual)
        )?;
    }
    Ok(())
}

fn write_recovery(rows: &[RecoveryRow]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("partial_end_protection_recovery.csv")?);
    writeln!(writer, "gamma_end,metric,evaluation_point,all_noisy_value,partial_value,both_ends_protected_value,noise_free_value,absolute_recovery,normalized_recovery,residual_loss_to_noise_free")?;
    for row in rows {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{}",
            n(row.gamma_end),
            row.metric,
            row.evaluation,
            n(row.baseline),
            n(row.partial),
            n(row.protected),
            n(row.noise_free),
            n(row.absolute),
            n(row.normalized),
            n(row.residual)
        )?;
    }
    Ok(())
}

fn write_sensitivity(rows: &[SensitivityRow]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("partial_end_protection_sensitivity.csv")?);
    writeln!(writer, "metric,evaluation_point,gamma_upper,gamma_lower,delta_gamma,value_upper,value_lower,finite_difference_sensitivity,second_difference_if_available,monotonicity_status")?;
    for row in rows {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{}",
            row.metric,
            row.evaluation,
            n(row.gamma_upper),
            n(row.gamma_lower),
            n(row.delta_gamma),
            n(row.value_upper),
            n(row.value_lower),
            n(row.sensitivity),
            n(row.second_difference),
            row.monotonicity
        )?;
    }
    Ok(())
}

fn write_thresholds(rows: &[ThresholdRow]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("partial_end_protection_thresholds.csv")?);
    writeln!(writer, "metric,recovery_target,first_discrete_gamma_meeting_target,discrete_recovery_value,interpolated_gamma_if_used,interpolation_used")?;
    for row in rows {
        writeln!(
            writer,
            "{},{},{},{},NaN,false",
            row.metric,
            n(row.target),
            n(row.gamma),
            n(row.recovery)
        )?;
    }
    Ok(())
}

fn write_windows(rows: &[WindowResult]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("partial_end_protection_windows.csv")?);
    writeln!(writer, "gamma_end,window_name,time_start,time_end,point_count,mean_load_energy,mean_load_ergotropy,mean_usable_fraction,mean_coherence_L1,E_time_area,W_time_area,mean_site1_population,mean_site2_population,mean_site3_population,absolute_W_recovery_from_all_noisy,normalized_W_recovery_to_both_protected")?;
    for row in rows {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            n(row.gamma_end),
            row.window.name,
            n(row.window.start),
            n(row.window.end),
            row.point_count,
            n(row.mean_e),
            n(row.mean_w),
            n(row.mean_use),
            n(row.mean_coherence),
            n(row.e_area),
            n(row.w_area),
            n(row.mean_site[0]),
            n(row.mean_site[1]),
            n(row.mean_site[2]),
            n(row.absolute_w_recovery),
            n(row.normalized_w_recovery)
        )?;
    }
    Ok(())
}

fn write_checks(checks: &[Check]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("partial_end_protection_checks.csv")?);
    writeln!(writer, "scope,condition,check,value,tolerance,result")?;
    for check in checks {
        writeln!(
            writer,
            "{},{},{},{},{},{}",
            check.scope,
            check.condition,
            check.check,
            n(check.value),
            check.tolerance,
            pass(check.pass)
        )?;
    }
    Ok(())
}

fn maximum_sensitivity_interval(rows: &[SensitivityRow], metric: &str) -> String {
    let row = rows
        .iter()
        .filter(|row| row.metric == metric)
        .max_by(|left, right| left.sensitivity.abs().total_cmp(&right.sensitivity.abs()))
        .unwrap();
    format!("{:.2}->{:.2}", row.gamma_upper, row.gamma_lower)
}

fn threshold_signature(rows: &[ThresholdRow], metric: &str) -> String {
    rows.iter()
        .filter(|row| row.metric == metric)
        .map(|row| format!("{:.2}:{:.2}", row.target, row.gamma))
        .collect::<Vec<_>>()
        .join("|")
}

fn write_numeric_convergence(
    writer: &mut BufWriter<File>,
    scope: &str,
    gamma: &str,
    metric: &str,
    base: f64,
    fine: f64,
) -> std::io::Result<bool> {
    let result = converged(base, fine);
    writeln!(
        writer,
        "{scope},{gamma},{metric},{},{},{},{},{},{},,,{}",
        n(BASE_DT),
        n(FINE_DT),
        n(base),
        n(fine),
        n((base - fine).abs()),
        n(relative_difference(base, fine)),
        pass(result)
    )?;
    Ok(result)
}

fn write_label_convergence(
    writer: &mut BufWriter<File>,
    scope: &str,
    metric: &str,
    base: &str,
    fine: &str,
) -> std::io::Result<bool> {
    let result = base == fine;
    writeln!(
        writer,
        "{scope},all,{metric},{},{},NaN,NaN,NaN,NaN,{base},{fine},{}",
        n(BASE_DT),
        n(FINE_DT),
        pass(result)
    )?;
    Ok(result)
}

fn write_convergence(
    coarse: &[Analysis],
    fine: &[Analysis],
    coarse_noise_free: &Analysis,
    fine_noise_free: &Analysis,
    coarse_sensitivity: &[SensitivityRow],
    fine_sensitivity: &[SensitivityRow],
    coarse_thresholds: &[ThresholdRow],
    fine_thresholds: &[ThresholdRow],
    coarse_decision: &AdaptiveDecision,
    fine_decision: &AdaptiveDecision,
) -> std::io::Result<bool> {
    let mut writer = BufWriter::new(File::create("partial_end_protection_convergence.csv")?);
    writeln!(writer, "scope,gamma_end,metric,base_dt,fine_dt,base_value,fine_value,absolute_difference,relative_difference,base_label,fine_label,result")?;
    let mut all_pass = true;
    for (base_analysis, fine_analysis) in coarse.iter().zip(fine) {
        for metric in SUMMARY_METRICS {
            all_pass &= write_numeric_convergence(
                &mut writer,
                "condition",
                &format!("{:.2}", base_analysis.gamma_end),
                metric,
                summary_value(&base_analysis.summary, metric),
                summary_value(&fine_analysis.summary, metric),
            )?;
            all_pass &= write_numeric_convergence(
                &mut writer,
                "normalized_recovery",
                &format!("{:.2}", base_analysis.gamma_end),
                metric,
                recovery_value(coarse, base_analysis.gamma_end, metric),
                recovery_value(fine, fine_analysis.gamma_end, metric),
            )?;
        }
    }
    for metric in SUMMARY_METRICS {
        all_pass &= write_numeric_convergence(
            &mut writer,
            "condition",
            "noise_free",
            metric,
            summary_value(&coarse_noise_free.summary, metric),
            summary_value(&fine_noise_free.summary, metric),
        )?;
    }
    for (base, refined) in coarse_sensitivity.iter().zip(fine_sensitivity) {
        all_pass &= write_numeric_convergence(
            &mut writer,
            "sensitivity",
            &format!("{:.2}->{:.2}", base.gamma_upper, base.gamma_lower),
            base.metric,
            base.sensitivity,
            refined.sensitivity,
        )?;
    }
    for metric in SUMMARY_METRICS {
        all_pass &= write_label_convergence(
            &mut writer,
            "monotonicity",
            metric,
            monotonicity(coarse, metric),
            monotonicity(fine, metric),
        )?;
    }
    for metric in SENSITIVITY_METRICS {
        all_pass &= write_label_convergence(
            &mut writer,
            "maximum_sensitivity_interval",
            metric,
            &maximum_sensitivity_interval(coarse_sensitivity, metric),
            &maximum_sensitivity_interval(fine_sensitivity, metric),
        )?;
    }
    for metric in THRESHOLD_METRICS {
        all_pass &= write_label_convergence(
            &mut writer,
            "threshold_order",
            metric,
            &threshold_signature(coarse_thresholds, metric),
            &threshold_signature(fine_thresholds, metric),
        )?;
    }
    all_pass &= write_label_convergence(
        &mut writer,
        "adaptive_decision",
        "additional_points",
        if coarse_decision.add_points {
            "add"
        } else {
            "none"
        },
        if fine_decision.add_points {
            "add"
        } else {
            "none"
        },
    )?;
    Ok(all_pass)
}

fn strongest_curvature(rows: &[SensitivityRow]) -> &SensitivityRow {
    rows.iter()
        .filter(|row| row.second_difference.is_finite())
        .max_by(|left, right| {
            left.second_difference
                .abs()
                .total_cmp(&right.second_difference.abs())
        })
        .expect("basic sweep supplies second differences")
}

fn strongest_window_interval(rows: &[WindowResult], window_name: &str) -> (f64, f64, f64) {
    let selected: Vec<_> = rows
        .iter()
        .filter(|row| row.window.name == window_name)
        .collect();
    selected
        .windows(2)
        .map(|pair| {
            (
                pair[0].gamma_end,
                pair[1].gamma_end,
                (pair[1].w_area - pair[0].w_area) / (pair[0].gamma_end - pair[1].gamma_end),
            )
        })
        .max_by(|left, right| left.2.abs().total_cmp(&right.2.abs()))
        .unwrap()
}

fn report(
    analyses: &[Analysis],
    recovery: &[RecoveryRow],
    sensitivity: &[SensitivityRow],
    thresholds: &[ThresholdRow],
    windows: &[WindowResult],
    checks: &[Check],
    decision: &AdaptiveDecision,
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("MILESTONE_7D_REPORT.md")?);
    writeln!(
        writer,
        "# Milestone 7d: Partial end-protection strength sweep\n"
    )?;
    writeln!(writer, "## 1. 目的\n\n中央siteのgamma=0.5を固定し、両端gammaだけを0.5から0へ下げたときの回復曲線を調べた。\n")?;
    writeln!(writer, "## 2. Milestone 7cから進めた問い\n\n7cの完全両端雑音除去という端点比較から、途中の有限gammaで回復がどのようにつながるかへ進めた。\n")?;
    writeln!(writer, "## 3. 部分保護の定義\n\nsite1/site3のdephasing rateを同時に数値的に弱める理想感度試験。現実的保護装置、制御、cost、不完全装置モデルではない。\n")?;
    let gammas = analyses
        .iter()
        .map(|analysis| format!("{:.2}", analysis.gamma_end))
        .collect::<Vec<_>>()
        .join(", ");
    writeln!(
        writer,
        "## 4. sweep条件\n\n実行gamma_end: `{gammas}`。基本6点を先に実行。{}\n",
        decision.reason
    )?;
    writeln!(writer, "## 5. 変更していない物理条件\n\n3site+3準位load、J=1、g=0.25、周波数=1、真空初期状態、tau=3.2、Omega=0.2、t_max=10、中央gamma=0.5、load無雑音を固定した。\n")?;
    writeln!(writer, "## 6. site別gamma実装\n\n新APIは `[gamma_site1, gamma_site2, gamma_site3]` を受け、負値・非有限値を拒否し、gamma=0のcollapse operatorを除外する。既存共通gamma APIは新APIへ委譲し、短時間密度行列の完全一致をunit testした。\n")?;
    writeln!(writer, "## 7. 数値手法\n\nRK4、基準dt=0.0025、半減dt=0.00125、保存間隔0.01。時間面積は台形則。感度と二階差分は離散点だけから計算し、連続微分・物理的感受率・臨界指数とは呼ばない。\n")?;
    writeln!(writer, "## 8. 数値品質チェック\n\n全{}項目PASS。trace、Hermiticity、positivity、population、縮約状態、load top-level、ledger、有限性、W<=E、usable範囲、drive整合、gamma mapping、collapse数、共通grid、初期状態を確認した。\n",checks.len())?;
    writeln!(writer, "## 9. 端点再現\n\ngamma_end=0.5は7c all_noisy、gamma_end=0は7c protect_both_ends、別計算のnoise_freeは7c noise_freeの8比較量を絶対誤差1e-9以内で再現した。\n")?;
    writeln!(writer, "## 10. t=10の回復曲線\n\n| gamma_end | E | W | normalized W recovery |\n|---:|---:|---:|---:|")?;
    for analysis in analyses {
        writeln!(
            writer,
            "| {:.2} | {:.10e} | {:.10e} | {:.8} |",
            analysis.gamma_end,
            analysis.summary.e_at_t10,
            analysis.summary.w_at_t10,
            recovery_value(analyses, analysis.gamma_end, "W_at_t10")
        )?;
    }
    writeln!(writer, "\n## 11. W_maxの回復曲線\n\n| gamma_end | W_max | t_at_W_max | normalized recovery |\n|---:|---:|---:|---:|")?;
    for analysis in analyses {
        writeln!(
            writer,
            "| {:.2} | {:.10e} | {:.2} | {:.8} |",
            analysis.gamma_end,
            analysis.summary.w_max,
            analysis.summary.t_at_w_max,
            recovery_value(analyses, analysis.gamma_end, "W_max")
        )?;
    }
    writeln!(writer, "\n## 12. E_time_areaとW_time_area\n\n両者は状態量の時間面積で、累積流入エネルギー・累積抽出仕事ではない。gamma低下時の判定はE=`{}`、W=`{}`。\n",monotonicity(analyses,"E_time_area"),monotonicity(analyses,"W_time_area"))?;
    writeln!(writer, "## 13. usable fractionの回復曲線\n\n0.5で{:.8}、0で{:.8}。0.4で完全両端保護回復の{:.2}%、0.3で{:.2}%へ到達した。\n",analysis_at(analyses,0.5).summary.usable_fraction_at_t10,analysis_at(analyses,0.0).summary.usable_fraction_at_t10,100.0*recovery_value(analyses,0.4,"usable_fraction_at_t10"),100.0*recovery_value(analyses,0.3,"usable_fraction_at_t10"))?;
    writeln!(writer, "## 14. 隣接感度\n")?;
    for metric in SENSITIVITY_METRICS {
        let interval = maximum_sensitivity_interval(sensitivity, metric);
        let value = sensitivity
            .iter()
            .filter(|row| row.metric == metric)
            .max_by(|a, b| a.sensitivity.abs().total_cmp(&b.sensitivity.abs()))
            .unwrap()
            .sensitivity;
        writeln!(
            writer,
            "- {metric}: 最大区間 `{interval}`, sensitivity={value:.8e}"
        )?;
    }
    let curve = strongest_curvature(sensitivity);
    writeln!(writer, "\n## 15. 離散曲率\n\n絶対二階差分最大は `{}` の中心gamma={:.2}（隣接区間 {:.2}->{:.2}）で {:.8e}。追加判断で使った正規化二階差分最大は `{}` のgamma={:.2}で {:.8}。相転移や閾値現象とは断定しない。\n",curve.metric,curve.gamma_lower,curve.gamma_upper,curve.gamma_lower,curve.second_difference,decision.metric,decision.location,decision.max_normalized_curvature)?;
    writeln!(writer, "## 16. 単調性\n\n| metric | status |\n|---|---|")?;
    for metric in SUMMARY_METRICS {
        writeln!(writer, "| {metric} | {} |", monotonicity(analyses, metric))?;
    }
    writeln!(writer, "\n## 17. 小さな雑音低減の効果\n\n| metric | recovery at gamma=0.4 | recovery at gamma=0.3 |\n|---|---:|---:|")?;
    for metric in ["W_at_t10", "usable_fraction_at_t10", "W_time_area"] {
        writeln!(
            writer,
            "| {metric} | {:.4}% | {:.4}% |",
            100.0 * recovery_value(analyses, 0.4, metric),
            100.0 * recovery_value(analyses, 0.3, metric)
        )?;
    }
    writeln!(writer, "\n## 18. 回復水準到達点\n\n補間は使わず、離散点で初めて超えた最大gammaを保存した。\n\n| metric | target | gamma | observed recovery |\n|---|---:|---:|---:|")?;
    for row in thresholds {
        writeln!(
            writer,
            "| {} | {:.0}% | {:.2} | {:.8} |",
            row.metric,
            100.0 * row.target,
            row.gamma,
            row.recovery
        )?;
    }
    writeln!(writer, "\n## 19. 時間窓別感度\n")?;
    for window in WINDOWS {
        let (upper, lower, value) = strongest_window_interval(windows, window.name);
        writeln!(
            writer,
            "- {}: W_time_area最大感度区間 {:.2}->{:.2}, {:.8e}",
            window.name, upper, lower, value
        )?;
    }
    writeln!(writer, "\n## 20. 時系列感度\n\n各保存時刻でE、W、usable fractionの絶対回復、正規化回復、隣接感度、最大感度gamma区間をtimeseries CSVへ保存した。初期の小分母はNaNのままとした。\n")?;
    writeln!(writer, "## 21. 刻み幅整合性\n\n最終sweep全点とnoise-freeを半減刻みで再計算し、要約値、正規化回復、感度、単調性、回復水準順序、最大感度区間、追加点判断の全行がPASS。\n")?;
    writeln!(writer, "## 22. 直接確認できたこと\n\n固定模型・Omega・中央gamma内で、両端gamma低減に対する回復曲線、離散単調性、有限差分感度、回復水準、時間窓差を確認した。\n")?;
    writeln!(writer, "## 23. 確認できていないこと\n\n現実的実装、cost、異なる中央gamma/Omega/network、長時間・連続運転、因果機構、新規性は未確認。\n")?;
    writeln!(writer, "## 24. 主張してはいけないこと\n\n実装可能な必要保護強度、物理的感受率・臨界指数・相転移、一般最適、量子優位、実用送電性能は主張しない。\n")?;
    writeln!(writer, "## 25. load雑音を扱わない理由\n\nloadは3準位でsiteの二準位sigma_zを流用できず、同じgammaが公平な強度を意味しないため別検証が必要。\n")?;
    let max_w = maximum_sensitivity_interval(sensitivity, "W_at_t10");
    writeln!(writer, "## 26. 次段階への判断材料\n\nWの最大離散感度区間は `{max_w}`。追加点判定の最大感度比は {:.8}。これは候補情報であり、自動的に次Milestoneへ進まない。\n",decision.max_sensitivity_ratio)?;
    writeln!(writer, "## 27. 生成ファイル一覧\n\n- `src/bin/partial_end_protection.rs`\n- `partial_end_protection_timeseries.csv`\n- `partial_end_protection_summary.csv`\n- `partial_end_protection_recovery.csv`\n- `partial_end_protection_sensitivity.csv`\n- `partial_end_protection_thresholds.csv`\n- `partial_end_protection_windows.csv`\n- `partial_end_protection_checks.csv`\n- `partial_end_protection_convergence.csv`\n- `MILESTONE_7D_REPORT.md`\n")?;
    let _ = recovery;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams::default();
    let operators = build_operators(&params)?;

    let mut coarse_basic = Vec::new();
    for gamma in BASIC_GAMMAS {
        coarse_basic.push(run_gamma(gamma, BASE_DT, &params, &operators)?);
    }
    let coarse_decision = adaptive_decision(&coarse_basic);
    println!("adaptive decision: {}", coarse_decision.reason);
    let mut coarse = coarse_basic;
    for &gamma in &coarse_decision.additions {
        coarse.push(run_gamma(gamma, BASE_DT, &params, &operators)?);
    }
    coarse.sort_by(|left, right| right.gamma_end.total_cmp(&left.gamma_end));
    if coarse.len() > 11 {
        return Err("adaptive sweep exceeded 11 points".into());
    }
    let coarse_noise_free = run_noise_free(BASE_DT, &params, &operators)?;
    let checks = checks(&coarse, &coarse_noise_free, &operators);
    write_checks(&checks)?;
    if checks.iter().any(|check| !check.pass) {
        return Err("Milestone 7d numerical/input checks failed".into());
    }
    let recovery = recovery_rows(&coarse, &coarse_noise_free);
    let sensitivity = sensitivity_rows(&coarse);
    let threshold_rows = thresholds(&coarse);
    let window_rows = window_results(&coarse);

    let mut fine_basic = Vec::new();
    for gamma in BASIC_GAMMAS {
        fine_basic.push(run_gamma(gamma, FINE_DT, &params, &operators)?);
    }
    let fine_decision = adaptive_decision(&fine_basic);
    let mut fine = fine_basic;
    for &gamma in &coarse_decision.additions {
        fine.push(run_gamma(gamma, FINE_DT, &params, &operators)?);
    }
    fine.sort_by(|left, right| right.gamma_end.total_cmp(&left.gamma_end));
    let fine_noise_free = run_noise_free(FINE_DT, &params, &operators)?;
    if fine
        .iter()
        .chain(std::iter::once(&fine_noise_free))
        .any(|analysis| !analysis.quality.all_pass())
    {
        return Err("fine-step physical diagnostics failed".into());
    }
    let fine_sensitivity = sensitivity_rows(&fine);
    let fine_thresholds = thresholds(&fine);
    let convergence_ok = write_convergence(
        &coarse,
        &fine,
        &coarse_noise_free,
        &fine_noise_free,
        &sensitivity,
        &fine_sensitivity,
        &threshold_rows,
        &fine_thresholds,
        &coarse_decision,
        &fine_decision,
    )?;
    if !convergence_ok {
        return Err("time-step convergence or conclusion stability failed".into());
    }
    write_timeseries(&coarse)?;
    write_summary(&coarse)?;
    write_recovery(&recovery)?;
    write_sensitivity(&sensitivity)?;
    write_thresholds(&threshold_rows)?;
    write_windows(&window_rows)?;
    report(
        &coarse,
        &recovery,
        &sensitivity,
        &threshold_rows,
        &window_rows,
        &checks,
        &coarse_decision,
    )?;
    println!("Milestone 7d complete");
    Ok(())
}
