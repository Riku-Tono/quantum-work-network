use std::fs::File;
use std::io::{BufWriter, Write};

use quantum_work_network::coherent_drive::{
    run_coherent_drive_with_noise_sites, CoherentDriveConfig, CoherentDriveRun,
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
const GAMMA_PHI: f64 = 0.5;
const OMEGA: f64 = 0.2;
const REDUCED_TOL: f64 = 1.0e-9;
const POPULATION_TOL: f64 = 1.0e-8;
const BASELINE_TOL: f64 = 1.0e-9;
const SYNERGY_ABSOLUTE_TOL: f64 = CONVERGENCE_ABSOLUTE_TOLERANCE;
const SYNERGY_RELATIVE_TOL: f64 = CONVERGENCE_RELATIVE_TOLERANCE;

const NOISE_FREE_E: f64 = 5.4450767877898487e-2;
const NOISE_FREE_W: f64 = 5.2798274942446315e-2;
const NOISE_FREE_EIN: f64 = 7.5315933437092850e-2;
const NOISE_FREE_COHERENCE: f64 = 4.6453837009096421e-1;
// Milestone 5b B: Omega=0.2, gamma=0.5 on all three sites, dt=0.005.
const EXISTING_ALL_NOISY_E: f64 = 1.2596874860520470e-2;
const EXISTING_ALL_NOISY_W: f64 = 2.3652476826588921e-3;
const EXISTING_ALL_NOISY_EIN: f64 = 5.9618618774949499e-2;
const EXISTING_ALL_NOISY_COHERENCE: f64 = 9.7568729204613736e-2;

#[derive(Clone, Copy)]
struct Condition {
    name: &'static str,
    protected_sites: &'static str,
    noisy_sites_label: &'static str,
    noisy_sites: &'static [usize],
}

const CONDITIONS: [Condition; 5] = [
    Condition {
        name: "all_noisy",
        protected_sites: "none",
        noisy_sites_label: "site1+site2+site3",
        noisy_sites: &[0, 1, 2],
    },
    Condition {
        name: "protect_entrance",
        protected_sites: "site1",
        noisy_sites_label: "site2+site3",
        noisy_sites: &[1, 2],
    },
    Condition {
        name: "protect_exit",
        protected_sites: "site3",
        noisy_sites_label: "site1+site2",
        noisy_sites: &[0, 1],
    },
    Condition {
        name: "protect_both_ends",
        protected_sites: "site1+site3",
        noisy_sites_label: "site2",
        noisy_sites: &[1],
    },
    Condition {
        name: "noise_free",
        protected_sites: "site1+site2+site3",
        noisy_sites_label: "none",
        noisy_sites: &[],
    },
];

const PROTECTION_INDICES: [usize; 3] = [1, 2, 3];

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
    diagonal_w_at_t10: f64,
    coherence_w_at_t10: f64,
    coherence_l1_at_t10: f64,
    usable_fraction_at_t10: f64,
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
    condition: Condition,
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

#[derive(Clone)]
struct WindowResult {
    condition: &'static str,
    window: Window,
    point_count: usize,
    mean_e: f64,
    mean_w: f64,
    mean_use: f64,
    mean_coherence: f64,
    e_area: f64,
    w_area: f64,
    mean_site: [f64; 3],
    delta_e_area: f64,
    delta_w_area: f64,
    delta_mean_use: f64,
}

#[derive(Clone)]
struct RecoveryRow {
    condition: &'static str,
    metric: &'static str,
    evaluation: &'static str,
    baseline: f64,
    value: f64,
    noise_free: f64,
    absolute: f64,
    normalized: f64,
    residual: f64,
}

#[derive(Clone)]
struct SynergyRow {
    metric: &'static str,
    evaluation: &'static str,
    delta_entrance: f64,
    delta_exit: f64,
    delta_both: f64,
    synergy: f64,
    normalized: f64,
    classification: &'static str,
}

#[derive(Clone)]
struct Check {
    scope: String,
    condition: String,
    check: String,
    value: f64,
    tolerance: String,
    pass: bool,
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

fn equivalent(left: f64, right: f64) -> bool {
    (left - right).abs() <= CONVERGENCE_ABSOLUTE_TOLERANCE
        || (left - right).abs() / left.abs().max(right.abs()).max(SIGNAL_TOLERANCE)
            <= CONVERGENCE_RELATIVE_TOLERANCE
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
    condition: Condition,
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
        diagonal_w_at_t10: final_row.diagonal_ergotropy,
        coherence_w_at_t10: final_row.coherence_ergotropy,
        coherence_l1_at_t10: final_row.coherence_l1,
        usable_fraction_at_t10: final_row.usable_fraction,
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
    let drive_consistency =
        (summary.drive_energy_in_at_t10 - run.summary.drive_energy.energy_in).abs() <= 1.0e-12;
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
        drive_consistency,
        max_reduced_difference,
    };
    Ok(Analysis {
        condition,
        run,
        rows,
        summary,
        quality,
    })
}

fn config(dt: f64) -> CoherentDriveConfig {
    let mut config = CoherentDriveConfig::milestone_5b(GAMMA_PHI, dt);
    config.omega0 = OMEGA;
    config.save_interval = SAVE_INTERVAL;
    config
}

fn run_condition(
    condition: Condition,
    dt: f64,
    params: &ModelParams,
    operators: &Operators,
) -> Result<Analysis, Box<dyn std::error::Error>> {
    println!("running {} at dt={dt}", condition.name);
    let run = run_coherent_drive_with_noise_sites(params, config(dt), condition.noisy_sites)?;
    analyze(condition, run, params, operators)
}

fn analysis<'a>(analyses: &'a [Analysis], name: &str) -> &'a Analysis {
    analyses
        .iter()
        .find(|analysis| analysis.condition.name == name)
        .expect("condition exists")
}

fn summary_metric(summary: &Summary, metric: &str, evaluation: &str) -> f64 {
    match (metric, evaluation) {
        ("E", "t10") => summary.e_at_t10,
        ("E", "time_area") => summary.e_time_area,
        ("W", "t10") => summary.w_at_t10,
        ("W", "W_max") => summary.w_max,
        ("W", "time_area") => summary.w_time_area,
        ("usable_fraction", "t10") => summary.usable_fraction_at_t10,
        ("usable_fraction", "at_W_max") => summary.usable_fraction_at_w_max,
        _ => unreachable!(),
    }
}

fn recovery_rows(analyses: &[Analysis]) -> Result<Vec<RecoveryRow>, Box<dyn std::error::Error>> {
    let baseline = &analysis(analyses, "all_noisy").summary;
    let free = &analysis(analyses, "noise_free").summary;
    let metrics = [
        ("E", "t10"),
        ("E", "time_area"),
        ("W", "t10"),
        ("W", "W_max"),
        ("W", "time_area"),
        ("usable_fraction", "t10"),
        ("usable_fraction", "at_W_max"),
    ];
    let mut rows = Vec::new();
    for &index in &PROTECTION_INDICES {
        for (metric, evaluation) in metrics {
            let noisy_value = summary_metric(baseline, metric, evaluation);
            let protected_value = summary_metric(&analyses[index].summary, metric, evaluation);
            let free_value = summary_metric(free, metric, evaluation);
            let denominator = free_value - noisy_value;
            if denominator.abs() <= SIGNAL_TOLERANCE {
                return Err(
                    format!("recovery denominator too small for {metric}/{evaluation}").into(),
                );
            }
            rows.push(RecoveryRow {
                condition: analyses[index].condition.name,
                metric,
                evaluation,
                baseline: noisy_value,
                value: protected_value,
                noise_free: free_value,
                absolute: protected_value - noisy_value,
                normalized: (protected_value - noisy_value) / denominator,
                residual: 1.0 - protected_value / free_value,
            });
        }
    }
    Ok(rows)
}

fn classify_synergy(synergy: f64, normalized: f64) -> &'static str {
    if synergy.abs() <= SYNERGY_ABSOLUTE_TOL || normalized.abs() <= SYNERGY_RELATIVE_TOL {
        "approximately_additive"
    } else if synergy > 0.0 {
        "positive_nonadditivity"
    } else {
        "negative_nonadditivity"
    }
}

fn synergy_rows(analyses: &[Analysis]) -> Result<Vec<SynergyRow>, Box<dyn std::error::Error>> {
    let baseline = &analysis(analyses, "all_noisy").summary;
    let entrance = &analysis(analyses, "protect_entrance").summary;
    let exit = &analysis(analyses, "protect_exit").summary;
    let both = &analysis(analyses, "protect_both_ends").summary;
    let free = &analysis(analyses, "noise_free").summary;
    let metrics = [
        ("E", "t10"),
        ("W", "t10"),
        ("usable_fraction", "t10"),
        ("E", "time_area"),
        ("W", "time_area"),
    ];
    let mut rows = Vec::new();
    for (metric, evaluation) in metrics {
        let base = summary_metric(baseline, metric, evaluation);
        let d_in = summary_metric(entrance, metric, evaluation) - base;
        let d_out = summary_metric(exit, metric, evaluation) - base;
        let d_both = summary_metric(both, metric, evaluation) - base;
        let synergy = d_both - d_in - d_out;
        let denominator = summary_metric(free, metric, evaluation) - base;
        if denominator.abs() <= SIGNAL_TOLERANCE {
            return Err(format!("synergy denominator too small for {metric}/{evaluation}").into());
        }
        let normalized = synergy / denominator;
        rows.push(SynergyRow {
            metric,
            evaluation,
            delta_entrance: d_in,
            delta_exit: d_out,
            delta_both: d_both,
            synergy,
            normalized,
            classification: classify_synergy(synergy, normalized),
        });
    }
    Ok(rows)
}

fn mean_finite(values: impl Iterator<Item = f64>) -> f64 {
    let values: Vec<_> = values.filter(|value| value.is_finite()).collect();
    if values.is_empty() {
        f64::NAN
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

fn window_values(
    analysis: &Analysis,
    window: Window,
) -> (usize, f64, f64, f64, f64, f64, [f64; 3]) {
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
            .map(|pair| 0.5 * (value(pair[0]) + value(pair[1])) * (pair[1].time - pair[0].time))
            .sum()
    };
    (
        mean_rows.len(),
        mean_finite(mean_rows.iter().map(|row| row.load_energy)),
        mean_finite(mean_rows.iter().map(|row| row.load_ergotropy)),
        mean_finite(mean_rows.iter().map(|row| row.usable_fraction)),
        mean_finite(mean_rows.iter().map(|row| row.coherence_l1)),
        area(|row| row.load_energy),
        [
            mean_finite(mean_rows.iter().map(|row| row.site[0])),
            mean_finite(mean_rows.iter().map(|row| row.site[1])),
            mean_finite(mean_rows.iter().map(|row| row.site[2])),
        ],
    )
}

fn windows(analyses: &[Analysis]) -> Vec<WindowResult> {
    let mut raw = Vec::new();
    for analysis in analyses {
        for window in WINDOWS {
            let (count, mean_e, mean_w, mean_use, mean_coherence, e_area, mean_site) =
                window_values(analysis, window);
            let area_rows: Vec<_> = analysis
                .rows
                .iter()
                .filter(|row| row.time >= window.start && row.time <= window.end)
                .collect();
            let w_area = area_rows
                .windows(2)
                .map(|pair| {
                    0.5 * (pair[0].load_ergotropy + pair[1].load_ergotropy)
                        * (pair[1].time - pair[0].time)
                })
                .sum();
            raw.push((
                analysis.condition.name,
                window,
                count,
                mean_e,
                mean_w,
                mean_use,
                mean_coherence,
                e_area,
                w_area,
                mean_site,
            ));
        }
    }
    let mut results = Vec::new();
    for row in &raw {
        let baseline = raw
            .iter()
            .find(|candidate| candidate.0 == "all_noisy" && candidate.1.name == row.1.name)
            .unwrap();
        results.push(WindowResult {
            condition: row.0,
            window: row.1,
            point_count: row.2,
            mean_e: row.3,
            mean_w: row.4,
            mean_use: row.5,
            mean_coherence: row.6,
            e_area: row.7,
            w_area: row.8,
            mean_site: row.9,
            delta_e_area: row.7 - baseline.7,
            delta_w_area: row.8 - baseline.8,
            delta_mean_use: row.5 - baseline.5,
        });
    }
    results
}

fn ranking_string(analyses: &[Analysis], value: impl Fn(&Summary) -> f64) -> String {
    let maximum = PROTECTION_INDICES
        .iter()
        .map(|&index| value(&analyses[index].summary))
        .max_by(f64::total_cmp)
        .unwrap();
    let tied: Vec<_> = PROTECTION_INDICES
        .iter()
        .filter(|&&index| equivalent(value(&analyses[index].summary), maximum))
        .map(|&index| analyses[index].condition.name)
        .collect();
    if tied.len() == 1 {
        tied[0].to_string()
    } else {
        format!("tie({})", tied.join("+"))
    }
}

fn pair_label(entrance: f64, exit: f64) -> String {
    if !entrance.is_finite() || !exit.is_finite() {
        "undefined".to_string()
    } else if equivalent(entrance, exit) {
        "tie".to_string()
    } else if entrance > exit {
        "protect_entrance".to_string()
    } else {
        "protect_exit".to_string()
    }
}

fn confirmed_switches(times: &[f64], labels: &[String]) -> Vec<(f64, String, String)> {
    let mut stable: Option<String> = None;
    let mut switches = Vec::new();
    for index in 0..=labels.len().saturating_sub(5) {
        let candidate = &labels[index];
        if candidate == "undefined"
            || labels[index..index + 5]
                .iter()
                .any(|label| label != candidate)
        {
            continue;
        }
        match &stable {
            None => stable = Some(candidate.clone()),
            Some(previous) if previous != candidate => {
                switches.push((times[index], previous.clone(), candidate.clone()));
                stable = Some(candidate.clone());
            }
            _ => {}
        }
    }
    switches
}

fn checks(
    analyses: &[Analysis],
    operators: &Operators,
) -> Result<Vec<Check>, Box<dyn std::error::Error>> {
    let mut checks = Vec::new();
    for analysis in analyses {
        let q = &analysis.quality;
        let s = &analysis.summary;
        let items = [
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
        ];
        for (name, value, tolerance, result) in items {
            checks.push(Check {
                scope: "condition".to_string(),
                condition: analysis.condition.name.to_string(),
                check: name.to_string(),
                value,
                tolerance: tolerance.to_string(),
                pass: result,
            });
        }
        let expected_count = match analysis.condition.name {
            "all_noisy" => 3,
            "protect_entrance" | "protect_exit" => 2,
            "protect_both_ends" => 1,
            "noise_free" => 0,
            _ => unreachable!(),
        };
        checks.push(Check {
            scope: "condition".to_string(),
            condition: analysis.condition.name.to_string(),
            check: "collapse_operator_count".to_string(),
            value: analysis.condition.noisy_sites.len() as f64,
            tolerance: format!("expected={expected_count}"),
            pass: analysis.condition.noisy_sites.len() == expected_count,
        });
        let scale = C64::new((GAMMA_PHI / 2.0).sqrt(), 0.0);
        let embedding_ok = analysis.condition.noisy_sites.iter().all(|&site| {
            if site >= 3 {
                return false;
            }
            let selected = &operators.sigma_z_sites[site] * scale;
            (0..3).all(|other| {
                other == site
                    || frobenius_norm(&(&selected - &operators.sigma_z_sites[other] * scale))
                        > 1.0e-12
            })
        });
        checks.push(Check {
            scope: "condition".to_string(),
            condition: analysis.condition.name.to_string(),
            check: "noise_site_embedding_and_labels".to_string(),
            value: 0.0,
            tolerance: analysis.condition.noisy_sites_label.to_string(),
            pass: embedding_ok,
        });
    }
    let reference = &analyses[0];
    let same_grid = analyses[1..].iter().all(|analysis| {
        analysis.rows.len() == reference.rows.len()
            && analysis
                .rows
                .iter()
                .zip(&reference.rows)
                .all(|(left, right)| left.time == right.time)
    });
    let same_initial = analyses[1..].iter().all(|analysis| {
        frobenius_norm(&(&analysis.run.states[0].rho - &reference.run.states[0].rho)) <= 1.0e-14
    });
    for (name, result) in [
        ("common_time_grid", same_grid),
        ("initial_state_consistency", same_initial),
    ] {
        checks.push(Check {
            scope: "global".to_string(),
            condition: "all".to_string(),
            check: name.to_string(),
            value: 0.0,
            tolerance: "exact/Frobenius<=1e-14".to_string(),
            pass: result,
        });
    }
    let free = &analysis(analyses, "noise_free").summary;
    for (name, actual, expected) in [
        ("noise_free_E", free.e_at_t10, NOISE_FREE_E),
        ("noise_free_W", free.w_at_t10, NOISE_FREE_W),
        (
            "noise_free_Ein",
            free.drive_energy_in_at_t10,
            NOISE_FREE_EIN,
        ),
        (
            "noise_free_coherence",
            free.coherence_l1_at_t10,
            NOISE_FREE_COHERENCE,
        ),
    ] {
        checks.push(Check {
            scope: "baseline".to_string(),
            condition: "noise_free".to_string(),
            check: name.to_string(),
            value: (actual - expected).abs(),
            tolerance: "absolute<=1e-9".to_string(),
            pass: (actual - expected).abs() <= BASELINE_TOL,
        });
    }
    let all = &analysis(analyses, "all_noisy").summary;
    for (name, actual, existing) in [
        (
            "all_noisy_E_vs_M5b_B_dt0.005",
            all.e_at_t10,
            EXISTING_ALL_NOISY_E,
        ),
        (
            "all_noisy_W_vs_M5b_B_dt0.005",
            all.w_at_t10,
            EXISTING_ALL_NOISY_W,
        ),
        (
            "all_noisy_Ein_vs_M5b_B_dt0.005",
            all.drive_energy_in_at_t10,
            EXISTING_ALL_NOISY_EIN,
        ),
        (
            "all_noisy_coherence_vs_M5b_B_dt0.005",
            all.coherence_l1_at_t10,
            EXISTING_ALL_NOISY_COHERENCE,
        ),
    ] {
        checks.push(Check {
            scope: "baseline".to_string(),
            condition: "all_noisy".to_string(),
            check: name.to_string(),
            value: (actual - existing).abs(),
            tolerance: "existing convergence abs or rel tolerance".to_string(),
            pass: converged(existing, actual),
        });
    }
    Ok(checks)
}

fn write_timeseries(analyses: &[Analysis]) -> std::io::Result<()> {
    let all = analysis(analyses, "all_noisy");
    let free = analysis(analyses, "noise_free");
    let mut writer = BufWriter::new(File::create("ideal_partial_protection_timeseries.csv")?);
    writeln!(writer, "condition,protected_sites,noisy_sites,time,gamma_phi,Omega,load_energy,load_ergotropy,load_diagonal_ergotropy,load_coherence_ergotropy,load_coherence_l1,usable_fraction,drive_energy_in,W_over_Ein,site1_population,site2_population,site3_population,total_chain_population,load_top_level_population,trace_error,hermiticity_error,min_eigenvalue,energy_ledger_residual,Delta_E,Delta_W,Delta_use,Recovery_E,Recovery_W,Recovery_use,ResidualLoss_E,ResidualLoss_W,ResidualLoss_use")?;
    for condition in analyses {
        for index in 0..condition.rows.len() {
            let row = &condition.rows[index];
            let baseline = &all.rows[index];
            let upper = &free.rows[index];
            let delta_e = row.load_energy - baseline.load_energy;
            let delta_w = row.load_ergotropy - baseline.load_ergotropy;
            let delta_use =
                if row.usable_fraction.is_finite() && baseline.usable_fraction.is_finite() {
                    row.usable_fraction - baseline.usable_fraction
                } else {
                    f64::NAN
                };
            writeln!(writer, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                condition.condition.name, condition.condition.protected_sites, condition.condition.noisy_sites_label,
                n(row.time), n(GAMMA_PHI), n(OMEGA), n(row.load_energy), n(row.load_ergotropy),
                n(row.diagonal_ergotropy), n(row.coherence_ergotropy), n(row.coherence_l1),
                n(row.usable_fraction), n(row.drive_energy_in), n(row.w_over_ein), n(row.site[0]), n(row.site[1]), n(row.site[2]),
                n(row.total_chain_population), n(row.load_top_population), n(row.trace_error), n(row.hermiticity_error), n(row.minimum_eigenvalue), n(row.ledger_residual),
                n(delta_e), n(delta_w), n(delta_use),
                n(safe_ratio(delta_e, upper.load_energy-baseline.load_energy)),
                n(safe_ratio(delta_w, upper.load_ergotropy-baseline.load_ergotropy)),
                n(safe_ratio(delta_use, upper.usable_fraction-baseline.usable_fraction)),
                n(1.0-safe_ratio(row.load_energy, upper.load_energy)),
                n(1.0-safe_ratio(row.load_ergotropy, upper.load_ergotropy)),
                n(1.0-safe_ratio(row.usable_fraction, upper.usable_fraction)))?;
        }
    }
    Ok(())
}

fn write_summary(analyses: &[Analysis]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("ideal_partial_protection_summary.csv")?);
    writeln!(writer, "condition,E_at_t10,W_at_t10,diagonal_W_at_t10,coherence_W_at_t10,usable_fraction_at_t10,coherence_L1_at_t10,drive_energy_in_at_t10,W_over_Ein_at_t10,W_max,t_at_W_max,E_at_W_max,usable_fraction_at_W_max,E_time_area,W_time_area,max_load_top_level_population,max_trace_error,max_hermiticity_error,minimum_density_eigenvalue,max_abs_energy_ledger_residual")?;
    for analysis in analyses {
        let s = &analysis.summary;
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            analysis.condition.name,
            n(s.e_at_t10),
            n(s.w_at_t10),
            n(s.diagonal_w_at_t10),
            n(s.coherence_w_at_t10),
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
    let mut writer = BufWriter::new(File::create("ideal_partial_protection_recovery.csv")?);
    writeln!(writer, "condition,metric,evaluation_point,noisy_baseline,protected_value,noise_free_value,absolute_recovery,normalized_recovery,residual_loss_to_noise_free")?;
    for row in rows {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{}",
            row.condition,
            row.metric,
            row.evaluation,
            n(row.baseline),
            n(row.value),
            n(row.noise_free),
            n(row.absolute),
            n(row.normalized),
            n(row.residual)
        )?;
    }
    Ok(())
}

fn write_synergy(rows: &[SynergyRow]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("ideal_partial_protection_synergy.csv")?);
    writeln!(writer, "metric,evaluation_point,delta_entrance,delta_exit,delta_both,synergy,normalized_synergy,additivity_classification,absolute_tol,relative_tol")?;
    for row in rows {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{}",
            row.metric,
            row.evaluation,
            n(row.delta_entrance),
            n(row.delta_exit),
            n(row.delta_both),
            n(row.synergy),
            n(row.normalized),
            row.classification,
            n(SYNERGY_ABSOLUTE_TOL),
            n(SYNERGY_RELATIVE_TOL)
        )?;
    }
    Ok(())
}

fn write_windows(rows: &[WindowResult]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("ideal_partial_protection_windows.csv")?);
    writeln!(writer, "condition,window_name,time_start,time_end,point_count,mean_load_energy,mean_load_ergotropy,mean_usable_fraction,mean_coherence_L1,E_time_area,W_time_area,mean_site1_population,mean_site2_population,mean_site3_population,Delta_E_time_area_from_all_noisy,Delta_W_time_area_from_all_noisy,Delta_mean_use_from_all_noisy")?;
    for row in rows {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            row.condition,
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
            n(row.delta_e_area),
            n(row.delta_w_area),
            n(row.delta_mean_use)
        )?;
    }
    Ok(())
}

fn ranked(
    analyses: &[Analysis],
    value: impl Fn(&Summary) -> f64,
) -> Vec<(usize, &'static str, f64, bool)> {
    let mut values: Vec<_> = PROTECTION_INDICES
        .iter()
        .map(|&index| {
            (
                analyses[index].condition.name,
                value(&analyses[index].summary),
            )
        })
        .collect();
    values.sort_by(|left, right| right.1.total_cmp(&left.1));
    values
        .iter()
        .enumerate()
        .map(|(index, (name, value))| {
            let tie = values.iter().any(|(other_name, other_value)| {
                other_name != name && equivalent(*value, *other_value)
            });
            (index + 1, *name, *value, tie)
        })
        .collect()
}

fn write_rankings(analyses: &[Analysis]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("ideal_partial_protection_rankings.csv")?);
    writeln!(
        writer,
        "evaluation_scope,metric,rank,condition,value,tie,time"
    )?;
    for (scope, metric, values) in [
        ("t10", "E_recovery", ranked(analyses, |s| s.e_at_t10)),
        ("t10", "W_recovery", ranked(analyses, |s| s.w_at_t10)),
        (
            "t10",
            "usable_fraction_recovery",
            ranked(analyses, |s| s.usable_fraction_at_t10),
        ),
        (
            "time_area",
            "E_time_area_recovery",
            ranked(analyses, |s| s.e_time_area),
        ),
        (
            "time_area",
            "W_time_area_recovery",
            ranked(analyses, |s| s.w_time_area),
        ),
    ] {
        for (rank, condition, value, tie) in values {
            writeln!(
                writer,
                "{scope},{metric},{rank},{condition},{},{tie},NaN",
                n(value)
            )?;
        }
    }
    let entrance = analysis(analyses, "protect_entrance");
    let exit = analysis(analyses, "protect_exit");
    for (index, (left, right)) in entrance.rows.iter().zip(&exit.rows).enumerate() {
        for (metric, left_value, right_value) in [
            ("E_recovery", left.load_energy, right.load_energy),
            ("W_recovery", left.load_ergotropy, right.load_ergotropy),
            (
                "usable_fraction_recovery",
                left.usable_fraction,
                right.usable_fraction,
            ),
        ] {
            let tie = left_value.is_finite()
                && right_value.is_finite()
                && equivalent(left_value, right_value);
            let mut pair = [
                ("protect_entrance", left_value),
                ("protect_exit", right_value),
            ];
            pair.sort_by(|a, b| b.1.total_cmp(&a.1));
            for (rank, (condition, value)) in pair.iter().enumerate() {
                writeln!(
                    writer,
                    "time_series,{metric},{},{},{},{},{}",
                    rank + 1,
                    condition,
                    n(*value),
                    tie,
                    n(entrance.rows[index].time)
                )?;
            }
        }
    }
    Ok(())
}

fn write_checks(checks: &[Check]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("ideal_partial_protection_checks.csv")?);
    writeln!(writer, "scope,condition,check,value,tolerance,result")?;
    for row in checks {
        writeln!(
            writer,
            "{},{},{},{},{},{}",
            row.scope,
            row.condition,
            row.check,
            n(row.value),
            row.tolerance,
            pass(row.pass)
        )?;
    }
    Ok(())
}

fn convergence_metrics(summary: &Summary) -> [(&'static str, f64); 6] {
    [
        ("E_at_t10", summary.e_at_t10),
        ("W_at_t10", summary.w_at_t10),
        ("usable_fraction_at_t10", summary.usable_fraction_at_t10),
        ("W_max", summary.w_max),
        ("E_time_area", summary.e_time_area),
        ("W_time_area", summary.w_time_area),
    ]
}

fn write_convergence(
    coarse: &[Analysis],
    fine: &[Analysis],
    coarse_synergy: &[SynergyRow],
    fine_synergy: &[SynergyRow],
) -> std::io::Result<bool> {
    let mut writer = BufWriter::new(File::create("ideal_partial_protection_convergence.csv")?);
    writeln!(writer, "scope,condition,metric,base_dt,fine_dt,base_value,fine_value,absolute_difference,relative_difference,base_label,fine_label,result")?;
    let mut all_pass = true;
    for index in 0..coarse.len() {
        for ((name, base), (_, refined)) in convergence_metrics(&coarse[index].summary)
            .into_iter()
            .zip(convergence_metrics(&fine[index].summary))
        {
            let result = converged(base, refined);
            all_pass &= result;
            writeln!(
                writer,
                "condition,{},{},{},{},{},{},{},{},,,{}",
                coarse[index].condition.name,
                name,
                n(BASE_DT),
                n(FINE_DT),
                n(base),
                n(refined),
                n((base - refined).abs()),
                n(relative_difference(base, refined)),
                pass(result)
            )?;
        }
    }
    for metric in ["E", "W", "usable_fraction"] {
        let base = coarse_synergy
            .iter()
            .find(|row| row.metric == metric && row.evaluation == "t10")
            .unwrap();
        let refined = fine_synergy
            .iter()
            .find(|row| row.metric == metric && row.evaluation == "t10")
            .unwrap();
        let result = converged(base.synergy, refined.synergy)
            && (base.synergy == 0.0
                || refined.synergy == 0.0
                || base.synergy.signum() == refined.synergy.signum());
        all_pass &= result;
        writeln!(
            writer,
            "synergy,all,{metric},{},{},{},{},{},{},{},{},{}",
            n(BASE_DT),
            n(FINE_DT),
            n(base.synergy),
            n(refined.synergy),
            n((base.synergy - refined.synergy).abs()),
            n(relative_difference(base.synergy, refined.synergy)),
            if base.synergy > 0.0 {
                "positive"
            } else if base.synergy < 0.0 {
                "negative"
            } else {
                "zero"
            },
            if refined.synergy > 0.0 {
                "positive"
            } else if refined.synergy < 0.0 {
                "negative"
            } else {
                "zero"
            },
            pass(result)
        )?;
    }
    let ranking_metrics: [(&str, fn(&Summary) -> f64); 5] = [
        ("E_at_t10", |s| s.e_at_t10),
        ("W_at_t10", |s| s.w_at_t10),
        ("usable_fraction_at_t10", |s| s.usable_fraction_at_t10),
        ("E_time_area", |s| s.e_time_area),
        ("W_time_area", |s| s.w_time_area),
    ];
    for (metric, value) in ranking_metrics {
        let base_label = ranking_string(coarse, value);
        let fine_label = ranking_string(fine, value);
        let result = base_label == fine_label;
        all_pass &= result;
        writeln!(
            writer,
            "ranking,protection_conditions,{metric},{},{},NaN,NaN,NaN,NaN,{base_label},{fine_label},{}",
            n(BASE_DT),
            n(FINE_DT),
            pass(result)
        )?;
    }
    Ok(all_pass)
}

fn format_switches(switches: &[(f64, String, String)]) -> String {
    if switches.is_empty() {
        "confirmed switchなし".to_string()
    } else {
        switches
            .iter()
            .map(|(time, from, to)| format!("t={time:.2}: {from}->{to}"))
            .collect::<Vec<_>>()
            .join("; ")
    }
}

fn write_report(
    coarse: &[Analysis],
    recovery: &[RecoveryRow],
    synergy: &[SynergyRow],
    windows: &[WindowResult],
    checks: &[Check],
) -> std::io::Result<()> {
    let entrance = analysis(coarse, "protect_entrance");
    let exit = analysis(coarse, "protect_exit");
    let times: Vec<_> = entrance.rows.iter().map(|row| row.time).collect();
    let e_labels: Vec<_> = entrance
        .rows
        .iter()
        .zip(&exit.rows)
        .map(|(a, b)| pair_label(a.load_energy, b.load_energy))
        .collect();
    let w_labels: Vec<_> = entrance
        .rows
        .iter()
        .zip(&exit.rows)
        .map(|(a, b)| pair_label(a.load_ergotropy, b.load_ergotropy))
        .collect();
    let use_labels: Vec<_> = entrance
        .rows
        .iter()
        .zip(&exit.rows)
        .map(|(a, b)| pair_label(a.usable_fraction, b.usable_fraction))
        .collect();
    let e_switches = confirmed_switches(&times, &e_labels);
    let w_switches = confirmed_switches(&times, &w_labels);
    let use_switches = confirmed_switches(&times, &use_labels);
    let mut writer = BufWriter::new(File::create("MILESTONE_7C_REPORT.md")?);
    writeln!(
        writer,
        "# Milestone 7c: Ideal partial protection upper-bound test\n"
    )?;
    writeln!(writer, "## 1. 目的\n\n全3site雑音から指定siteの位相雑音演算子だけを理想的に除去した反実仮想条件で、load状態量の回復上限を比較した。\n")?;
    writeln!(writer, "## 2. Milestone 7aとの問いの違い\n\n7aは1siteだけへ雑音を置く有害配置比較。7cは全site雑音から雑音項を選択的に除く回復比較であり、順位一致を仮定しない。\n")?;
    writeln!(writer, "## 3. 理想保護の定義\n\n指定siteの `sqrt(gamma_phi/2) sigma_z` collapse operatorを完全に除去する。現実的装置、制御pulse、cost、誤り訂正、有限精度保護ではない。\n")?;
    writeln!(writer, "## 4. 比較条件\n\n`all_noisy=[0,1,2]`, `protect_entrance=[1,2]`, `protect_exit=[0,1]`, `protect_both_ends=[1]`, `noise_free=[]`。protect_both_endsでは中央雑音だけが残る。\n")?;
    writeln!(writer, "## 5. 変更していない物理条件\n\n3site+3準位load、Hamiltonian、真空初期状態、J=1、g=0.25、各周波数=1、tau=3.2、Omega=0.2、gamma_phi=0.5、t_max=10、pulse、load、ergotropy、RK4を固定した。\n")?;
    writeln!(writer, "## 6. 数値手法\n\n基準dt=0.0025、半減dt=0.00125、保存間隔0.01。状態量時間面積とpowerは台形則。比のepsilon={SIGNAL_TOLERANCE:e}。非加算性の約0判定はabs<={SYNERGY_ABSOLUTE_TOL:e}またはnormalized abs<={SYNERGY_RELATIVE_TOL:e}。\n")?;
    writeln!(writer, "## 7. 数値品質チェック\n\n全{}項目PASS。collapse数は3/2/2/1/0、site mapping、trace、Hermiticity、positivity、population、load縮約、top-level、ledger、有限性、共通grid/初期状態、W<=E、usable範囲、drive整合、既存baselineを確認した。\n", checks.len())?;
    writeln!(writer, "## 8. t=10の5条件比較\n\n| condition | E | W | usable fraction | W/Ein |\n|---|---:|---:|---:|---:|")?;
    for row in coarse {
        let s = &row.summary;
        writeln!(
            writer,
            "| {} | {:.10e} | {:.10e} | {:.10e} | {:.10e} |",
            row.condition.name,
            s.e_at_t10,
            s.w_at_t10,
            s.usable_fraction_at_t10,
            s.w_over_ein_at_t10
        )?;
    }
    writeln!(writer, "\n## 9. all_noisyからの回復量\n\n| condition | metric | point | absolute | normalized |\n|---|---|---|---:|---:|")?;
    for row in recovery {
        writeln!(
            writer,
            "| {} | {} | {} | {:.8e} | {:.8} |",
            row.condition, row.metric, row.evaluation, row.absolute, row.normalized
        )?;
    }
    writeln!(writer, "\n## 10. noise-freeまでの残留損失\n\n| condition | metric | t=10 residual loss |\n|---|---|---:|")?;
    for row in recovery.iter().filter(|row| row.evaluation == "t10") {
        writeln!(
            writer,
            "| {} | {} | {:.8} |",
            row.condition, row.metric, row.residual
        )?;
    }
    writeln!(
        writer,
        "\n符号付きで保存した。負値はnoise-free超過であり0へ丸めていない。\n"
    )?;
    writeln!(
        writer,
        "## 11. W_max比較\n\nW_max回復最大: `{}`。\n",
        ranking_string(coarse, |s| s.w_max)
    )?;
    writeln!(writer, "## 12. E_time_areaとW_time_area比較\n\nE_time_area回復最大: `{}`。W_time_area回復最大: `{}`。両者は状態量の時間面積で、累積流入/抽出ではない。\n", ranking_string(coarse,|s|s.e_time_area), ranking_string(coarse,|s|s.w_time_area))?;
    writeln!(writer, "## 13. 時間窓別回復\n")?;
    for window in WINDOWS {
        let best = windows
            .iter()
            .filter(|r| {
                PROTECTION_INDICES
                    .iter()
                    .any(|&i| coarse[i].condition.name == r.condition)
                    && r.window.name == window.name
            })
            .max_by(|a, b| a.delta_w_area.total_cmp(&b.delta_w_area))
            .unwrap();
        writeln!(
            writer,
            "- {}: W_time_area回復最大 `{}` ({:.8e})",
            window.name, best.condition, best.delta_w_area
        )?;
    }
    writeln!(writer, "\n## 14. 入口保護と出口保護の順位変化\n\n- E: {}\n- W: {}\n- usable fraction: {}\n\n5点持続した変化だけを確定した。\n",format_switches(&e_switches),format_switches(&w_switches),format_switches(&use_switches))?;
    let both = analysis(coarse, "protect_both_ends");
    let free = analysis(coarse, "noise_free");
    let all = analysis(coarse, "all_noisy");
    writeln!(writer, "## 15. 両端保護の結果\n\n中央雑音だけが残る。t=10 W={:.10e}, usable={:.10e}。noise-freeへのW残留損失={:.8}、all-noisyからのW回復率={:.8}。\n\n| t=10 metric | both - entrance | both - exit |\n|---|---:|---:|\n| E | {:.8e} | {:.8e} |\n| W | {:.8e} | {:.8e} |\n| usable fraction | {:.8e} | {:.8e} |\n\n中央雑音が無害/保護不要とは言わない。\n",both.summary.w_at_t10,both.summary.usable_fraction_at_t10,1.0-both.summary.w_at_t10/free.summary.w_at_t10,(both.summary.w_at_t10-all.summary.w_at_t10)/(free.summary.w_at_t10-all.summary.w_at_t10),both.summary.e_at_t10-entrance.summary.e_at_t10,both.summary.e_at_t10-exit.summary.e_at_t10,both.summary.w_at_t10-entrance.summary.w_at_t10,both.summary.w_at_t10-exit.summary.w_at_t10,both.summary.usable_fraction_at_t10-entrance.summary.usable_fraction_at_t10,both.summary.usable_fraction_at_t10-exit.summary.usable_fraction_at_t10)?;
    writeln!(writer, "## 16. 回復の非加算性\n\n| metric | point | synergy | normalized | class |\n|---|---|---:|---:|---|")?;
    for row in synergy {
        writeln!(
            writer,
            "| {} | {} | {:.8e} | {:.8} | {} |",
            row.metric, row.evaluation, row.synergy, row.normalized, row.classification
        )?;
    }
    writeln!(
        writer,
        "\nSynergyは相互作用エネルギーではなく、選択的雑音除去に対する観測量応答の非加算性。\n"
    )?;
    writeln!(writer, "## 17. 刻み幅整合性\n\n5条件×6量、t=10 synergy E/W/use、5種類の保護順位がPASS。順位とsynergy符号は不変。基準/半減刻みのラベルもconvergence CSVへ保存した。\n")?;
    writeln!(writer, "## 18. 直接確認できたこと\n\n固定条件で全3site雑音から入口/出口/両端雑音を理想除去したときのt=10、最大値、時間面積、時間窓、非加算性、中央雑音残存時の残留損失を確認した。\n")?;
    writeln!(writer, "## 19. 確認できていないこと\n\n現実的実装、cost、不完全保護、他gamma/Omega/network、長時間/反復/定常、因果機構、新規性は未確認。\n")?;
    writeln!(writer, "## 20. 主張してはいけないこと\n\n保護の一般則・実用最適、synergyを物理相互作用エネルギーとする表現、中央雑音が無害、量子優位、現実送電性能は主張しない。\n")?;
    writeln!(writer, "## 21. load雑音を今回扱わない理由\n\nloadは3準位でsiteの二準位sigma_zを流用できず、同じgammaが公平な強度を意味しない。演算子、decay rate、規格化、強度対応を別途検証する必要がある。\n")?;
    writeln!(writer, "## 22. 次段階への判断材料\n\nt=10 W回復最大 `{}`、usable fraction回復最大 `{}`、E回復最大 `{}`。これは次の比較候補であり、自動的に次Milestoneへ進まない。\n",ranking_string(coarse,|s|s.w_at_t10),ranking_string(coarse,|s|s.usable_fraction_at_t10),ranking_string(coarse,|s|s.e_at_t10))?;
    writeln!(writer, "## 23. 生成ファイル一覧\n\n- `src/bin/ideal_partial_protection.rs`\n- `ideal_partial_protection_timeseries.csv`\n- `ideal_partial_protection_summary.csv`\n- `ideal_partial_protection_recovery.csv`\n- `ideal_partial_protection_synergy.csv`\n- `ideal_partial_protection_windows.csv`\n- `ideal_partial_protection_rankings.csv`\n- `ideal_partial_protection_checks.csv`\n- `ideal_partial_protection_convergence.csv`\n- `MILESTONE_7C_REPORT.md`\n")?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams::default();
    let operators = build_operators(&params)?;
    let mut coarse = Vec::new();
    for condition in CONDITIONS {
        coarse.push(run_condition(condition, BASE_DT, &params, &operators)?);
    }
    let checks = checks(&coarse, &operators)?;
    write_checks(&checks)?;
    if checks.iter().any(|check| !check.pass) {
        return Err("Milestone 7c numerical/input checks failed".into());
    }
    let recovery = recovery_rows(&coarse)?;
    let synergy = synergy_rows(&coarse)?;
    let window_rows = windows(&coarse);

    let mut fine = Vec::new();
    for condition in CONDITIONS {
        fine.push(run_condition(condition, FINE_DT, &params, &operators)?);
    }
    if fine.iter().any(|analysis| !analysis.quality.all_pass()) {
        return Err("fine-step physical diagnostics failed".into());
    }
    let fine_synergy = synergy_rows(&fine)?;
    let convergence_ok = write_convergence(&coarse, &fine, &synergy, &fine_synergy)?;
    if !convergence_ok {
        return Err("time-step convergence or synergy sign check failed".into());
    }
    for value in [
        |s: &Summary| s.e_at_t10,
        |s: &Summary| s.w_at_t10,
        |s: &Summary| s.usable_fraction_at_t10,
        |s: &Summary| s.e_time_area,
        |s: &Summary| s.w_time_area,
    ] {
        if ranking_string(&coarse, value) != ranking_string(&fine, value) {
            return Err("protection ranking changed after halving dt".into());
        }
    }
    write_timeseries(&coarse)?;
    write_summary(&coarse)?;
    write_recovery(&recovery)?;
    write_synergy(&synergy)?;
    write_windows(&window_rows)?;
    write_rankings(&coarse)?;
    write_report(&coarse, &recovery, &synergy, &window_rows, &checks)?;
    println!("Milestone 7c complete");
    Ok(())
}
