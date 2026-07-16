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
const RELATIVE_TOL: f64 = CONVERGENCE_RELATIVE_TOLERANCE;
const BASELINE_TOLERANCE: f64 = 1.0e-9;
const REDUCED_STATE_TOLERANCE: f64 = 1.0e-9;
const POPULATION_TOLERANCE: f64 = 1.0e-8;

const BASELINE_E: f64 = 5.4450767877898487e-2;
const BASELINE_W: f64 = 5.2798274942446315e-2;
const BASELINE_EIN: f64 = 7.5315933437092850e-2;
const BASELINE_COHERENCE_L1: f64 = 4.6453837009096421e-1;

#[derive(Clone, Copy)]
struct Condition {
    name: &'static str,
    noise_site: &'static str,
    gamma_phi: f64,
    sites: &'static [usize],
}

const CONDITIONS: [Condition; 4] = [
    Condition {
        name: "noise_free",
        noise_site: "none",
        gamma_phi: 0.0,
        sites: &[],
    },
    Condition {
        name: "noise_entrance",
        noise_site: "site1",
        gamma_phi: 0.5,
        sites: &[0],
    },
    Condition {
        name: "noise_middle",
        noise_site: "site2",
        gamma_phi: 0.5,
        sites: &[1],
    },
    Condition {
        name: "noise_exit",
        noise_site: "site3",
        gamma_phi: 0.5,
        sites: &[2],
    },
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
    site_populations: [f64; 3],
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
    w_time_area: f64,
    e_time_area: f64,
    w_time_mean: f64,
    e_time_mean: f64,
    max_load_top_population: f64,
    max_trace_error: f64,
    max_hermiticity_error: f64,
    minimum_density_eigenvalue: f64,
    max_abs_energy_ledger_residual: f64,
}

#[derive(Clone)]
struct Quality {
    trace: bool,
    hermiticity: bool,
    positivity: bool,
    reduced_state: bool,
    population_bounds: bool,
    top_level: bool,
    ledger: bool,
    finite: bool,
    max_reduced_difference: f64,
}

impl Quality {
    fn all_pass(&self) -> bool {
        self.trace
            && self.hermiticity
            && self.positivity
            && self.reduced_state
            && self.population_bounds
            && self.top_level
            && self.ledger
            && self.finite
    }
}

struct Analysis {
    condition: Condition,
    run: CoherentDriveRun,
    rows: Vec<Row>,
    summary: Summary,
    quality: Quality,
}

#[derive(Clone)]
struct RatioRow {
    condition: &'static str,
    noise_site: &'static str,
    r_e: f64,
    r_w: f64,
    r_use: f64,
    r_wmax: f64,
    r_warea: f64,
    r_wmean: f64,
    classification: &'static str,
}

fn n(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else {
        format!("{value:.16e}")
    }
}

fn b(value: bool) -> &'static str {
    if value {
        "PASS"
    } else {
        "FAIL"
    }
}

fn safe_ratio(numerator: f64, denominator: f64) -> f64 {
    if denominator.abs() <= SIGNAL_TOLERANCE {
        f64::NAN
    } else {
        numerator / denominator
    }
}

fn trapezoid(rows: &[Row], value: impl Fn(&Row) -> f64) -> f64 {
    rows.windows(2)
        .map(|window| {
            0.5 * (value(&window[0]) + value(&window[1])) * (window[1].time - window[0].time)
        })
        .sum()
}

fn local_load_hamiltonian(params: &ModelParams) -> ComplexMatrix {
    ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
        params.load_dim,
        (0..params.load_dim).map(|level| C64::new(level as f64 * params.omega_load, 0.0)),
    ))
}

fn analyze(
    condition: Condition,
    run: CoherentDriveRun,
    params: &ModelParams,
    operators: &Operators,
) -> Result<Analysis, Box<dyn std::error::Error>> {
    if run.samples.len() != run.states.len() {
        return Err("saved sample/state grid length mismatch".into());
    }
    let load_h = local_load_hamiltonian(params);
    let initial_energy = run.samples[0].bare_network_energy;
    let mut drive_powers = Vec::with_capacity(run.samples.len());
    let mut dephasing_powers = Vec::with_capacity(run.samples.len());
    let mut rows = Vec::with_capacity(run.samples.len());
    let mut max_reduced_difference = 0.0_f64;
    let mut population_bounds = true;
    let mut ledger_ok = true;

    for (index, (sample, state)) in run.samples.iter().zip(&run.states).enumerate() {
        if sample.time != state.time {
            return Err("saved sample/state time mismatch".into());
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
        let ledger_scale = (sample.bare_network_energy - initial_energy)
            .abs()
            .max(drive_net.abs())
            .max(dephasing_net.abs());
        ledger_ok &= ledger_residual.abs()
            <= LEDGER_ABSOLUTE_TOLERANCE + LEDGER_RELATIVE_TOLERANCE * ledger_scale;

        let site_populations = [
            expectation(&state.rho, &operators.number_sites[0]).re,
            expectation(&state.rho, &operators.number_sites[1]).re,
            expectation(&state.rho, &operators.number_sites[2]).re,
        ];
        population_bounds &= site_populations
            .iter()
            .chain(sample.load_populations.iter())
            .all(|value| *value >= -POPULATION_TOLERANCE && *value <= 1.0 + POPULATION_TOLERANCE);

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

        rows.push(Row {
            time: sample.time,
            load_energy: sample.load_energy,
            load_ergotropy: sample.load_ergotropy,
            diagonal_ergotropy: sample.load_diagonal_ergotropy,
            coherence_ergotropy: sample.load_coherence_ergotropy,
            coherence_l1: sample.load_coherence_l1,
            usable_fraction: safe_ratio(sample.load_ergotropy, sample.load_energy),
            drive_energy_in: drive_in,
            w_over_ein: safe_ratio(sample.load_ergotropy, drive_in),
            site_populations,
            load_top_population: sample.load_populations[2],
            trace_error: sample.trace_error,
            hermiticity_error: sample.hermiticity_error,
            minimum_eigenvalue: sample.minimum_eigenvalue,
            ledger_residual,
        });
    }

    let final_row = rows.last().expect("run has final sample");
    let (w_max_index, _) = rows
        .iter()
        .enumerate()
        .max_by(|(_, left), (_, right)| left.load_ergotropy.total_cmp(&right.load_ergotropy))
        .expect("run has samples");
    let w_max_row = &rows[w_max_index];
    let w_time_area = trapezoid(&rows, |row| row.load_ergotropy);
    let e_time_area = trapezoid(&rows, |row| row.load_energy);
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
        w_time_area,
        e_time_area,
        w_time_mean: w_time_area / 10.0,
        e_time_mean: e_time_area / 10.0,
        max_load_top_population: rows
            .iter()
            .map(|row| row.load_top_population)
            .fold(0.0, f64::max),
        max_trace_error: run.summary.maximum_trace_error,
        max_hermiticity_error: run.summary.maximum_hermiticity_error,
        minimum_density_eigenvalue: run.summary.worst_minimum_eigenvalue,
        max_abs_energy_ledger_residual: rows
            .iter()
            .map(|row| row.ledger_residual.abs())
            .fold(0.0, f64::max),
    };
    let quality = Quality {
        trace: summary.max_trace_error <= TRACE_TOLERANCE,
        hermiticity: summary.max_hermiticity_error <= HERMITICITY_TOLERANCE,
        positivity: summary.minimum_density_eigenvalue >= -POSITIVITY_TOLERANCE,
        reduced_state: max_reduced_difference <= REDUCED_STATE_TOLERANCE,
        population_bounds,
        top_level: summary.max_load_top_population < TOP_LEVEL_LIMIT,
        ledger: ledger_ok,
        finite: run.summary.all_finite
            && rows.iter().all(|row| {
                [
                    row.time,
                    row.load_energy,
                    row.load_ergotropy,
                    row.diagonal_ergotropy,
                    row.coherence_ergotropy,
                    row.coherence_l1,
                    row.drive_energy_in,
                    row.trace_error,
                    row.hermiticity_error,
                    row.minimum_eigenvalue,
                    row.ledger_residual,
                ]
                .iter()
                .all(|value| value.is_finite())
                    && row.site_populations.iter().all(|value| value.is_finite())
            }),
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

fn config(gamma_phi: f64, dt: f64) -> CoherentDriveConfig {
    let mut config = CoherentDriveConfig::milestone_5b(gamma_phi, dt);
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
    let run = run_coherent_drive_with_noise_sites(
        params,
        config(condition.gamma_phi, dt),
        condition.sites,
    )?;
    analyze(condition, run, params, operators)
}

fn classification(r_e: f64, r_use: f64) -> &'static str {
    if r_e < 1.0 - RELATIVE_TOL && r_use >= 1.0 - RELATIVE_TOL {
        "transport_loss"
    } else if (r_e - 1.0).abs() <= RELATIVE_TOL && r_use < 1.0 - RELATIVE_TOL {
        "quality_loss"
    } else if r_e < 1.0 - RELATIVE_TOL && r_use < 1.0 - RELATIVE_TOL {
        "transport_and_quality_loss"
    } else {
        "no_clear_loss"
    }
}

fn ratios(analyses: &[Analysis]) -> Vec<RatioRow> {
    let reference = &analyses[0].summary;
    analyses[1..]
        .iter()
        .map(|analysis| {
            let summary = &analysis.summary;
            let r_e = safe_ratio(summary.e_at_t10, reference.e_at_t10);
            let r_w = safe_ratio(summary.w_at_t10, reference.w_at_t10);
            let r_use = safe_ratio(
                summary.usable_fraction_at_t10,
                reference.usable_fraction_at_t10,
            );
            RatioRow {
                condition: analysis.condition.name,
                noise_site: analysis.condition.noise_site,
                r_e,
                r_w,
                r_use,
                r_wmax: safe_ratio(summary.w_max, reference.w_max),
                r_warea: safe_ratio(summary.w_time_area, reference.w_time_area),
                r_wmean: safe_ratio(summary.w_time_mean, reference.w_time_mean),
                classification: classification(r_e, r_use),
            }
        })
        .collect()
}

fn equivalent(left: f64, right: f64) -> bool {
    (left - right).abs() <= CONVERGENCE_ABSOLUTE_TOLERANCE
        || (left - right).abs() / left.abs().max(right.abs()).max(SIGNAL_TOLERANCE)
            <= CONVERGENCE_RELATIVE_TOLERANCE
}

fn ranking(analyses: &[Analysis], value: impl Fn(&Summary) -> f64) -> String {
    let minimum = analyses[1..]
        .iter()
        .map(|analysis| value(&analysis.summary))
        .min_by(f64::total_cmp)
        .expect("three noisy conditions");
    let tied: Vec<_> = analyses[1..]
        .iter()
        .filter(|analysis| equivalent(value(&analysis.summary), minimum))
        .map(|analysis| analysis.condition.noise_site)
        .collect();
    if tied.len() == 1 {
        tied[0].to_string()
    } else {
        format!("tie({})", tied.join("+"))
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

fn convergence_metrics(summary: &Summary) -> [(&'static str, f64); 5] {
    [
        ("E_at_t10", summary.e_at_t10),
        ("W_at_t10", summary.w_at_t10),
        ("usable_fraction_at_t10", summary.usable_fraction_at_t10),
        ("W_max", summary.w_max),
        ("W_time_area", summary.w_time_area),
    ]
}

fn write_timeseries(analyses: &[Analysis]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("local_noise_placement_timeseries.csv")?);
    writeln!(writer, "condition,noise_site,gamma_phi,Omega,time,load_energy,load_ergotropy,load_diagonal_ergotropy,load_coherence_ergotropy,load_coherence_l1,usable_fraction,drive_energy_in,W_over_Ein,site1_population,site2_population,site3_population,load_top_level_population,trace_error,hermiticity_error,min_eigenvalue,energy_ledger_residual")?;
    for analysis in analyses {
        for row in &analysis.rows {
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                analysis.condition.name,
                analysis.condition.noise_site,
                n(analysis.condition.gamma_phi),
                n(0.2),
                n(row.time),
                n(row.load_energy),
                n(row.load_ergotropy),
                n(row.diagonal_ergotropy),
                n(row.coherence_ergotropy),
                n(row.coherence_l1),
                n(row.usable_fraction),
                n(row.drive_energy_in),
                n(row.w_over_ein),
                n(row.site_populations[0]),
                n(row.site_populations[1]),
                n(row.site_populations[2]),
                n(row.load_top_population),
                n(row.trace_error),
                n(row.hermiticity_error),
                n(row.minimum_eigenvalue),
                n(row.ledger_residual)
            )?;
        }
    }
    Ok(())
}

fn write_summary(analyses: &[Analysis]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("local_noise_placement_summary.csv")?);
    writeln!(writer, "condition,noise_site,gamma_phi,Omega,E_at_t10,W_at_t10,diagonal_W_at_t10,coherence_W_at_t10,coherence_L1_at_t10,usable_fraction_at_t10,drive_energy_in_at_t10,W_over_Ein_at_t10,W_max,t_at_W_max,E_at_W_max,usable_fraction_at_W_max,W_time_area,E_time_area,W_time_mean,E_time_mean,max_load_top_level_population,max_trace_error,max_hermiticity_error,minimum_density_eigenvalue,max_abs_energy_ledger_residual,all_checks")?;
    for analysis in analyses {
        let x = &analysis.summary;
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            analysis.condition.name,
            analysis.condition.noise_site,
            n(analysis.condition.gamma_phi),
            n(0.2),
            n(x.e_at_t10),
            n(x.w_at_t10),
            n(x.diagonal_w_at_t10),
            n(x.coherence_w_at_t10),
            n(x.coherence_l1_at_t10),
            n(x.usable_fraction_at_t10),
            n(x.drive_energy_in_at_t10),
            n(x.w_over_ein_at_t10),
            n(x.w_max),
            n(x.t_at_w_max),
            n(x.e_at_w_max),
            n(x.usable_fraction_at_w_max),
            n(x.w_time_area),
            n(x.e_time_area),
            n(x.w_time_mean),
            n(x.e_time_mean),
            n(x.max_load_top_population),
            n(x.max_trace_error),
            n(x.max_hermiticity_error),
            n(x.minimum_density_eigenvalue),
            n(x.max_abs_energy_ledger_residual),
            b(analysis.quality.all_pass())
        )?;
    }
    Ok(())
}

fn write_ratios(rows: &[RatioRow]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("local_noise_placement_ratios.csv")?);
    writeln!(
        writer,
        "condition,noise_site,R_E,R_W,R_use,R_Wmax,R_Warea,R_Wmean,relative_tol,classification"
    )?;
    for row in rows {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{}",
            row.condition,
            row.noise_site,
            n(row.r_e),
            n(row.r_w),
            n(row.r_use),
            n(row.r_wmax),
            n(row.r_warea),
            n(row.r_wmean),
            n(RELATIVE_TOL),
            row.classification
        )?;
    }
    Ok(())
}

fn write_check(
    writer: &mut impl Write,
    scope: &str,
    condition: &str,
    check: &str,
    value: f64,
    tolerance: &str,
    pass: bool,
) -> std::io::Result<()> {
    writeln!(
        writer,
        "{scope},{condition},{check},{},{tolerance},{}",
        n(value),
        b(pass)
    )
}

fn write_checks(
    analyses: &[Analysis],
    same_grid: bool,
    same_initial: bool,
    baseline: &[(&str, f64, f64, f64, bool)],
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("local_noise_placement_checks.csv")?);
    writeln!(writer, "scope,condition,check,value,tolerance,result")?;
    for analysis in analyses {
        let x = &analysis.summary;
        let q = &analysis.quality;
        write_check(
            &mut writer,
            "condition",
            analysis.condition.name,
            "max_trace_error",
            x.max_trace_error,
            "<=1e-8",
            q.trace,
        )?;
        write_check(
            &mut writer,
            "condition",
            analysis.condition.name,
            "max_hermiticity_error",
            x.max_hermiticity_error,
            "<=1e-8",
            q.hermiticity,
        )?;
        write_check(
            &mut writer,
            "condition",
            analysis.condition.name,
            "minimum_density_eigenvalue",
            x.minimum_density_eigenvalue,
            ">=-1e-8",
            q.positivity,
        )?;
        write_check(
            &mut writer,
            "condition",
            analysis.condition.name,
            "load_reduced_state_max_difference",
            q.max_reduced_difference,
            "<=1e-9",
            q.reduced_state,
        )?;
        write_check(
            &mut writer,
            "condition",
            analysis.condition.name,
            "population_bounds",
            0.0,
            "-1e-8<=p<=1+1e-8",
            q.population_bounds,
        )?;
        write_check(
            &mut writer,
            "condition",
            analysis.condition.name,
            "max_load_top_level_population",
            x.max_load_top_population,
            "<0.05",
            q.top_level,
        )?;
        write_check(
            &mut writer,
            "condition",
            analysis.condition.name,
            "max_abs_energy_ledger_residual",
            x.max_abs_energy_ledger_residual,
            "<=5e-5+5e-4*scale",
            q.ledger,
        )?;
        write_check(
            &mut writer,
            "condition",
            analysis.condition.name,
            "unexpected_nonfinite",
            0.0,
            "none",
            q.finite,
        )?;
        write_check(
            &mut writer,
            "condition",
            analysis.condition.name,
            "all_checks",
            0.0,
            "all above",
            q.all_pass(),
        )?;
    }
    write_check(
        &mut writer,
        "global",
        "all",
        "same_time_grid",
        0.0,
        "exact",
        same_grid,
    )?;
    write_check(
        &mut writer,
        "global",
        "all",
        "same_initial_state",
        0.0,
        "Frobenius<=1e-14",
        same_initial,
    )?;
    for (name, _actual, _expected, difference, pass) in baseline {
        write_check(
            &mut writer,
            "baseline",
            "noise_free",
            name,
            *difference,
            "absolute<=1e-9",
            *pass,
        )?;
    }
    Ok(())
}

fn write_convergence(pairs: &[(&Analysis, &Analysis)]) -> Result<bool, std::io::Error> {
    let mut writer = BufWriter::new(File::create("local_noise_placement_convergence.csv")?);
    writeln!(writer, "condition,noise_site,metric,base_dt,fine_dt,base_value,fine_value,absolute_difference,relative_difference,result")?;
    let mut all_pass = true;
    for (coarse, fine) in pairs {
        for ((name_c, value_c), (name_f, value_f)) in convergence_metrics(&coarse.summary)
            .into_iter()
            .zip(convergence_metrics(&fine.summary))
        {
            assert_eq!(name_c, name_f);
            let pass = converged(value_c, value_f);
            all_pass &= pass;
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{},{}",
                coarse.condition.name,
                coarse.condition.noise_site,
                name_c,
                n(BASE_DT),
                n(FINE_DT),
                n(value_c),
                n(value_f),
                n((value_c - value_f).abs()),
                n(relative_difference(value_c, value_f)),
                b(pass)
            )?;
        }
    }
    Ok(all_pass)
}

fn report_table(analyses: &[Analysis]) -> String {
    let mut text = String::from("| condition | noise site | E(t=10) | W(t=10) | usable fraction | W_max | W_time_area |\n|---|---|---:|---:|---:|---:|---:|\n");
    for analysis in analyses {
        let x = &analysis.summary;
        text.push_str(&format!(
            "| {} | {} | {:.10e} | {:.10e} | {:.10e} | {:.10e} | {:.10e} |\n",
            analysis.condition.name,
            analysis.condition.noise_site,
            x.e_at_t10,
            x.w_at_t10,
            x.usable_fraction_at_t10,
            x.w_max,
            x.w_time_area
        ));
    }
    text
}

fn write_report(
    analyses: &[Analysis],
    ratio_rows: &[RatioRow],
    convergence_pairs: &[(&Analysis, &Analysis)],
    rankings: &[(&str, String)],
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create("MILESTONE_7A_REPORT.md")?);
    writeln!(writer, "# Milestone 7a: 雑音位置の比較\n")?;
    writeln!(writer, "## 1. 目的\n\n固定した有限3サイト模型で、位相雑音を1サイトだけに置き、最後のloadへの有限時間内の影響を比較した。原因機構の断定は行わない。\n")?;
    writeln!(writer, "## 2. 既存模型から変更した点\n\n位相雑音を入れる0始まりsite集合を指定できる入口を追加した。既存の全3サイト雑音APIは全site `[0,1,2]` を渡すラッパーとして保持した。Milestone 7a専用binで比較・CSV・診断を実装した。\n")?;
    writeln!(writer, "## 3. 変更していない条件\n\nHamiltonian、3つの二準位site、3準位load、全系真空、J=1、g=0.25、各角周波数=1、tau=3.2、Omega=0.2、t_max=10、pulse、load定義、ergotropy、RK4を変更していない。\n")?;
    writeln!(writer, "## 4. 雑音位置の定義\n\nsite 1/2/3は内部index 0/1/2で、それぞれ入口・中央・出口側と呼ぶ。雑音あり条件は `sqrt(gamma_phi/2) sigma_z` を指定した1サイトにだけ置き、gamma_phi=0.5とした。noise_freeはcollapse operatorを作らない。\n")?;
    writeln!(writer, "## 5. Milestone 5cとの違い\n\nMilestone 5cの雑音ありBは3サイトすべてに位相雑音があった。今回は1サイトだけであり、数値を直接比較しない。\n")?;
    writeln!(writer, "## 6. 数値手法\n\nMilestone 5c本計算と同じ基準刻み `dt={BASE_DT}`、保存間隔 `0.01`、固定刻みRK4を使用した。半減確認は `dt={FINE_DT}`。積分は共通保存グリッド上の台形則。比の分母許容値は既存 `SIGNAL_TOLERANCE={SIGNAL_TOLERANCE:e}`。分類の `relative_tol={RELATIVE_TOL:e}` は既存の収束相対許容値を採用した。\n")?;
    writeln!(writer, "## 7. 数値品質チェック\n\n4条件のtrace、Hermiticity、positivity、load縮約整合、population bounds、top-level、energy ledger、予期しない非有限値、共通時間グリッド、共通初期状態はすべてPASS。ledger基準は `|r| <= {LEDGER_ABSOLUTE_TOLERANCE:e} + {LEDGER_RELATIVE_TOLERANCE:e}*scale`。詳細は `local_noise_placement_checks.csv`。\n")?;
    writeln!(writer, "## 8. t=10の比較\n\n{}", report_table(analyses))?;
    writeln!(writer, "## 9. 最大ergotropyの比較\n\n各条件の `W_max`, `t_at_W_max`, `E_at_W_max` はsummary CSVに保存した。最小順位は `{}`。\n", rankings[2].1)?;
    writeln!(writer, "## 10. W_time_areaとW_time_meanの比較\n\n`W_time_area`は0〜10にergotropyがどの程度存在したかの補助値であり、累積仕事・総仕事・供給された仕事ではない。最小順位は `{}`。\n", rankings[3].1)?;
    writeln!(writer, "## 11. 雑音なし条件に対する比\n\n| condition | R_E | R_W | R_use | R_Wmax | R_Warea | R_Wmean |\n|---|---:|---:|---:|---:|---:|---:|")?;
    for row in ratio_rows {
        writeln!(
            writer,
            "| {} | {:.8} | {:.8} | {:.8} | {:.8} | {:.8} | {:.8} |",
            row.condition, row.r_e, row.r_w, row.r_use, row.r_wmax, row.r_warea, row.r_wmean
        )?;
    }
    writeln!(writer)?;
    writeln!(writer, "## 12. 壊れ方の便宜的分類\n")?;
    for row in ratio_rows {
        writeln!(writer, "- {}: `{}`", row.condition, row.classification)?;
    }
    writeln!(
        writer,
        "\nこれは今回だけの整理で、一般法則や物理定理ではない。\n"
    )?;
    writeln!(writer, "## 13. 雑音位置の順位\n")?;
    for (label, result) in rankings {
        writeln!(writer, "- {label}: `{result}`")?;
    }
    writeln!(writer)?;
    writeln!(writer, "## 14. 時間刻み半減の結果\n\n既存基準に従い、絶対差 `<= {CONVERGENCE_ABSOLUTE_TOLERANCE:e}` または相対差 `<= {CONVERGENCE_RELATIVE_TOLERANCE:e}` をPASSとした。\n\n| condition | metric | abs diff | rel diff |\n|---|---|---:|---:|")?;
    for (coarse, fine) in convergence_pairs {
        for ((name, base), (_, refined)) in convergence_metrics(&coarse.summary)
            .into_iter()
            .zip(convergence_metrics(&fine.summary))
        {
            writeln!(
                writer,
                "| {} | {} | {:.6e} | {:.6e} |",
                coarse.condition.name,
                name,
                (base - refined).abs(),
                relative_difference(base, refined)
            )?;
        }
    }
    writeln!(
        writer,
        "\nW_at_t10が最小の雑音位置は刻み半減後も `{}` で変わらず、全10比較量がPASSした。\n",
        rankings[0].1
    )?;
    writeln!(writer, "## 15. 直接確認できたこと\n\n固定した有限3サイト模型、gamma_phi=0.5、Omega=0.2、固定単発pulse、0<=t<=10で、1サイトだけに雑音を置いた場合のload energy・ergotropy等の位置別差を直接確認した。t=10のW最小は `{}`、usable fraction最小は `{}` だった。\n", rankings[0].1, rankings[1].1)?;
    writeln!(writer, "## 16. 確認できていないこと\n\n他のgamma_phi、Omega、pulse、長いnetwork、長時間、連続運転、放電、抽出操作、現実装置、古典比較、文献上の新規性は確認していない。\n")?;
    writeln!(writer, "## 17. 主張してはいけないこと\n\n雑音位置の普遍則、任意パラメータや長いnetworkへの一般化、量子優位、現実の送電・装置効率、雑音が常にergotropyを減らす一般論、各site雑音が注入・輸送・受け渡しだけを壊したという因果断定はできない。\n")?;
    writeln!(writer, "## 18. 作成ファイル一覧\n\n- `src/bin/local_noise_placement.rs`\n- `local_noise_placement_timeseries.csv`\n- `local_noise_placement_summary.csv`\n- `local_noise_placement_ratios.csv`\n- `local_noise_placement_checks.csv`\n- `local_noise_placement_convergence.csv`\n- `MILESTONE_7A_REPORT.md`\n")?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams::default();
    let operators = build_operators(&params)?;
    let mut analyses = Vec::with_capacity(CONDITIONS.len());
    for condition in CONDITIONS {
        analyses.push(run_condition(condition, BASE_DT, &params, &operators)?);
    }

    let same_grid = analyses[1..].iter().all(|analysis| {
        analysis.rows.len() == analyses[0].rows.len()
            && analysis
                .rows
                .iter()
                .zip(&analyses[0].rows)
                .all(|(left, right)| left.time == right.time)
    });
    let same_initial = analyses[1..].iter().all(|analysis| {
        frobenius_norm(&(&analysis.run.states[0].rho - &analyses[0].run.states[0].rho)) <= 1.0e-14
    });
    if !same_grid || !same_initial {
        return Err("conditions do not share the exact time grid and initial state".into());
    }
    if let Some(failed) = analyses
        .iter()
        .find(|analysis| !analysis.quality.all_pass())
    {
        return Err(format!("numerical diagnostics failed for {}", failed.condition.name).into());
    }

    let baseline_actual = [
        ("load_energy_t10", analyses[0].summary.e_at_t10, BASELINE_E),
        (
            "load_ergotropy_t10",
            analyses[0].summary.w_at_t10,
            BASELINE_W,
        ),
        (
            "drive_energy_in_t10",
            analyses[0].summary.drive_energy_in_at_t10,
            BASELINE_EIN,
        ),
        (
            "load_coherence_l1_t10",
            analyses[0].summary.coherence_l1_at_t10,
            BASELINE_COHERENCE_L1,
        ),
    ];
    let baseline: Vec<_> = baseline_actual
        .iter()
        .map(|(name, actual, expected)| {
            let difference = (actual - expected).abs();
            (
                *name,
                *actual,
                *expected,
                difference,
                difference <= BASELINE_TOLERANCE,
            )
        })
        .collect();
    if baseline.iter().any(|row| !row.4) {
        return Err("noise_free did not reproduce the Milestone 5c baseline".into());
    }

    let worst_index = (1..analyses.len())
        .min_by(|&left, &right| {
            analyses[left]
                .summary
                .w_at_t10
                .total_cmp(&analyses[right].summary.w_at_t10)
        })
        .expect("three noisy conditions");
    let mut fine_analyses = Vec::with_capacity(CONDITIONS.len());
    for condition in CONDITIONS {
        fine_analyses.push(run_condition(condition, FINE_DT, &params, &operators)?);
    }
    if fine_analyses
        .iter()
        .any(|analysis| !analysis.quality.all_pass())
    {
        return Err("fine-step numerical diagnostics failed".into());
    }
    let convergence_pairs = [
        (&analyses[0], &fine_analyses[0]),
        (&analyses[worst_index], &fine_analyses[worst_index]),
    ];
    let all_converged = convergence_pairs.iter().all(|(coarse, fine)| {
        convergence_metrics(&coarse.summary)
            .into_iter()
            .zip(convergence_metrics(&fine.summary))
            .all(|((_, base), (_, refined))| converged(base, refined))
    });
    if !all_converged {
        return Err("time-step halving check failed".into());
    }
    let fine_worst_name = ranking(&fine_analyses, |summary| summary.w_at_t10);
    if fine_worst_name != analyses[worst_index].condition.noise_site {
        return Err("worst-position ranking changed after halving dt".into());
    }

    let ratio_rows = ratios(&analyses);
    if ratio_rows.iter().any(|row| {
        [
            row.r_e,
            row.r_w,
            row.r_use,
            row.r_wmax,
            row.r_warea,
            row.r_wmean,
        ]
        .iter()
        .any(|value| value.is_nan())
    }) {
        return Err("ratio denominator was too small".into());
    }
    let rankings = [
        (
            "t=10でload_ergotropyを最も小さくした位置",
            ranking(&analyses, |x| x.w_at_t10),
        ),
        (
            "t=10でusable_fractionを最も小さくした位置",
            ranking(&analyses, |x| x.usable_fraction_at_t10),
        ),
        ("W_maxを最も小さくした位置", ranking(&analyses, |x| x.w_max)),
        (
            "W_time_areaを最も小さくした位置",
            ranking(&analyses, |x| x.w_time_area),
        ),
    ];
    if rankings[0].1 != analyses[worst_index].condition.noise_site {
        return Err("base ranking unexpectedly resolved as a tie".into());
    }

    write_timeseries(&analyses)?;
    write_summary(&analyses)?;
    write_ratios(&ratio_rows)?;
    write_checks(&analyses, same_grid, same_initial, &baseline)?;
    if !write_convergence(&convergence_pairs)? {
        return Err("convergence CSV contains a failed comparison".into());
    }
    write_report(&analyses, &ratio_rows, &convergence_pairs, &rankings)?;
    println!(
        "Milestone 7a complete; worst W(t=10) position={}",
        rankings[0].1
    );
    Ok(())
}
