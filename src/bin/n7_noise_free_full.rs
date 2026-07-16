use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;
use std::time::Instant;

use nalgebra::linalg::SymmetricEigen;
use quantum_work_network::coherent_drive::{
    drive_envelope, drive_hamiltonian, CoherentDriveConfig,
};
use quantum_work_network::ergotropy::ergotropy;
use quantum_work_network::matrix::{
    commutator, expectation, frobenius_norm, hermiticity_error, ComplexMatrix, C64,
};
use quantum_work_network::operators::{build_operators_for_chain, ModelParams, Operators};
use quantum_work_network::partial_trace::partial_trace;

const N: usize = 7;
const DIM: usize = 384;
const DT: f64 = 0.0025;
const T_END: f64 = 10.0;
const SAVE_STEPS: usize = 4;
const SAVE_INTERVAL: f64 = 0.01;
const OMEGA: f64 = 0.2;
const EPS: f64 = 1.0e-14;
const TRACE_TOL: f64 = 1.0e-8;
const HERM_TOL: f64 = 1.0e-8;
const POS_TOL: f64 = 1.0e-8;
const LEDGER_TOL: f64 = 5.0e-5;
const TOP_LEVEL_LIMIT: f64 = 0.05;
const OLD_DIR: &str =
    r"C:\Users\yauki\Documents\Codex\2026-07-14\codex-codex-a-b-2-a\work\quantum_work_network";

const OUTPUTS: [&str; 9] = [
    "n7_noise_free_timeseries.csv",
    "n7_noise_free_site_populations.csv",
    "n7_noise_free_summary.csv",
    "n7_noise_free_arrivals.csv",
    "n7_noise_free_windows.csv",
    "n7_noise_free_length_comparison.csv",
    "n7_noise_free_checks.csv",
    "n7_noise_free_performance.csv",
    "MILESTONE_9A_REPORT.md",
];

#[derive(Clone)]
struct Row {
    time: f64,
    envelope: f64,
    energy: f64,
    work: f64,
    diagonal_work: f64,
    coherence_work: f64,
    coherence_l1: f64,
    usable: f64,
    drive_in: f64,
    drive_net: f64,
    w_over_ein: f64,
    sites: Vec<f64>,
    chain_population: f64,
    load_populations: [f64; 3],
    bare_energy: f64,
    drive_power: f64,
    drive_power_imag: f64,
    trace_error: f64,
    herm_error: f64,
    min_eigenvalue: f64,
    ledger: f64,
    finite: bool,
    reduced_trace_error: f64,
}

struct RunResult {
    rows: Vec<Row>,
    construction_seconds: f64,
    propagation_seconds: f64,
    diagnostics_seconds: f64,
    total_seconds: f64,
    working_set_before: u64,
    working_set_after: u64,
    peak_working_set: u64,
}

#[derive(Clone)]
struct Arrival {
    name: &'static str,
    threshold: f64,
    consecutive: usize,
    time: f64,
    value: f64,
    w_reference: f64,
}

struct Summary {
    endpoint: Row,
    e_max: Row,
    w_max: Row,
    coherence_max: Row,
    arrivals: Vec<Arrival>,
    e_area: f64,
    w_area: f64,
    max_top: f64,
    max_trace: f64,
    max_herm: f64,
    min_eig: f64,
    max_ledger: f64,
    max_drive_power_imag: f64,
    all_finite: bool,
    w_delta_99_100: f64,
    e_delta_99_100: f64,
    w_final_slope: f64,
    e_final_slope: f64,
    peak_class: String,
}

fn config() -> CoherentDriveConfig {
    CoherentDriveConfig {
        omega0: OMEGA,
        omega_drive: 1.0,
        tau: 3.2,
        t_end: T_END,
        dt: DT,
        save_interval: SAVE_INTERVAL,
        gamma_phi: 0.0,
    }
}

fn fmt(x: f64) -> String {
    if x.is_finite() {
        format!("{x:.16e}")
    } else {
        "NaN".to_string()
    }
}

fn ratio(a: f64, b: f64) -> f64 {
    if b.abs() <= EPS {
        f64::NAN
    } else {
        a / b
    }
}

fn process_memory() -> (u64, u64) {
    let pid = std::process::id();
    let script = format!(
        "$p=Get-Process -Id {pid}; Write-Output ($p.WorkingSet64.ToString()+','+$p.PeakWorkingSet64.ToString())"
    );
    if let Ok(output) = Command::new("powershell")
        .args(["-NoProfile", "-Command", &script])
        .output()
    {
        if let Ok(text) = String::from_utf8(output.stdout) {
            let fields: Vec<_> = text.trim().split(',').collect();
            if fields.len() == 2 {
                return (
                    fields[0].parse().unwrap_or(0),
                    fields[1].parse().unwrap_or(0),
                );
            }
        }
    }
    (0, 0)
}

fn ensure_new_outputs() -> Result<(), Box<dyn std::error::Error>> {
    for output in OUTPUTS {
        if std::path::Path::new(output).exists() {
            return Err(format!("refusing to overwrite existing output {output}").into());
        }
    }
    Ok(())
}

fn ensure_references() -> Result<(), Box<dyn std::error::Error>> {
    for name in [
        "chain_length_reachability_summary.csv",
        "chain_length_reachability_arrivals.csv",
        "chain_length_reachability_windows.csv",
        "chain_length_reachability_timeseries.csv",
        "chain_length_reachability_checks.csv",
        "MILESTONE_8A_REPORT.md",
    ] {
        let path = std::path::Path::new(OLD_DIR).join(name);
        if !path.exists() {
            return Err(format!("missing read-only comparison file {}", path.display()).into());
        }
    }
    let checks = std::fs::read_to_string("dephasing_kernel_checks.csv")?;
    if checks.lines().skip(1).any(|line| line.ends_with(",FAIL")) {
        return Err("Milestone 8c contains failed checks".into());
    }
    Ok(())
}

fn construction_checks(
    ops: &Operators,
    rho0: &ComplexMatrix,
) -> Result<Vec<(String, String, String, bool)>, Box<dyn std::error::Error>> {
    let cfg = config();
    let load = partial_trace(rho0, &ops.dims, &[N])?;
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
    let drive0 = drive_hamiltonian(0.0, &cfg, &ops.sigma_1_plus);
    let drive_tau = drive_hamiltonian(cfg.tau, &cfg, &ops.sigma_1_plus);
    let drive_after = drive_hamiltonian(cfg.tau + 0.1, &cfg, &ops.sigma_1_plus);
    let mut checks = vec![
        ("chain_length".into(), N.to_string(), "7".into(), N == 7),
        (
            "hilbert_dimension".into(),
            ops.h_total.nrows().to_string(),
            "384".into(),
            ops.h_total.shape() == (DIM, DIM),
        ),
        (
            "density_matrix_shape".into(),
            format!("{}x{}", rho0.nrows(), rho0.ncols()),
            "384x384".into(),
            rho0.shape() == (DIM, DIM),
        ),
        (
            "density_matrix_elements".into(),
            (rho0.nrows() * rho0.ncols()).to_string(),
            "147456".into(),
            rho0.len() == 147456,
        ),
        (
            "bond_count".into(),
            (N - 1).to_string(),
            "6".into(),
            N - 1 == 6,
        ),
        (
            "drive_site".into(),
            format!("{:?}", driven_sites),
            "[0]".into(),
            driven_sites == vec![0],
        ),
        (
            "load_coupling_site".into(),
            format!("{:?}", load_sites),
            "[6]".into(),
            load_sites == vec![6],
        ),
        (
            "collapse_operator_count".into(),
            "0".into(),
            "0".into(),
            true,
        ),
        (
            "load_reduced_shape".into(),
            format!("{}x{}", load.nrows(), load.ncols()),
            "3x3".into(),
            load.shape() == (3, 3),
        ),
        (
            "initial_trace".into(),
            fmt(rho0.trace().re),
            "1".into(),
            (rho0.trace().re - 1.0).abs() < 1.0e-12,
        ),
        (
            "initial_vacuum".into(),
            fmt(expectation(rho0, &ops.number_total).re),
            "0".into(),
            expectation(rho0, &ops.number_total).re.abs() < 1.0e-12,
        ),
        (
            "operator_dimensions".into(),
            "384x384".into(),
            "all full operators 384x384".into(),
            ops.number_sites
                .iter()
                .chain(ops.sigma_z_sites.iter())
                .all(|x| x.shape() == (DIM, DIM)),
        ),
        (
            "H_drive_0".into(),
            fmt(frobenius_norm(&drive0)),
            "<=1e-12".into(),
            frobenius_norm(&drive0) <= 1.0e-12,
        ),
        (
            "H_drive_tau".into(),
            fmt(frobenius_norm(&drive_tau)),
            "<=1e-12".into(),
            frobenius_norm(&drive_tau) <= 1.0e-12,
        ),
        (
            "H_drive_after_tau".into(),
            fmt(frobenius_norm(&drive_after)),
            "0".into(),
            frobenius_norm(&drive_after) == 0.0,
        ),
    ];
    for time in [0.0, 0.7, cfg.tau, cfg.tau + 0.1] {
        let h = &ops.h_total + drive_hamiltonian(time, &cfg, &ops.sigma_1_plus);
        let error = hermiticity_error(&h);
        checks.push((
            format!("Hamiltonian_hermitian_t{time}"),
            fmt(error),
            "<=1e-12".into(),
            error <= 1.0e-12,
        ));
    }
    Ok(checks)
}

fn rhs(rho: &ComplexMatrix, h: &ComplexMatrix) -> ComplexMatrix {
    (h * rho - rho * h) * C64::new(0.0, -1.0)
}

fn rk4_step(rho: &ComplexMatrix, time: f64, ops: &Operators) -> ComplexMatrix {
    let cfg = config();
    let h = |t| &ops.h_total + drive_hamiltonian(t, &cfg, &ops.sigma_1_plus);
    let half = C64::new(0.5 * DT, 0.0);
    let full = C64::new(DT, 0.0);
    let k1 = rhs(rho, &h(time));
    let k2 = rhs(&(rho + &k1 * half), &h(time + 0.5 * DT));
    let k3 = rhs(&(rho + &k2 * half), &h(time + 0.5 * DT));
    let k4 = rhs(&(rho + &k3 * full), &h(time + DT));
    rho + (k1 + k2 * C64::new(2.0, 0.0) + k3 * C64::new(2.0, 0.0) + k4) * C64::new(DT / 6.0, 0.0)
}

fn minimum_eigenvalue(rho: &ComplexMatrix) -> f64 {
    SymmetricEigen::new(rho.clone())
        .eigenvalues
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min)
}

fn diagnose(
    rho: &ComplexMatrix,
    time: f64,
    ops: &Operators,
    params: &ModelParams,
    drive_in: f64,
    drive_net: f64,
    bare0: f64,
) -> Result<Row, Box<dyn std::error::Error>> {
    let load = partial_trace(rho, &ops.dims, &[N])?;
    let h_load = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
        params.load_dim,
        (0..params.load_dim).map(|level| C64::new(level as f64 * params.omega_load, 0.0)),
    ));
    let work = ergotropy(&load, &h_load, 1.0e-9)?;
    let mut diagonal = ComplexMatrix::zeros(params.load_dim, params.load_dim);
    for level in 0..params.load_dim {
        diagonal[(level, level)] = load[(level, level)];
    }
    let diagonal_work = ergotropy(&diagonal, &h_load, 1.0e-9)?.ergotropy;
    let coherence_l1: f64 = (0..params.load_dim)
        .flat_map(|i| (0..params.load_dim).map(move |j| (i, j)))
        .filter(|(i, j)| i != j)
        .map(|(i, j)| load[(i, j)].norm())
        .sum();
    let sites: Vec<f64> = ops
        .number_sites
        .iter()
        .map(|number| expectation(rho, number).re)
        .collect();
    let chain_population = sites.iter().sum();
    let bare_energy = expectation(rho, &ops.h_total).re;
    let drive = drive_hamiltonian(time, &config(), &ops.sigma_1_plus);
    let drive_power = expectation(rho, &commutator(&drive, &ops.h_total)) * C64::new(0.0, 1.0);
    let trace_error = (rho.trace() - C64::new(1.0, 0.0)).norm();
    let herm_error = hermiticity_error(rho);
    let min_eigenvalue = minimum_eigenvalue(rho);
    let load_populations = [load[(0, 0)].re, load[(1, 1)].re, load[(2, 2)].re];
    let values = [
        work.energy,
        work.ergotropy,
        diagonal_work,
        coherence_l1,
        drive_in,
        drive_net,
        chain_population,
        bare_energy,
        drive_power.re,
        drive_power.im,
        trace_error,
        herm_error,
        min_eigenvalue,
    ];
    Ok(Row {
        time,
        envelope: drive_envelope(time, config().tau),
        energy: work.energy,
        work: work.ergotropy,
        diagonal_work,
        coherence_work: work.ergotropy - diagonal_work,
        coherence_l1,
        usable: ratio(work.ergotropy, work.energy),
        drive_in,
        drive_net,
        w_over_ein: ratio(work.ergotropy, drive_in),
        sites,
        chain_population,
        load_populations,
        bare_energy,
        drive_power: drive_power.re,
        drive_power_imag: drive_power.im,
        trace_error,
        herm_error,
        min_eigenvalue,
        ledger: bare_energy - bare0 - drive_net,
        finite: values.iter().all(|x| x.is_finite())
            && rho.iter().all(|z| z.re.is_finite() && z.im.is_finite()),
        reduced_trace_error: (load.trace() - C64::new(1.0, 0.0)).norm(),
    })
}

fn run_full() -> Result<(RunResult, Vec<(String, String, String, bool)>), Box<dyn std::error::Error>>
{
    let total_start = Instant::now();
    let construction_start = Instant::now();
    let params = ModelParams::default();
    let ops = build_operators_for_chain(&params, N)?;
    let mut rho = ComplexMatrix::zeros(DIM, DIM);
    rho[(0, 0)] = C64::new(1.0, 0.0);
    let checks = construction_checks(&ops, &rho)?;
    if let Some(failed) = checks.iter().find(|x| !x.3) {
        return Err(format!(
            "construction check failed: {} observed {}",
            failed.0, failed.1
        )
        .into());
    }
    let construction_seconds = construction_start.elapsed().as_secs_f64();
    let (working_set_before, _) = process_memory();
    let bare0 = expectation(&rho, &ops.h_total).re;
    let mut drive_in = 0.0;
    let mut drive_net = 0.0;
    let mut previous_power = 0.0;
    let mut rows = Vec::with_capacity(1001);
    let diag0 = Instant::now();
    rows.push(diagnose(
        &rho, 0.0, &ops, &params, drive_in, drive_net, bare0,
    )?);
    let mut diagnostics_seconds = diag0.elapsed().as_secs_f64();
    let mut propagation_seconds = 0.0;
    let steps = (T_END / DT).round() as usize;
    for step in 0..steps {
        let time = step as f64 * DT;
        let propagation_start = Instant::now();
        rho = rk4_step(&rho, time, &ops);
        propagation_seconds += propagation_start.elapsed().as_secs_f64();
        if (step + 1) % SAVE_STEPS == 0 {
            let now = (step + 1) as f64 * DT;
            let diagnostic_start = Instant::now();
            let drive = drive_hamiltonian(now, &config(), &ops.sigma_1_plus);
            let current_power =
                (expectation(&rho, &commutator(&drive, &ops.h_total)) * C64::new(0.0, 1.0)).re;
            drive_net += 0.5 * SAVE_INTERVAL * (previous_power + current_power);
            drive_in += 0.5 * SAVE_INTERVAL * (previous_power.max(0.0) + current_power.max(0.0));
            previous_power = current_power;
            rows.push(diagnose(
                &rho, now, &ops, &params, drive_in, drive_net, bare0,
            )?);
            diagnostics_seconds += diagnostic_start.elapsed().as_secs_f64();
            if rows.len() % 100 == 1 {
                println!(
                    "progress t={now:.2} saved={} propagation={:.1}s diagnostics={:.1}s",
                    rows.len(),
                    propagation_seconds,
                    diagnostics_seconds
                );
            }
        }
    }
    let total_seconds = total_start.elapsed().as_secs_f64();
    let (working_set_after, peak_working_set) = process_memory();
    Ok((
        RunResult {
            rows,
            construction_seconds,
            propagation_seconds,
            diagnostics_seconds,
            total_seconds,
            working_set_before,
            working_set_after,
            peak_working_set,
        },
        checks,
    ))
}

fn area(rows: &[Row], value: impl Fn(&Row) -> f64) -> f64 {
    rows.windows(2)
        .map(|pair| 0.5 * (pair[1].time - pair[0].time) * (value(&pair[0]) + value(&pair[1])))
        .sum()
}

fn sustained_arrival(
    rows: &[Row],
    name: &'static str,
    threshold: f64,
    value: impl Fn(&Row) -> f64,
    w_reference: f64,
) -> Arrival {
    let found = rows
        .windows(5)
        .find(|window| window.iter().all(|row| value(row) >= threshold))
        .map(|window| (window[0].time, value(&window[0])));
    Arrival {
        name,
        threshold,
        consecutive: 5,
        time: found.map(|x| x.0).unwrap_or(f64::NAN),
        value: found.map(|x| x.1).unwrap_or(f64::NAN),
        w_reference,
    }
}

fn linear_slope(rows: &[Row], value: impl Fn(&Row) -> f64) -> f64 {
    let n = rows.len() as f64;
    let mean_t = rows.iter().map(|r| r.time).sum::<f64>() / n;
    let mean_y = rows.iter().map(&value).sum::<f64>() / n;
    let numerator: f64 = rows
        .iter()
        .map(|r| (r.time - mean_t) * (value(r) - mean_y))
        .sum();
    let denominator: f64 = rows.iter().map(|r| (r.time - mean_t).powi(2)).sum();
    numerator / denominator
}

fn summarize(rows: &[Row]) -> Summary {
    let endpoint = rows.last().unwrap().clone();
    let e_max = rows
        .iter()
        .max_by(|a, b| a.energy.total_cmp(&b.energy))
        .unwrap()
        .clone();
    let w_max = rows
        .iter()
        .max_by(|a, b| a.work.total_cmp(&b.work))
        .unwrap()
        .clone();
    let coherence_max = rows
        .iter()
        .max_by(|a, b| a.coherence_l1.total_cmp(&b.coherence_l1))
        .unwrap()
        .clone();
    let arrivals = vec![
        sustained_arrival(rows, "energy_ge_1e-4", 1.0e-4, |r| r.energy, f64::NAN),
        sustained_arrival(rows, "ergotropy_ge_1e-5", 1.0e-5, |r| r.work, f64::NAN),
        sustained_arrival(
            rows,
            "W_ge_10pct_Wmax",
            0.1 * w_max.work,
            |r| r.work,
            w_max.work,
        ),
        sustained_arrival(
            rows,
            "W_ge_50pct_Wmax",
            0.5 * w_max.work,
            |r| r.work,
            w_max.work,
        ),
    ];
    let r99 = rows
        .iter()
        .find(|r| (r.time - 9.9).abs() < 1.0e-12)
        .unwrap();
    let final_ten = &rows[rows.len() - 10..];
    let w_final_slope = linear_slope(final_ten, |r| r.work);
    let e_final_slope = linear_slope(final_ten, |r| r.energy);
    let peak_class = if (w_max.time - T_END).abs() < 1.0e-12 && w_final_slope > 0.0 {
        "still_rising_at_boundary"
    } else if w_max.time > 9.0 {
        "peak_near_boundary"
    } else if w_max.time <= 9.0 && w_final_slope <= 1.0e-6 {
        "peak_resolved"
    } else {
        "ambiguous_oscillatory"
    }
    .to_string();
    let w_delta_99_100 = endpoint.work - r99.work;
    let e_delta_99_100 = endpoint.energy - r99.energy;
    Summary {
        endpoint,
        e_max,
        w_max,
        coherence_max,
        arrivals,
        e_area: area(rows, |r| r.energy),
        w_area: area(rows, |r| r.work),
        max_top: rows
            .iter()
            .map(|r| r.load_populations[2])
            .fold(0.0, f64::max),
        max_trace: rows.iter().map(|r| r.trace_error).fold(0.0, f64::max),
        max_herm: rows.iter().map(|r| r.herm_error).fold(0.0, f64::max),
        min_eig: rows
            .iter()
            .map(|r| r.min_eigenvalue)
            .fold(f64::INFINITY, f64::min),
        max_ledger: rows.iter().map(|r| r.ledger.abs()).fold(0.0, f64::max),
        max_drive_power_imag: rows
            .iter()
            .map(|r| r.drive_power_imag.abs())
            .fold(0.0, f64::max),
        all_finite: rows.iter().all(|r| r.finite),
        w_delta_99_100,
        e_delta_99_100,
        w_final_slope,
        e_final_slope,
        peak_class,
    }
}

fn write_timeseries(rows: &[Row]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("n7_noise_free_timeseries.csv")?);
    writeln!(w, "condition,chain_length,noise_condition,noisy_site_count,gamma_phi_per_noisy_site,hilbert_dimension,time,Omega,drive_envelope,load_energy,load_ergotropy,load_diagonal_ergotropy,load_coherence_ergotropy,load_coherence_l1,usable_fraction,drive_energy_in,drive_energy_net,W_over_Ein,total_chain_population,load_population_0,load_population_1,load_population_2,load_top_level_population,bare_network_energy,drive_power,trace_error,hermiticity_error,min_eigenvalue,energy_ledger_residual")?;
    for r in rows {
        let values = [
            r.time,
            OMEGA,
            r.envelope,
            r.energy,
            r.work,
            r.diagonal_work,
            r.coherence_work,
            r.coherence_l1,
            r.usable,
            r.drive_in,
            r.drive_net,
            r.w_over_ein,
            r.chain_population,
            r.load_populations[0],
            r.load_populations[1],
            r.load_populations[2],
            r.load_populations[2],
            r.bare_energy,
            r.drive_power,
            r.trace_error,
            r.herm_error,
            r.min_eigenvalue,
            r.ledger,
        ];
        writeln!(
            w,
            "N7_noise_free,7,noise_free,0,0,384,{}",
            values.iter().map(|x| fmt(*x)).collect::<Vec<_>>().join(",")
        )?;
    }
    Ok(())
}

fn write_sites(rows: &[Row]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("n7_noise_free_site_populations.csv")?);
    writeln!(
        w,
        "condition,chain_length,time,site_index,site_label,population"
    )?;
    for r in rows {
        for (site, population) in r.sites.iter().enumerate() {
            writeln!(
                w,
                "N7_noise_free,7,{},{},site{},{}",
                fmt(r.time),
                site,
                site + 1,
                fmt(*population)
            )?;
        }
    }
    Ok(())
}

fn write_arrivals(summary: &Summary) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("n7_noise_free_arrivals.csv")?);
    writeln!(w, "condition,arrival_definition,threshold,consecutive_points,arrival_time,value_at_arrival,W_max_reference_if_used")?;
    for a in &summary.arrivals {
        writeln!(
            w,
            "N7_noise_free,{},{},{},{},{},{}",
            a.name,
            fmt(a.threshold),
            a.consecutive,
            fmt(a.time),
            fmt(a.value),
            fmt(a.w_reference)
        )?;
    }
    Ok(())
}

fn write_summary(summary: &Summary) -> std::io::Result<()> {
    let e = &summary.endpoint;
    let w = &summary.w_max;
    let a = &summary.arrivals;
    let mut out = BufWriter::new(File::create("n7_noise_free_summary.csv")?);
    writeln!(out, "condition,chain_length,noise_condition,E_at_t10,W_at_t10,usable_fraction_at_t10,coherence_L1_at_t10,drive_energy_in_at_t10,drive_energy_net_at_t10,W_over_Ein_at_t10,E_max,t_at_E_max,coherence_L1_max,t_at_coherence_L1_max,W_max,t_at_W_max,E_at_W_max,usable_fraction_at_W_max,coherence_L1_at_W_max,drive_energy_in_at_W_max,W_over_Ein_at_W_max,energy_arrival_time,ergotropy_arrival_time,W_10_percent_arrival_time,W_50_percent_arrival_time,E_time_area_0_to_t10,W_time_area_0_to_t10,max_load_top_level_population,max_trace_error,max_hermiticity_error,minimum_density_eigenvalue,max_abs_energy_ledger_residual,maximum_drive_power_imaginary_part,finite_values_pass,endpoint_peak_classification")?;
    let values = [
        e.energy,
        e.work,
        e.usable,
        e.coherence_l1,
        e.drive_in,
        e.drive_net,
        e.w_over_ein,
        summary.e_max.energy,
        summary.e_max.time,
        summary.coherence_max.coherence_l1,
        summary.coherence_max.time,
        w.work,
        w.time,
        w.energy,
        w.usable,
        w.coherence_l1,
        w.drive_in,
        w.w_over_ein,
        a[0].time,
        a[1].time,
        a[2].time,
        a[3].time,
        summary.e_area,
        summary.w_area,
        summary.max_top,
        summary.max_trace,
        summary.max_herm,
        summary.min_eig,
        summary.max_ledger,
        summary.max_drive_power_imag,
    ];
    writeln!(
        out,
        "N7_noise_free,7,noise_free,{},{},{}",
        values.iter().map(|x| fmt(*x)).collect::<Vec<_>>().join(","),
        summary.all_finite,
        summary.peak_class
    )?;
    Ok(())
}

fn window_rows<'a>(rows: &'a [Row], start: f64, end: f64, include_start: bool) -> Vec<&'a Row> {
    rows.iter()
        .filter(|r| {
            if include_start {
                r.time >= start - 1.0e-12
            } else {
                r.time > start + 1.0e-12
            }
        })
        .filter(|r| r.time <= end + 1.0e-12)
        .collect()
}

fn subset_area(rows: &[&Row], value: impl Fn(&Row) -> f64) -> f64 {
    rows.windows(2)
        .map(|p| 0.5 * (p[1].time - p[0].time) * (value(p[0]) + value(p[1])))
        .sum()
}

fn write_windows(rows: &[Row]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("n7_noise_free_windows.csv")?);
    writeln!(w, "condition,window_name,time_start,time_end,point_count,mean_load_energy,mean_load_ergotropy,mean_usable_fraction,mean_coherence_L1,E_time_area,W_time_area,mean_total_chain_population,maximum_load_ergotropy,time_of_window_max_W")?;
    for (name, start, end, include_start) in [
        ("pulse_interval", 0.0, 3.2, true),
        ("early_post_pulse", 3.2, 5.0, false),
        ("middle_interval", 5.0, 7.5, false),
        ("late_interval", 7.5, 10.0, false),
    ] {
        let subset = window_rows(rows, start, end, include_start);
        let n = subset.len() as f64;
        let finite_usable: Vec<f64> = subset
            .iter()
            .map(|r| r.usable)
            .filter(|x| x.is_finite())
            .collect();
        let peak = subset
            .iter()
            .max_by(|a, b| a.work.total_cmp(&b.work))
            .unwrap();
        let values = [
            subset.iter().map(|r| r.energy).sum::<f64>() / n,
            subset.iter().map(|r| r.work).sum::<f64>() / n,
            if finite_usable.is_empty() {
                f64::NAN
            } else {
                finite_usable.iter().sum::<f64>() / finite_usable.len() as f64
            },
            subset.iter().map(|r| r.coherence_l1).sum::<f64>() / n,
            subset_area(&subset, |r| r.energy),
            subset_area(&subset, |r| r.work),
            subset.iter().map(|r| r.chain_population).sum::<f64>() / n,
            peak.work,
            peak.time,
        ];
        writeln!(
            w,
            "N7_noise_free,{name},{},{},{},{}",
            fmt(start),
            fmt(end),
            subset.len(),
            values.iter().map(|x| fmt(*x)).collect::<Vec<_>>().join(",")
        )?;
    }
    Ok(())
}

struct CsvTable {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl CsvTable {
    fn read(path: &std::path::Path) -> Result<Self, Box<dyn std::error::Error>> {
        let mut lines = BufReader::new(File::open(path)?).lines();
        let header = lines
            .next()
            .ok_or("missing CSV header")??
            .split(',')
            .map(str::to_string)
            .collect();
        let rows = lines
            .map(|line| line.map(|x| x.split(',').map(str::to_string).collect()))
            .collect::<Result<_, _>>()?;
        Ok(Self { header, rows })
    }
    fn col(&self, name: &str) -> usize {
        self.header
            .iter()
            .position(|x| x == name)
            .unwrap_or_else(|| panic!("missing column {name}"))
    }
    fn value(&self, condition: &str, metric: &str) -> f64 {
        let condition_col = self.col("condition");
        let metric_col = self.col(metric);
        self.rows
            .iter()
            .find(|r| r[condition_col] == condition)
            .unwrap()[metric_col]
            .parse()
            .unwrap()
    }
}

fn old_arrival(table: &CsvTable, condition: &str, definition: &str) -> f64 {
    let c = table.col("condition");
    let d = table.col("arrival_definition");
    let t = table.col("arrival_time");
    table
        .rows
        .iter()
        .find(|r| r[c] == condition && r[d] == definition)
        .unwrap()[t]
        .parse()
        .unwrap()
}

fn write_comparison(summary: &Summary) -> Result<(), Box<dyn std::error::Error>> {
    let old_summary = CsvTable::read(
        &std::path::Path::new(OLD_DIR).join("chain_length_reachability_summary.csv"),
    )?;
    let old_arrivals = CsvTable::read(
        &std::path::Path::new(OLD_DIR).join("chain_length_reachability_arrivals.csv"),
    )?;
    let e = &summary.endpoint;
    let wmax = &summary.w_max;
    let metrics = vec![
        (
            "E",
            "t10",
            old_summary.value("N3_noise_free", "E_at_t10"),
            old_summary.value("N5_noise_free", "E_at_t10"),
            e.energy,
        ),
        (
            "W",
            "t10",
            old_summary.value("N3_noise_free", "W_at_t10"),
            old_summary.value("N5_noise_free", "W_at_t10"),
            e.work,
        ),
        (
            "usable_fraction",
            "t10",
            old_summary.value("N3_noise_free", "usable_fraction_at_t10"),
            old_summary.value("N5_noise_free", "usable_fraction_at_t10"),
            e.usable,
        ),
        (
            "drive_energy_in",
            "t10",
            old_summary.value("N3_noise_free", "drive_energy_in_at_t10"),
            old_summary.value("N5_noise_free", "drive_energy_in_at_t10"),
            e.drive_in,
        ),
        (
            "W_over_Ein",
            "t10",
            old_summary.value("N3_noise_free", "W_over_Ein_at_t10"),
            old_summary.value("N5_noise_free", "W_over_Ein_at_t10"),
            e.w_over_ein,
        ),
        (
            "W",
            "individual_peak",
            old_summary.value("N3_noise_free", "W_max"),
            old_summary.value("N5_noise_free", "W_max"),
            wmax.work,
        ),
        (
            "t_at_W_max",
            "individual_peak",
            old_summary.value("N3_noise_free", "t_at_W_max"),
            old_summary.value("N5_noise_free", "t_at_W_max"),
            wmax.time,
        ),
        (
            "E_at_W_max",
            "individual_peak",
            old_summary.value("N3_noise_free", "E_at_W_max"),
            old_summary.value("N5_noise_free", "E_at_W_max"),
            wmax.energy,
        ),
        (
            "usable_fraction_at_W_max",
            "individual_peak",
            old_summary.value("N3_noise_free", "usable_fraction_at_W_max"),
            old_summary.value("N5_noise_free", "usable_fraction_at_W_max"),
            wmax.usable,
        ),
        (
            "energy_arrival",
            "sustained_threshold",
            old_arrival(&old_arrivals, "N3_noise_free", "energy_ge_1e-4"),
            old_arrival(&old_arrivals, "N5_noise_free", "energy_ge_1e-4"),
            summary.arrivals[0].time,
        ),
        (
            "ergotropy_arrival",
            "sustained_threshold",
            old_arrival(&old_arrivals, "N3_noise_free", "ergotropy_ge_1e-5"),
            old_arrival(&old_arrivals, "N5_noise_free", "ergotropy_ge_1e-5"),
            summary.arrivals[1].time,
        ),
        (
            "W_10pct_arrival",
            "relative_sustained",
            old_arrival(&old_arrivals, "N3_noise_free", "W_ge_10pct_Wmax"),
            old_arrival(&old_arrivals, "N5_noise_free", "W_ge_10pct_Wmax"),
            summary.arrivals[2].time,
        ),
        (
            "W_50pct_arrival",
            "relative_sustained",
            old_arrival(&old_arrivals, "N3_noise_free", "W_ge_50pct_Wmax"),
            old_arrival(&old_arrivals, "N5_noise_free", "W_ge_50pct_Wmax"),
            summary.arrivals[3].time,
        ),
        (
            "E_time_area",
            "0_to_t10",
            old_summary.value("N3_noise_free", "E_time_area_0_to_t10"),
            old_summary.value("N5_noise_free", "E_time_area_0_to_t10"),
            summary.e_area,
        ),
        (
            "W_time_area",
            "0_to_t10",
            old_summary.value("N3_noise_free", "W_time_area_0_to_t10"),
            old_summary.value("N5_noise_free", "W_time_area_0_to_t10"),
            summary.w_area,
        ),
    ];
    let mut out = BufWriter::new(File::create("n7_noise_free_length_comparison.csv")?);
    writeln!(out, "metric,evaluation_point,N3_value,N5_value,N7_value,N5_over_N3,N7_over_N5,N7_over_N3,N5_minus_N3,N7_minus_N5,N7_minus_N3")?;
    for (metric, point, n3, n5, n7) in metrics {
        writeln!(
            out,
            "{metric},{point},{},{},{},{},{},{},{},{},{}",
            fmt(n3),
            fmt(n5),
            fmt(n7),
            fmt(ratio(n5, n3)),
            fmt(ratio(n7, n5)),
            fmt(ratio(n7, n3)),
            fmt(n5 - n3),
            fmt(n7 - n5),
            fmt(n7 - n3)
        )?;
    }
    Ok(())
}

fn append_check(
    w: &mut impl Write,
    stage: &str,
    check: &str,
    observed: &str,
    expected: &str,
    pass: bool,
) -> std::io::Result<()> {
    writeln!(
        w,
        "{stage},N7_noise_free,{check},{observed},{expected},{}",
        if pass { "PASS" } else { "FAIL" }
    )
}

fn write_checks(
    rows: &[Row],
    summary: &Summary,
    construction: &[(String, String, String, bool)],
) -> Result<bool, Box<dyn std::error::Error>> {
    let mut w = BufWriter::new(File::create("n7_noise_free_checks.csv")?);
    writeln!(w, "stage,condition,check,observed,expected,status")?;
    let mut all = true;
    for (name, observed, expected, pass) in construction {
        append_check(&mut w, "construction", name, observed, expected, *pass)?;
        all &= *pass;
    }
    let time_grid = rows
        .windows(2)
        .all(|p| ((p[1].time - p[0].time) - SAVE_INTERVAL).abs() < 1.0e-12);
    let unique = rows.windows(2).all(|p| p[1].time > p[0].time);
    let site_bounds = rows.iter().all(|r| {
        r.sites
            .iter()
            .all(|x| *x >= -1.0e-10 && *x <= 1.0 + 1.0e-10)
    });
    let load_bounds = rows.iter().all(|r| {
        r.load_populations
            .iter()
            .all(|x| *x >= -1.0e-10 && *x <= 1.0 + 1.0e-10)
    });
    let load_sum = rows
        .iter()
        .map(|r| (r.load_populations.iter().sum::<f64>() - 1.0).abs())
        .fold(0.0, f64::max);
    let chain_sum = rows
        .iter()
        .map(|r| (r.chain_population - r.sites.iter().sum::<f64>()).abs())
        .fold(0.0, f64::max);
    let w_bound = rows
        .iter()
        .all(|r| r.work >= -1.0e-10 && r.work <= r.energy + 1.0e-10);
    let usable_bound = rows
        .iter()
        .all(|r| !r.usable.is_finite() || (r.usable >= -1.0e-9 && r.usable <= 1.0 + 1.0e-9));
    let checks: Vec<(&str, String, String, bool)> = vec![
        (
            "trace_preservation",
            fmt(summary.max_trace),
            "<=1e-8".into(),
            summary.max_trace <= TRACE_TOL,
        ),
        (
            "hermiticity",
            fmt(summary.max_herm),
            "<=1e-8".into(),
            summary.max_herm <= HERM_TOL,
        ),
        (
            "positivity",
            fmt(summary.min_eig),
            ">=-1e-8".into(),
            summary.min_eig >= -POS_TOL,
        ),
        (
            "finite_values",
            summary.all_finite.to_string(),
            "true".into(),
            summary.all_finite,
        ),
        (
            "load_reduced_trace",
            fmt(rows
                .iter()
                .map(|r| r.reduced_trace_error)
                .fold(0.0, f64::max)),
            "<=1e-8".into(),
            rows.iter().all(|r| r.reduced_trace_error <= TRACE_TOL),
        ),
        (
            "site_population_bounds",
            site_bounds.to_string(),
            "true".into(),
            site_bounds,
        ),
        (
            "load_population_bounds",
            load_bounds.to_string(),
            "true".into(),
            load_bounds,
        ),
        (
            "load_population_sum",
            fmt(load_sum),
            "<=1e-8".into(),
            load_sum <= TRACE_TOL,
        ),
        (
            "chain_population_consistency",
            fmt(chain_sum),
            "<=1e-12".into(),
            chain_sum <= 1.0e-12,
        ),
        (
            "load_top_level",
            fmt(summary.max_top),
            "<0.05".into(),
            summary.max_top < TOP_LEVEL_LIMIT,
        ),
        ("W_le_E", w_bound.to_string(), "true".into(), w_bound),
        (
            "usable_fraction_range",
            usable_bound.to_string(),
            "NaN before signal or 0..1".into(),
            usable_bound,
        ),
        (
            "drive_energy_consistency",
            format!(
                "in={} net={}",
                fmt(summary.endpoint.drive_in),
                fmt(summary.endpoint.drive_net)
            ),
            "energy_in nonnegative and >= net".into(),
            summary.endpoint.drive_in >= -1.0e-12
                && summary.endpoint.drive_in + 1.0e-12 >= summary.endpoint.drive_net,
        ),
        (
            "energy_ledger",
            fmt(summary.max_ledger),
            "<=5e-5".into(),
            summary.max_ledger <= LEDGER_TOL,
        ),
        (
            "common_time_grid",
            time_grid.to_string(),
            "true".into(),
            time_grid,
        ),
        (
            "time_monotonic_unique",
            unique.to_string(),
            "true".into(),
            unique,
        ),
        (
            "saved_point_count",
            rows.len().to_string(),
            "1001".into(),
            rows.len() == 1001,
        ),
        (
            "timeseries_row_count",
            rows.len().to_string(),
            "1001".into(),
            rows.len() == 1001,
        ),
        (
            "site_population_row_count",
            (rows.len() * N).to_string(),
            "7007".into(),
            rows.len() * N == 7007,
        ),
        (
            "initial_load_zero",
            fmt(rows[0].energy + rows[0].work + rows[0].coherence_l1),
            "0".into(),
            (rows[0].energy + rows[0].work + rows[0].coherence_l1).abs() < 1.0e-12,
        ),
        (
            "drive_envelope_endpoints",
            format!(
                "{} {}",
                fmt(rows[0].envelope),
                fmt(rows.last().unwrap().envelope)
            ),
            "0 0".into(),
            rows[0].envelope == 0.0 && rows.last().unwrap().envelope == 0.0,
        ),
        (
            "no_allocation_failure",
            "completed".into(),
            "completed".into(),
            true,
        ),
    ];
    for (name, observed, expected, pass) in checks {
        append_check(&mut w, "trajectory", name, &observed, &expected, pass)?;
        all &= pass;
    }
    Ok(all)
}

fn read_estimate(
    path: &str,
    condition: &str,
    value_col: &str,
) -> Result<f64, Box<dyn std::error::Error>> {
    let table = CsvTable::read(std::path::Path::new(path))?;
    Ok(table.value(condition, value_col))
}

fn write_performance(run: &RunResult) -> Result<(), Box<dyn std::error::Error>> {
    let estimated = read_estimate(
        "n7_feasibility_estimates.csv",
        "N7_noise_free",
        "estimated_seconds",
    )?;
    let measured = run.propagation_seconds + run.diagnostics_seconds;
    let mut w = BufWriter::new(File::create("n7_noise_free_performance.csv")?);
    writeln!(w, "condition,chain_length,hilbert_dimension,density_matrix_element_count,step_count,saved_time_points,construction_seconds,propagation_seconds,diagnostics_seconds,total_seconds,seconds_per_step,estimated_seconds_from_8b,measured_to_estimated_ratio,peak_memory_bytes_if_available,process_working_set_before,process_working_set_after,timeseries_rows,site_population_rows")?;
    let fields = vec![
        "N7_noise_free".to_string(),
        "7".to_string(),
        "384".to_string(),
        "147456".to_string(),
        "4000".to_string(),
        run.rows.len().to_string(),
        fmt(run.construction_seconds),
        fmt(run.propagation_seconds),
        fmt(run.diagnostics_seconds),
        fmt(run.total_seconds),
        fmt(run.propagation_seconds / 4000.0),
        fmt(estimated),
        fmt(ratio(measured, estimated)),
        run.peak_working_set.to_string(),
        run.working_set_before.to_string(),
        run.working_set_after.to_string(),
        run.rows.len().to_string(),
        (run.rows.len() * N).to_string(),
    ];
    writeln!(w, "{}", fields.join(","))?;
    Ok(())
}

fn write_report(
    run: &RunResult,
    summary: &Summary,
    checks_pass: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let n5_wmax = CsvTable::read(
        &std::path::Path::new(OLD_DIR).join("chain_length_reachability_summary.csv"),
    )?
    .value("N5_noise_free", "W_max");
    let n3_wmax = CsvTable::read(
        &std::path::Path::new(OLD_DIR).join("chain_length_reachability_summary.csv"),
    )?
    .value("N3_noise_free", "W_max");
    let noisy_hours = read_estimate(
        "dephasing_kernel_estimates.csv",
        "N7_all_site_noisy",
        "estimated_t10_hours",
    )?;
    let reached_e = summary.arrivals[0].time.is_finite();
    let reached_w = summary.arrivals[1].time.is_finite();
    let decision = if !checks_pass {
        "numerical_issue_stop"
    } else if !reached_e || !reached_w {
        "insufficient_reachability"
    } else if summary.peak_class != "peak_resolved" {
        "extend_time_before_noisy_candidate"
    } else if noisy_hours <= 6.0 {
        "proceed_candidate"
    } else {
        "performance_issue_stop"
    };
    let e = &summary.endpoint;
    let p = &summary.w_max;
    let mut w = BufWriter::new(File::create("MILESTONE_9A_REPORT.md")?);
    writeln!(w, "# Milestone 9a: N=7 noise-free full reachability run\n")?;
    let sections = [
        ("1. 目的", "N=7 noise-freeをdt=0.0025でt=10まで本計算し、load energy・ergotropy到達とN=3/N=5との差を確認した。".to_string()),
        ("2. Milestone 8aから8cまでの位置づけ", "8aのN=3/N=5結果、8bのN=7 feasibility、8cの厳密dephasing kernel検証を前提にした。".to_string()),
        ("3. 今回はN=7 noise-freeだけであること", "N=7 all-site noisy、細刻み、t>10、追加最適化は実行していない。".to_string()),
        ("4. 変更していない物理模型", "Hamiltonian、drive、RK4、dt、基底、load、初期真空、観測量を既存実装から変更していない。".to_string()),
        ("5. N=7模型構成", "7二準位site、3準位load、dim=384、bond=6、drive site=0、load coupling site=6、collapse=0。".to_string()),
        ("6. 数値手法", "dense complex density matrixのtime-dependent RK4。4000 step、保存間隔0.01、1001点。各保存点で縮約・ergotropy・全系最小固有値を診断した。".to_string()),
        ("7. 構成検査", "次元、mapping、真空、drive端点、Hamiltonian Hermiticityを時間発展前に検査し全合格した。".to_string()),
        ("8. 実行時間とメモリ", format!("construction {:.3}s、propagation {:.3}s、diagnostics {:.3}s、total {:.3}s。peak working set {} bytes。", run.construction_seconds, run.propagation_seconds, run.diagnostics_seconds, run.total_seconds, run.peak_working_set)),
        ("9. 数値品質チェック", format!("全チェック{}。max trace={:.3e}、max Hermiticity={:.3e}、min eigenvalue={:.3e}、max ledger={:.3e}。", if checks_pass { "PASS" } else { "FAIL" }, summary.max_trace, summary.max_herm, summary.min_eig, summary.max_ledger)),
        ("10. load energy到達", format!("持続閾値1e-4の到達時刻は {}。", fmt(summary.arrivals[0].time))),
        ("11. load ergotropy到達", format!("持続閾値1e-5の到達時刻は {}。", fmt(summary.arrivals[1].time))),
        ("12. coherence生成", format!("coherence L1最大値は {:.10e}、時刻 {:.2}。", summary.coherence_max.coherence_l1, summary.coherence_max.time)),
        ("13. t=10結果", format!("E={:.10e}、W={:.10e}、usable={:.10e}、W/Ein={:.10e}。", e.energy, e.work, e.usable, e.w_over_ein)),
        ("14. W最大値と時刻", format!("W_max={:.10e}、t={:.2}、その時のE={:.10e}。", p.work, p.time, p.energy)),
        ("15. 終端挙動とピーク判定", format!("分類 `{}`。W(10)-W(9.9)={:.3e}、E(10)-E(9.9)={:.3e}、最終10点W slope={:.3e}、E slope={:.3e}。", summary.peak_class, summary.w_delta_99_100, summary.e_delta_99_100, summary.w_final_slope, summary.e_final_slope)),
        ("16. usable fraction", format!("t=10で {:.10e}、W最大時で {:.10e}。", e.usable, p.usable)),
        ("17. W/Ein", format!("t=10で {:.10e}、W最大時で {:.10e}。制御費用を含む総合効率ではない。", e.w_over_ein, p.w_over_ein)),
        ("18. site population時間発展", "1001時刻×7site=7007行をlong形式で保存した。site_indexは0始まり、site_labelは1始まり。".to_string()),
        ("19. 時間窓解析", "pulse、early post-pulse、middle、lateの4窓をwindows CSVへ保存した。時間面積は状態量の面積であり累積仕事ではない。".to_string()),
        ("20. N=3、N=5、N=7比較", format!("WmaxはN3 {:.10e}、N5 {:.10e}、N7 {:.10e}。N7/N5={:.10e}、N7/N3={:.10e}。", n3_wmax, n5_wmax, p.work, ratio(p.work, n5_wmax), ratio(p.work, n3_wmax))),
        ("21. N=3からN=7への有限差", "3点の有限長比較のみであり、指数・べき則・漸近scalingは推定しない。Nとともに距離、次元、bond数、干渉構造が同時に変わる。".to_string()),
        ("22. 直接確認できたこと", "固定模型でN=7のt=10までのenergy/W到達、W最大値、usable、W/Ein、有限長3点差、実測時間を確認した。".to_string()),
        ("23. 確認できていないこと", "N=7 noisy、細刻み、t>10、最終到達上限、N>7、scaling、実機性能は未確認。".to_string()),
        ("24. 主張してはいけないこと", "指数/べき減衰、熱力学極限、N>7外挿、物理的輸送限界、noisy予測、量子優位、新規性、距離だけの純粋因果。".to_string()),
        ("25. N=7 noisy本計算へ進む判断", format!("判定 **{}**。8c noisy推定は約{:.3}時間。ただしnoisy本計算は自動実行していない。", decision, noisy_hours)),
        ("26. 生成ファイル一覧", "`src/bin/n7_noise_free_full.rs` と指定9成果物を新規作成した。".to_string()),
    ];
    for (title, body) in sections {
        writeln!(w, "## {title}\n\n{body}\n")?;
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    ensure_new_outputs()?;
    ensure_references()?;
    if std::env::args().nth(1).as_deref() == Some("--preflight") {
        let params = ModelParams::default();
        let ops = build_operators_for_chain(&params, N)?;
        let mut rho = ComplexMatrix::zeros(DIM, DIM);
        rho[(0, 0)] = C64::new(1.0, 0.0);
        let checks = construction_checks(&ops, &rho)?;
        for (name, observed, expected, pass) in &checks {
            println!("{name}: observed={observed} expected={expected} pass={pass}");
        }
        if checks.iter().any(|x| !x.3) {
            return Err("preflight construction checks failed".into());
        }
        println!("preflight PASS: {} construction checks", checks.len());
        return Ok(());
    }
    let (run, construction_checks) = run_full()?;
    let summary = summarize(&run.rows);
    write_timeseries(&run.rows)?;
    write_sites(&run.rows)?;
    write_summary(&summary)?;
    write_arrivals(&summary)?;
    write_windows(&run.rows)?;
    write_comparison(&summary)?;
    let checks_pass = write_checks(&run.rows, &summary, &construction_checks)?;
    write_performance(&run)?;
    write_report(&run, &summary, checks_pass)?;
    println!(
        "completed checks={} E10={:.10e} W10={:.10e} Wmax={:.10e} tWmax={:.2} peak={}",
        checks_pass,
        summary.endpoint.energy,
        summary.endpoint.work,
        summary.w_max.work,
        summary.w_max.time,
        summary.peak_class
    );
    if !checks_pass {
        return Err("numerical quality checks failed".into());
    }
    Ok(())
}
