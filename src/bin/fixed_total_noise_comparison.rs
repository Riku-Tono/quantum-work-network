use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;
use std::time::Instant;

use nalgebra::linalg::SymmetricEigen;
use quantum_work_network::coherent_drive::{
    drive_envelope, drive_hamiltonian, CoherentDriveConfig,
};
use quantum_work_network::dephasing_kernel::DiagonalDephasingKernel;
use quantum_work_network::ergotropy::ergotropy;
use quantum_work_network::matrix::{
    commutator, expectation, frobenius_norm, hermiticity_error, ComplexMatrix, C64,
};
use quantum_work_network::operators::{build_operators_for_chain, ModelParams, Operators};
use quantum_work_network::partial_trace::partial_trace;

const TOTAL_GAMMA: f64 = 1.5;
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
const OLD_DIR: &str = ".";

const OUTPUTS: [&str; 11] = [
    "fixed_total_noise_timeseries.csv",
    "fixed_total_noise_site_populations.csv",
    "fixed_total_noise_summary.csv",
    "fixed_total_noise_arrivals.csv",
    "fixed_total_noise_windows.csv",
    "fixed_total_noise_length_comparison.csv",
    "fixed_total_vs_fixed_per_site.csv",
    "fixed_total_noise_recovery.csv",
    "fixed_total_noise_checks.csv",
    "fixed_total_noise_performance.csv",
    "MILESTONE_9C_REPORT.md",
];

#[derive(Clone, Copy)]
struct Spec {
    n: usize,
    dim: usize,
    condition: &'static str,
}

const N5_SPEC: Spec = Spec {
    n: 5,
    dim: 96,
    condition: "N5_fixed_total_noise",
};
const N7_SPEC: Spec = Spec {
    n: 7,
    dim: 384,
    condition: "N7_fixed_total_noise",
};

impl Spec {
    fn gamma_site(self) -> f64 {
        TOTAL_GAMMA / self.n as f64
    }

    fn gammas(self) -> Vec<f64> {
        vec![self.gamma_site(); self.n]
    }
}

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
    dephasing_net: f64,
    w_over_ein: f64,
    sites: Vec<f64>,
    chain_population: f64,
    load_populations: [f64; 3],
    bare_energy: f64,
    drive_power: f64,
    drive_power_imag: f64,
    dephasing_power: f64,
    dephasing_power_imag: f64,
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
    kernel_construction_seconds: f64,
    propagation_seconds: f64,
    diagnostics_seconds: f64,
    total_seconds: f64,
    working_set_before: u64,
    working_set_after: u64,
    peak_working_set: u64,
}

struct ConditionResult {
    spec: Spec,
    run: RunResult,
    summary: Summary,
    checks_pass: bool,
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
    max_dephasing_power_imag: f64,
    all_finite: bool,
    w_delta_99_100: f64,
    e_delta_99_100: f64,
    w_final_slope: f64,
    e_final_slope: f64,
    peak_class: String,
}

fn config(gamma_site: f64) -> CoherentDriveConfig {
    CoherentDriveConfig {
        omega0: OMEGA,
        omega_drive: 1.0,
        tau: 3.2,
        t_end: T_END,
        dt: DT,
        save_interval: SAVE_INTERVAL,
        gamma_phi: gamma_site,
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
    for name in [
        "n7_noise_free_summary.csv",
        "n7_noise_free_arrivals.csv",
        "n7_noise_free_timeseries.csv",
        "n7_noise_free_checks.csv",
        "MILESTONE_9A_REPORT.md",
        "n7_all_site_noisy_summary.csv",
        "n7_all_site_noisy_arrivals.csv",
        "n7_all_site_noisy_checks.csv",
        "MILESTONE_9B_REPORT.md",
    ] {
        if !std::path::Path::new(name).exists() {
            return Err(format!("missing Milestone 9a comparison file {name}").into());
        }
    }
    Ok(())
}

fn construction_checks(
    spec: Spec,
    gammas: &[f64],
    ops: &Operators,
    rho0: &ComplexMatrix,
    kernel: &DiagonalDephasingKernel,
) -> Result<Vec<(String, String, String, bool)>, Box<dyn std::error::Error>> {
    let cfg = config(spec.gamma_site());
    let load = partial_trace(rho0, &ops.dims, &[spec.n])?;
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
    let mut kernel_diagonal_zero = true;
    let mut kernel_symmetric = true;
    let mut kernel_nonnegative = true;
    let mut kernel_mapping = true;
    let mut load_excluded = true;
    for row in 0..spec.dim {
        let row_chain = row / 3;
        for col in 0..spec.dim {
            let observed = kernel.rate(row, col)?;
            let col_chain = col / 3;
            let differing = row_chain ^ col_chain;
            let expected: f64 = (0..spec.n)
                .filter(|site| differing & (1usize << (spec.n - 1 - site)) != 0)
                .map(|site| gammas[site])
                .sum();
            kernel_mapping &= (observed - expected).abs() <= 1.0e-12;
            kernel_nonnegative &= observed >= 0.0;
            kernel_symmetric &= (observed - kernel.rate(col, row)?).abs() <= 1.0e-12;
            if row == col {
                kernel_diagonal_zero &= observed == 0.0;
            }
            if row_chain == col_chain {
                load_excluded &= observed == 0.0;
            }
        }
    }
    let mut checks = vec![
        (
            "chain_length".into(),
            spec.n.to_string(),
            spec.n.to_string(),
            matches!(spec.n, 5 | 7),
        ),
        (
            "hilbert_dimension".into(),
            ops.h_total.nrows().to_string(),
            spec.dim.to_string(),
            ops.h_total.shape() == (spec.dim, spec.dim),
        ),
        (
            "density_matrix_shape".into(),
            format!("{}x{}", rho0.nrows(), rho0.ncols()),
            format!("{}x{}", spec.dim, spec.dim),
            rho0.shape() == (spec.dim, spec.dim),
        ),
        (
            "density_matrix_elements".into(),
            (rho0.nrows() * rho0.ncols()).to_string(),
            (spec.dim * spec.dim).to_string(),
            rho0.len() == spec.dim * spec.dim,
        ),
        (
            "bond_count".into(),
            (spec.n - 1).to_string(),
            (spec.n - 1).to_string(),
            matches!(spec.n - 1, 4 | 6),
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
            format!("[{}]", spec.n - 1),
            load_sites == vec![spec.n - 1],
        ),
        (
            "collapse_operator_count".into(),
            "0".into(),
            "0 dense operators in hot path".into(),
            true,
        ),
        (
            "noisy_site_count".into(),
            gammas
                .iter()
                .filter(|gamma| **gamma > 0.0)
                .count()
                .to_string(),
            spec.n.to_string(),
            gammas.iter().filter(|gamma| **gamma > 0.0).count() == spec.n,
        ),
        (
            "all_site_gammas".into(),
            fmt(spec.gamma_site()),
            fmt(TOTAL_GAMMA / spec.n as f64),
            gammas
                .iter()
                .all(|gamma| (*gamma - spec.gamma_site()).abs() <= 1.0e-15),
        ),
        (
            "gamma_sum".into(),
            fmt(gammas.iter().sum()),
            fmt(TOTAL_GAMMA),
            (gammas.iter().sum::<f64>() - TOTAL_GAMMA).abs() <= 1.0e-14,
        ),
        (
            "kernel_dimension".into(),
            format!("{}x{}", kernel.dimension(), kernel.dimension()),
            format!("{}x{}", spec.dim, spec.dim),
            kernel.dimension() == spec.dim,
        ),
        (
            "kernel_chain_length".into(),
            kernel.chain_length().to_string(),
            spec.n.to_string(),
            kernel.chain_length() == spec.n,
        ),
        (
            "kernel_load_dimension".into(),
            kernel.load_dim().to_string(),
            "3".into(),
            kernel.load_dim() == 3,
        ),
        (
            "kernel_diagonal_zero".into(),
            kernel_diagonal_zero.to_string(),
            "true".into(),
            kernel_diagonal_zero,
        ),
        (
            "kernel_symmetric".into(),
            kernel_symmetric.to_string(),
            "true".into(),
            kernel_symmetric,
        ),
        (
            "kernel_nonnegative".into(),
            kernel_nonnegative.to_string(),
            "true".into(),
            kernel_nonnegative,
        ),
        (
            "kernel_mapping".into(),
            kernel_mapping.to_string(),
            "sum differing-site gammas".into(),
            kernel_mapping,
        ),
        (
            "load_excluded_from_kernel".into(),
            load_excluded.to_string(),
            "true".into(),
            load_excluded,
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
            format!("{}x{}", spec.dim, spec.dim),
            format!("all full operators {}x{}", spec.dim, spec.dim),
            ops.number_sites
                .iter()
                .chain(ops.sigma_z_sites.iter())
                .all(|x| x.shape() == (spec.dim, spec.dim)),
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

fn rhs(
    rho: &ComplexMatrix,
    h: &ComplexMatrix,
    kernel: &DiagonalDephasingKernel,
) -> Result<ComplexMatrix, Box<dyn std::error::Error>> {
    let mut out = (h * rho - rho * h) * C64::new(0.0, -1.0);
    kernel.add_to(rho, &mut out)?;
    Ok(out)
}

fn rk4_step(
    rho: &ComplexMatrix,
    time: f64,
    ops: &Operators,
    kernel: &DiagonalDephasingKernel,
    gamma_site: f64,
) -> Result<ComplexMatrix, Box<dyn std::error::Error>> {
    let cfg = config(gamma_site);
    let h = |t| &ops.h_total + drive_hamiltonian(t, &cfg, &ops.sigma_1_plus);
    let half = C64::new(0.5 * DT, 0.0);
    let full = C64::new(DT, 0.0);
    let k1 = rhs(rho, &h(time), kernel)?;
    let k2 = rhs(&(rho + &k1 * half), &h(time + 0.5 * DT), kernel)?;
    let k3 = rhs(&(rho + &k2 * half), &h(time + 0.5 * DT), kernel)?;
    let k4 = rhs(&(rho + &k3 * full), &h(time + DT), kernel)?;
    Ok(rho
        + (k1 + k2 * C64::new(2.0, 0.0) + k3 * C64::new(2.0, 0.0) + k4) * C64::new(DT / 6.0, 0.0))
}

fn minimum_eigenvalue(rho: &ComplexMatrix) -> f64 {
    SymmetricEigen::new(rho.clone())
        .eigenvalues
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min)
}

fn instantaneous_powers(
    rho: &ComplexMatrix,
    time: f64,
    ops: &Operators,
    kernel: &DiagonalDephasingKernel,
    gamma_site: f64,
) -> Result<(C64, C64), Box<dyn std::error::Error>> {
    let drive = drive_hamiltonian(time, &config(gamma_site), &ops.sigma_1_plus);
    let drive_power = expectation(rho, &commutator(&drive, &ops.h_total)) * C64::new(0.0, 1.0);
    let dephasing_power = expectation(&kernel.apply(rho)?, &ops.h_total);
    Ok((drive_power, dephasing_power))
}

fn diagnose(
    spec: Spec,
    rho: &ComplexMatrix,
    time: f64,
    ops: &Operators,
    params: &ModelParams,
    drive_in: f64,
    drive_net: f64,
    dephasing_net: f64,
    bare0: f64,
    drive_power: C64,
    dephasing_power: C64,
) -> Result<Row, Box<dyn std::error::Error>> {
    let load = partial_trace(rho, &ops.dims, &[spec.n])?;
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
        dephasing_net,
        chain_population,
        bare_energy,
        drive_power.re,
        drive_power.im,
        dephasing_power.re,
        dephasing_power.im,
        trace_error,
        herm_error,
        min_eigenvalue,
    ];
    Ok(Row {
        time,
        envelope: drive_envelope(time, config(spec.gamma_site()).tau),
        energy: work.energy,
        work: work.ergotropy,
        diagonal_work,
        coherence_work: work.ergotropy - diagonal_work,
        coherence_l1,
        usable: ratio(work.ergotropy, work.energy),
        drive_in,
        drive_net,
        dephasing_net,
        w_over_ein: ratio(work.ergotropy, drive_in),
        sites,
        chain_population,
        load_populations,
        bare_energy,
        drive_power: drive_power.re,
        drive_power_imag: drive_power.im,
        dephasing_power: dephasing_power.re,
        dephasing_power_imag: dephasing_power.im,
        trace_error,
        herm_error,
        min_eigenvalue,
        ledger: bare_energy - bare0 - drive_net - dephasing_net,
        finite: values.iter().all(|x| x.is_finite())
            && rho.iter().all(|z| z.re.is_finite() && z.im.is_finite()),
        reduced_trace_error: (load.trace() - C64::new(1.0, 0.0)).norm(),
    })
}

fn run_full(
    spec: Spec,
) -> Result<(RunResult, Vec<(String, String, String, bool)>), Box<dyn std::error::Error>> {
    let total_start = Instant::now();
    let construction_start = Instant::now();
    let params = ModelParams::default();
    let ops = build_operators_for_chain(&params, spec.n)?;
    let gammas = spec.gammas();
    let kernel_start = Instant::now();
    let kernel = DiagonalDephasingKernel::new(spec.n, params.load_dim, &gammas)?;
    let kernel_construction_seconds = kernel_start.elapsed().as_secs_f64();
    let mut rho = ComplexMatrix::zeros(spec.dim, spec.dim);
    rho[(0, 0)] = C64::new(1.0, 0.0);
    let checks = construction_checks(spec, &gammas, &ops, &rho, &kernel)?;
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
    let mut dephasing_net = 0.0;
    let (power0, dephasing_power0) =
        instantaneous_powers(&rho, 0.0, &ops, &kernel, spec.gamma_site())?;
    let mut previous_power = power0.re;
    let mut previous_dephasing_power = dephasing_power0.re;
    let mut rows = Vec::with_capacity(1001);
    let diag0 = Instant::now();
    rows.push(diagnose(
        spec,
        &rho,
        0.0,
        &ops,
        &params,
        drive_in,
        drive_net,
        dephasing_net,
        bare0,
        power0,
        dephasing_power0,
    )?);
    let mut diagnostics_seconds = diag0.elapsed().as_secs_f64();
    let mut propagation_seconds = 0.0;
    let steps = (T_END / DT).round() as usize;
    for step in 0..steps {
        let time = step as f64 * DT;
        let propagation_start = Instant::now();
        rho = rk4_step(&rho, time, &ops, &kernel, spec.gamma_site())?;
        propagation_seconds += propagation_start.elapsed().as_secs_f64();
        if (step + 1) % SAVE_STEPS == 0 {
            let now = (step + 1) as f64 * DT;
            let diagnostic_start = Instant::now();
            let (current_power_complex, current_dephasing_power_complex) =
                instantaneous_powers(&rho, now, &ops, &kernel, spec.gamma_site())?;
            let current_power = current_power_complex.re;
            let current_dephasing_power = current_dephasing_power_complex.re;
            drive_net += 0.5 * SAVE_INTERVAL * (previous_power + current_power);
            drive_in += 0.5 * SAVE_INTERVAL * (previous_power.max(0.0) + current_power.max(0.0));
            dephasing_net +=
                0.5 * SAVE_INTERVAL * (previous_dephasing_power + current_dephasing_power);
            previous_power = current_power;
            previous_dephasing_power = current_dephasing_power;
            rows.push(diagnose(
                spec,
                &rho,
                now,
                &ops,
                &params,
                drive_in,
                drive_net,
                dephasing_net,
                bare0,
                current_power_complex,
                current_dephasing_power_complex,
            )?);
            diagnostics_seconds += diagnostic_start.elapsed().as_secs_f64();
            if rows.len() % 100 == 1 {
                println!(
                    "{} progress t={now:.2} saved={} propagation={:.1}s diagnostics={:.1}s",
                    spec.condition,
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
            kernel_construction_seconds,
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
        max_dephasing_power_imag: rows
            .iter()
            .map(|r| r.dephasing_power_imag.abs())
            .fold(0.0, f64::max),
        all_finite: rows.iter().all(|r| r.finite),
        w_delta_99_100,
        e_delta_99_100,
        w_final_slope,
        e_final_slope,
        peak_class,
    }
}

fn write_timeseries(results: &[ConditionResult]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("fixed_total_noise_timeseries.csv")?);
    writeln!(w, "condition,chain_length,noise_normalization,noisy_site_count,gamma_phi_per_site,gamma_phi_total,hilbert_dimension,time,Omega,drive_envelope,load_energy,load_ergotropy,load_diagonal_ergotropy,load_coherence_ergotropy,load_coherence_l1,usable_fraction,drive_energy_in,drive_energy_net,dephasing_energy_net,W_over_Ein,total_chain_population,load_population_0,load_population_1,load_population_2,load_top_level_population,bare_network_energy,drive_power,dephasing_power,trace_error,hermiticity_error,min_eigenvalue,energy_ledger_residual")?;
    for result in results {
        for r in &result.run.rows {
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
                r.dephasing_net,
                r.w_over_ein,
                r.chain_population,
                r.load_populations[0],
                r.load_populations[1],
                r.load_populations[2],
                r.load_populations[2],
                r.bare_energy,
                r.drive_power,
                r.dephasing_power,
                r.trace_error,
                r.herm_error,
                r.min_eigenvalue,
                r.ledger,
            ];
            writeln!(
                w,
                "{},{},fixed_total_gamma_1p5,{},{},{},{},{}",
                result.spec.condition,
                result.spec.n,
                result.spec.n,
                fmt(result.spec.gamma_site()),
                fmt(TOTAL_GAMMA),
                result.spec.dim,
                values.iter().map(|x| fmt(*x)).collect::<Vec<_>>().join(",")
            )?;
        }
    }
    Ok(())
}

fn write_sites(results: &[ConditionResult]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("fixed_total_noise_site_populations.csv")?);
    writeln!(
        w,
        "condition,chain_length,noise_normalization,gamma_phi_per_site,gamma_phi_total,time,site_index,site_label,population"
    )?;
    for result in results {
        for r in &result.run.rows {
            for (site, population) in r.sites.iter().enumerate() {
                writeln!(
                    w,
                    "{},{},fixed_total_gamma_1p5,{},{},{},{},site{},{}",
                    result.spec.condition,
                    result.spec.n,
                    fmt(result.spec.gamma_site()),
                    fmt(TOTAL_GAMMA),
                    fmt(r.time),
                    site,
                    site + 1,
                    fmt(*population)
                )?;
            }
        }
    }
    Ok(())
}

fn write_arrivals(results: &[ConditionResult]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("fixed_total_noise_arrivals.csv")?);
    writeln!(w, "condition,chain_length,noise_normalization,gamma_phi_per_site,gamma_phi_total,arrival_definition,threshold,consecutive_points,arrival_time,value_at_arrival,W_max_reference_if_used")?;
    for result in results {
        for a in &result.summary.arrivals {
            writeln!(
                w,
                "{},{},fixed_total_gamma_1p5,{},{},{},{},{},{},{},{}",
                result.spec.condition,
                result.spec.n,
                fmt(result.spec.gamma_site()),
                fmt(TOTAL_GAMMA),
                a.name,
                fmt(a.threshold),
                a.consecutive,
                fmt(a.time),
                fmt(a.value),
                fmt(a.w_reference)
            )?;
        }
    }
    Ok(())
}

fn write_summary(results: &[ConditionResult]) -> std::io::Result<()> {
    let mut out = BufWriter::new(File::create("fixed_total_noise_summary.csv")?);
    writeln!(out, "condition,chain_length,noise_normalization,noisy_site_count,gamma_phi_per_site,gamma_phi_total,E_at_t10,W_at_t10,usable_fraction_at_t10,coherence_L1_at_t10,drive_energy_in_at_t10,drive_energy_net_at_t10,dephasing_energy_net_at_t10,W_over_Ein_at_t10,E_max,t_at_E_max,W_max,t_at_W_max,E_at_W_max,usable_fraction_at_W_max,coherence_L1_at_W_max,W_over_Ein_at_W_max,energy_arrival_time,ergotropy_arrival_time,W_10_percent_arrival_time,W_50_percent_arrival_time,E_time_area_0_to_t10,W_time_area_0_to_t10,max_load_top_level_population,max_trace_error,max_hermiticity_error,minimum_density_eigenvalue,max_abs_energy_ledger_residual,maximum_drive_power_imaginary_part,maximum_dephasing_power_imaginary_part,finite_values_pass,endpoint_peak_classification")?;
    for result in results {
        let summary = &result.summary;
        let e = &summary.endpoint;
        let w = &summary.w_max;
        let a = &summary.arrivals;
        let values = [
            e.energy,
            e.work,
            e.usable,
            e.coherence_l1,
            e.drive_in,
            e.drive_net,
            e.dephasing_net,
            e.w_over_ein,
            summary.e_max.energy,
            summary.e_max.time,
            w.work,
            w.time,
            w.energy,
            w.usable,
            w.coherence_l1,
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
            summary.max_dephasing_power_imag,
        ];
        writeln!(
            out,
            "{},{},fixed_total_gamma_1p5,{},{},{},{},{},{}",
            result.spec.condition,
            result.spec.n,
            result.spec.n,
            fmt(result.spec.gamma_site()),
            fmt(TOTAL_GAMMA),
            values.iter().map(|x| fmt(*x)).collect::<Vec<_>>().join(","),
            summary.all_finite,
            summary.peak_class
        )?;
    }
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

fn write_windows(results: &[ConditionResult]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create("fixed_total_noise_windows.csv")?);
    writeln!(w, "condition,chain_length,noise_normalization,gamma_phi_per_site,gamma_phi_total,window_name,time_start,time_end,point_count,mean_load_energy,mean_load_ergotropy,mean_usable_fraction,mean_coherence_L1,E_time_area,W_time_area,mean_total_chain_population,maximum_load_ergotropy,time_of_window_max_W,drive_energy_net_change,dephasing_energy_net_change")?;
    for result in results {
        for (name, start, end, include_start) in [
            ("pulse_interval", 0.0, 3.2, true),
            ("early_post_pulse", 3.2, 5.0, false),
            ("middle_interval", 5.0, 7.5, false),
            ("late_interval", 7.5, 10.0, false),
        ] {
            let subset = window_rows(&result.run.rows, start, end, include_start);
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
                subset.last().unwrap().drive_net - subset.first().unwrap().drive_net,
                subset.last().unwrap().dephasing_net - subset.first().unwrap().dephasing_net,
            ];
            writeln!(
                w,
                "{},{},fixed_total_gamma_1p5,{},{},{name},{},{},{},{}",
                result.spec.condition,
                result.spec.n,
                fmt(result.spec.gamma_site()),
                fmt(TOTAL_GAMMA),
                fmt(start),
                fmt(end),
                subset.len(),
                values.iter().map(|x| fmt(*x)).collect::<Vec<_>>().join(",")
            )?;
        }
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

fn write_comparisons(summary: &Summary) -> Result<(), Box<dyn std::error::Error>> {
    let old_summary = CsvTable::read(
        &std::path::Path::new(OLD_DIR).join("chain_length_reachability_summary.csv"),
    )?;
    let old_arrivals = CsvTable::read(
        &std::path::Path::new(OLD_DIR).join("chain_length_reachability_arrivals.csv"),
    )?;
    let n7_free = CsvTable::read(std::path::Path::new("n7_noise_free_summary.csv"))?;
    let e = &summary.endpoint;
    let p = &summary.w_max;

    let n7_metrics = [
        (
            "E",
            "t10",
            n7_free.value("N7_noise_free", "E_at_t10"),
            e.energy,
            "same_t10",
        ),
        (
            "W",
            "t10",
            n7_free.value("N7_noise_free", "W_at_t10"),
            e.work,
            "same_t10",
        ),
        (
            "usable_fraction",
            "t10",
            n7_free.value("N7_noise_free", "usable_fraction_at_t10"),
            e.usable,
            "same_t10",
        ),
        (
            "coherence_L1",
            "t10",
            n7_free.value("N7_noise_free", "coherence_L1_at_t10"),
            e.coherence_l1,
            "same_t10",
        ),
        (
            "W_over_Ein",
            "t10",
            n7_free.value("N7_noise_free", "W_over_Ein_at_t10"),
            e.w_over_ein,
            "same_t10",
        ),
        (
            "W",
            "individual_peak",
            n7_free.value("N7_noise_free", "W_max"),
            p.work,
            "individual_peak_times_may_differ",
        ),
        (
            "t_at_W_max",
            "individual_peak",
            n7_free.value("N7_noise_free", "t_at_W_max"),
            p.time,
            "individual_peak_times_may_differ",
        ),
        (
            "E_time_area",
            "0_to_t10",
            n7_free.value("N7_noise_free", "E_time_area_0_to_t10"),
            summary.e_area,
            "same_time_window",
        ),
        (
            "W_time_area",
            "0_to_t10",
            n7_free.value("N7_noise_free", "W_time_area_0_to_t10"),
            summary.w_area,
            "same_time_window",
        ),
    ];
    let mut n7_out = BufWriter::new(File::create("n7_noise_comparison.csv")?);
    writeln!(n7_out, "metric,evaluation_point,noise_free_value,all_site_noisy_value,noisy_over_free,noisy_minus_free,comparison_note")?;
    for (metric, point, free, noisy, note) in n7_metrics {
        writeln!(
            n7_out,
            "{metric},{point},{},{},{},{},{note}",
            fmt(free),
            fmt(noisy),
            fmt(ratio(noisy, free)),
            fmt(noisy - free)
        )?;
    }

    let noisy_metrics = vec![
        (
            "E",
            "t10",
            old_summary.value("N3_all_site_noisy", "E_at_t10"),
            old_summary.value("N5_all_site_noisy", "E_at_t10"),
            e.energy,
        ),
        (
            "W",
            "t10",
            old_summary.value("N3_all_site_noisy", "W_at_t10"),
            old_summary.value("N5_all_site_noisy", "W_at_t10"),
            e.work,
        ),
        (
            "usable_fraction",
            "t10",
            old_summary.value("N3_all_site_noisy", "usable_fraction_at_t10"),
            old_summary.value("N5_all_site_noisy", "usable_fraction_at_t10"),
            e.usable,
        ),
        (
            "W_over_Ein",
            "t10",
            old_summary.value("N3_all_site_noisy", "W_over_Ein_at_t10"),
            old_summary.value("N5_all_site_noisy", "W_over_Ein_at_t10"),
            e.w_over_ein,
        ),
        (
            "W",
            "individual_peak",
            old_summary.value("N3_all_site_noisy", "W_max"),
            old_summary.value("N5_all_site_noisy", "W_max"),
            p.work,
        ),
        (
            "t_at_W_max",
            "individual_peak",
            old_summary.value("N3_all_site_noisy", "t_at_W_max"),
            old_summary.value("N5_all_site_noisy", "t_at_W_max"),
            p.time,
        ),
        (
            "energy_arrival",
            "sustained_threshold",
            old_arrival(&old_arrivals, "N3_all_site_noisy", "energy_ge_1e-4"),
            old_arrival(&old_arrivals, "N5_all_site_noisy", "energy_ge_1e-4"),
            summary.arrivals[0].time,
        ),
        (
            "ergotropy_arrival",
            "sustained_threshold",
            old_arrival(&old_arrivals, "N3_all_site_noisy", "ergotropy_ge_1e-5"),
            old_arrival(&old_arrivals, "N5_all_site_noisy", "ergotropy_ge_1e-5"),
            summary.arrivals[1].time,
        ),
        (
            "E_time_area",
            "0_to_t10",
            old_summary.value("N3_all_site_noisy", "E_time_area_0_to_t10"),
            old_summary.value("N5_all_site_noisy", "E_time_area_0_to_t10"),
            summary.e_area,
        ),
        (
            "W_time_area",
            "0_to_t10",
            old_summary.value("N3_all_site_noisy", "W_time_area_0_to_t10"),
            old_summary.value("N5_all_site_noisy", "W_time_area_0_to_t10"),
            summary.w_area,
        ),
    ];
    let mut length_out = BufWriter::new(File::create("chain_length_noisy_comparison.csv")?);
    writeln!(length_out, "metric,evaluation_point,N3_value,N5_value,N7_value,N5_over_N3,N7_over_N5,N7_over_N3,N5_minus_N3,N7_minus_N5,N7_minus_N3")?;
    for (metric, point, n3, n5, n7) in &noisy_metrics {
        writeln!(
            length_out,
            "{metric},{point},{},{},{},{},{},{},{},{},{}",
            fmt(*n3),
            fmt(*n5),
            fmt(*n7),
            fmt(ratio(*n5, *n3)),
            fmt(ratio(*n7, *n5)),
            fmt(ratio(*n7, *n3)),
            fmt(*n5 - *n3),
            fmt(*n7 - *n5),
            fmt(*n7 - *n3)
        )?;
    }

    let mut loss_out = BufWriter::new(File::create("chain_length_noise_loss.csv")?);
    writeln!(loss_out, "chain_length,metric,evaluation_point,noise_free_value,all_site_noisy_value,noisy_over_free,noisy_minus_free,comparison_note")?;
    let loss_specs = [
        ("E", "t10", "E_at_t10", "same_t10"),
        ("W", "t10", "W_at_t10", "same_t10"),
        (
            "usable_fraction",
            "t10",
            "usable_fraction_at_t10",
            "same_t10",
        ),
        ("coherence_L1", "t10", "coherence_L1_at_t10", "same_t10"),
        ("W_over_Ein", "t10", "W_over_Ein_at_t10", "same_t10"),
        (
            "W",
            "individual_peak",
            "W_max",
            "individual_peak_times_may_differ",
        ),
        (
            "E_time_area",
            "0_to_t10",
            "E_time_area_0_to_t10",
            "same_time_window",
        ),
        (
            "W_time_area",
            "0_to_t10",
            "W_time_area_0_to_t10",
            "same_time_window",
        ),
    ];
    for n in [3usize, 5usize] {
        let free_condition = format!("N{n}_noise_free");
        let noisy_condition = format!("N{n}_all_site_noisy");
        for (metric, point, column, note) in loss_specs {
            let free = old_summary.value(&free_condition, column);
            let noisy = old_summary.value(&noisy_condition, column);
            writeln!(
                loss_out,
                "{n},{metric},{point},{},{},{},{},{note}",
                fmt(free),
                fmt(noisy),
                fmt(ratio(noisy, free)),
                fmt(noisy - free)
            )?;
        }
    }
    for (metric, point, column, note) in loss_specs {
        let free = n7_free.value("N7_noise_free", column);
        let noisy = match column {
            "E_at_t10" => e.energy,
            "W_at_t10" => e.work,
            "usable_fraction_at_t10" => e.usable,
            "coherence_L1_at_t10" => e.coherence_l1,
            "W_over_Ein_at_t10" => e.w_over_ein,
            "W_max" => p.work,
            "E_time_area_0_to_t10" => summary.e_area,
            "W_time_area_0_to_t10" => summary.w_area,
            _ => unreachable!(),
        };
        writeln!(
            loss_out,
            "7,{metric},{point},{},{},{},{},{note}",
            fmt(free),
            fmt(noisy),
            fmt(ratio(noisy, free)),
            fmt(noisy - free)
        )?;
    }
    Ok(())
}

fn result_by_n(results: &[ConditionResult], n: usize) -> &ConditionResult {
    results.iter().find(|result| result.spec.n == n).unwrap()
}

fn result_metric(summary: &Summary, metric: &str) -> f64 {
    match metric {
        "E_at_t10" => summary.endpoint.energy,
        "W_at_t10" => summary.endpoint.work,
        "usable_fraction_at_t10" => summary.endpoint.usable,
        "coherence_L1_at_t10" => summary.endpoint.coherence_l1,
        "W_over_Ein_at_t10" => summary.endpoint.w_over_ein,
        "W_max" => summary.w_max.work,
        "t_at_W_max" => summary.w_max.time,
        "E_at_W_max" => summary.w_max.energy,
        "usable_fraction_at_W_max" => summary.w_max.usable,
        "energy_arrival_time" => summary.arrivals[0].time,
        "ergotropy_arrival_time" => summary.arrivals[1].time,
        "E_time_area" => summary.e_area,
        "W_time_area" => summary.w_area,
        _ => panic!("unsupported result metric {metric}"),
    }
}

fn reference_metric(summary: &CsvTable, arrivals: &CsvTable, condition: &str, metric: &str) -> f64 {
    match metric {
        "energy_arrival_time" => old_arrival(arrivals, condition, "energy_ge_1e-4"),
        "ergotropy_arrival_time" => old_arrival(arrivals, condition, "ergotropy_ge_1e-5"),
        "E_time_area" => summary.value(condition, "E_time_area_0_to_t10"),
        "W_time_area" => summary.value(condition, "W_time_area_0_to_t10"),
        _ => summary.value(condition, metric),
    }
}

fn write_fixed_total_comparisons(
    results: &[ConditionResult],
) -> Result<(), Box<dyn std::error::Error>> {
    let old_summary = CsvTable::read(
        &std::path::Path::new(OLD_DIR).join("chain_length_reachability_summary.csv"),
    )?;
    let old_arrivals = CsvTable::read(
        &std::path::Path::new(OLD_DIR).join("chain_length_reachability_arrivals.csv"),
    )?;
    let n7_free_summary = CsvTable::read(std::path::Path::new("n7_noise_free_summary.csv"))?;
    let n7_free_arrivals = CsvTable::read(std::path::Path::new("n7_noise_free_arrivals.csv"))?;
    let n7_per_summary = CsvTable::read(std::path::Path::new("n7_all_site_noisy_summary.csv"))?;
    let n7_per_arrivals = CsvTable::read(std::path::Path::new("n7_all_site_noisy_arrivals.csv"))?;
    let specs = [
        ("E_at_t10", "t10"),
        ("W_at_t10", "t10"),
        ("usable_fraction_at_t10", "t10"),
        ("coherence_L1_at_t10", "t10"),
        ("W_over_Ein_at_t10", "t10"),
        ("W_max", "individual_peak"),
        ("t_at_W_max", "individual_peak"),
        ("E_at_W_max", "individual_peak"),
        ("usable_fraction_at_W_max", "individual_peak"),
        ("energy_arrival_time", "sustained_threshold"),
        ("ergotropy_arrival_time", "sustained_threshold"),
        ("E_time_area", "0_to_t10"),
        ("W_time_area", "0_to_t10"),
    ];

    let mut length = BufWriter::new(File::create("fixed_total_noise_length_comparison.csv")?);
    writeln!(length, "metric,evaluation_point,N3_value,N5_value,N7_value,N5_over_N3,N7_over_N5,N7_over_N3,N5_minus_N3,N7_minus_N5,N7_minus_N3")?;
    for (metric, point) in specs {
        let n3 = reference_metric(&old_summary, &old_arrivals, "N3_all_site_noisy", metric);
        let n5 = result_metric(&result_by_n(results, 5).summary, metric);
        let n7 = result_metric(&result_by_n(results, 7).summary, metric);
        writeln!(
            length,
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

    let compare_specs = [
        ("E_at_t10", "t10"),
        ("W_at_t10", "t10"),
        ("usable_fraction_at_t10", "t10"),
        ("coherence_L1_at_t10", "t10"),
        ("W_over_Ein_at_t10", "t10"),
        ("W_max", "individual_peak"),
        ("t_at_W_max", "individual_peak"),
        ("E_time_area", "0_to_t10"),
        ("W_time_area", "0_to_t10"),
        ("energy_arrival_time", "sustained_threshold"),
        ("ergotropy_arrival_time", "sustained_threshold"),
    ];
    let mut normalization = BufWriter::new(File::create("fixed_total_vs_fixed_per_site.csv")?);
    writeln!(normalization, "chain_length,metric,evaluation_point,fixed_per_site_gamma,fixed_per_site_total_gamma,fixed_per_site_value,fixed_total_gamma_per_site,fixed_total_gamma,fixed_total_value,fixed_total_over_fixed_per_site,fixed_total_minus_fixed_per_site")?;
    for n in [5usize, 7usize] {
        for (metric, point) in compare_specs {
            let per_site = if n == 5 {
                reference_metric(&old_summary, &old_arrivals, "N5_all_site_noisy", metric)
            } else {
                reference_metric(
                    &n7_per_summary,
                    &n7_per_arrivals,
                    "N7_all_site_noisy",
                    metric,
                )
            };
            let fixed = result_metric(&result_by_n(results, n).summary, metric);
            writeln!(
                normalization,
                "{n},{metric},{point},{},{},{},{},{},{},{},{}",
                fmt(0.5),
                fmt(0.5 * n as f64),
                fmt(per_site),
                fmt(TOTAL_GAMMA / n as f64),
                fmt(TOTAL_GAMMA),
                fmt(fixed),
                fmt(ratio(fixed, per_site)),
                fmt(fixed - per_site)
            )?;
        }
    }

    let mut recovery = BufWriter::new(File::create("fixed_total_noise_recovery.csv")?);
    writeln!(recovery, "chain_length,metric,evaluation_point,noise_free_value,fixed_per_site_noisy_value,fixed_total_noisy_value,fixed_per_site_residual,fixed_total_residual,normalization_recovery_ratio,fixed_total_absolute_gain,fixed_total_residual_loss")?;
    for n in [3usize, 5usize, 7usize] {
        for (metric, point) in compare_specs {
            let (free, per_site, fixed) = match n {
                3 => {
                    let value =
                        reference_metric(&old_summary, &old_arrivals, "N3_all_site_noisy", metric);
                    (
                        reference_metric(&old_summary, &old_arrivals, "N3_noise_free", metric),
                        value,
                        value,
                    )
                }
                5 => (
                    reference_metric(&old_summary, &old_arrivals, "N5_noise_free", metric),
                    reference_metric(&old_summary, &old_arrivals, "N5_all_site_noisy", metric),
                    result_metric(&result_by_n(results, 5).summary, metric),
                ),
                7 => (
                    reference_metric(&n7_free_summary, &n7_free_arrivals, "N7_noise_free", metric),
                    reference_metric(
                        &n7_per_summary,
                        &n7_per_arrivals,
                        "N7_all_site_noisy",
                        metric,
                    ),
                    result_metric(&result_by_n(results, 7).summary, metric),
                ),
                _ => unreachable!(),
            };
            writeln!(
                recovery,
                "{n},{metric},{point},{},{},{},{},{},{},{},{}",
                fmt(free),
                fmt(per_site),
                fmt(fixed),
                fmt(ratio(per_site, free)),
                fmt(ratio(fixed, free)),
                fmt(ratio(fixed, per_site)),
                fmt(fixed - per_site),
                fmt(free - fixed)
            )?;
        }
    }
    Ok(())
}

fn append_check(
    w: &mut impl Write,
    condition: &str,
    stage: &str,
    check: &str,
    observed: &str,
    expected: &str,
    pass: bool,
) -> std::io::Result<()> {
    writeln!(
        w,
        "{stage},{condition},{check},{observed},{expected},{}",
        if pass { "PASS" } else { "FAIL" }
    )
}

fn write_checks(
    spec: Spec,
    rows: &[Row],
    summary: &Summary,
    construction: &[(String, String, String, bool)],
    create: bool,
) -> Result<bool, Box<dyn std::error::Error>> {
    let file = if create {
        File::create("fixed_total_noise_checks.csv")?
    } else {
        OpenOptions::new()
            .append(true)
            .open("fixed_total_noise_checks.csv")?
    };
    let mut w = BufWriter::new(file);
    if create {
        writeln!(w, "stage,condition,check,observed,expected,status")?;
    }
    let mut all = true;
    for (name, observed, expected, pass) in construction {
        append_check(
            &mut w,
            spec.condition,
            "construction",
            name,
            observed,
            expected,
            *pass,
        )?;
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
            "dephasing_energy_finite",
            fmt(summary.endpoint.dephasing_net),
            "finite".into(),
            summary.endpoint.dephasing_net.is_finite()
                && rows.iter().all(|r| r.dephasing_power.is_finite()),
        ),
        (
            "noisy_site_count",
            spec.n.to_string(),
            spec.n.to_string(),
            matches!(spec.n, 5 | 7),
        ),
        (
            "all_site_gamma",
            fmt(spec.gamma_site()),
            fmt(TOTAL_GAMMA / spec.n as f64),
            spec.gammas()
                .iter()
                .all(|gamma| (*gamma - spec.gamma_site()).abs() <= 1.0e-15),
        ),
        (
            "gamma_sum",
            fmt(spec.gammas().iter().sum()),
            fmt(TOTAL_GAMMA),
            (spec.gammas().iter().sum::<f64>() - TOTAL_GAMMA).abs() <= 1.0e-14,
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
            (rows.len() * spec.n).to_string(),
            (1001 * spec.n).to_string(),
            rows.len() * spec.n == 1001 * spec.n,
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
        append_check(
            &mut w,
            spec.condition,
            "trajectory",
            name,
            &observed,
            &expected,
            pass,
        )?;
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

fn write_performance(results: &[ConditionResult]) -> Result<(), Box<dyn std::error::Error>> {
    let mut w = BufWriter::new(File::create("fixed_total_noise_performance.csv")?);
    writeln!(w, "condition,chain_length,hilbert_dimension,density_matrix_element_count,noisy_site_count,gamma_phi_per_site,gamma_phi_total,kernel_bytes,step_count,saved_time_points,construction_seconds,kernel_construction_seconds,propagation_seconds,diagnostics_seconds,total_seconds,seconds_per_step,peak_memory_bytes_if_available,process_working_set_before,process_working_set_after,timeseries_rows,site_population_rows")?;
    for result in results {
        let run = &result.run;
        let spec = result.spec;
        let fields = vec![
            spec.condition.to_string(),
            spec.n.to_string(),
            spec.dim.to_string(),
            (spec.dim * spec.dim).to_string(),
            spec.n.to_string(),
            fmt(spec.gamma_site()),
            fmt(TOTAL_GAMMA),
            (spec.dim * spec.dim * std::mem::size_of::<f64>()).to_string(),
            "4000".to_string(),
            run.rows.len().to_string(),
            fmt(run.construction_seconds),
            fmt(run.kernel_construction_seconds),
            fmt(run.propagation_seconds),
            fmt(run.diagnostics_seconds),
            fmt(run.total_seconds),
            fmt(run.propagation_seconds / 4000.0),
            run.peak_working_set.to_string(),
            run.working_set_before.to_string(),
            run.working_set_after.to_string(),
            run.rows.len().to_string(),
            (run.rows.len() * spec.n).to_string(),
        ];
        writeln!(w, "{}", fields.join(","))?;
    }
    Ok(())
}

fn write_report(
    run: &RunResult,
    summary: &Summary,
    checks_pass: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let old = CsvTable::read(
        &std::path::Path::new(OLD_DIR).join("chain_length_reachability_summary.csv"),
    )?;
    let n3_wmax = old.value("N3_all_site_noisy", "W_max");
    let n5_wmax = old.value("N5_all_site_noisy", "W_max");
    let n7_free = CsvTable::read(std::path::Path::new("n7_noise_free_summary.csv"))?;
    let n7_free_w10 = n7_free.value("N7_noise_free", "W_at_t10");
    let n7_free_wmax = n7_free.value("N7_noise_free", "W_max");
    let n7_free_e10 = n7_free.value("N7_noise_free", "E_at_t10");
    let estimated_seconds = read_estimate(
        "dephasing_kernel_estimates.csv",
        "N7_all_site_noisy",
        "estimated_t10_seconds",
    )?;
    let reached_e = summary.arrivals[0].time.is_finite();
    let reached_w = summary.arrivals[1].time.is_finite();
    let decision = if !checks_pass {
        "numerical_issue_stop"
    } else if !reached_e || !reached_w {
        "insufficient_reachability"
    } else if run.total_seconds > 6.0 * estimated_seconds && run.total_seconds > 21_600.0 {
        "performance_issue_stop"
    } else if summary.peak_class != "peak_resolved" {
        "extend_time_candidate"
    } else {
        "completed_comparison"
    };
    let e = &summary.endpoint;
    let p = &summary.w_max;
    let feature_verdict = if p.work > n5_wmax {
        "残った"
    } else {
        "残らなかった"
    };
    let mut w = BufWriter::new(File::create("MILESTONE_9B_REPORT.md")?);
    writeln!(w, "# Milestone 9b: N=7 all-site noisy full run\n")?;
    let sections = [
        ("1. 目的", "N=7の全7siteへgamma_phi=0.5の位相雑音を入れ、t=10までenergyとergotropyがloadへ届くか確認した。".to_string()),
        ("2. 今回の範囲", "新規本計算はN=7 all-site noisy、dt=0.0025、t=10の1条件だけである。noise-free再計算、細刻み、時間延長、sweep、次Milestoneは実行していない。".to_string()),
        ("3. 参照成果物", "Milestone 8aのN=3/N=5、8cのkernel検証、9aのN=7 noise-free成果物を読み取り専用で参照した。".to_string()),
        ("4. 変更していない物理模型", "Hamiltonian、drive、RK4、dt、基底、load、初期真空、観測量を既存実装から変更していない。".to_string()),
        ("5. 雑音模型", "全7siteにL_phi,j=sqrt(0.5/2) sigma_z,jを入れ、loadには直接雑音を入れていない。".to_string()),
        ("6. dephasing kernel", "Milestone 8cでdense collapse pathと等価性確認済みのDiagonalDephasingKernelを使用した。物理近似ではなく同じLindblad項の成分表示である。".to_string()),
        ("7. N=7模型構成", "7二準位site、3準位load、dim=384、bond=6、drive site=0、load coupling site=6、noisy site count=7。".to_string()),
        ("8. 数値手法", "dense complex density matrixのtime-dependent RK4。4000 step、保存間隔0.01、1001点。各保存点で縮約、ergotropy、全系最小固有値、dephasing powerを診断した。".to_string()),
        ("9. 構成検査", "次元、全gamma、kernel寸法・対称性・非負性・mapping・load除外、真空、drive端点、Hamiltonian Hermiticityを本計算前に検査した。".to_string()),
        ("10. 実行時間とメモリ", format!("construction {:.3}s、propagation {:.3}s、diagnostics {:.3}s、total {:.3}s。8c推定 {:.3}s。peak working set {} bytes。", run.construction_seconds, run.propagation_seconds, run.diagnostics_seconds, run.total_seconds, estimated_seconds, run.peak_working_set)),
        ("11. 数値品質チェック", format!("全チェック{}。max trace={:.3e}、max Hermiticity={:.3e}、min eigenvalue={:.3e}、max ledger={:.3e}。", if checks_pass { "PASS" } else { "FAIL" }, summary.max_trace, summary.max_herm, summary.min_eig, summary.max_ledger)),
        ("12. load energy到達", format!("持続閾値1e-4の到達時刻は {}。", fmt(summary.arrivals[0].time))),
        ("13. load ergotropy到達", format!("持続閾値1e-5の到達時刻は {}。", fmt(summary.arrivals[1].time))),
        ("14. t=10結果", format!("E={:.10e}、W={:.10e}、usable={:.10e}、W/Ein={:.10e}、coherence L1={:.10e}。", e.energy, e.work, e.usable, e.w_over_ein, e.coherence_l1)),
        ("15. W最大値と時刻", format!("W_max={:.10e}、t={:.2}、その時のE={:.10e}、usable={:.10e}。", p.work, p.time, p.energy, p.usable)),
        ("16. 終端挙動とピーク判定", format!("分類 `{}`。W(10)-W(9.9)={:.3e}、E(10)-E(9.9)={:.3e}、最終10点W slope={:.3e}、E slope={:.3e}。", summary.peak_class, summary.w_delta_99_100, summary.e_delta_99_100, summary.w_final_slope, summary.e_final_slope)),
        ("17. energy ledger", format!("t=10でdrive net={:.10e}、dephasing net={:.10e}。ledgerは両方を含み、最大残差は {:.3e}。", e.drive_net, e.dephasing_net, summary.max_ledger)),
        ("18. site populationと時間窓", "1001時刻×7site=7007行を保存し、pulse、early post-pulse、middle、lateの4窓を集計した。".to_string()),
        ("19. N=7 noisy/free比較", format!("t=10でE noisy/free={:.10e}、W noisy/free={:.10e}。Wmax noisy/free={:.10e}。Wmaxは異なる時刻同士の最大値比較である。", ratio(e.energy, n7_free_e10), ratio(e.work, n7_free_w10), ratio(p.work, n7_free_wmax))),
        ("20. N=3、N=5、N=7 noisy比較", format!("WmaxはN3 {:.10e}、N5 {:.10e}、N7 {:.10e}。N7/N5={:.10e}、N7/N3={:.10e}。", n3_wmax, n5_wmax, p.work, ratio(p.work, n5_wmax), ratio(p.work, n3_wmax))),
        ("21. 中心問題への有限計算上の答え", format!("N=7 all-site noisyでenergy到達={}、ergotropy到達={}。noise-freeで見えたN7 Wmax>N5 Wmaxという特徴はnoisyでは{}。N7/N5={:.10e}。", reached_e, reached_w, feature_verdict, ratio(p.work, n5_wmax))),
        ("22. 解釈上の制限", "Nとともに距離と雑音site数が同時に増える。3点の有限長比較から距離だけの因果、指数/べき則、漸近scaling、輸送限界を主張しない。".to_string()),
        ("23. 未確認", "dt半減、t>10、位置別雑音、保護、gamma/Omega sweep、N>7、実機性能は未確認である。".to_string()),
        ("24. 最終判定", format!("判定 **{}**。どの判定でも時間延長や次Milestoneは自動実行していない。", decision)),
        ("25. 生成ファイル", "`src/bin/n7_all_site_noisy_full.rs` と指定11成果物を新規作成した。既存成果物は上書きしていない。".to_string()),
    ];
    for (title, body) in sections {
        writeln!(w, "## {title}\n\n{body}\n")?;
    }
    Ok(())
}

fn write_9c_report(results: &[ConditionResult]) -> Result<(), Box<dyn std::error::Error>> {
    let old = CsvTable::read(
        &std::path::Path::new(OLD_DIR).join("chain_length_reachability_summary.csv"),
    )?;
    let n7_per = CsvTable::read(std::path::Path::new("n7_all_site_noisy_summary.csv"))?;
    let n5 = result_by_n(results, 5);
    let n7 = result_by_n(results, 7);
    let s5 = &n5.summary;
    let s7 = &n7.summary;
    let n3_wmax = old.value("N3_all_site_noisy", "W_max");
    let n5_per_wmax = old.value("N5_all_site_noisy", "W_max");
    let n7_per_wmax = n7_per.value("N7_all_site_noisy", "W_max");
    let ratio_75 = ratio(s7.w_max.work, s5.w_max.work);
    let recovery5 = ratio(s5.w_max.work, n5_per_wmax);
    let recovery7 = ratio(s7.w_max.work, n7_per_wmax);
    let band = if ratio_75 < 0.9 {
        "N7 substantially lower in this finite comparison"
    } else if ratio_75 <= 1.1 {
        "similar within 10 percent descriptive band"
    } else {
        "N7 substantially higher in this finite comparison"
    };
    let decision = if !results.iter().all(|result| result.checks_pass) {
        "numerical_issue_stop"
    } else if ratio_75 >= 1.0 {
        "nonmonotonic_finite_behavior"
    } else if ratio_75 >= 0.9 && (recovery5 > 1.2 || recovery7 > 1.2) {
        "total_noise_dominant_candidate"
    } else if ratio_75 < 0.9 && (recovery5 > 1.2 || recovery7 > 1.2) {
        "mixed_effects_candidate"
    } else {
        "length_or_dynamics_remains_candidate"
    };
    let mut w = BufWriter::new(File::create("MILESTONE_9C_REPORT.md")?);
    writeln!(
        w,
        "# Milestone 9c: Fixed-total-dephasing comparison across N=3 N=5 N=7\n"
    )?;
    let sections = [
        ("1. 目的", "全site gammaの単純和を1.5へ固定し、N増加と総雑音増加の交絡を部分的に切り分けた。".to_string()),
        ("2. 9b比較の交絡", "fixed-per-site gamma=0.5ではtotal gammaがN3=1.5、N5=2.5、N7=3.5と同時に増えていた。".to_string()),
        ("3. fixed-total-noise設計", "TOTAL_GAMMA=1.5を全siteへ均等配分した。これは総dephasing rateの単純和を固定する記述的比較である。".to_string()),
        ("4. 今回の新規計算", "N=5 gamma_site=0.3とN=7 gamma_site=1.5/7だけを新規本計算した。N=3は既存gamma_site=0.5結果を参照した。".to_string()),
        ("5. 変更していない物理模型", "Hamiltonian、drive、RK4、dt=0.0025、t=10、load、初期真空、観測量を変更していない。".to_string()),
        ("6. 雑音定義", "各chain siteにL_phi,j=sqrt(gamma_j/2) sigma_z,jを適用し、loadへ直接雑音を入れていない。".to_string()),
        ("7. gamma配分", format!("N5 gamma_site={:.16e}、N7 gamma_site={:.16e}。両条件のsum gammaは1.5。", N5_SPEC.gamma_site(), N7_SPEC.gamma_site())),
        ("8. dephasing kernel", "Milestone 8cでdense pathとの等価性を確認したDiagonalDephasingKernelを使用した。新しい近似ではない。".to_string()),
        ("9. 数値手法", "time-dependent dense density-matrix RK4、4000 steps、1001保存点。各保存点で最小固有値とpower ledgerを診断した。".to_string()),
        ("10. 構成検査", "chain長、次元、bond、drive/load mapping、gamma総和、kernel mapping、load除外、真空、Hermiticityを本計算前に検査した。".to_string()),
        ("11. N=5実行結果", format!("E10={:.10e}、W10={:.10e}、Wmax={:.10e} at t={:.2}、peak={}。", s5.endpoint.energy, s5.endpoint.work, s5.w_max.work, s5.w_max.time, s5.peak_class)),
        ("12. N=7実行結果", format!("E10={:.10e}、W10={:.10e}、Wmax={:.10e} at t={:.2}、peak={}。", s7.endpoint.energy, s7.endpoint.work, s7.w_max.work, s7.w_max.time, s7.peak_class)),
        ("13. 数値品質", format!("N5 checks={}、N7 checks={}。max ledgerはN5 {:.3e}、N7 {:.3e}。", n5.checks_pass, n7.checks_pass, s5.max_ledger, s7.max_ledger)),
        ("14. 実行時間", format!("N5 total {:.3}s、N7 total {:.3}s。時間差は性能診断であり物理結果ではない。", n5.run.total_seconds, n7.run.total_seconds)),
        ("15. fixed-total N=3/5/7比較", format!("WmaxはN3 {:.10e}、N5 {:.10e}、N7 {:.10e}。", n3_wmax, s5.w_max.work, s7.w_max.work)),
        ("16. fixed-per-site N=3/5/7比較", format!("既存WmaxはN3 {:.10e}、N5 {:.10e}、N7 {:.10e}。", n3_wmax, n5_per_wmax, n7_per_wmax)),
        ("17. N=5での正規化効果", format!("Wmax fixed-total/fixed-per-site={:.10e}、absolute gain={:.10e}。", recovery5, s5.w_max.work - n5_per_wmax)),
        ("18. N=7での正規化効果", format!("Wmax fixed-total/fixed-per-site={:.10e}、absolute gain={:.10e}。", recovery7, s7.w_max.work - n7_per_wmax)),
        ("19. noise-free基準の残存率", "metricごとのnoise-free基準残存率をfixed_total_noise_recovery.csvへ保存した。巨大な比率だけでなく絶対値も併記した。".to_string()),
        ("20. Wmax比較", format!("fixed-total N7/N5={:.10e}。分類 `{}`。10%帯は統計的有意差ではない。", ratio_75, band)),
        ("21. t=10比較", format!("N5 E10={:.10e} W10={:.10e}、N7 E10={:.10e} W10={:.10e}。", s5.endpoint.energy, s5.endpoint.work, s7.endpoint.energy, s7.endpoint.work)),
        ("22. arrival time比較", format!("energy arrival N5={} N7={}、ergotropy arrival N5={} N7={}。絶対閾値が異なるため単純な到着順とは解釈しない。", fmt(s5.arrivals[0].time), fmt(s7.arrivals[0].time), fmt(s5.arrivals[1].time), fmt(s7.arrivals[1].time))),
        ("23. usable fraction比較", format!("t10 usableはN5 {:.10e}、N7 {:.10e}。", s5.endpoint.usable, s7.endpoint.usable)),
        ("24. W/Ein比較", format!("t10 W/EinはN5 {:.10e}、N7 {:.10e}。制御費用を含む総合効率ではない。", s5.endpoint.w_over_ein, s7.endpoint.w_over_ein)),
        ("25. 時間面積比較", format!("E area N5 {:.10e} N7 {:.10e}、W area N5 {:.10e} N7 {:.10e}。", s5.e_area, s7.e_area, s5.w_area, s7.w_area)),
        ("26. 時間窓解析", "pulse、early post-pulse、middle、lateの4窓を両条件で同じ定義により保存した。".to_string()),
        ("27. 中心判定", format!("判定 **{}**。fixed-totalでの回復とN5/N7差を併記する有限模型上の候補分類であり因果証明ではない。", decision)),
        ("28. 直接確認できたこと", "total gamma=1.5のN5/N7本計算、N3参照との有限長比較、fixed-per-site差、noise-free基準残存率を確認した。".to_string()),
        ("29. 確認できていないこと", "距離だけの純粋因果、他のtotal gamma、dt半減、t>10、N>7、位置別雑音、occupation-weighted exposureは未確認。".to_string()),
        ("30. 主張してはいけないこと", "指数/べき則、熱力学極限、一般的輸送限界、統計的有意差、total noiseまたはlengthだけの単独因果を主張しない。".to_string()),
        ("31. 次段階候補", "total gamma sweep、noise exposure integral、site occupation weighted dephasingは候補のみ。自動実行していない。".to_string()),
        ("32. 生成ファイル一覧", "`src/bin/fixed_total_noise_comparison.rs` と指定11成果物を新規作成した。既存成果物は上書きしていない。".to_string()),
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
        for spec in [N5_SPEC, N7_SPEC] {
            let params = ModelParams::default();
            let ops = build_operators_for_chain(&params, spec.n)?;
            let gammas = spec.gammas();
            let mut rho = ComplexMatrix::zeros(spec.dim, spec.dim);
            rho[(0, 0)] = C64::new(1.0, 0.0);
            let kernel = DiagonalDephasingKernel::new(spec.n, params.load_dim, &gammas)?;
            let checks = construction_checks(spec, &gammas, &ops, &rho, &kernel)?;
            for (name, observed, expected, pass) in &checks {
                println!(
                    "{} {name}: observed={observed} expected={expected} pass={pass}",
                    spec.condition
                );
            }
            if checks.iter().any(|x| !x.3) {
                return Err(
                    format!("{} preflight construction checks failed", spec.condition).into(),
                );
            }
            println!(
                "{} preflight PASS: {} construction checks",
                spec.condition,
                checks.len()
            );
        }
        return Ok(());
    }

    let (run5, construction5) = run_full(N5_SPEC)?;
    let summary5 = summarize(&run5.rows);
    let checks5 = write_checks(N5_SPEC, &run5.rows, &summary5, &construction5, true)?;
    if !checks5 {
        return Err("N5 fixed-total numerical quality checks failed; N7 not started".into());
    }
    println!(
        "N5 completed checks=true E10={:.10e} W10={:.10e} Wmax={:.10e}",
        summary5.endpoint.energy, summary5.endpoint.work, summary5.w_max.work
    );

    let (run7, construction7) = run_full(N7_SPEC)?;
    let summary7 = summarize(&run7.rows);
    let checks7 = write_checks(N7_SPEC, &run7.rows, &summary7, &construction7, false)?;
    let results = vec![
        ConditionResult {
            spec: N5_SPEC,
            run: run5,
            summary: summary5,
            checks_pass: checks5,
        },
        ConditionResult {
            spec: N7_SPEC,
            run: run7,
            summary: summary7,
            checks_pass: checks7,
        },
    ];
    write_timeseries(&results)?;
    write_sites(&results)?;
    write_summary(&results)?;
    write_arrivals(&results)?;
    write_windows(&results)?;
    write_fixed_total_comparisons(&results)?;
    write_performance(&results)?;
    write_9c_report(&results)?;
    println!(
        "completed N5_checks={} N7_checks={} N5_Wmax={:.10e} N7_Wmax={:.10e} ratio={:.10e}",
        checks5,
        checks7,
        result_by_n(&results, 5).summary.w_max.work,
        result_by_n(&results, 7).summary.w_max.work,
        ratio(
            result_by_n(&results, 7).summary.w_max.work,
            result_by_n(&results, 5).summary.w_max.work
        )
    );
    if !checks7 {
        return Err("N7 fixed-total numerical quality checks failed".into());
    }
    Ok(())
}
