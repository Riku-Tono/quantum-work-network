use std::collections::HashMap;
use std::env;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

use nalgebra::linalg::Schur;
use quantum_work_network::coherent_drive::{drive_hamiltonian, CoherentDriveConfig};
use quantum_work_network::dephasing_kernel::DiagonalDephasingKernel;
use quantum_work_network::diagnostics::lindblad_action;
use quantum_work_network::ergotropy::ergotropy;
use quantum_work_network::matrix::{
    commutator, expectation, frobenius_norm, hermiticity_error, ComplexMatrix, C64,
};
use quantum_work_network::operators::{build_operators_for_chain, ModelParams, Operators};
use quantum_work_network::partial_trace::partial_trace;
use quantum_work_network::time_dependent::lindblad_rhs;

const DT: f64 = 0.0025;
const SAVE_STEPS: usize = 4;
const GAMMA: f64 = 0.5;
const TRAJECTORY_TOL: f64 = 2.0e-9;
const RHS_MAX_TOL: f64 = 1.0e-12;
const TRACE_TOL: f64 = 1.0e-8;
const HERM_TOL: f64 = 1.0e-8;
const POS_TOL: f64 = 1.0e-8;
const LEDGER_TOL: f64 = 5.0e-5;

#[derive(Clone, Copy, Debug)]
enum Mode {
    Dense,
    Kernel,
}

struct Problem {
    params: ModelParams,
    ops: Operators,
    collapses: Vec<ComplexMatrix>,
    kernel: DiagonalDephasingKernel,
}

#[derive(Clone)]
struct Sample {
    time: f64,
    energy: f64,
    work: f64,
    usable: f64,
    coherence: f64,
    drive_in: f64,
}

struct FinalMetrics {
    energy: f64,
    work: f64,
    usable: f64,
    coherence: f64,
    drive_in: f64,
    w_over_ein: f64,
    trace_error: f64,
    herm_error: f64,
    min_eigenvalue: f64,
    ledger: f64,
    finite: bool,
}

struct Run {
    captures: Vec<(f64, ComplexMatrix)>,
    samples: Vec<Sample>,
    metrics: FinalMetrics,
    wall_seconds: f64,
    steps: usize,
    kernel_bytes: usize,
}

fn cfg() -> CoherentDriveConfig {
    CoherentDriveConfig {
        omega0: 0.2,
        omega_drive: 1.0,
        tau: 3.2,
        t_end: 10.0,
        dt: DT,
        save_interval: 0.01,
        gamma_phi: GAMMA,
    }
}

fn problem(n: usize, gammas: &[f64], mode: Mode) -> Result<Problem, Box<dyn std::error::Error>> {
    let params = ModelParams::default();
    let ops = build_operators_for_chain(&params, n)?;
    let collapses = if matches!(mode, Mode::Dense) {
        ops.sigma_z_sites
            .iter()
            .zip(gammas)
            .filter(|(_, gamma)| **gamma > 0.0)
            .map(|(z, gamma)| z * C64::new((*gamma / 2.0).sqrt(), 0.0))
            .collect()
    } else {
        Vec::new()
    };
    let kernel = DiagonalDephasingKernel::new(n, params.load_dim, gammas)?;
    Ok(Problem {
        params,
        ops,
        collapses,
        kernel,
    })
}

fn initial_rho(dim: usize) -> ComplexMatrix {
    let mut rho = ComplexMatrix::zeros(dim, dim);
    rho[(0, 0)] = C64::new(1.0, 0.0);
    rho
}

fn rhs(
    rho: &ComplexMatrix,
    h: &ComplexMatrix,
    p: &Problem,
    mode: Mode,
) -> Result<ComplexMatrix, Box<dyn std::error::Error>> {
    match mode {
        Mode::Dense => Ok(lindblad_rhs(rho, h, &p.collapses)?),
        Mode::Kernel => {
            let mut out = (h * rho - rho * h) * C64::new(0.0, -1.0);
            p.kernel.add_to(rho, &mut out)?;
            Ok(out)
        }
    }
}

fn rk4(
    rho: &ComplexMatrix,
    time: f64,
    p: &Problem,
    mode: Mode,
) -> Result<ComplexMatrix, Box<dyn std::error::Error>> {
    let h = |t| &p.ops.h_total + drive_hamiltonian(t, &cfg(), &p.ops.sigma_1_plus);
    let half = C64::new(0.5 * DT, 0.0);
    let full = C64::new(DT, 0.0);
    let k1 = rhs(rho, &h(time), p, mode)?;
    let k2 = rhs(&(rho + &k1 * half), &h(time + 0.5 * DT), p, mode)?;
    let k3 = rhs(&(rho + &k2 * half), &h(time + 0.5 * DT), p, mode)?;
    let k4 = rhs(&(rho + &k3 * full), &h(time + DT), p, mode)?;
    Ok(rho
        + (k1 + k2 * C64::new(2.0, 0.0) + k3 * C64::new(2.0, 0.0) + k4) * C64::new(DT / 6.0, 0.0))
}

fn dephasing_action(
    rho: &ComplexMatrix,
    p: &Problem,
    mode: Mode,
) -> Result<ComplexMatrix, Box<dyn std::error::Error>> {
    match mode {
        Mode::Kernel => Ok(p.kernel.apply(rho)?),
        Mode::Dense => {
            let mut out = ComplexMatrix::zeros(rho.nrows(), rho.ncols());
            for collapse in &p.collapses {
                out += lindblad_action(collapse, rho)?;
            }
            Ok(out)
        }
    }
}

fn powers(
    rho: &ComplexMatrix,
    time: f64,
    p: &Problem,
    mode: Mode,
) -> Result<(f64, f64), Box<dyn std::error::Error>> {
    let drive = drive_hamiltonian(time, &cfg(), &p.ops.sigma_1_plus);
    let drive_power =
        (expectation(rho, &commutator(&drive, &p.ops.h_total)) * C64::new(0.0, 1.0)).re;
    let dephasing_power = expectation(&dephasing_action(rho, p, mode)?, &p.ops.h_total).re;
    Ok((drive_power, dephasing_power))
}

fn load_values(
    rho: &ComplexMatrix,
    p: &Problem,
) -> Result<(f64, f64, f64, f64), Box<dyn std::error::Error>> {
    let load = partial_trace(rho, &p.ops.dims, &[p.ops.dims.len() - 1])?;
    let h = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
        p.params.load_dim,
        (0..p.params.load_dim).map(|i| C64::new(i as f64 * p.params.omega_load, 0.0)),
    ));
    let er = ergotropy(&load, &h, 1.0e-9)?;
    let coherence: f64 = (0..p.params.load_dim)
        .flat_map(|i| (0..p.params.load_dim).map(move |j| (i, j)))
        .filter(|(i, j)| i != j)
        .map(|(i, j)| load[(i, j)].norm())
        .sum();
    let usable = if er.energy.abs() > 1.0e-14 {
        er.ergotropy / er.energy
    } else {
        f64::NAN
    };
    Ok((er.energy, er.ergotropy, usable, coherence))
}

fn minimum_eigenvalue(rho: &ComplexMatrix) -> f64 {
    let (_, t) = Schur::new(rho.clone()).unpack();
    (0..t.nrows())
        .map(|i| t[(i, i)].re)
        .fold(f64::INFINITY, f64::min)
}

fn simulate(
    n: usize,
    gammas: &[f64],
    t_end: f64,
    mode: Mode,
    captures: &[f64],
    collect_series: bool,
) -> Result<Run, Box<dyn std::error::Error>> {
    let p = problem(n, gammas, mode)?;
    let dim = p.ops.h_total.nrows();
    let rho0 = initial_rho(dim);
    let mut rho = rho0.clone();
    let steps = (t_end / DT).round() as usize;
    let mut captured = Vec::new();
    let mut samples = Vec::new();
    let (e0, w0, u0, c0) = load_values(&rho, &p)?;
    if collect_series {
        samples.push(Sample {
            time: 0.0,
            energy: e0,
            work: w0,
            usable: u0,
            coherence: c0,
            drive_in: 0.0,
        });
    }
    let (mut prev_drive, mut prev_deph) = powers(&rho, 0.0, &p, mode)?;
    let mut last_power_time = 0.0;
    let mut drive_in = 0.0;
    let mut drive_net = 0.0;
    let mut deph_net = 0.0;
    let start = Instant::now();
    for step in 0..steps {
        let time = step as f64 * DT;
        rho = rk4(&rho, time, &p, mode)?;
        let now = (step + 1) as f64 * DT;
        if captures.iter().any(|target| (now - target).abs() < 1.0e-12) {
            captured.push((now, rho.clone()));
        }
        if (step + 1) % SAVE_STEPS == 0 || step + 1 == steps {
            let sample_dt = now - last_power_time;
            let (drive, deph) = powers(&rho, now, &p, mode)?;
            drive_in += 0.5 * sample_dt * (prev_drive.max(0.0) + drive.max(0.0));
            drive_net += 0.5 * sample_dt * (prev_drive + drive);
            deph_net += 0.5 * sample_dt * (prev_deph + deph);
            prev_drive = drive;
            prev_deph = deph;
            last_power_time = now;
            if collect_series {
                let (energy, work, usable, coherence) = load_values(&rho, &p)?;
                samples.push(Sample {
                    time: now,
                    energy,
                    work,
                    usable,
                    coherence,
                    drive_in,
                });
            }
        }
    }
    let wall_seconds = start.elapsed().as_secs_f64();
    let (energy, work, usable, coherence) = load_values(&rho, &p)?;
    let ledger = expectation(&rho, &p.ops.h_total).re
        - expectation(&rho0, &p.ops.h_total).re
        - drive_net
        - deph_net;
    let scalar_finite = [energy, work, coherence, drive_in, ledger]
        .iter()
        .all(|x| x.is_finite());
    let metrics = FinalMetrics {
        energy,
        work,
        usable,
        coherence,
        drive_in,
        w_over_ein: ratio(work, drive_in),
        trace_error: (rho.trace() - C64::new(1.0, 0.0)).norm(),
        herm_error: hermiticity_error(&rho),
        min_eigenvalue: minimum_eigenvalue(&rho),
        ledger,
        finite: scalar_finite && rho.iter().all(|z| z.re.is_finite() && z.im.is_finite()),
    };
    Ok(Run {
        captures: captured,
        samples,
        metrics,
        wall_seconds,
        steps,
        kernel_bytes: p.kernel.estimated_bytes(),
    })
}

fn ratio(a: f64, b: f64) -> f64 {
    if b.abs() <= 1.0e-14 {
        f64::NAN
    } else {
        a / b
    }
}

fn fmt(x: f64) -> String {
    if x.is_finite() {
        format!("{x:.16e}")
    } else {
        "NaN".to_string()
    }
}

fn append_check(
    stage: &str,
    condition: &str,
    check: &str,
    observed: &str,
    expected: &str,
    pass: bool,
) -> std::io::Result<()> {
    let exists = std::path::Path::new("dephasing_kernel_checks.csv").exists();
    let mut w = BufWriter::new(
        OpenOptions::new()
            .create(true)
            .append(true)
            .open("dephasing_kernel_checks.csv")?,
    );
    if !exists {
        writeln!(w, "stage,condition,check,observed,expected,status")?;
    }
    writeln!(
        w,
        "{stage},{condition},{check},{observed},{expected},{}",
        if pass { "PASS" } else { "FAIL" }
    )
}

fn fixed_hermitian(dim: usize) -> ComplexMatrix {
    let mut out = ComplexMatrix::zeros(dim, dim);
    for col in 0..dim {
        for row in 0..=col {
            let v = C64::new(
                ((row + 2 * col + 1) as f64).sin() / dim as f64,
                if row == col {
                    0.0
                } else {
                    ((3 * row + col + 2) as f64).cos() / dim as f64
                },
            );
            out[(row, col)] = v;
            out[(col, row)] = v.conj();
        }
    }
    out
}

fn write_unit_checks() -> Result<(), Box<dyn std::error::Error>> {
    let mut w = BufWriter::new(File::create("dephasing_kernel_unit_checks.csv")?);
    writeln!(w, "check,observed,expected,pass")?;
    let k = DiagonalDephasingKernel::new(3, 3, &[0.2, 0.3, 0.5])?;
    let idx = |chain: usize, load: usize| chain * 3 + load;
    let checks = [
        ("dimension", k.dimension() as f64, 24.0),
        (
            "load_only_difference_rate",
            k.rate(idx(0, 0), idx(0, 2))?,
            0.0,
        ),
        ("site0_bit_rate", k.rate(idx(0, 0), idx(4, 0))?, 0.2),
        ("site1_bit_rate", k.rate(idx(0, 0), idx(2, 0))?, 0.3),
        ("site2_bit_rate", k.rate(idx(0, 0), idx(1, 0))?, 0.5),
        ("three_bit_rate", k.rate(idx(5, 1), idx(2, 2))?, 1.0),
        ("diagonal_rate", k.rate(7, 7)?, 0.0),
    ];
    for (name, observed, expected) in checks {
        let pass = (observed - expected).abs() <= 1.0e-15;
        writeln!(w, "{name},{},{},{}", fmt(observed), fmt(expected), pass)?;
        append_check("unit", "N3", name, &fmt(observed), &fmt(expected), pass)?;
        if !pass {
            return Err(format!("unit check failed: {name}").into());
        }
    }
    let zero = DiagonalDephasingKernel::new(3, 3, &[0.0; 3])?.apply(&fixed_hermitian(24))?;
    let zero_ok = zero.iter().all(|z| z.norm() == 0.0);
    writeln!(w, "gamma_zero,{},true,{}", zero_ok, zero_ok)?;
    append_check(
        "unit",
        "N3",
        "gamma_zero",
        &zero_ok.to_string(),
        "true",
        zero_ok,
    )?;
    Ok(())
}

fn rhs_equivalence() -> Result<(), Box<dyn std::error::Error>> {
    let mut w = BufWriter::new(File::create("dephasing_kernel_rhs_equivalence.csv")?);
    writeln!(w, "chain_length,gamma_configuration,dimension,max_abs_difference,frobenius_difference,reference_frobenius_norm,relative_difference,trace_real,trace_imag,hermiticity_error,pass")?;
    for (n, label, gammas) in [
        (3usize, "all_site", vec![0.5; 3]),
        (3, "site_dependent_partial", vec![0.2, 0.0, 0.7]),
        (5, "all_site", vec![0.5; 5]),
        (5, "site_dependent_partial", vec![0.1, 0.0, 0.3, 0.0, 0.8]),
    ] {
        let p = problem(n, &gammas, Mode::Dense)?;
        let rho = fixed_hermitian(p.ops.h_total.nrows());
        let dense = dephasing_action(&rho, &p, Mode::Dense)?;
        let kernel = dephasing_action(&rho, &p, Mode::Kernel)?;
        let diff = &kernel - &dense;
        let max_abs = diff.iter().map(|z| z.norm()).fold(0.0, f64::max);
        let frob = frobenius_norm(&diff);
        let reference = frobenius_norm(&dense);
        let relative = ratio(frob, reference);
        let pass = max_abs <= RHS_MAX_TOL
            && diff.trace().norm() <= RHS_MAX_TOL
            && hermiticity_error(&kernel) <= 1.0e-12;
        writeln!(
            w,
            "{n},{label},{},{},{},{},{},{},{},{},{}",
            p.ops.h_total.nrows(),
            fmt(max_abs),
            fmt(frob),
            fmt(reference),
            fmt(relative),
            fmt(kernel.trace().re),
            fmt(kernel.trace().im),
            fmt(hermiticity_error(&kernel)),
            pass
        )?;
        append_check(
            "rhs",
            &format!("N{n}_{label}"),
            "dense_kernel_equivalence",
            &fmt(max_abs),
            "<=1e-12",
            pass,
        )?;
        if !pass {
            return Err(format!("RHS equivalence failed for N{n} {label}").into());
        }
    }
    Ok(())
}

fn sample_at(run: &Run, time: f64) -> &Sample {
    run.samples
        .iter()
        .find(|s| (s.time - time).abs() < 1.0e-12)
        .expect("requested saved sample")
}

fn capture_at(run: &Run, time: f64) -> &ComplexMatrix {
    &run.captures
        .iter()
        .find(|(t, _)| (*t - time).abs() < 1.0e-12)
        .expect("requested capture")
        .1
}

fn comparison_pass(a: f64, b: f64, tolerance: f64) -> (f64, f64, bool) {
    if !a.is_finite() && !b.is_finite() {
        (0.0, 0.0, true)
    } else {
        let abs = (a - b).abs();
        let rel = ratio(abs, a.abs().max(b.abs()));
        (abs, rel, abs <= tolerance)
    }
}

fn write_trajectory_metric(
    w: &mut impl Write,
    n: usize,
    condition: &str,
    time: f64,
    source: &str,
    metric: &str,
    dense: f64,
    kernel: f64,
    tolerance: f64,
) -> Result<bool, Box<dyn std::error::Error>> {
    let (abs, rel, pass) = comparison_pass(dense, kernel, tolerance);
    writeln!(
        w,
        "{n},{condition},{},{},{source},{metric},{},{},{},{},{},{}",
        fmt(time),
        fmt(DT),
        fmt(dense),
        fmt(kernel),
        fmt(abs),
        fmt(rel),
        fmt(tolerance),
        pass
    )?;
    append_check(
        "trajectory",
        &format!("N{n}_{condition}_t{time}"),
        metric,
        &fmt(abs),
        &format!("<={tolerance:e}"),
        pass,
    )?;
    Ok(pass)
}

#[derive(Clone)]
struct Reference {
    values: HashMap<String, f64>,
}

fn references() -> Result<HashMap<String, Reference>, Box<dyn std::error::Error>> {
    let file = BufReader::new(File::open("milestone_8a_reference_summary.csv")?);
    let mut lines = file.lines();
    let header: Vec<String> = lines
        .next()
        .ok_or("missing reference header")??
        .split(',')
        .map(str::to_string)
        .collect();
    let mut out = HashMap::new();
    for line in lines {
        let fields: Vec<String> = line?.split(',').map(str::to_string).collect();
        let mut values = HashMap::new();
        for (name, value) in header.iter().zip(&fields) {
            if let Ok(number) = value.parse::<f64>() {
                values.insert(name.clone(), number);
            }
        }
        out.insert(fields[0].clone(), Reference { values });
    }
    Ok(out)
}

fn summary_value(run: &Run, metric: &str) -> f64 {
    let end = run.samples.last().expect("series present");
    match metric {
        "E_at_t10" => end.energy,
        "W_at_t10" => end.work,
        "usable_fraction_at_t10" => end.usable,
        "coherence_L1_at_t10" => end.coherence,
        "drive_energy_in_at_t10" => end.drive_in,
        "W_over_Ein_at_t10" => ratio(end.work, end.drive_in),
        "W_max" => run
            .samples
            .iter()
            .map(|s| s.work)
            .fold(f64::NEG_INFINITY, f64::max),
        "t_at_W_max" => {
            run.samples
                .iter()
                .max_by(|a, b| a.work.total_cmp(&b.work))
                .unwrap()
                .time
        }
        "E_time_area_full" => run
            .samples
            .windows(2)
            .map(|p| 0.5 * (p[1].time - p[0].time) * (p[0].energy + p[1].energy))
            .sum(),
        "W_time_area_full" => run
            .samples
            .windows(2)
            .map(|p| 0.5 * (p[1].time - p[0].time) * (p[0].work + p[1].work))
            .sum(),
        _ => panic!("unknown summary metric {metric}"),
    }
}

fn validate() -> Result<(), Box<dyn std::error::Error>> {
    if std::path::Path::new("dephasing_kernel_checks.csv").exists() {
        return Err("validation outputs already exist; refusing to overwrite".into());
    }
    write_unit_checks()?;
    rhs_equivalence()?;
    let capture_n3 = [0.1, 1.0, 10.0];
    let capture_n5_dense = [0.1, 1.0];
    let n3_all_dense = simulate(3, &[0.5; 3], 10.0, Mode::Dense, &capture_n3, true)?;
    let n3_all_kernel = simulate(3, &[0.5; 3], 10.0, Mode::Kernel, &capture_n3, true)?;
    let partial = [0.2, 0.0, 0.7];
    let n3_partial_dense = simulate(3, &partial, 10.0, Mode::Dense, &capture_n3, true)?;
    let n3_partial_kernel = simulate(3, &partial, 10.0, Mode::Kernel, &capture_n3, true)?;
    let n5_dense = simulate(5, &[0.5; 5], 1.0, Mode::Dense, &capture_n5_dense, true)?;
    let n5_kernel = simulate(5, &[0.5; 5], 10.0, Mode::Kernel, &capture_n3, true)?;
    let n3_free = simulate(3, &[0.0; 3], 10.0, Mode::Kernel, &capture_n3, true)?;
    let mut tw = BufWriter::new(File::create("dephasing_kernel_trajectory_equivalence.csv")?);
    writeln!(tw, "chain_length,condition,final_time,dt,reference_source,metric,dense_value,kernel_value,absolute_difference,relative_difference,tolerance,pass")?;
    let mut all = true;
    for (n, condition, dense, kernel, times) in [
        (
            3usize,
            "all_site_noisy",
            &n3_all_dense,
            &n3_all_kernel,
            &capture_n3[..],
        ),
        (
            3,
            "partial_site_dependent",
            &n3_partial_dense,
            &n3_partial_kernel,
            &capture_n3[..],
        ),
        (
            5,
            "all_site_noisy",
            &n5_dense,
            &n5_kernel,
            &capture_n5_dense[..],
        ),
    ] {
        for &time in times {
            let drho = capture_at(dense, time);
            let krho = capture_at(kernel, time);
            let ds = sample_at(dense, time);
            let ks = sample_at(kernel, time);
            let max_density = (drho - krho).iter().map(|z| z.norm()).fold(0.0, f64::max);
            all &= write_trajectory_metric(
                &mut tw,
                n,
                condition,
                time,
                "dense_recomputed",
                "full_density_max_abs",
                0.0,
                max_density,
                TRAJECTORY_TOL,
            )?;
            let dload = partial_trace(
                drho,
                &build_operators_for_chain(&ModelParams::default(), n)?.dims,
                &[n],
            )?;
            let kload = partial_trace(
                krho,
                &build_operators_for_chain(&ModelParams::default(), n)?.dims,
                &[n],
            )?;
            let load_diff = (&dload - &kload)
                .iter()
                .map(|z| z.norm())
                .fold(0.0, f64::max);
            all &= write_trajectory_metric(
                &mut tw,
                n,
                condition,
                time,
                "dense_recomputed",
                "load_reduced_max_abs",
                0.0,
                load_diff,
                TRAJECTORY_TOL,
            )?;
            for (metric, a, b) in [
                ("load_energy", ds.energy, ks.energy),
                ("load_ergotropy", ds.work, ks.work),
                ("usable_fraction", ds.usable, ks.usable),
                ("coherence_L1", ds.coherence, ks.coherence),
                ("drive_energy_in", ds.drive_in, ks.drive_in),
                (
                    "W_over_Ein",
                    ratio(ds.work, ds.drive_in),
                    ratio(ks.work, ks.drive_in),
                ),
            ] {
                all &= write_trajectory_metric(
                    &mut tw,
                    n,
                    condition,
                    time,
                    "dense_recomputed",
                    metric,
                    a,
                    b,
                    TRAJECTORY_TOL,
                )?;
            }
        }
    }
    let refs = references()?;
    let n5_ref = refs
        .get("N5_all_site_noisy")
        .ok_or("missing N5 reference")?;
    let n5_end = sample_at(&n5_kernel, 10.0);
    for (metric, reference_key, kernel_value) in [
        ("load_energy", "E_at_t10", n5_end.energy),
        ("load_ergotropy", "W_at_t10", n5_end.work),
        ("usable_fraction", "usable_fraction_at_t10", n5_end.usable),
        ("coherence_L1", "coherence_L1_at_t10", n5_end.coherence),
        ("drive_energy_in", "drive_energy_in_at_t10", n5_end.drive_in),
        (
            "W_over_Ein",
            "W_over_Ein_at_t10",
            ratio(n5_end.work, n5_end.drive_in),
        ),
    ] {
        all &= write_trajectory_metric(
            &mut tw,
            5,
            "all_site_noisy",
            10.0,
            "milestone_8a_existing",
            metric,
            *n5_ref.values.get(reference_key).unwrap(),
            kernel_value,
            TRAJECTORY_TOL,
        )?;
    }
    let mut rw = BufWriter::new(File::create("dephasing_kernel_regression.csv")?);
    writeln!(rw, "chain_length,condition,metric,existing_value,kernel_value,absolute_difference,relative_difference,tolerance,pass")?;
    for (condition, n, run) in [
        ("N3_noise_free", 3usize, &n3_free),
        ("N3_all_site_noisy", 3, &n3_all_kernel),
        ("N5_all_site_noisy", 5, &n5_kernel),
    ] {
        let reference = refs.get(condition).ok_or("missing regression reference")?;
        for metric in [
            "E_at_t10",
            "W_at_t10",
            "usable_fraction_at_t10",
            "coherence_L1_at_t10",
            "drive_energy_in_at_t10",
            "W_over_Ein_at_t10",
            "W_max",
            "t_at_W_max",
            "E_time_area_full",
            "W_time_area_full",
        ] {
            let existing = *reference.values.get(metric).unwrap();
            let kernel = summary_value(run, metric);
            let tolerance = if metric == "t_at_W_max" {
                0.011
            } else {
                TRAJECTORY_TOL
            };
            let (abs, rel, pass) = comparison_pass(existing, kernel, tolerance);
            writeln!(
                rw,
                "{n},{condition},{metric},{},{},{},{},{},{}",
                fmt(existing),
                fmt(kernel),
                fmt(abs),
                fmt(rel),
                fmt(tolerance),
                pass
            )?;
            append_check(
                "regression",
                condition,
                metric,
                &fmt(abs),
                &format!("<={tolerance:e}"),
                pass,
            )?;
            all &= pass;
        }
    }
    for (condition, run) in [
        ("N3_noise_free", &n3_free),
        ("N3_all_site_noisy", &n3_all_kernel),
        ("N5_all_site_noisy", &n5_kernel),
    ] {
        for (name, value, pass) in [
            (
                "trace",
                run.metrics.trace_error,
                run.metrics.trace_error <= TRACE_TOL,
            ),
            (
                "hermiticity",
                run.metrics.herm_error,
                run.metrics.herm_error <= HERM_TOL,
            ),
            (
                "positivity",
                run.metrics.min_eigenvalue,
                run.metrics.min_eigenvalue >= -POS_TOL,
            ),
            (
                "ledger",
                run.metrics.ledger.abs(),
                run.metrics.ledger.abs() <= LEDGER_TOL,
            ),
        ] {
            append_check(
                "quality",
                condition,
                name,
                &fmt(value),
                "within tolerance",
                pass,
            )?;
            all &= pass;
        }
        append_check(
            "quality",
            condition,
            "finite",
            &run.metrics.finite.to_string(),
            "true",
            run.metrics.finite,
        )?;
        all &= run.metrics.finite;
        let w_bound = run.metrics.work <= run.metrics.energy + 1.0e-10;
        append_check(
            "quality",
            condition,
            "W_le_E",
            &w_bound.to_string(),
            "true",
            w_bound,
        )?;
        all &= w_bound;
        let usable_bound = !run.metrics.usable.is_finite()
            || (run.metrics.usable >= -1.0e-9 && run.metrics.usable <= 1.0 + 1.0e-9);
        append_check(
            "quality",
            condition,
            "usable_fraction_range",
            &fmt(run.metrics.usable),
            "NaN before signal or 0..1",
            usable_bound,
        )?;
        all &= usable_bound;
        let scalar_finite = [
            run.metrics.coherence,
            run.metrics.drive_in,
            run.metrics.w_over_ein,
        ]
        .iter()
        .all(|x| x.is_finite() || x.is_nan());
        append_check(
            "quality",
            condition,
            "derived_scalar_validity",
            &scalar_finite.to_string(),
            "true",
            scalar_finite,
        )?;
        all &= scalar_finite;
    }
    if !all {
        return Err("one or more equivalence/regression checks failed".into());
    }
    println!("validation PASS: N3/N5 RHS, trajectories, regressions");
    Ok(())
}

fn benchmark(duration: f64) -> Result<(), Box<dyn std::error::Error>> {
    if ![0.1, 0.5].iter().any(|x| (duration - x).abs() < 1.0e-12) {
        return Err("benchmark duration must be 0.1 or 0.5".into());
    }
    let exists = std::path::Path::new("dephasing_kernel_benchmarks.csv").exists();
    let mut w = BufWriter::new(
        OpenOptions::new()
            .create(true)
            .append(true)
            .open("dephasing_kernel_benchmarks.csv")?,
    );
    if !exists {
        writeln!(w, "chain_length,condition,implementation,probe_duration,dt,step_count,run_index,wall_clock_seconds,seconds_per_step,min_seconds_per_step,median_seconds_per_step,max_seconds_per_step")?;
        for (condition, gammas) in [
            ("N7_noise_free", vec![0.0; 7]),
            ("N7_all_site_noisy", vec![0.5; 7]),
        ] {
            let one = simulate(7, &gammas, DT, Mode::Kernel, &[DT], false)?;
            append_check(
                "N7_one_step",
                condition,
                "finite",
                &one.metrics.finite.to_string(),
                "true",
                one.metrics.finite,
            )?;
            append_check(
                "N7_one_step",
                condition,
                "trace",
                &fmt(one.metrics.trace_error),
                "<=1e-8",
                one.metrics.trace_error <= TRACE_TOL,
            )?;
            append_check(
                "N7_construction",
                condition,
                "kernel_bytes",
                &one.kernel_bytes.to_string(),
                "1179648",
                one.kernel_bytes == 384 * 384 * 8,
            )?;
        }
    }
    for (condition, gammas) in [
        ("N7_noise_free", vec![0.0; 7]),
        ("N7_all_site_noisy", vec![0.5; 7]),
    ] {
        let mut runs = Vec::new();
        for run_index in 0..4 {
            let run = simulate(7, &gammas, duration, Mode::Kernel, &[], false)?;
            let sec = run.wall_seconds / run.steps as f64;
            println!(
                "{condition} duration={duration} run={run_index} wall={:.3}s sec/step={sec:.6}",
                run.wall_seconds
            );
            if !run.metrics.finite
                || run.metrics.trace_error > TRACE_TOL
                || run.metrics.herm_error > HERM_TOL
                || run.metrics.min_eigenvalue < -POS_TOL
                || run.metrics.ledger.abs() > LEDGER_TOL
            {
                return Err(format!("N7 quality failure: {condition} run {run_index}").into());
            }
            runs.push((run_index, run.wall_seconds, sec));
        }
        let mut measured: Vec<f64> = runs
            .iter()
            .filter(|(i, _, _)| *i > 0)
            .map(|(_, _, s)| *s)
            .collect();
        measured.sort_by(f64::total_cmp);
        let min = measured[0];
        let median = measured[1];
        let max = measured[2];
        for (index, wall, sec) in runs {
            writeln!(
                w,
                "7,{condition},DiagonalDephasingKernel,{},{},{},{index},{},{},{},{},{}",
                fmt(duration),
                fmt(DT),
                (duration / DT).round() as usize,
                fmt(wall),
                fmt(sec),
                fmt(min),
                fmt(median),
                fmt(max)
            )?;
        }
        append_check(
            "N7_benchmark",
            condition,
            &format!("t{duration}_quality"),
            "all four runs stable",
            "PASS",
            true,
        )?;
    }
    Ok(())
}

fn finalize() -> Result<(), Box<dyn std::error::Error>> {
    let lines: Vec<String> = BufReader::new(File::open("dephasing_kernel_benchmarks.csv")?)
        .lines()
        .skip(1)
        .collect::<Result<_, _>>()?;
    let mut longest: HashMap<String, (f64, f64, f64, f64)> = HashMap::new();
    for line in &lines {
        let v: Vec<_> = line.split(',').collect();
        let condition = v[1].to_string();
        let duration: f64 = v[3].parse()?;
        let min: f64 = v[9].parse()?;
        let median: f64 = v[10].parse()?;
        let max: f64 = v[11].parse()?;
        if longest
            .get(&condition)
            .map(|x| duration > x.0)
            .unwrap_or(true)
        {
            longest.insert(condition, (duration, min, median, max));
        }
    }
    let old_free = 1.1001046649999999;
    let old_noisy = 21.317895377500001;
    let mut ew = BufWriter::new(File::create("dephasing_kernel_estimates.csv")?);
    writeln!(ew, "chain_length,condition,old_seconds_per_step,new_seconds_per_step,speedup_factor,estimated_t10_seconds,estimated_t10_hours,feasibility_class")?;
    for (condition, old) in [
        ("N7_noise_free", old_free),
        ("N7_all_site_noisy", old_noisy),
    ] {
        let (_, _, median, _) = longest
            .get(condition)
            .ok_or("missing benchmark condition")?;
        let speedup = old / median;
        let seconds = median * 4000.0;
        let hours = seconds / 3600.0;
        let class = if condition.contains("noisy") && speedup < 5.0 {
            "insufficient_speedup"
        } else if hours <= 6.0 {
            "feasible_candidate"
        } else {
            "borderline"
        };
        writeln!(
            ew,
            "7,{condition},{},{},{},{},{},{class}",
            fmt(old),
            fmt(*median),
            fmt(speedup),
            fmt(seconds),
            fmt(hours)
        )?;
    }
    ew.flush()?;
    drop(ew);
    let checks = std::fs::read_to_string("dephasing_kernel_checks.csv")?;
    let pass_count = checks
        .lines()
        .skip(1)
        .filter(|x| x.ends_with(",PASS"))
        .count();
    let fail_count = checks
        .lines()
        .skip(1)
        .filter(|x| x.ends_with(",FAIL"))
        .count();
    if fail_count > 0 {
        return Err("cannot finalize with failed checks".into());
    }
    let estimates = std::fs::read_to_string("dephasing_kernel_estimates.csv")?;
    let noisy_line = estimates
        .lines()
        .find(|x| x.starts_with("7,N7_all_site_noisy,"))
        .ok_or("missing noisy estimate")?;
    let nv: Vec<_> = noisy_line.split(',').collect();
    let new_noisy: f64 = nv[3].parse()?;
    let speedup: f64 = nv[4].parse()?;
    let hours: f64 = nv[6].parse()?;
    let decision = nv[7];
    let mut r = BufWriter::new(File::create("MILESTONE_8C_REPORT.md")?);
    writeln!(
        r,
        "# Milestone 8c: Exact dephasing-kernel optimization and equivalence validation\n"
    )?;
    let sections = [
        ("1. 目的", "局所sigma_z位相雑音のLindblad項だけを、物理模型を変えずに厳密高速化した。"),
        ("2. Milestone 8bでのボトルネック", "N=7 all-site noisyは旧dense collapse実装で1 step約21.317895秒、t=10推定約23.69時間だった。"),
        ("3. 物理模型を変更していないこと", "Hamiltonian、drive、RK4、dt、load、初期状態、gamma、ergotropy、Hilbert空間、dense density matrixを変更していない。"),
        ("4. dephasing kernelの数式", "各要素へ -Gamma[a,b] rho[a,b] を加え、Gammaはchain site bitが異なるsiteのgamma和として一度だけ構築した。これは同じLindblad dissipatorの厳密な成分表示である。"),
        ("5. basis mapping", "既存tensor順序 |q1,...,qN,load> をoperator対角要素と全basis pairで照合した。load levelはGammaへ含めない。"),
        ("6. 実装", "旧Dense pathを残し、DiagonalDephasingKernelを独立moduleとして追加した。kernelはf64連続配列で、N=7では1,179,648 bytes。"),
        ("7. unit checks", "Hamming rate、load除外、site順序、diagonal zero、gamma zero、trace、Hermiticity、入力検証を追加した。"),
        ("8. dissipator単体一致", "N=3/N=5のall-siteとsite-dependent partial条件でdense dissipatorと比較した。詳細はrhs equivalence CSV。"),
        ("9. trajectory一致", "N=3 all-site/partialはt=0.1,1,10、N=5 all-siteはt=0.1,1をdense再計算し、N=5 t=10はMilestone 8a既存値と照合した。"),
        ("10. N=3回帰", "noise-freeとall-site noisyのt=10主要量、W peak、時間面積をMilestone 8aへ回帰した。"),
        ("11. N=5回帰", "all-site noisyのE、W、usable fraction、W peak、時間面積、drive inputをMilestone 8aへ回帰した。"),
    ];
    for (title, body) in sections {
        writeln!(r, "## {title}\n\n{body}\n")?;
    }
    writeln!(r, "## 12. N=7 benchmark\n\n最長probeについてwarmup 1回、measurement 3回をrelease buildで実行した。詳細はbenchmarks CSV。\n")?;
    writeln!(r, "## 13. speedup\n\nN=7 noisyのmedianは `{:.6}` s/step、旧値に対するspeedupは **{:.2}x**。\n", new_noisy, speedup)?;
    writeln!(r, "## 14. t=10実行時間再推定\n\nN=7 noisyはmedianから約 **{:.3}時間**。分類は **{}**。今回はt=10本計算を実行していない。\n", hours, decision)?;
    writeln!(r, "## 15. memory\n\nN=7 Gamma kernelは1,179,648 bytes（1.125 MiB）。dephasing評価用の7個のcollapse行列はkernel hot pathで使用しない。既存operator bundleは診断・reference互換性のため維持した。\n")?;
    writeln!(r, "## 16. 数値品質チェック\n\nchecks CSVは **{pass_count} PASS / {fail_count} FAIL**。許容値はRHS max abs 1e-12、trajectory主要量2e-9、trace/Hermiticity/positivity 1e-8、ledger 5e-5。\n")?;
    writeln!(r, "## 17. 直接確認できたこと\n\n- sigma_z dephasingをelementwise kernelで再現した。\n- N=3/N=5で旧dense path・既存値と許容内一致した。\n- N=7 noisyの短時間速度とt=10再推定を測定した。\n")?;
    writeln!(r, "## 18. 確認できていないこと\n\nN=7 t=10の最終物理量、長時間の実測性能、N scaling、他のLindblad operatorへの適用。\n")?;
    writeln!(r, "## 19. 主張してはいけないこと\n\n物理近似、新現象、一般Lindbladへの普遍適用、N=7最終結果、scaling則、GPU/sparse法より優れるという主張。\n")?;
    writeln!(r, "## 20. N=7本計算へ進む判断\n\n数値等価性、speedup>=5、推定<=6時間を満たす場合だけ **次段階候補** とする。実行承認ではなく、今回自動実行していない。現在判定: **{}**。\n", decision)?;
    writeln!(r, "## 21. 生成ファイル一覧\n\n- `src/dephasing_kernel.rs`\n- `src/bin/dephasing_kernel_benchmark.rs`\n- `dephasing_kernel_unit_checks.csv`\n- `dephasing_kernel_rhs_equivalence.csv`\n- `dephasing_kernel_trajectory_equivalence.csv`\n- `dephasing_kernel_regression.csv`\n- `dephasing_kernel_benchmarks.csv`\n- `dephasing_kernel_estimates.csv`\n- `dephasing_kernel_checks.csv`\n- `MILESTONE_8C_REPORT.md`\n")?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    match env::args().nth(1).as_deref() {
        Some("validate") => validate(),
        Some("benchmark") => benchmark(env::args().nth(2).ok_or("missing duration")?.parse()?),
        Some("finalize") => finalize(),
        _ => Err(
            "usage: dephasing_kernel_benchmark validate | benchmark <0.1|0.5> | finalize".into(),
        ),
    }
}
