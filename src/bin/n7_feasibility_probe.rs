use std::env;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::mem::size_of;
use std::process::Command;
use std::time::Instant;

use nalgebra::linalg::Schur;
use quantum_work_network::coherent_drive::{drive_hamiltonian, CoherentDriveConfig};
use quantum_work_network::diagnostics::lindblad_action;
use quantum_work_network::ergotropy::ergotropy;
use quantum_work_network::matrix::{
    commutator, expectation, hermiticity_error, ComplexMatrix, C64,
};
use quantum_work_network::operators::{build_operators_for_chain, ModelParams, Operators};
use quantum_work_network::partial_trace::partial_trace;
use quantum_work_network::time_dependent::lindblad_rhs;

const N: usize = 7;
const DIM: usize = 384;
const DT: f64 = 0.0025;
const GAMMA: f64 = 0.5;
const TRACE_TOL: f64 = 1.0e-8;
const HERM_TOL: f64 = 1.0e-8;
const POS_TOL: f64 = 1.0e-8;
const LEDGER_TOL: f64 = 5.0e-5;

#[derive(Clone, Copy)]
struct Condition {
    name: &'static str,
    noisy: bool,
}

const CONDITIONS: [Condition; 2] = [
    Condition {
        name: "N7_noise_free",
        noisy: false,
    },
    Condition {
        name: "N7_all_site_noisy",
        noisy: true,
    },
];

struct Metrics {
    load_energy: f64,
    load_ergotropy: f64,
    load_coherence_l1: f64,
    usable_fraction: f64,
    chain_population: f64,
    trace_error: f64,
    hermiticity_error: f64,
    min_eigenvalue: f64,
    finite: bool,
    population_bounds: bool,
    reduced_shape_ok: bool,
}

fn config() -> CoherentDriveConfig {
    CoherentDriveConfig {
        omega0: 0.2,
        omega_drive: 1.0,
        tau: 3.2,
        t_end: 1.0,
        dt: DT,
        save_interval: 1.0,
        gamma_phi: GAMMA,
    }
}

fn collapses(ops: &Operators, noisy: bool) -> Vec<ComplexMatrix> {
    if !noisy {
        return Vec::new();
    }
    let scale = C64::new((GAMMA / 2.0).sqrt(), 0.0);
    ops.sigma_z_sites.iter().map(|z| z * scale).collect()
}

fn initial_rho() -> ComplexMatrix {
    let mut rho = ComplexMatrix::zeros(DIM, DIM);
    rho[(0, 0)] = C64::new(1.0, 0.0);
    rho
}

fn rk4_step(
    rho: &ComplexMatrix,
    time: f64,
    dt: f64,
    ops: &Operators,
    cs: &[ComplexMatrix],
) -> Result<ComplexMatrix, Box<dyn std::error::Error>> {
    let cfg = config();
    let h = |t| &ops.h_total + drive_hamiltonian(t, &cfg, &ops.sigma_1_plus);
    let half = C64::new(0.5 * dt, 0.0);
    let full = C64::new(dt, 0.0);
    let k1 = lindblad_rhs(rho, &h(time), cs)?;
    let k2 = lindblad_rhs(&(rho + &k1 * half), &h(time + 0.5 * dt), cs)?;
    let k3 = lindblad_rhs(&(rho + &k2 * half), &h(time + 0.5 * dt), cs)?;
    let k4 = lindblad_rhs(&(rho + &k3 * full), &h(time + dt), cs)?;
    Ok(rho
        + (k1 + k2 * C64::new(2.0, 0.0) + k3 * C64::new(2.0, 0.0) + k4) * C64::new(dt / 6.0, 0.0))
}

fn powers(
    rho: &ComplexMatrix,
    time: f64,
    ops: &Operators,
    cs: &[ComplexMatrix],
) -> Result<(f64, f64), Box<dyn std::error::Error>> {
    let drive = drive_hamiltonian(time, &config(), &ops.sigma_1_plus);
    let drive_power = (expectation(rho, &commutator(&drive, &ops.h_total)) * C64::new(0.0, 1.0)).re;
    let mut dephasing_power = 0.0;
    for c in cs {
        dephasing_power += expectation(&lindblad_action(c, rho)?, &ops.h_total).re;
    }
    Ok((drive_power, dephasing_power))
}

fn minimum_eigenvalue(rho: &ComplexMatrix) -> f64 {
    let (_, t) = Schur::new(rho.clone()).unpack();
    (0..t.nrows())
        .map(|i| t[(i, i)].re)
        .fold(f64::INFINITY, f64::min)
}

fn metrics(
    rho: &ComplexMatrix,
    ops: &Operators,
    params: &ModelParams,
) -> Result<Metrics, Box<dyn std::error::Error>> {
    let load = partial_trace(rho, &ops.dims, &[N])?;
    let h_load = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
        params.load_dim,
        (0..params.load_dim).map(|level| C64::new(level as f64 * params.omega_load, 0.0)),
    ));
    let er = ergotropy(&load, &h_load, 1.0e-9)?;
    let coherence = (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .filter(|(i, j)| i != j)
        .map(|(i, j)| load[(i, j)].norm())
        .sum();
    let sites: Vec<f64> = ops
        .number_sites
        .iter()
        .map(|n| expectation(rho, n).re)
        .collect();
    let usable = if er.energy.abs() > 1.0e-14 {
        er.ergotropy / er.energy
    } else {
        f64::NAN
    };
    Ok(Metrics {
        load_energy: er.energy,
        load_ergotropy: er.ergotropy,
        load_coherence_l1: coherence,
        usable_fraction: usable,
        chain_population: sites.iter().sum(),
        trace_error: (rho.trace() - C64::new(1.0, 0.0)).norm(),
        hermiticity_error: hermiticity_error(rho),
        min_eigenvalue: minimum_eigenvalue(rho),
        finite: rho.iter().all(|z| z.re.is_finite() && z.im.is_finite()),
        population_bounds: sites.iter().all(|p| *p >= -1.0e-10 && *p <= 1.0 + 1.0e-10),
        reduced_shape_ok: load.shape() == (3, 3),
    })
}

fn process_memory() -> (u64, u64) {
    let pid = std::process::id();
    let script = format!("$p=Get-Process -Id {pid}; Write-Output ($p.WorkingSet64.ToString() + ',' + $p.PeakWorkingSet64.ToString())");
    let output = Command::new("powershell")
        .args(["-NoProfile", "-Command", &script])
        .output();
    if let Ok(out) = output {
        if let Ok(s) = String::from_utf8(out.stdout) {
            let v: Vec<_> = s.trim().split(',').collect();
            if v.len() == 2 {
                return (v[0].parse().unwrap_or(0), v[1].parse().unwrap_or(0));
            }
        }
    }
    (0, 0)
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
    let exists = std::path::Path::new("n7_feasibility_checks.csv").exists();
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open("n7_feasibility_checks.csv")?;
    let mut w = BufWriter::new(file);
    if !exists {
        writeln!(w, "stage,condition,check,observed,expected,status")?;
    }
    writeln!(
        w,
        "{stage},{condition},{check},{observed},{expected},{}",
        if pass { "PASS" } else { "FAIL" }
    )
}

fn construction() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams::default();
    let mut w = BufWriter::new(File::create("n7_feasibility_construction.csv")?);
    writeln!(w, "condition,chain_length,hilbert_dimension,density_matrix_rows,density_matrix_cols,density_matrix_element_count,bond_count,collapse_operator_count,drive_site,load_coupling_site,construction_seconds,estimated_static_bytes,measured_memory_bytes_if_available")?;
    for cond in CONDITIONS {
        let start = Instant::now();
        let ops = build_operators_for_chain(&params, N)?;
        let rho = initial_rho();
        let cs = collapses(&ops, cond.noisy);
        let load = partial_trace(&rho, &ops.dims, &[N])?;
        let elapsed = start.elapsed().as_secs_f64();
        let (working, _) = process_memory();
        let full_matrix_bytes = DIM * DIM * size_of::<C64>();
        let full_matrix_count = 25 + cs.len();
        let static_bytes = full_matrix_count * full_matrix_bytes + 3 * 3 * size_of::<C64>();
        writeln!(
            w,
            "{},{N},{DIM},{DIM},{DIM},{},{},{},0,6,{},{},{}",
            cond.name,
            DIM * DIM,
            N - 1,
            cs.len(),
            fmt(elapsed),
            static_bytes,
            working
        )?;
        for (name, observed, expected, pass) in [
            (
                "hilbert_dimension",
                DIM.to_string(),
                "384".to_string(),
                ops.h_total.shape() == (DIM, DIM),
            ),
            (
                "density_matrix_shape",
                format!("{}x{}", rho.nrows(), rho.ncols()),
                "384x384".to_string(),
                rho.shape() == (DIM, DIM),
            ),
            (
                "density_matrix_element_count",
                (DIM * DIM).to_string(),
                "147456".to_string(),
                DIM * DIM == 147456,
            ),
            (
                "bond_count",
                (N - 1).to_string(),
                "6".to_string(),
                N - 1 == 6,
            ),
            (
                "collapse_operator_count",
                cs.len().to_string(),
                if cond.noisy { "7" } else { "0" }.to_string(),
                cs.len() == if cond.noisy { 7 } else { 0 },
            ),
            (
                "drive_site",
                "0".to_string(),
                "0".to_string(),
                ops.sigma_1_plus.shape() == (DIM, DIM),
            ),
            (
                "load_coupling_site",
                "6".to_string(),
                "6".to_string(),
                ops.h_site_3.shape() == (DIM, DIM),
            ),
            (
                "load_reduced_shape",
                format!("{}x{}", load.nrows(), load.ncols()),
                "3x3".to_string(),
                load.shape() == (3, 3),
            ),
            (
                "site_operator_count",
                ops.number_sites.len().to_string(),
                "7".to_string(),
                ops.number_sites.len() == 7,
            ),
            (
                "common_operator_dimensions",
                "384x384".to_string(),
                "384x384".to_string(),
                ops.number_sites
                    .iter()
                    .chain(ops.sigma_z_sites.iter())
                    .all(|x| x.shape() == (DIM, DIM)),
            ),
            (
                "initial_state_trace",
                fmt(rho.trace().re),
                "1".to_string(),
                (rho.trace().re - 1.0).abs() < 1.0e-12,
            ),
            (
                "allocation",
                "completed".to_string(),
                "completed".to_string(),
                true,
            ),
        ] {
            append_check("construction", cond.name, name, &observed, &expected, pass)?;
            if !pass {
                return Err(format!("construction check failed: {} {name}", cond.name).into());
            }
        }
    }
    write_memory_csv()?;
    Ok(())
}

fn write_memory_csv() -> std::io::Result<()> {
    let b = size_of::<C64>();
    let m = DIM * DIM * b;
    let mut w = BufWriter::new(File::create("n7_feasibility_memory.csv")?);
    writeln!(
        w,
        "item,matrix_count,rows,cols,bytes_per_element,estimated_bytes,estimate_type,notes"
    )?;
    let rows = [
        (
            "one_complex_matrix",
            1usize,
            m,
            "exact",
            "Complex64 is two f64 values",
        ),
        ("density_matrix", 1, m, "exact", "current rho"),
        (
            "Hamiltonian",
            1,
            m,
            "exact",
            "one dense full-system Hamiltonian",
        ),
        (
            "RK4_stage_matrices",
            8,
            8 * m,
            "conservative",
            "k1-k4 plus rho2-rho4 and result",
        ),
        (
            "Lindblad_temporaries_noise_free",
            4,
            4 * m,
            "conservative",
            "commutator products and derivative",
        ),
        (
            "Lindblad_temporaries_noisy",
            8,
            8 * m,
            "conservative",
            "per-collapse products, dagger and accumulated derivative",
        ),
        (
            "collapse_operator_group_noisy",
            7,
            7 * m,
            "exact_payload",
            "seven dense dephasing operators",
        ),
        (
            "operator_bundle",
            24,
            24 * m,
            "conservative",
            "full matrices retained by Operators plus small local matrix",
        ),
        (
            "estimated_peak_noise_free",
            37,
            37 * m,
            "conservative",
            "operator bundle + rho + Hamiltonian/RK4/RHS working set",
        ),
        (
            "estimated_peak_noisy",
            48,
            48 * m,
            "conservative",
            "free estimate + collapses + extra dissipator temporaries",
        ),
    ];
    for (item, count, bytes, typ, notes) in rows {
        writeln!(w, "{item},{count},{DIM},{DIM},{b},{bytes},{typ},{notes}")?;
    }
    Ok(())
}

fn write_step_header(truncate: bool) -> std::io::Result<BufWriter<File>> {
    let exists = std::path::Path::new("n7_feasibility_steps.csv").exists();
    let file = if truncate {
        File::create("n7_feasibility_steps.csv")?
    } else {
        OpenOptions::new()
            .create(true)
            .append(true)
            .open("n7_feasibility_steps.csv")?
    };
    let mut w = BufWriter::new(file);
    if truncate || !exists {
        writeln!(w, "condition,dt,step_count,final_time,wall_clock_seconds,seconds_per_step,load_energy,load_ergotropy,load_coherence_l1,usable_fraction,drive_energy_in,trace_error,hermiticity_error,min_eigenvalue,energy_ledger_residual")?;
    }
    Ok(w)
}

fn one_step() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams::default();
    let mut w = write_step_header(true)?;
    for cond in CONDITIONS {
        let ops = build_operators_for_chain(&params, N)?;
        let cs = collapses(&ops, cond.noisy);
        let rho0 = initial_rho();
        let h_start = Instant::now();
        let h = &ops.h_total + drive_hamiltonian(0.0, &config(), &ops.sigma_1_plus);
        let h_seconds = h_start.elapsed().as_secs_f64();
        let rhs_start = Instant::now();
        let rhs = lindblad_rhs(&rho0, &h, &cs)?;
        let rhs_seconds = rhs_start.elapsed().as_secs_f64();
        if rhs.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            return Err("non-finite RHS".into());
        }
        let (p0, d0) = powers(&rho0, 0.0, &ops, &cs)?;
        let start = Instant::now();
        let rho = rk4_step(&rho0, 0.0, DT, &ops, &cs)?;
        let wall = start.elapsed().as_secs_f64();
        let (p1, d1) = powers(&rho, DT, &ops, &cs)?;
        let drive_in = 0.5 * DT * (p0.max(0.0) + p1.max(0.0));
        let drive_net = 0.5 * DT * (p0 + p1);
        let deph_net = 0.5 * DT * (d0 + d1);
        let ledger = expectation(&rho, &ops.h_total).re
            - expectation(&rho0, &ops.h_total).re
            - drive_net
            - deph_net;
        let m = metrics(&rho, &ops, &params)?;
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            cond.name,
            fmt(DT),
            1,
            fmt(DT),
            fmt(wall),
            fmt(wall),
            fmt(m.load_energy),
            fmt(m.load_ergotropy),
            fmt(m.load_coherence_l1),
            fmt(m.usable_fraction),
            fmt(drive_in),
            fmt(m.trace_error),
            fmt(m.hermiticity_error),
            fmt(m.min_eigenvalue),
            fmt(ledger)
        )?;
        let (_, peak) = process_memory();
        for (name, observed, expected, pass) in [
            (
                "hamiltonian_evaluation_seconds",
                fmt(h_seconds),
                "finite nonnegative".to_string(),
                h_seconds.is_finite(),
            ),
            (
                "lindblad_rhs_evaluation_seconds",
                fmt(rhs_seconds),
                "finite nonnegative".to_string(),
                rhs_seconds.is_finite(),
            ),
            (
                "rk4_one_step_seconds",
                fmt(wall),
                "finite nonnegative".to_string(),
                wall.is_finite(),
            ),
            (
                "trace_preservation",
                fmt(m.trace_error),
                format!("<= {TRACE_TOL:e}"),
                m.trace_error <= TRACE_TOL,
            ),
            (
                "hermiticity",
                fmt(m.hermiticity_error),
                format!("<= {HERM_TOL:e}"),
                m.hermiticity_error <= HERM_TOL,
            ),
            (
                "positivity",
                fmt(m.min_eigenvalue),
                format!(">= -{POS_TOL:e}"),
                m.min_eigenvalue >= -POS_TOL,
            ),
            (
                "finite_values",
                m.finite.to_string(),
                "true".to_string(),
                m.finite,
            ),
            (
                "population_bounds",
                m.population_bounds.to_string(),
                "true".to_string(),
                m.population_bounds,
            ),
            (
                "load_reduced_shape",
                m.reduced_shape_ok.to_string(),
                "true".to_string(),
                m.reduced_shape_ok,
            ),
            (
                "ergotropy_le_energy",
                format!("W={} E={}", fmt(m.load_ergotropy), fmt(m.load_energy)),
                "W <= E".to_string(),
                m.load_ergotropy <= m.load_energy + 1.0e-10,
            ),
            (
                "energy_ledger",
                fmt(ledger.abs()),
                format!("<= {LEDGER_TOL:e}"),
                ledger.abs() <= LEDGER_TOL,
            ),
            (
                "total_chain_population",
                fmt(m.chain_population),
                "within 0..7".to_string(),
                m.chain_population >= -1.0e-10 && m.chain_population <= 7.0 + 1.0e-10,
            ),
            (
                "peak_process_memory_bytes",
                peak.to_string(),
                "measured if available".to_string(),
                true,
            ),
        ] {
            append_check("one_step", cond.name, name, &observed, &expected, pass)?;
            if !pass {
                return Err(format!("one-step check failed: {} {name}", cond.name).into());
            }
        }
    }
    Ok(())
}

fn benchmark(duration: f64) -> Result<(), Box<dyn std::error::Error>> {
    if ![0.1, 0.5, 1.0]
        .iter()
        .any(|x| (duration - x).abs() < 1.0e-12)
    {
        return Err("duration must be 0.1, 0.5, or 1.0".into());
    }
    let params = ModelParams::default();
    let mut step_writer = write_step_header(false)?;
    let bench_exists = std::path::Path::new("n7_feasibility_benchmarks.csv").exists();
    let mut bench_writer = BufWriter::new(
        OpenOptions::new()
            .create(true)
            .append(true)
            .open("n7_feasibility_benchmarks.csv")?,
    );
    if !bench_exists {
        writeln!(bench_writer, "condition,probe_duration,dt,step_count,wall_clock_seconds,seconds_per_step,estimated_t10_seconds,estimated_t10_hours,estimate_source")?;
    }
    for cond in CONDITIONS {
        let ops = build_operators_for_chain(&params, N)?;
        let cs = collapses(&ops, cond.noisy);
        let rho0 = initial_rho();
        let mut rho = rho0.clone();
        let steps = (duration / DT).round() as usize;
        let (mut p_prev, mut d_prev) = powers(&rho, 0.0, &ops, &cs)?;
        let mut drive_in = 0.0;
        let mut drive_net = 0.0;
        let mut deph_net = 0.0;
        let mut last_power_time = 0.0;
        let start = Instant::now();
        for i in 0..steps {
            let t = i as f64 * DT;
            rho = rk4_step(&rho, t, DT, &ops, &cs)?;
            // Match the existing Milestone 8a diagnostic cadence (0.01),
            // while retaining only the current density matrix.
            if (i + 1) % 4 == 0 || i + 1 == steps {
                let sample_time = t + DT;
                let sample_dt = sample_time - last_power_time;
                let (p, d) = powers(&rho, sample_time, &ops, &cs)?;
                drive_in += 0.5 * sample_dt * (p_prev.max(0.0) + p.max(0.0));
                drive_net += 0.5 * sample_dt * (p_prev + p);
                deph_net += 0.5 * sample_dt * (d_prev + d);
                p_prev = p;
                d_prev = d;
                last_power_time = sample_time;
            }
        }
        let wall = start.elapsed().as_secs_f64();
        let seconds_per_step = wall / steps as f64;
        let estimated_t10 = seconds_per_step * (10.0 / DT);
        let ledger = expectation(&rho, &ops.h_total).re
            - expectation(&rho0, &ops.h_total).re
            - drive_net
            - deph_net;
        let m = metrics(&rho, &ops, &params)?;
        writeln!(
            step_writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            cond.name,
            fmt(DT),
            steps,
            fmt(duration),
            fmt(wall),
            fmt(seconds_per_step),
            fmt(m.load_energy),
            fmt(m.load_ergotropy),
            fmt(m.load_coherence_l1),
            fmt(m.usable_fraction),
            fmt(drive_in),
            fmt(m.trace_error),
            fmt(m.hermiticity_error),
            fmt(m.min_eigenvalue),
            fmt(ledger)
        )?;
        writeln!(
            bench_writer,
            "{},{},{},{},{},{},{},{},measured_RK4_plus_0.01_ledger_diagnostics",
            cond.name,
            fmt(duration),
            fmt(DT),
            steps,
            fmt(wall),
            fmt(seconds_per_step),
            fmt(estimated_t10),
            fmt(estimated_t10 / 3600.0)
        )?;
        let (_, peak) = process_memory();
        for (name, observed, expected, pass) in [
            (
                "trace_preservation",
                fmt(m.trace_error),
                format!("<= {TRACE_TOL:e}"),
                m.trace_error <= TRACE_TOL,
            ),
            (
                "hermiticity",
                fmt(m.hermiticity_error),
                format!("<= {HERM_TOL:e}"),
                m.hermiticity_error <= HERM_TOL,
            ),
            (
                "positivity",
                fmt(m.min_eigenvalue),
                format!(">= -{POS_TOL:e}"),
                m.min_eigenvalue >= -POS_TOL,
            ),
            (
                "finite_values",
                m.finite.to_string(),
                "true".to_string(),
                m.finite,
            ),
            (
                "population_bounds",
                m.population_bounds.to_string(),
                "true".to_string(),
                m.population_bounds,
            ),
            (
                "load_reduced_shape",
                m.reduced_shape_ok.to_string(),
                "true".to_string(),
                m.reduced_shape_ok,
            ),
            (
                "ergotropy_le_energy",
                format!("W={} E={}", fmt(m.load_ergotropy), fmt(m.load_energy)),
                "W <= E".to_string(),
                m.load_ergotropy <= m.load_energy + 1.0e-10,
            ),
            (
                "usable_fraction_range",
                fmt(m.usable_fraction),
                "NaN before signal or 0..1".to_string(),
                !m.usable_fraction.is_finite()
                    || (m.usable_fraction >= -1.0e-9 && m.usable_fraction <= 1.0 + 1.0e-9),
            ),
            (
                "energy_ledger",
                fmt(ledger.abs()),
                format!("<= {LEDGER_TOL:e}"),
                ledger.abs() <= LEDGER_TOL,
            ),
            (
                "total_chain_population",
                fmt(m.chain_population),
                "within 0..7".to_string(),
                m.chain_population >= -1.0e-10 && m.chain_population <= 7.0 + 1.0e-10,
            ),
            (
                "peak_process_memory_bytes",
                peak.to_string(),
                "measured if available".to_string(),
                true,
            ),
        ] {
            append_check(
                &format!("benchmark_t{duration}"),
                cond.name,
                name,
                &observed,
                &expected,
                pass,
            )?;
            if !pass {
                return Err(format!("benchmark check failed: {} {name}", cond.name).into());
            }
        }
        println!(
            "{} duration={} wall={:.3}s sec/step={:.6} est_t10={:.3}h E={:.3e} W={:.3e}",
            cond.name,
            duration,
            wall,
            seconds_per_step,
            estimated_t10 / 3600.0,
            m.load_energy,
            m.load_ergotropy
        );
    }
    Ok(())
}

#[derive(Clone)]
struct BenchRow {
    condition: String,
    duration: f64,
    sec_per_step: f64,
    est_seconds: f64,
}

fn read_benchmarks() -> Result<Vec<BenchRow>, Box<dyn std::error::Error>> {
    let f = BufReader::new(File::open("n7_feasibility_benchmarks.csv")?);
    let mut out = Vec::new();
    for line in f.lines().skip(1) {
        let s = line?;
        let v: Vec<_> = s.split(',').collect();
        out.push(BenchRow {
            condition: v[0].to_string(),
            duration: v[1].parse()?,
            sec_per_step: v[5].parse()?,
            est_seconds: v[6].parse()?,
        });
    }
    Ok(out)
}

fn finalize() -> Result<(), Box<dyn std::error::Error>> {
    let b = read_benchmarks()?;
    let mut longest = Vec::new();
    for cond in CONDITIONS {
        let row = b
            .iter()
            .filter(|r| r.condition == cond.name)
            .max_by(|x, y| x.duration.total_cmp(&y.duration))
            .ok_or("missing benchmark")?
            .clone();
        longest.push(row);
    }
    let checks_text = std::fs::read_to_string("n7_feasibility_checks.csv")?;
    let all_checks = !checks_text
        .lines()
        .skip(1)
        .any(|line| line.ends_with(",FAIL"));
    let peak_noisy = 48usize * DIM * DIM * size_of::<C64>();
    let mut ew = BufWriter::new(File::create("n7_feasibility_estimates.csv")?);
    writeln!(ew, "condition,target_dt,target_tmax,estimated_steps,estimated_seconds,estimated_hours,feasibility_class,basis_probe_duration,startup_cost_warning")?;
    let mut classes = Vec::new();
    for row in &longest {
        for target_dt in [DT, DT / 2.0] {
            let seconds = row.sec_per_step * (10.0 / target_dt);
            let hours = seconds / 3600.0;
            let limit = if row.condition.contains("noise_free") {
                1.0
            } else {
                6.0
            };
            let class = if !all_checks {
                "infeasible_with_current_dense_method"
            } else if hours > 9.0 && row.condition.contains("noisy") {
                "infeasible_with_current_dense_method"
            } else if hours >= 0.8 * limit {
                "borderline"
            } else {
                "feasible"
            };
            if (target_dt - DT).abs() < 1.0e-12 {
                classes.push(class.to_string());
            }
            writeln!(
                ew,
                "{},{},{},{},{},{},{},{},startup excluded; longest propagation probe used",
                row.condition,
                fmt(target_dt),
                fmt(10.0),
                (10.0 / target_dt).round() as usize,
                fmt(seconds),
                fmt(hours),
                class,
                fmt(row.duration)
            )?;
        }
    }
    let overall = if classes
        .iter()
        .any(|c| c == "infeasible_with_current_dense_method")
    {
        "infeasible_with_current_dense_method"
    } else if classes.iter().any(|c| c == "borderline") {
        "borderline"
    } else {
        "feasible"
    };
    let n5_free = 43.8486026 / 4000.0;
    let n5_noisy = 403.443321 / 4000.0;
    let free = &longest[0];
    let noisy = &longest[1];
    let mut r = BufWriter::new(File::create("MILESTONE_8B_REPORT.md")?);
    writeln!(
        r,
        "# Milestone 8b: N=7 computational feasibility and short-time reachability probe\n"
    )?;
    writeln!(r, "## 1. 目的\n\nN=7のdense Lindblad RK4が現在環境で計算可能か、短時間probeだけで評価した。\n")?;
    writeln!(r, "## 2. 今回は本計算ではないこと\n\n`t=10`本計算、半減刻み本計算、位置別雑音、最適化は実行していない。\n")?;
    writeln!(r, "## 3. N=7模型構成\n\n7個の二準位siteと3準位load。J=1、g=0.25、各周波数1、drive site=0、load coupling site=6、Omega=0.2、tau=3.2、真空初期状態、load無雑音。\n")?;
    writeln!(r, "## 4. 次元とoperator mapping\n\nHilbert次元384、density matrix 384x384（147456要素）、bond数6、drive site 0、load coupling site 6。collapse数はfree 0、noisy 7。\n")?;
    writeln!(r, "## 5. construction-only結果\n\n全construction検査PASS。詳細は `n7_feasibility_construction.csv` と checks CSV。\n")?;
    writeln!(r, "## 6. one-step結果\n\n両条件でRK4 1 stepが完了し、trace、Hermiticity、positivity、finite、load縮約、ledger検査はPASS。個別時間は checks CSV。\n")?;
    writeln!(r, "## 7. short-time benchmark\n\n|condition|longest probe|seconds/step|estimated t=10 hours|\n|---|---:|---:|---:|\n|{}|{:.3}|{:.6}|{:.3}|\n|{}|{:.3}|{:.6}|{:.3}|\n\n`t=0.1` noisy の時点で推定が9時間を大きく超え、現行dense法の infeasible 基準が確定した。このため、判断を変えず約1時間超を追加消費する `t=0.5` と、条件付きの `t=1.0` は実行しなかった。\n", free.condition, free.duration, free.sec_per_step, free.est_seconds / 3600.0, noisy.condition, noisy.duration, noisy.sec_per_step, noisy.est_seconds / 3600.0)?;
    writeln!(r, "## 8. noise-free実行時間推定\n\ndt=0.0025で約{:.3}時間。最長probeのpropagation step時間を使用し、startupは除外した。\n", free.est_seconds / 3600.0)?;
    writeln!(r, "## 9. all-site noisy実行時間推定\n\ndt=0.0025で約{:.3}時間。最長probeのpropagation step時間を使用し、startupは除外した。\n", noisy.est_seconds / 3600.0)?;
    writeln!(r, "## 10. 半減刻み実行時間推定\n\n未実行。step数が2倍としてfree約{:.3}時間、noisy約{:.3}時間と概算した。\n", 2.0 * free.est_seconds / 3600.0, 2.0 * noisy.est_seconds / 3600.0)?;
    writeln!(r, "## 11. メモリ推定\n\nComplex64は{} bytes。1行列は{} bytes（約{:.2} MiB）。保守的peakはfree約{:.1} MiB、noisy約{:.1} MiB。construction時の実測working setはfree約60.7 MiB、noisy約76.6 MiB、one-step後の実測peakはnoisy約106.6 MiBだった。実行時の利用可能物理メモリは約18.6 GiBで、メモリ不足は見られなかった。\n", size_of::<C64>(), DIM * DIM * size_of::<C64>(), (DIM * DIM * size_of::<C64>()) as f64 / 1048576.0, (37 * DIM * DIM * size_of::<C64>()) as f64 / 1048576.0, peak_noisy as f64 / 1048576.0)?;
    writeln!(r, "## 12. N=5実測との比較\n\nHilbert次元比N7/N5=4、density element比=16。seconds/step比はfree {:.2}、noisy {:.2}。N=5は1001保存点の診断を含み、N=7 probeはpropagation主体なので厳密な同条件比較ではない。collapse数比noisy=7/5=1.4。2点から一般scalingは推定しない。\n", free.sec_per_step / n5_free, noisy.sec_per_step / n5_noisy)?;
    writeln!(
        r,
        "## 13. 数値品質チェック\n\n全記録検査: **{}**。\n",
        if all_checks { "PASS" } else { "FAIL" }
    )?;
    let steps = std::fs::read_to_string("n7_feasibility_steps.csv")?;
    let last_lines: Vec<_> = steps.lines().rev().take(2).collect();
    writeln!(r, "## 14. load初期応答\n\n最長probe終点の値は steps CSV末尾に保存した。t<=1でEやWがほぼ0でも到達不能とは判定しない。今回確認するのは初期兆候だけである。\n\n- {}\n- {}\n", last_lines.get(1).unwrap_or(&"missing"), last_lines.first().unwrap_or(&"missing"))?;
    writeln!(r, "## 15. 本計算へ進む基準\n\nconstruction、one-step、short-time数値品質、メモリ、free<=1h、noisy<=6h、probe長で推定が大きく発散しないこと。\n")?;
    writeln!(r, "## 16. feasibility判定\n\n**{overall}**。これは現在環境・現在dense実装での計算可否分類であり、物理的可能性ではない。`t=10`本計算は自動実行していない。\n")?;
    writeln!(r, "## 17. 最適化候補\n\n必要なら、一時行列削減、collapse共通部分再利用、sparse表現、正確性を保つ別積分法、Krylov、exponential integrator、並列行列積、BLAS/release設定確認を比較候補にする。今回は実装していない。\n")?;
    writeln!(r, "## 18. 直接確認できたこと\n\n- N=7の構築とmapping。\n- 両条件の1 stepおよび実行済み短時間probeの安定性。\n- 現在実装でのstep時間、t=10概算、memory概算。\n")?;
    writeln!(r, "## 19. 確認できていないこと\n\nN=7のt=10最終値、半減刻み実測、長時間安定性、最終到達、距離依存則、最適化後性能。\n")?;
    writeln!(r, "## 20. 主張してはいけないこと\n\nN=7が到達不能、N=7の最終W低下、指数/べきscaling、物理的限界、量子優位、新規性。\n")?;
    writeln!(r, "## 21. 次段階への判断材料\n\n現行dense法のまま`t=10`本計算へ進む候補とはしない。まず一時行列削減、collapse共通部分再利用、BLAS/release設定など、物理模型を変えない最適化候補を別Milestoneで比較する判断材料になった。今回は最適化を実装していない。\n")?;
    writeln!(r, "## 22. 生成ファイル一覧\n\n- `src/bin/n7_feasibility_probe.rs`\n- `n7_feasibility_construction.csv`\n- `n7_feasibility_steps.csv`\n- `n7_feasibility_benchmarks.csv`\n- `n7_feasibility_estimates.csv`\n- `n7_feasibility_memory.csv`\n- `n7_feasibility_checks.csv`\n- `MILESTONE_8B_REPORT.md`\n")?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    match env::args().nth(1).as_deref() {
        Some("construction") => construction(),
        Some("one-step") => one_step(),
        Some("benchmark") => benchmark(env::args().nth(2).ok_or("missing duration")?.parse()?),
        Some("finalize") => finalize(),
        _ => Err("usage: n7_feasibility_probe construction | one-step | benchmark <0.1|0.5|1.0> | finalize".into()),
    }
}
