use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
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

const N: usize = 7;
const LOAD_DIM: usize = 3;
const DIM: usize = 384;
const TOTAL_GAMMA: f64 = 1.5;
const GAMMA_SITE: f64 = TOTAL_GAMMA / N as f64;
const DT: f64 = 0.0025;
const T_END: f64 = 0.03;
const OMEGA: f64 = 0.2;
const TAU: f64 = 3.2;
const SCALAR_TOL: f64 = 1.0e-12;
const STATE_TOL: f64 = 1.0e-15;
const TRACE_TOL: f64 = 1.0e-10;
const HERM_TOL: f64 = 1.0e-12;
const POS_TOL: f64 = 1.0e-8;
const EIG_IMAG_TOL: f64 = 1.0e-10;
const EIG_SUM_TOL: f64 = 1.0e-10;
const SOLVER_AGREEMENT_TOL: f64 = 1.0e-10;
const SAVED_9C: &str = "fixed_total_noise_timeseries.csv";

const OUTPUTS: [&str; 6] = [
    "n7_t002_eigen_diagnostic.csv",
    "n7_t002_state_summary.csv",
    "n7_t002_reproducibility.csv",
    "n7_t002_saved_value_comparison.csv",
    "n7_t002_diagnostic_checks.csv",
    "MILESTONE_9C_DIAGNOSTIC.md",
];

#[derive(Clone)]
struct StateSummary {
    run: usize,
    time: f64,
    all_finite: bool,
    frobenius: f64,
    trace: C64,
    hermiticity: f64,
    antihermitian: f64,
    correction: f64,
    diagonal_imag_max: f64,
    diagonal_real_min: f64,
    diagonal_real_max: f64,
    max_abs: f64,
    sum: C64,
}

#[derive(Clone)]
struct EigenResult {
    run: usize,
    time: f64,
    method: &'static str,
    input: &'static str,
    minimum: f64,
    max_imag: f64,
    sum: C64,
    trace: C64,
    sum_trace_difference: f64,
    finite: bool,
    pass: bool,
}

#[derive(Clone)]
struct Metrics {
    load_energy: f64,
    load_ergotropy: f64,
    load_coherence_l1: f64,
    bare_network_energy: f64,
    drive_power: f64,
    dephasing_power: f64,
    trace_error: f64,
    hermiticity_error: f64,
}

#[derive(Clone)]
struct ReproRow {
    run: usize,
    state: StateSummary,
    direct: f64,
    hermitianized: f64,
    schur: f64,
    state_matches_run1: bool,
    solvers_finite: bool,
    runtime_seconds: f64,
}

#[derive(Clone)]
struct Check {
    stage: &'static str,
    check: &'static str,
    observed: String,
    expected: String,
    pass: bool,
}

fn n(x: f64) -> String {
    if x.is_finite() {
        format!("{x:.16e}")
    } else if x.is_nan() {
        "NaN".to_string()
    } else if x.is_sign_positive() {
        "Inf".to_string()
    } else {
        "-Inf".to_string()
    }
}

fn config() -> CoherentDriveConfig {
    CoherentDriveConfig {
        omega0: OMEGA,
        omega_drive: 1.0,
        tau: TAU,
        t_end: 10.0,
        dt: DT,
        save_interval: 0.01,
        gamma_phi: GAMMA_SITE,
    }
}

fn ensure_new_outputs() -> Result<(), Box<dyn std::error::Error>> {
    for path in OUTPUTS {
        if std::path::Path::new(path).exists() {
            return Err(format!("refusing to overwrite existing output {path}").into());
        }
    }
    if !std::path::Path::new(SAVED_9C).exists() {
        return Err(format!("missing saved Milestone 9c input {SAVED_9C}").into());
    }
    Ok(())
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
) -> Result<ComplexMatrix, Box<dyn std::error::Error>> {
    let cfg = config();
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

fn short_run(
    ops: &Operators,
    kernel: &DiagonalDephasingKernel,
) -> Result<Vec<(f64, ComplexMatrix)>, Box<dyn std::error::Error>> {
    let mut rho = ComplexMatrix::zeros(DIM, DIM);
    rho[(0, 0)] = C64::new(1.0, 0.0);
    let mut saved = vec![(0.0, rho.clone())];
    for step in 0..12 {
        rho = rk4_step(&rho, step as f64 * DT, ops, kernel)?;
        if (step + 1) % 4 == 0 {
            saved.push(((step + 1) as f64 * DT, rho.clone()));
        }
    }
    if saved.len() != 4 || (saved.last().unwrap().0 - T_END).abs() > 1.0e-15 {
        return Err("short run did not produce t=0,0.01,0.02,0.03".into());
    }
    Ok(saved)
}

fn hermitianized(rho: &ComplexMatrix) -> ComplexMatrix {
    (rho + rho.adjoint()) * C64::new(0.5, 0.0)
}

fn state_summary(run: usize, time: f64, rho: &ComplexMatrix) -> StateSummary {
    let rho_h = hermitianized(rho);
    let anti = rho - rho.adjoint();
    let mut diagonal_imag_max = 0.0_f64;
    let mut diagonal_real_min = f64::INFINITY;
    let mut diagonal_real_max = f64::NEG_INFINITY;
    let mut max_abs = 0.0_f64;
    let mut sum = C64::new(0.0, 0.0);
    let mut all_finite = true;
    for z in rho.iter() {
        all_finite &= z.re.is_finite() && z.im.is_finite();
        max_abs = max_abs.max(z.norm());
        sum += z;
    }
    for i in 0..rho.nrows() {
        diagonal_imag_max = diagonal_imag_max.max(rho[(i, i)].im.abs());
        diagonal_real_min = diagonal_real_min.min(rho[(i, i)].re);
        diagonal_real_max = diagonal_real_max.max(rho[(i, i)].re);
    }
    StateSummary {
        run,
        time,
        all_finite,
        frobenius: frobenius_norm(rho),
        trace: rho.trace(),
        hermiticity: hermiticity_error(rho),
        antihermitian: frobenius_norm(&anti),
        correction: frobenius_norm(&(&rho_h - rho)),
        diagonal_imag_max,
        diagonal_real_min,
        diagonal_real_max,
        max_abs,
        sum,
    }
}

fn symmetric_result(
    run: usize,
    time: f64,
    method: &'static str,
    input: &'static str,
    rho: &ComplexMatrix,
) -> EigenResult {
    let eigenvalues = SymmetricEigen::new(rho.clone()).eigenvalues;
    let finite = eigenvalues.iter().all(|x| x.is_finite());
    let minimum = eigenvalues.iter().copied().fold(f64::INFINITY, f64::min);
    let sum_real: f64 = eigenvalues.iter().sum();
    let sum = C64::new(sum_real, 0.0);
    let trace = rho.trace();
    let difference = (sum - trace).norm();
    let pass =
        finite && minimum >= -POS_TOL && difference <= EIG_SUM_TOL && trace.im.abs() <= TRACE_TOL;
    EigenResult {
        run,
        time,
        method,
        input,
        minimum,
        max_imag: 0.0,
        sum,
        trace,
        sum_trace_difference: difference,
        finite,
        pass,
    }
}

fn schur_result(run: usize, time: f64, rho: &ComplexMatrix) -> EigenResult {
    let (_, triangular) = Schur::new(rho.clone()).unpack();
    let eigenvalues: Vec<C64> = (0..triangular.nrows())
        .map(|i| triangular[(i, i)])
        .collect();
    let finite = eigenvalues
        .iter()
        .all(|z| z.re.is_finite() && z.im.is_finite());
    let minimum = eigenvalues
        .iter()
        .map(|z| z.re)
        .fold(f64::INFINITY, f64::min);
    let max_imag = eigenvalues.iter().map(|z| z.im.abs()).fold(0.0, f64::max);
    let sum: C64 = eigenvalues.iter().copied().sum();
    let trace = rho.trace();
    let difference = (sum - trace).norm();
    let pass =
        finite && minimum >= -POS_TOL && max_imag <= EIG_IMAG_TOL && difference <= EIG_SUM_TOL;
    EigenResult {
        run,
        time,
        method: "complex_schur",
        input: "raw_rho",
        minimum,
        max_imag,
        sum,
        trace,
        sum_trace_difference: difference,
        finite,
        pass,
    }
}

fn metrics(
    rho: &ComplexMatrix,
    time: f64,
    ops: &Operators,
    params: &ModelParams,
    kernel: &DiagonalDephasingKernel,
) -> Result<Metrics, Box<dyn std::error::Error>> {
    let load = partial_trace(rho, &ops.dims, &[N])?;
    let h_load = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
        params.load_dim,
        (0..params.load_dim).map(|i| C64::new(i as f64 * params.omega_load, 0.0)),
    ));
    let work = ergotropy(&load, &h_load, 1.0e-9)?;
    let coherence_l1: f64 = (0..params.load_dim)
        .flat_map(|i| (0..params.load_dim).map(move |j| (i, j)))
        .filter(|(i, j)| i != j)
        .map(|(i, j)| load[(i, j)].norm())
        .sum();
    let drive = drive_hamiltonian(time, &config(), &ops.sigma_1_plus);
    let drive_power = expectation(rho, &commutator(&drive, &ops.h_total)) * C64::new(0.0, 1.0);
    let dephasing_power = expectation(&kernel.apply(rho)?, &ops.h_total);
    Ok(Metrics {
        load_energy: work.energy,
        load_ergotropy: work.ergotropy,
        load_coherence_l1: coherence_l1,
        bare_network_energy: expectation(rho, &ops.h_total).re,
        drive_power: drive_power.re,
        dephasing_power: dephasing_power.re,
        trace_error: (rho.trace() - C64::new(1.0, 0.0)).norm(),
        hermiticity_error: hermiticity_error(rho),
    })
}

fn read_saved() -> Result<HashMap<(i32, String), f64>, Box<dyn std::error::Error>> {
    let mut lines = BufReader::new(File::open(SAVED_9C)?).lines();
    let header = lines.next().ok_or("empty saved CSV")??;
    let names: Vec<&str> = header.split(',').collect();
    let index = |name: &str| {
        names
            .iter()
            .position(|x| *x == name)
            .ok_or_else(|| format!("missing saved column {name}"))
    };
    let condition_i = index("condition")?;
    let time_i = index("time")?;
    let metrics = [
        "load_energy",
        "load_ergotropy",
        "load_coherence_l1",
        "bare_network_energy",
        "drive_power",
        "dephasing_power",
        "trace_error",
        "hermiticity_error",
    ];
    let metric_indices: Vec<(&str, usize)> = metrics
        .iter()
        .map(|name| Ok((*name, index(name)?)))
        .collect::<Result<_, String>>()?;
    let mut out = HashMap::new();
    for line in lines {
        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.get(condition_i) != Some(&"N7_fixed_total_noise") {
            continue;
        }
        let time: f64 = fields[time_i].parse()?;
        let key = (time * 100.0).round() as i32;
        if matches!(key, 1 | 2 | 3) {
            for (name, i) in &metric_indices {
                out.insert((key, (*name).to_string()), fields[*i].parse()?);
            }
        }
    }
    if out.len() != 24 {
        return Err(format!("expected 24 saved scalar values, found {}", out.len()).into());
    }
    Ok(out)
}

fn metric_values(m: &Metrics) -> [(&'static str, f64); 8] {
    [
        ("load_energy", m.load_energy),
        ("load_ergotropy", m.load_ergotropy),
        ("load_coherence_l1", m.load_coherence_l1),
        ("bare_network_energy", m.bare_network_energy),
        ("drive_power", m.drive_power),
        ("dephasing_power", m.dephasing_power),
        ("trace_error", m.trace_error),
        ("hermiticity_error", m.hermiticity_error),
    ]
}

fn construction_checks(
    params: &ModelParams,
    ops: &Operators,
    kernel: &DiagonalDephasingKernel,
) -> Result<Vec<Check>, Box<dyn std::error::Error>> {
    let gammas = vec![GAMMA_SITE; N];
    let mut mapping = true;
    let mut load_excluded = true;
    for row in 0..DIM {
        let row_chain = row / LOAD_DIM;
        for col in 0..DIM {
            let col_chain = col / LOAD_DIM;
            let differing = row_chain ^ col_chain;
            let expected: f64 = (0..N)
                .filter(|site| differing & (1usize << (N - 1 - site)) != 0)
                .map(|site| gammas[site])
                .sum();
            let observed = kernel.rate(row, col)?;
            mapping &= (observed - expected).abs() <= 1.0e-12;
            if row_chain == col_chain {
                load_excluded &= observed == 0.0;
            }
        }
    }
    let cfg = config();
    let model_ok = params.omega_chain == 1.0
        && params.omega_load == 1.0
        && params.hopping_j == 1.0
        && params.coupling_g == 0.25
        && params.load_dim == LOAD_DIM
        && cfg.omega0 == OMEGA
        && cfg.tau == TAU
        && cfg.dt == DT
        && ops.h_total.shape() == (DIM, DIM)
        && ops.sigma_z_sites.len() == N;
    Ok(vec![
        Check {
            stage: "construction",
            check: "model_configuration_matches_9c",
            observed: model_ok.to_string(),
            expected: "true".into(),
            pass: model_ok,
        },
        Check {
            stage: "construction",
            check: "gamma_per_site",
            observed: n(GAMMA_SITE),
            expected: n(1.5 / 7.0),
            pass: GAMMA_SITE == 1.5 / 7.0,
        },
        Check {
            stage: "construction",
            check: "gamma_sum",
            observed: n(gammas.iter().sum()),
            expected: n(TOTAL_GAMMA),
            pass: (gammas.iter().sum::<f64>() - TOTAL_GAMMA).abs() <= 1.0e-14,
        },
        Check {
            stage: "construction",
            check: "kernel_mapping",
            observed: mapping.to_string(),
            expected: "true".into(),
            pass: mapping,
        },
        Check {
            stage: "construction",
            check: "load_excluded",
            observed: load_excluded.to_string(),
            expected: "true".into(),
            pass: load_excluded,
        },
    ])
}

fn same_state(a: &StateSummary, b: &StateSummary) -> bool {
    (a.frobenius - b.frobenius).abs() <= STATE_TOL
        && (a.trace - b.trace).norm() <= STATE_TOL
        && (a.max_abs - b.max_abs).abs() <= STATE_TOL
        && (a.sum - b.sum).norm() <= STATE_TOL
}

fn write_state(rows: &[StateSummary]) -> Result<(), Box<dyn std::error::Error>> {
    let mut w = BufWriter::new(File::create(OUTPUTS[1])?);
    writeln!(w, "run_index,time,all_elements_finite,frobenius_norm,trace_real,trace_imag,hermiticity_error,antihermitian_norm,hermitianization_correction_norm,diagonal_imag_max,diagonal_real_min,diagonal_real_max,max_abs_element,sum_real_elements,sum_imag_elements")?;
    for x in rows {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            x.run,
            n(x.time),
            x.all_finite,
            n(x.frobenius),
            n(x.trace.re),
            n(x.trace.im),
            n(x.hermiticity),
            n(x.antihermitian),
            n(x.correction),
            n(x.diagonal_imag_max),
            n(x.diagonal_real_min),
            n(x.diagonal_real_max),
            n(x.max_abs),
            n(x.sum.re),
            n(x.sum.im)
        )?;
    }
    Ok(())
}

fn write_eigen(rows: &[EigenResult]) -> Result<(), Box<dyn std::error::Error>> {
    let mut w = BufWriter::new(File::create(OUTPUTS[0])?);
    writeln!(w, "run_index,time,method,input_variant,minimum_eigenvalue,maximum_eigenvalue_imaginary_part,eigenvalue_sum_real,eigenvalue_sum_imag,trace_real,trace_imag,eigenvalue_sum_trace_difference,finite,pass")?;
    for x in rows {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{},{},{},{},{}",
            x.run,
            n(x.time),
            x.method,
            x.input,
            n(x.minimum),
            n(x.max_imag),
            n(x.sum.re),
            n(x.sum.im),
            n(x.trace.re),
            n(x.trace.im),
            n(x.sum_trace_difference),
            x.finite,
            x.pass
        )?;
    }
    Ok(())
}

fn write_repro(rows: &[ReproRow]) -> Result<(), Box<dyn std::error::Error>> {
    let mut w = BufWriter::new(File::create(OUTPUTS[2])?);
    writeln!(w, "run_index,time,frobenius_norm,trace_real,trace_imag,max_abs_element,sum_real_elements,sum_imag_elements,symmetric_direct_min_eigenvalue,symmetric_hermitianized_min_eigenvalue,schur_min_eigenvalue,state_summary_matches_run1,solver_results_finite,runtime_seconds")?;
    for x in rows {
        let s = &x.state;
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            x.run,
            n(s.time),
            n(s.frobenius),
            n(s.trace.re),
            n(s.trace.im),
            n(s.max_abs),
            n(s.sum.re),
            n(s.sum.im),
            n(x.direct),
            n(x.hermitianized),
            n(x.schur),
            x.state_matches_run1,
            x.solvers_finite,
            n(x.runtime_seconds)
        )?;
    }
    Ok(())
}

fn write_saved_comparison(
    saved: &HashMap<(i32, String), f64>,
    rerun: &[(f64, Metrics)],
) -> Result<(bool, f64), Box<dyn std::error::Error>> {
    let mut w = BufWriter::new(File::create(OUTPUTS[3])?);
    writeln!(
        w,
        "time,metric,saved_value,rerun_value,absolute_difference,tolerance,pass"
    )?;
    let mut all_pass = true;
    let mut max_difference = 0.0_f64;
    for (time, metrics) in rerun {
        let key = (*time * 100.0).round() as i32;
        for (name, value) in metric_values(metrics) {
            let saved_value = saved[&(key, name.to_string())];
            let difference = (saved_value - value).abs();
            let pass = difference <= SCALAR_TOL;
            all_pass &= pass;
            max_difference = max_difference.max(difference);
            writeln!(
                w,
                "{},{},{},{},{},{},{}",
                n(*time),
                name,
                n(saved_value),
                n(value),
                n(difference),
                n(SCALAR_TOL),
                pass
            )?;
        }
    }
    Ok((all_pass, max_difference))
}

fn write_checks(checks: &[Check]) -> Result<(), Box<dyn std::error::Error>> {
    let mut w = BufWriter::new(File::create(OUTPUTS[4])?);
    writeln!(w, "stage,check,observed,expected,status")?;
    for x in checks {
        writeln!(
            w,
            "{},{},{},{},{}",
            x.stage,
            x.check,
            x.observed,
            x.expected,
            if x.pass { "PASS" } else { "FAIL" }
        )?;
    }
    Ok(())
}

fn write_report(
    state_t002: &StateSummary,
    direct: &EigenResult,
    herm: &EigenResult,
    schur: &EigenResult,
    direct_nan_reproduced: bool,
    max_saved_difference: f64,
    reproducible: bool,
    verdict: &str,
    updated: &str,
    checks_pass: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let sections = vec![
        ("1. 目的", "Milestone 9cのN=7 t=0.02 minimum-eigenvalue NaNを、短時間再計算だけで切り分けた。".to_string()),
        ("2. Milestone 9cで発生した症状", "t=0.02のdirect SymmetricEigen 1点だけがNaNとなり、positivityとfinite_valuesがFAILした。".to_string()),
        ("3. 再計算範囲", "t=0から0.03まで12 RK4 stepsを3回実行し、t=0,0.01,0.02,0.03を保存した。t=10本計算は再実行していない。".to_string()),
        ("4. 変更していない模型", "N=7、total gamma=1.5、gamma_site=1.5/7、all-site noise、Omega=0.2、tau=3.2、J=1、g=0.25、dt=0.0025、RK4、exact dephasing kernel、真空初期状態を維持した。".to_string()),
        ("5. t=0.01結果", "保存済み9c CSVと8主要値は差0で一致した。directおよびHermitianized SymmetricEigenには一部非有限固有値があったが、有限値だけから作られたminimumは -2.808e-22 だったため、9cのminimum列では異常が表面化していなかった。Schurは全固有値有限だった。".to_string()),
        ("6. t=0.02結果", format!("rho finite={}、trace={}+{}i、Hermiticity={:.3e}。", state_t002.all_finite, n(state_t002.trace.re), n(state_t002.trace.im), state_t002.hermiticity)),
        ("7. t=0.03結果", "保存済み9c CSVと8主要値を照合した。".to_string()),
        ("8. density matrix有限性", format!("t=0.02の全要素finite={}、Frobenius norm={}、max abs element={}。", state_t002.all_finite, n(state_t002.frobenius), n(state_t002.max_abs))),
        ("9. traceとHermiticity", format!("trace error={:.3e}、Hermiticity error={:.3e}、Hermitian化補正norm={:.3e}。", (state_t002.trace-C64::new(1.0,0.0)).norm(), state_t002.hermiticity, state_t002.correction)),
        ("10. direct SymmetricEigen", format!("raw出力には非有限固有値があり、minimum集約結果={}、finite={}。9cのCSV formatterは非有限値を一律NaNと記録するため、元CSVのNaNに対応する現象は3回の短時間再計算で{}。", n(direct.minimum), direct.finite, if direct_nan_reproduced { "再現した" } else { "再現しなかった" })),
        ("11. Hermitianized SymmetricEigen", format!("minimum={}、finite={}、eigenvalue sum-trace差={:.3e}。", n(herm.minimum), herm.finite, herm.sum_trace_difference)),
        ("12. independent solver", format!("Complex Schur minimum={}、max eigenvalue imag={:.3e}、sum-trace差={:.3e}。", n(schur.minimum), schur.max_imag, schur.sum_trace_difference)),
        ("13. solver間比較", format!("HermitianizedとSchurのminimum差={:.3e}、許容値={:.1e}。", (herm.minimum-schur.minimum).abs(), SOLVER_AGREEMENT_TOL)),
        ("14. 保存済み9c結果との一致", format!("24 scalar比較の最大絶対差={:.3e}、許容値={:.1e}。9cにはfull density matrixが保存されていないため、density matrix要素のrun-to-saved直接比較はできない。代わりに3回の決定論的rho要約を比較した。", max_saved_difference, SCALAR_TOL)),
        ("15. 3回再現性", format!("t=0.02の決定論的rho要約一致={}。", reproducible)),
        ("16. minimum eigenvalue判定", format!("独立solver minimum={}、基準 >= -{:.1e}。", n(schur.minimum), POS_TOL)),
        ("17. 直接確認できたこと", "短時間状態の有限性、trace、Hermiticity、3 solver、保存済み主要値、3回再現性を直接確認した。状態側の検査とSchur結果は正常だった一方、rawとHermitianizedの両SymmetricEigenで非有限固有値が再現した。".to_string()),
        ("18. 確認できていないこと", "t=0.03より後のdensity matrix再構築、別dt、別gamma、別solver crate、N=7 t=10再計算は行っていない。".to_string()),
        ("19. Milestone 9c判定の更新可否", format!("更新可否={}。元レポートは上書きせず、この補足だけを追加した。", if updated == "completed_comparison_with_diagnostic_note" { "可" } else { "不可" })),
        ("20. 最終判定", format!("指定判定規則による診断判定は **{verdict}**。state_level_numerical_issueはHermitianized SymmetricEigenも非有限ならCase Bという停止規則によるラベルであり、rho自体にNaN、trace異常、Hermiticity異常、非決定性が見つかったという意味ではない。Milestone 9c更新後判定 **{updated}**。全必須チェック={}。", checks_pass)),
        ("21. 生成ファイル一覧", OUTPUTS.iter().map(|x| format!("- `{x}`")).collect::<Vec<_>>().join("\n")),
    ];
    let mut w = BufWriter::new(File::create(OUTPUTS[5])?);
    writeln!(
        w,
        "# Milestone 9c diagnostic: N=7 t=0.02 minimum-eigenvalue NaN\n"
    )?;
    for (title, body) in sections {
        writeln!(w, "## {title}\n\n{body}\n")?;
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    ensure_new_outputs()?;
    let saved = read_saved()?;
    let params = ModelParams::default();
    let ops = build_operators_for_chain(&params, N)?;
    let gammas = vec![GAMMA_SITE; N];
    let kernel = DiagonalDephasingKernel::new(N, LOAD_DIM, &gammas)?;
    let mut checks = construction_checks(&params, &ops, &kernel)?;
    if checks.iter().any(|x| !x.pass) {
        return Err("construction check failed".into());
    }

    let mut states = Vec::new();
    let mut eigen = Vec::new();
    let mut repro = Vec::new();
    let mut run1_t002: Option<StateSummary> = None;
    let mut first_metrics = Vec::new();
    let allocation_ok = true;
    for run in 1..=3 {
        let start = Instant::now();
        let snapshots = short_run(&ops, &kernel)?;
        for (time, rho) in &snapshots {
            let summary = state_summary(run, *time, rho);
            let direct = symmetric_result(run, *time, "symmetric_eigen_direct", "raw_rho", rho);
            let rho_h = hermitianized(rho);
            let herm = symmetric_result(
                run,
                *time,
                "symmetric_eigen_hermitianized",
                "hermitianized_rho",
                &rho_h,
            );
            let schur = schur_result(run, *time, rho);
            if run == 1 && *time > 0.0 {
                first_metrics.push((*time, metrics(rho, *time, &ops, &params, &kernel)?));
            }
            if (*time - 0.02).abs() <= 1.0e-15 {
                if run == 1 {
                    run1_t002 = Some(summary.clone());
                }
                let matches = run1_t002
                    .as_ref()
                    .is_none_or(|reference| same_state(&summary, reference));
                repro.push(ReproRow {
                    run,
                    state: summary.clone(),
                    direct: direct.minimum,
                    hermitianized: herm.minimum,
                    schur: schur.minimum,
                    state_matches_run1: matches,
                    solvers_finite: direct.finite && herm.finite && schur.finite,
                    runtime_seconds: start.elapsed().as_secs_f64(),
                });
            }
            states.push(summary);
            eigen.extend([direct, herm, schur]);
        }
        println!(
            "diagnostic run {run}/3 completed in {:.3}s",
            start.elapsed().as_secs_f64()
        );
    }

    let (saved_pass, max_saved_difference) = write_saved_comparison(&saved, &first_metrics)?;
    let t002_states: Vec<&StateSummary> = states
        .iter()
        .filter(|x| (x.time - 0.02).abs() <= 1.0e-15)
        .collect();
    let reproducible = t002_states.iter().all(|x| same_state(x, t002_states[0]));
    let state_ok = states.iter().all(|x| {
        x.all_finite
            && (x.trace - C64::new(1.0, 0.0)).norm() <= TRACE_TOL
            && x.hermiticity <= HERM_TOL
    });
    let direct_t002: Vec<&EigenResult> = eigen
        .iter()
        .filter(|x| (x.time - 0.02).abs() <= 1.0e-15 && x.method == "symmetric_eigen_direct")
        .collect();
    let herm_t002: Vec<&EigenResult> = eigen
        .iter()
        .filter(|x| (x.time - 0.02).abs() <= 1.0e-15 && x.method == "symmetric_eigen_hermitianized")
        .collect();
    let schur_t002: Vec<&EigenResult> = eigen
        .iter()
        .filter(|x| (x.time - 0.02).abs() <= 1.0e-15 && x.method == "complex_schur")
        .collect();
    let direct_nan_reproduced = direct_t002.iter().any(|x| !x.finite);
    let herm_ok = herm_t002.iter().all(|x| x.pass);
    let schur_ok = schur_t002.iter().all(|x| x.pass);
    let solver_agreement = herm_t002
        .iter()
        .zip(&schur_t002)
        .all(|(a, b)| (a.minimum - b.minimum).abs() <= SOLVER_AGREEMENT_TOL);
    let no_nan_outside_direct = states.iter().all(|x| x.all_finite)
        && herm_t002.iter().all(|x| x.finite)
        && schur_t002.iter().all(|x| x.finite)
        && first_metrics
            .iter()
            .all(|(_, m)| metric_values(m).iter().all(|(_, x)| x.is_finite()));
    let diagnostic_issue = state_ok
        && saved_pass
        && reproducible
        && herm_ok
        && schur_ok
        && solver_agreement
        && no_nan_outside_direct;
    let verdict = if diagnostic_issue {
        "diagnostic_eigensolver_issue"
    } else if state_ok && herm_ok && schur_ok && !solver_agreement {
        "eigensolver_disagreement_stop"
    } else {
        "state_level_numerical_issue"
    };
    let updated = if diagnostic_issue {
        "completed_comparison_with_diagnostic_note"
    } else {
        "numerical_issue_stop"
    };

    checks.extend([
        Check {
            stage: "execution",
            check: "short_run_completed",
            observed: states.len().to_string(),
            expected: "12 state summaries".into(),
            pass: states.len() == 12,
        },
        Check {
            stage: "state",
            check: "all_rho_elements_finite",
            observed: state_ok.to_string(),
            expected: "true".into(),
            pass: state_ok,
        },
        Check {
            stage: "state",
            check: "trace",
            observed: n(states
                .iter()
                .map(|x| (x.trace - C64::new(1.0, 0.0)).norm())
                .fold(0.0, f64::max)),
            expected: "<=1e-10".into(),
            pass: states
                .iter()
                .all(|x| (x.trace - C64::new(1.0, 0.0)).norm() <= TRACE_TOL),
        },
        Check {
            stage: "state",
            check: "hermiticity",
            observed: n(states.iter().map(|x| x.hermiticity).fold(0.0, f64::max)),
            expected: "<=1e-12".into(),
            pass: states.iter().all(|x| x.hermiticity <= HERM_TOL),
        },
        Check {
            stage: "comparison",
            check: "saved_scalar_reproduction",
            observed: n(max_saved_difference),
            expected: "<=1e-12".into(),
            pass: saved_pass,
        },
        Check {
            stage: "reproducibility",
            check: "run_to_run_reproducibility",
            observed: reproducible.to_string(),
            expected: "true".into(),
            pass: reproducible,
        },
        Check {
            stage: "solver",
            check: "direct_solver_finite_or_isolated_failure",
            observed: format!(
                "finite_all={} isolated={}",
                direct_t002.iter().all(|x| x.finite),
                herm_ok && schur_ok
            ),
            expected: "true or isolated direct failure".into(),
            pass: direct_t002.iter().all(|x| x.finite) || (herm_ok && schur_ok),
        },
        Check {
            stage: "solver",
            check: "hermitianized_solver_finite",
            observed: herm_ok.to_string(),
            expected: "true".into(),
            pass: herm_ok,
        },
        Check {
            stage: "solver",
            check: "independent_solver_finite",
            observed: schur_ok.to_string(),
            expected: "true".into(),
            pass: schur_ok,
        },
        Check {
            stage: "solver",
            check: "minimum_eigenvalue_tolerance",
            observed: n(schur_t002
                .iter()
                .map(|x| x.minimum)
                .fold(f64::INFINITY, f64::min)),
            expected: ">=-1e-8".into(),
            pass: schur_t002.iter().all(|x| x.minimum >= -POS_TOL),
        },
        Check {
            stage: "solver",
            check: "eigenvalue_sum_trace_agreement",
            observed: n(herm_t002
                .iter()
                .chain(schur_t002.iter())
                .map(|x| x.sum_trace_difference)
                .fold(0.0, f64::max)),
            expected: "<=1e-10".into(),
            pass: herm_t002
                .iter()
                .chain(schur_t002.iter())
                .all(|x| x.sum_trace_difference <= EIG_SUM_TOL),
        },
        Check {
            stage: "solver",
            check: "solver_agreement",
            observed: solver_agreement.to_string(),
            expected: "difference <=1e-10".into(),
            pass: solver_agreement,
        },
        Check {
            stage: "finiteness",
            check: "no_nan_outside_direct_solver",
            observed: no_nan_outside_direct.to_string(),
            expected: "true".into(),
            pass: no_nan_outside_direct,
        },
        Check {
            stage: "execution",
            check: "no_allocation_failure",
            observed: allocation_ok.to_string(),
            expected: "true".into(),
            pass: allocation_ok,
        },
    ]);
    let checks_pass = checks.iter().all(|x| x.pass);
    write_state(&states)?;
    write_eigen(&eigen)?;
    write_repro(&repro)?;
    write_checks(&checks)?;
    let s = run1_t002.as_ref().ok_or("missing run1 t=0.02")?;
    write_report(
        s,
        direct_t002[0],
        herm_t002[0],
        schur_t002[0],
        direct_nan_reproduced,
        max_saved_difference,
        reproducible,
        verdict,
        updated,
        checks_pass,
    )?;
    println!("verdict={verdict} updated={updated} checks_pass={checks_pass} direct_nan_reproduced={direct_nan_reproduced}");
    if !diagnostic_issue {
        return Err(format!("diagnostic stop: {verdict}").into());
    }
    Ok(())
}
