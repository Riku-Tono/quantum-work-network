use std::fs::File;
use std::io::{BufWriter, Write};

use quantum_work_network::coherent_drive::{
    drive_hamiltonian, run_coherent_drive, CoherentDriveConfig, CoherentDriveRun,
};
use quantum_work_network::ergotropy::ergotropy;
use quantum_work_network::load_extraction::{
    build_passive_extraction, mutual_information, sorted_hermitian_eigenvalues,
    DEGENERACY_TOLERANCE, ENTROPY_EIGENVALUE_TOLERANCE,
};
use quantum_work_network::matrix::{
    eye, frobenius_norm, hermiticity_error, kron, ComplexMatrix, C64,
};
use quantum_work_network::operators::{build_operators, ModelParams};
use quantum_work_network::partial_trace::partial_trace;

const DT: f64 = 0.0025;
const OMEGA_A: f64 = 0.2;
const OMEGA_B: f64 = 0.431953125;
const RECONSTRUCTION_TOLERANCE: f64 = 1.0e-9;
const UNITARY_TOLERANCE: f64 = 1.0e-10;
const TRACE_TOLERANCE: f64 = 1.0e-10;
const HERMITICITY_TOLERANCE: f64 = 1.0e-10;
const POSITIVITY_TOLERANCE: f64 = 1.0e-8;
const STATE_AGREEMENT_TOLERANCE: f64 = 1.0e-9;
const ENERGY_AGREEMENT_TOLERANCE: f64 = 1.0e-9;
const POST_ERGOTROPY_TOLERANCE: f64 = 1.0e-9;
const MUTUAL_INFORMATION_TOLERANCE: f64 = 1.0e-9;

#[derive(Clone)]
struct CheckRow {
    condition: &'static str,
    check: &'static str,
    error: f64,
    tolerance: f64,
    pass: bool,
}

struct ExtractionResult {
    condition: &'static str,
    omega: f64,
    gamma: f64,
    load_energy_before: f64,
    load_ergotropy_before: f64,
    load_passive_energy_before: f64,
    load_purity_before: f64,
    load_coherence_before: f64,
    populations_before: Vec<f64>,
    state_eigenvalues: Vec<f64>,
    load_energies: Vec<f64>,
    switch_work: f64,
    switch_work_imaginary: f64,
    gross_work: f64,
    load_energy_after: f64,
    load_ergotropy_after: f64,
    load_coherence_after: f64,
    signed_net_work: f64,
    conservative_net_work: f64,
    drive_energy_in: f64,
    full_energy_on_before: f64,
    full_energy_off_after_switch: f64,
    full_energy_off_before_extraction: f64,
    full_energy_off_after_extraction: f64,
    trace_error_after: f64,
    hermiticity_error_after: f64,
    minimum_eigenvalue_after: f64,
    mutual_information_before: f64,
    mutual_information_after: f64,
    finite: bool,
    physical: bool,
    checks: Vec<CheckRow>,
}

fn n(value: f64) -> String {
    format!("{value:.16e}")
}
fn b(value: bool) -> &'static str {
    if value {
        "true"
    } else {
        "false"
    }
}
fn ratio(value: f64, denominator: f64) -> String {
    if denominator.abs() <= 1.0e-12 {
        "undefined".to_string()
    } else {
        n(value / denominator)
    }
}

fn config(omega: f64, gamma: f64) -> CoherentDriveConfig {
    let mut config = CoherentDriveConfig::milestone_5b(gamma, DT);
    config.omega0 = omega;
    config
}

fn load_hamiltonian(params: &ModelParams) -> ComplexMatrix {
    ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
        params.load_dim,
        (0..params.load_dim).map(|level| C64::new(level as f64 * params.omega_load, 0.0)),
    ))
}

fn coherence_l1(rho: &ComplexMatrix) -> f64 {
    let dim = rho.nrows();
    (0..dim)
        .flat_map(|row| (0..dim).map(move |col| (row, col)))
        .filter(|(row, col)| row != col)
        .map(|(row, col)| rho[(row, col)].norm())
        .sum()
}

fn max_difference(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

fn add_check(
    checks: &mut Vec<CheckRow>,
    condition: &'static str,
    check: &'static str,
    error: f64,
    tolerance: f64,
) {
    checks.push(CheckRow {
        condition,
        check,
        error,
        tolerance,
        pass: error <= tolerance,
    });
}

fn extract(
    condition: &'static str,
    omega: f64,
    gamma: f64,
    run: &CoherentDriveRun,
    params: &ModelParams,
) -> Result<ExtractionResult, Box<dyn std::error::Error>> {
    let operators = build_operators(params)?;
    let rho_before = &run.final_state;
    let h_on = &operators.h_total;
    let h_off = &operators.h_chain + &operators.h_load;
    let drive_at_end = drive_hamiltonian(10.0, &config(omega, gamma), &operators.sigma_1_plus);
    let drive_end_error = frobenius_norm(&drive_at_end);
    let switch_complex = (rho_before * (&h_off - h_on)).trace();
    let interaction_switch = -(rho_before * &operators.h_interaction).trace();
    let full_energy_on_before = (rho_before * h_on).trace().re;
    let full_energy_off_after_switch = (rho_before * &h_off).trace().re;
    let switch_identity_error =
        (full_energy_off_after_switch - full_energy_on_before - switch_complex.re)
            .abs()
            .max((switch_complex - interaction_switch).norm());

    let rho_load = partial_trace(rho_before, &operators.dims, &[3])?;
    let h_load = load_hamiltonian(params);
    let load_before = ergotropy(&rho_load, &h_load, 1.0e-9)?;
    let passive = build_passive_extraction(&rho_load, &h_load)?;
    let u_full = kron(&eye(8), &passive.unitary);
    let rho_after = &u_full * rho_before * u_full.adjoint();
    let rho_load_after = partial_trace(&rho_after, &operators.dims, &[3])?;
    let rho_chain_before = partial_trace(rho_before, &operators.dims, &[0, 1, 2])?;
    let rho_chain_after = partial_trace(&rho_after, &operators.dims, &[0, 1, 2])?;
    let load_after = ergotropy(&rho_load_after, &h_load, 1.0e-9)?;
    let full_energy_off_before_extraction = (rho_before * &h_off).trace().re;
    let full_energy_off_after_extraction = (&rho_after * &h_off).trace().re;
    let gross_work = full_energy_off_before_extraction - full_energy_off_after_extraction;
    let load_energy_drop = load_before.energy - load_after.energy;
    let switch_work = switch_complex.re;
    let signed_net_work = gross_work - switch_work;
    let conservative_net_work = gross_work - switch_work.max(0.0);
    let trace_error_after = (rho_after.trace() - C64::new(1.0, 0.0)).norm();
    let hermiticity_error_after = hermiticity_error(&rho_after);
    let full_eigen_before = sorted_hermitian_eigenvalues(rho_before)?;
    let full_eigen_after = sorted_hermitian_eigenvalues(&rho_after)?;
    let load_eigen_after = sorted_hermitian_eigenvalues(&rho_load_after)?;
    let minimum_eigenvalue_after = full_eigen_after
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let mutual_information_before =
        mutual_information(rho_before, &operators.dims, &[0, 1, 2], &[3])?;
    let mutual_information_after =
        mutual_information(&rho_after, &operators.dims, &[0, 1, 2], &[3])?;
    let full_unitary_error = frobenius_norm(&(u_full.adjoint() * &u_full - eye(24)));
    let finite = rho_after
        .iter()
        .all(|z| z.re.is_finite() && z.im.is_finite())
        && [
            gross_work,
            switch_work,
            signed_net_work,
            conservative_net_work,
            mutual_information_before,
            mutual_information_after,
        ]
        .iter()
        .all(|v| v.is_finite());
    let mut checks = Vec::new();
    add_check(
        &mut checks,
        condition,
        "drive_at_t10_is_zero",
        drive_end_error,
        UNITARY_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "U_ext_unitarity",
        passive.unitary_error,
        UNITARY_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "U_full_unitarity",
        full_unitary_error,
        UNITARY_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "passive_mapping",
        passive.passive_mapping_error,
        STATE_AGREEMENT_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "trace_after",
        trace_error_after,
        TRACE_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "hermiticity_after",
        hermiticity_error_after,
        HERMITICITY_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "positivity_after",
        (-minimum_eigenvalue_after).max(0.0),
        POSITIVITY_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "full_eigenvalues_preserved",
        max_difference(&full_eigen_before, &full_eigen_after),
        STATE_AGREEMENT_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "load_eigenvalues_preserved",
        max_difference(
            &passive
                .state_eigenvalues_descending
                .iter()
                .rev()
                .copied()
                .collect::<Vec<_>>(),
            &load_eigen_after,
        ),
        STATE_AGREEMENT_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "chain_reduced_state_unchanged",
        frobenius_norm(&(rho_chain_before - rho_chain_after)),
        STATE_AGREEMENT_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "load_full_vs_direct_mapping",
        frobenius_norm(&(&rho_load_after - &passive.transformed_state)),
        STATE_AGREEMENT_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "post_load_ergotropy",
        load_after.ergotropy.abs(),
        POST_ERGOTROPY_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "gross_vs_prior_ergotropy",
        (gross_work - load_before.ergotropy).abs(),
        ENERGY_AGREEMENT_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "H_off_drop_vs_gross",
        ((full_energy_off_before_extraction - full_energy_off_after_extraction) - gross_work).abs(),
        ENERGY_AGREEMENT_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "load_drop_vs_gross",
        (load_energy_drop - gross_work).abs(),
        ENERGY_AGREEMENT_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "switch_work_identity",
        switch_identity_error,
        ENERGY_AGREEMENT_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "mutual_information_invariant",
        (mutual_information_before - mutual_information_after).abs(),
        MUTUAL_INFORMATION_TOLERANCE,
    );
    add_check(
        &mut checks,
        condition,
        "finite",
        if finite { 0.0 } else { f64::INFINITY },
        0.0,
    );
    let physical = trace_error_after <= TRACE_TOLERANCE
        && hermiticity_error_after <= HERMITICITY_TOLERANCE
        && minimum_eigenvalue_after >= -POSITIVITY_TOLERANCE
        && finite;
    Ok(ExtractionResult {
        condition,
        omega,
        gamma,
        load_energy_before: load_before.energy,
        load_ergotropy_before: load_before.ergotropy,
        load_passive_energy_before: load_before.passive_energy,
        load_purity_before: (&rho_load * &rho_load).trace().re,
        load_coherence_before: coherence_l1(&rho_load),
        populations_before: (0..3).map(|i| rho_load[(i, i)].re).collect(),
        state_eigenvalues: passive.state_eigenvalues_descending,
        load_energies: passive.hamiltonian_energies_ascending,
        switch_work,
        switch_work_imaginary: switch_complex.im.abs(),
        gross_work,
        load_energy_after: load_after.energy,
        load_ergotropy_after: load_after.ergotropy,
        load_coherence_after: coherence_l1(&rho_load_after),
        signed_net_work,
        conservative_net_work,
        drive_energy_in: run.summary.drive_energy.energy_in,
        full_energy_on_before,
        full_energy_off_after_switch,
        full_energy_off_before_extraction,
        full_energy_off_after_extraction,
        trace_error_after,
        hermiticity_error_after,
        minimum_eigenvalue_after,
        mutual_information_before,
        mutual_information_after,
        finite,
        physical,
        checks,
    })
}

fn write_results(path: &str, results: &[ExtractionResult]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w,"condition,Omega,gamma_phi,load_energy_before,load_ergotropy_before,load_passive_energy_before,load_purity_before,load_coherence_l1_before,load_population_0_before,load_population_1_before,load_population_2_before,switch_work,switch_work_imaginary_part,gross_extracted_work,load_energy_after,load_ergotropy_after,load_coherence_l1_after,signed_net_work,conservative_net_work,drive_energy_in,gross_work_recovery_fraction,signed_net_work_recovery_fraction,conservative_net_work_recovery_fraction,full_energy_on_before_switch,full_energy_off_after_switch,full_energy_off_before_extraction,full_energy_off_after_extraction,trace_error_after,hermiticity_error_after,minimum_eigenvalue_after,mutual_information_before,mutual_information_after,finite,physical")?;
    for r in results {
        writeln!(w,"{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",r.condition,n(r.omega),n(r.gamma),n(r.load_energy_before),n(r.load_ergotropy_before),n(r.load_passive_energy_before),n(r.load_purity_before),n(r.load_coherence_before),n(r.populations_before[0]),n(r.populations_before[1]),n(r.populations_before[2]),n(r.switch_work),n(r.switch_work_imaginary),n(r.gross_work),n(r.load_energy_after),n(r.load_ergotropy_after),n(r.load_coherence_after),n(r.signed_net_work),n(r.conservative_net_work),n(r.drive_energy_in),ratio(r.gross_work,r.drive_energy_in),ratio(r.signed_net_work,r.drive_energy_in),ratio(r.conservative_net_work,r.drive_energy_in),n(r.full_energy_on_before),n(r.full_energy_off_after_switch),n(r.full_energy_off_before_extraction),n(r.full_energy_off_after_extraction),n(r.trace_error_after),n(r.hermiticity_error_after),n(r.minimum_eigenvalue_after),n(r.mutual_information_before),n(r.mutual_information_after),b(r.finite),b(r.physical))?;
    }
    Ok(())
}

fn write_checks(path: &str, results: &[ExtractionResult]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "condition,check,error,tolerance,pass")?;
    for r in results {
        for c in &r.checks {
            writeln!(
                w,
                "{},{},{},{},{}",
                c.condition,
                c.check,
                n(c.error),
                n(c.tolerance),
                b(c.pass)
            )?;
        }
    }
    Ok(())
}

fn write_mapping(path: &str, results: &[ExtractionResult]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w,"condition,rank,state_eigenvalue_descending,load_energy_ascending,passive_population_mapping")?;
    for r in results {
        for i in 0..r.state_eigenvalues.len() {
            writeln!(
                w,
                "{},{},{},{},r_{} -> epsilon_{}",
                r.condition,
                i,
                n(r.state_eigenvalues[i]),
                n(r.load_energies[i]),
                i,
                i
            )?;
        }
    }
    Ok(())
}

fn write_report(
    path: &str,
    a: &ExtractionResult,
    bres: &ExtractionResult,
    recon_errors: [f64; 4],
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w,"# Milestone 6a explicit load work extraction\n\n## 実装済み\n\nMilestone 5cではfull density matrixを保存していなかったため、同一設定と`dt=0.0025`で確定状態を決定論的に再構成した。新しい探索・matching・実験ではない。`H_drive(10)=0`を明示検算し、突然切断後に理想的瞬時load-local unitaryを適用した。\n")?;
    writeln!(w,"再構成差 A energy/ergotropy `{:.3e}`/`{:.3e}`, B `{:.3e}`/`{:.3e}`。許容値 `{RECONSTRUCTION_TOLERANCE:e}`。\n",recon_errors[0],recon_errors[1],recon_errors[2],recon_errors[3])?;
    writeln!(w, "| quantity | A | B |\n|---|---:|---:|")?;
    for (name, left, right) in [
        (
            "load energy before",
            a.load_energy_before,
            bres.load_energy_before,
        ),
        (
            "ergotropy before",
            a.load_ergotropy_before,
            bres.load_ergotropy_before,
        ),
        ("switch work", a.switch_work, bres.switch_work),
        ("gross extracted work", a.gross_work, bres.gross_work),
        ("signed net work", a.signed_net_work, bres.signed_net_work),
        (
            "conservative net work",
            a.conservative_net_work,
            bres.conservative_net_work,
        ),
        (
            "load energy after",
            a.load_energy_after,
            bres.load_energy_after,
        ),
        (
            "ergotropy after",
            a.load_ergotropy_after,
            bres.load_ergotropy_after,
        ),
    ] {
        writeln!(w, "| {name} | {left:.10e} | {right:.10e} |")?;
    }
    writeln!(w,"\n- gross work A-B `{:.10e}`, ratio `{:.8}`\n- switch work A-B `{:.10e}`\n- signed net work A-B `{:.10e}`, ratio `{:.8}`\n- conservative net work A-B `{:.10e}`, ratio `{:.8}`\n- gross work recovery fraction A/B `{:.8}` / `{:.8}`\n- conservative work recovery fraction A/B `{:.8}` / `{:.8}`\n",a.gross_work-bres.gross_work,a.gross_work/bres.gross_work,a.switch_work-bres.switch_work,a.signed_net_work-bres.signed_net_work,a.signed_net_work/bres.signed_net_work,a.conservative_net_work-bres.conservative_net_work,a.conservative_net_work/bres.conservative_net_work,a.gross_work/a.drive_energy_in,bres.gross_work/bres.drive_energy_in,a.conservative_net_work/a.drive_energy_in,bres.conservative_net_work/bres.drive_energy_in)?;
    writeln!(
        w,
        "全検算: A `{}/{} PASS`, B `{}/{} PASS`。相互情報量変化 A `{:.3e}`, B `{:.3e}`。\n",
        a.checks.iter().filter(|c| c.pass).count(),
        a.checks.len(),
        bres.checks.iter().filter(|c| c.pass).count(),
        bres.checks.len(),
        (a.mutual_information_before - a.mutual_information_after).abs(),
        (bres.mutual_information_before - bres.mutual_information_after).abs()
    )?;
    writeln!(w,"閾値: unitarity/trace/Hermiticity `1e-10`; positivity `-1e-8`; state/eigenvalue and energy agreement `1e-9`; post ergotropy `1e-9`; mutual information `1e-9`; degeneracy `{DEGENERACY_TOLERANCE:e}`; entropy negative-eigenvalue clip `{ENTROPY_EIGENVALUE_TOLERANCE:e}`。\n\n## 未確認\n\n抽出unitaryの制御費用、有限時間切断・抽出、繰り返し放電、連続運転、正味周期仕事、相関からのglobal extraction、古典比較、量子優位は未確認。switch workのsigned/conservative定義は理想化であり、装置費用を完全評価していない。")?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams::default();
    println!("reconstructing Milestone 5c A state");
    let a_run = run_coherent_drive(&params, config(OMEGA_A, 0.0))?;
    println!("reconstructing Milestone 5c B state");
    let b_run = run_coherent_drive(&params, config(OMEGA_B, 0.5))?;
    let refs = [
        0.054450767878,
        0.052798274942,
        0.054452946589,
        0.008284636248,
    ];
    let errors = [
        (a_run.summary.at_end.load_energy - refs[0]).abs(),
        (a_run.summary.at_end.load_ergotropy - refs[1]).abs(),
        (b_run.summary.at_end.load_energy - refs[2]).abs(),
        (b_run.summary.at_end.load_ergotropy - refs[3]).abs(),
    ];
    if errors.iter().any(|&e| e > RECONSTRUCTION_TOLERANCE) {
        return Err("Milestone 5c reconstruction mismatch".into());
    }
    let results = vec![
        extract("A", OMEGA_A, 0.0, &a_run, &params)?,
        extract("B", OMEGA_B, 0.5, &b_run, &params)?,
    ];
    write_results("explicit_load_extraction_results.csv", &results)?;
    write_checks("explicit_load_extraction_checks.csv", &results)?;
    write_mapping("explicit_load_extraction_mapping.csv", &results)?;
    let a = &results[0];
    let bres = &results[1];
    write_report("MILESTONE_6A_REPORT.md", &a, &bres, errors)?;
    println!(
        "A gross={:.10e}, switch={:.10e}; B gross={:.10e}, switch={:.10e}",
        a.gross_work, a.switch_work, bres.gross_work, bres.switch_work
    );
    Ok(())
}
