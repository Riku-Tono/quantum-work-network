use std::fs::File;
use std::io::{BufWriter, Write};

use quantum_work_network::coherent_drive::{
    run_coherent_drive, sample_at, CoherentDriveRun, CoherentDriveSample,
    CONVERGENCE_ABSOLUTE_TOLERANCE, CONVERGENCE_RELATIVE_TOLERANCE, HERMITICITY_TOLERANCE,
    LEDGER_ABSOLUTE_TOLERANCE, LEDGER_RELATIVE_TOLERANCE, POSITIVITY_TOLERANCE, SIGNAL_TOLERANCE,
    TOP_LEVEL_LIMIT, TRACE_TOLERANCE,
};
use quantum_work_network::operators::ModelParams;

const MAIN_DT: f64 = 0.005;
const COARSE_DT: f64 = 0.01;

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

fn write_timeseries(path: &str, runs: &[(&str, &CoherentDriveRun)]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "condition,time,drive_envelope,drive_amplitude,load_energy,load_ergotropy,load_diagonal_ergotropy,load_coherence_ergotropy,load_coherence_l1,load_population_0,load_population_1,load_population_2,bare_network_energy,drive_power,dephasing_power,trace_error,hermiticity_error,minimum_eigenvalue")?;
    for (condition, run) in runs {
        for s in &run.samples {
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                condition,
                n(s.time),
                n(s.drive_envelope),
                n(s.drive_amplitude),
                n(s.load_energy),
                n(s.load_ergotropy),
                n(s.load_diagonal_ergotropy),
                n(s.load_coherence_ergotropy),
                n(s.load_coherence_l1),
                n(s.load_populations[0]),
                n(s.load_populations[1]),
                n(s.load_populations[2]),
                n(s.bare_network_energy),
                n(s.drive_power),
                n(s.dephasing_power),
                n(s.trace_error),
                n(s.hermiticity_error),
                n(s.minimum_eigenvalue)
            )?;
        }
    }
    Ok(())
}

fn write_summary(path: &str, runs: &[(&str, &CoherentDriveRun)]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "condition,dt,gamma_phi,max_load_energy,max_load_energy_time,max_load_ergotropy,max_load_ergotropy_time,max_load_coherence_ergotropy,max_load_coherence_ergotropy_time,max_load_coherence_l1,max_load_coherence_l1_time,max_top_level_population,max_top_level_population_time,tau_load_energy,tau_load_ergotropy,tau_load_coherence_ergotropy,tau_load_coherence_l1,end_load_energy,end_load_ergotropy,end_load_coherence_ergotropy,end_load_coherence_l1,drive_energy_net,drive_energy_in,drive_energy_out,dephasing_energy_net,dephasing_energy_in,dephasing_energy_out,delta_bare_network_energy,ledger_residual,max_trace_error,max_hermiticity_error,worst_minimum_eigenvalue,max_drive_power_imaginary,max_dephasing_power_imaginary,all_finite,physical_checks_pass,ledger_check_pass,top_level_check_pass,all_checks_pass")?;
    for (condition, run) in runs {
        let s = &run.summary;
        writeln!(writer, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            condition, n(run.config.dt), n(run.config.gamma_phi),
            n(s.maximum_load_energy.0), n(s.maximum_load_energy.1),
            n(s.maximum_load_ergotropy.0), n(s.maximum_load_ergotropy.1),
            n(s.maximum_load_coherence_ergotropy.0), n(s.maximum_load_coherence_ergotropy.1),
            n(s.maximum_load_coherence_l1.0), n(s.maximum_load_coherence_l1.1),
            n(s.maximum_top_level_population.0), n(s.maximum_top_level_population.1),
            n(s.at_tau.load_energy), n(s.at_tau.load_ergotropy), n(s.at_tau.load_coherence_ergotropy), n(s.at_tau.load_coherence_l1),
            n(s.at_end.load_energy), n(s.at_end.load_ergotropy), n(s.at_end.load_coherence_ergotropy), n(s.at_end.load_coherence_l1),
            n(s.drive_energy.energy_net), n(s.drive_energy.energy_in), n(s.drive_energy.energy_out),
            n(s.dephasing_energy.energy_net), n(s.dephasing_energy.energy_in), n(s.dephasing_energy.energy_out),
            n(s.delta_bare_network_energy), n(s.ledger_residual), n(s.maximum_trace_error),
            n(s.maximum_hermiticity_error), n(s.worst_minimum_eigenvalue),
            n(s.maximum_drive_power_imaginary), n(s.maximum_dephasing_power_imaginary),
            b(s.all_finite), b(s.physical_checks_pass), b(s.ledger_check_pass), b(s.top_level_check_pass),
            b(s.physical_checks_pass && s.ledger_check_pass && s.top_level_check_pass))?;
    }
    Ok(())
}

fn metric_values(run: &CoherentDriveRun) -> [(&'static str, f64); 8] {
    let s = &run.summary;
    [
        ("tau_load_energy", s.at_tau.load_energy),
        ("tau_load_ergotropy", s.at_tau.load_ergotropy),
        ("end_load_energy", s.at_end.load_energy),
        ("end_load_ergotropy", s.at_end.load_ergotropy),
        ("maximum_load_ergotropy", s.maximum_load_ergotropy.0),
        ("drive_energy_net", s.drive_energy.energy_net),
        ("dephasing_energy_net", s.dephasing_energy.energy_net),
        ("ledger_residual", s.ledger_residual),
    ]
}

fn relative_difference(coarse: f64, fine: f64) -> Option<f64> {
    (fine.abs() > 1.0e-12).then_some((coarse - fine).abs() / fine.abs())
}

fn converged(coarse: f64, fine: f64) -> bool {
    let absolute = (coarse - fine).abs();
    absolute <= CONVERGENCE_ABSOLUTE_TOLERANCE
        || relative_difference(coarse, fine)
            .is_some_and(|relative| relative <= CONVERGENCE_RELATIVE_TOLERANCE)
}

fn write_convergence(
    path: &str,
    pairs: &[(&str, &CoherentDriveRun, &CoherentDriveRun)],
) -> std::io::Result<bool> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "condition,metric,dt_coarse,dt_fine,coarse_value,fine_value,absolute_difference,relative_difference,converged")?;
    let mut all_converged = true;
    for (condition, coarse, fine) in pairs {
        for ((name_c, value_c), (name_f, value_f)) in
            metric_values(coarse).into_iter().zip(metric_values(fine))
        {
            assert_eq!(name_c, name_f);
            let absolute = (value_c - value_f).abs();
            let relative = relative_difference(value_c, value_f);
            let pass = converged(value_c, value_f);
            all_converged &= pass;
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{}",
                condition,
                name_c,
                n(coarse.config.dt),
                n(fine.config.dt),
                n(value_c),
                n(value_f),
                n(absolute),
                relative.map(n).unwrap_or_default(),
                b(pass)
            )?;
        }
    }
    Ok(all_converged)
}

fn sample_comparison(
    label: &str,
    a: &CoherentDriveSample,
    b_sample: &CoherentDriveSample,
) -> String {
    format!(
        "- {label}: load energy A `{:.10e}`, B `{:.10e}`; ergotropy A `{:.10e}`, B `{:.10e}`; coherence-derived ergotropy A `{:.10e}`, B `{:.10e}`; coherence L1 A `{:.10e}`, B `{:.10e}`\n",
        a.load_energy,
        b_sample.load_energy,
        a.load_ergotropy,
        b_sample.load_ergotropy,
        a.load_coherence_ergotropy,
        b_sample.load_coherence_ergotropy,
        a.load_coherence_l1,
        b_sample.load_coherence_l1
    )
}

fn write_report(
    path: &str,
    a: &CoherentDriveRun,
    b_run: &CoherentDriveRun,
    convergence_ok: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let a_max_time = a.summary.maximum_load_ergotropy.1;
    let a_at_max = sample_at(&a.samples, a_max_time)?;
    let b_at_a_max = sample_at(&b_run.samples, a_max_time)?;
    let reduction = [
        (&a.summary.at_tau, &b_run.summary.at_tau),
        (&a.summary.at_end, &b_run.summary.at_end),
        (a_at_max, b_at_a_max),
    ]
    .iter()
    .any(|(left, right)| {
        right.load_ergotropy + SIGNAL_TOLERANCE < left.load_ergotropy
            || right.load_coherence_l1 + SIGNAL_TOLERANCE < left.load_coherence_l1
    });
    let checks = [
        a.summary.maximum_load_coherence_l1.0 > SIGNAL_TOLERANCE,
        a.summary.maximum_load_ergotropy.0 > SIGNAL_TOLERANCE,
        reduction,
        a.summary.top_level_check_pass && b_run.summary.top_level_check_pass,
        a.summary.physical_checks_pass && b_run.summary.physical_checks_pass,
        a.summary.ledger_check_pass && b_run.summary.ledger_check_pass,
        convergence_ok,
    ];
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "# Milestone 5b coherent-drive sanity check\n")?;
    writeln!(writer, "## 実装済み\n")?;
    writeln!(writer, "- 真空初期状態、既存24次元bare network `H0`。")?;
    writeln!(
        writer,
        "- `Omega0=0.2`, `omega_drive=1`, `tau=3.2`, `t_end=10` の単一 `sin^2` パルス。"
    )?;
    writeln!(writer, "- A: `gamma_phi=0`; B: `gamma_phi=0.5`。Bは3サイトすべてへ `sqrt(gamma_phi/2) sigma_z,j`。loadへの直接雑音なし。")?;
    writeln!(
        writer,
        "- Lindblad励起注入、energy matching、探索、最適化は未使用。\n"
    )?;
    writeln!(writer, "## 実行確認済み\n")?;
    writeln!(writer, "主計算は最大RK4刻み `0.005`、CSV間隔 `0.01`。\n")?;
    for (label, run) in [("A", a), ("B", b_run)] {
        let s = &run.summary;
        writeln!(writer, "### {label}\n")?;
        writeln!(
            writer,
            "- 最大load energy: `{:.10e}` at `t={:.2}`",
            s.maximum_load_energy.0, s.maximum_load_energy.1
        )?;
        writeln!(
            writer,
            "- 最大load ergotropy: `{:.10e}` at `t={:.2}`",
            s.maximum_load_ergotropy.0, s.maximum_load_ergotropy.1
        )?;
        writeln!(
            writer,
            "- 最大coherence-derived ergotropy: `{:.10e}` at `t={:.2}`",
            s.maximum_load_coherence_ergotropy.0, s.maximum_load_coherence_ergotropy.1
        )?;
        writeln!(
            writer,
            "- 最大coherence L1: `{:.10e}` at `t={:.2}`",
            s.maximum_load_coherence_l1.0, s.maximum_load_coherence_l1.1
        )?;
        writeln!(
            writer,
            "- 最大最上段占有率: `{:.10e}` at `t={:.2}`",
            s.maximum_top_level_population.0, s.maximum_top_level_population.1
        )?;
        writeln!(
            writer,
            "- drive energy net/in/out: `{:.10e}` / `{:.10e}` / `{:.10e}`",
            s.drive_energy.energy_net, s.drive_energy.energy_in, s.drive_energy.energy_out
        )?;
        writeln!(
            writer,
            "- dephasing energy net/in/out: `{:.10e}` / `{:.10e}` / `{:.10e}`",
            s.dephasing_energy.energy_net,
            s.dephasing_energy.energy_in,
            s.dephasing_energy.energy_out
        )?;
        writeln!(
            writer,
            "- ledger absolute residual: `{:.10e}`",
            s.ledger_residual.abs()
        )?;
        writeln!(
            writer,
            "- max trace/Hermiticity errors: `{:.3e}` / `{:.3e}`",
            s.maximum_trace_error, s.maximum_hermiticity_error
        )?;
        writeln!(
            writer,
            "- worst minimum eigenvalue: `{:.3e}`",
            s.worst_minimum_eigenvalue
        )?;
        writeln!(
            writer,
            "- max power imaginary parts drive/dephasing: `{:.3e}` / `{:.3e}`",
            s.maximum_drive_power_imaginary, s.maximum_dephasing_power_imaginary
        )?;
        writeln!(
            writer,
            "- finite={}, physical={}, ledger={}, top-level={}\n",
            s.all_finite, s.physical_checks_pass, s.ledger_check_pass, s.top_level_check_pass
        )?;
    }
    writeln!(writer, "## 同時刻A/B確認\n")?;
    writeln!(
        writer,
        "- 最大load ergotropy A/B: `{:.10e}` / `{:.10e}` (比 `{:.6}`)",
        a.summary.maximum_load_ergotropy.0,
        b_run.summary.maximum_load_ergotropy.0,
        a.summary.maximum_load_ergotropy.0 / b_run.summary.maximum_load_ergotropy.0
    )?;
    writeln!(
        writer,
        "- 最大load coherence L1 A/B: `{:.10e}` / `{:.10e}` (比 `{:.6}`)",
        a.summary.maximum_load_coherence_l1.0,
        b_run.summary.maximum_load_coherence_l1.0,
        a.summary.maximum_load_coherence_l1.0 / b_run.summary.maximum_load_coherence_l1.0
    )?;
    write!(
        writer,
        "{}",
        sample_comparison("t=tau", &a.summary.at_tau, &b_run.summary.at_tau)
    )?;
    write!(
        writer,
        "{}",
        sample_comparison("t=10", &a.summary.at_end, &b_run.summary.at_end)
    )?;
    write!(
        writer,
        "{}",
        sample_comparison(
            &format!("A最大ergotropy時刻 t={a_max_time:.2}"),
            a_at_max,
            b_at_a_max
        )
    )?;
    writeln!(writer, "\nA/Bそれぞれの最大値は時刻が異なる可能性があり、その比は公平な同時刻比較ではなくsanity check用。今回は同一load energyへのmatchingもしていない。\n")?;
    writeln!(writer, "刻み収束の全16指標はCSVに絶対差と相対差を保存した。ledger residualは分母がほぼゼロなので、相対差ではなく絶対差を主判定に使用した。\n")?;
    writeln!(writer, "## 成功確認\n")?;
    for (index, (description, pass)) in [
        ("Aの非ゼロload coherence", checks[0]),
        ("Aの非ゼロload ergotropy", checks[1]),
        ("同時刻でBのcoherenceまたはergotropy低下", checks[2]),
        ("A/Bの最上段占有率5%未満", checks[3]),
        ("A/Bの物理チェック", checks[4]),
        ("A/Bの固定H0エネルギー台帳", checks[5]),
        ("dt半減による主要量収束", checks[6]),
    ]
    .iter()
    .enumerate()
    {
        writeln!(
            writer,
            "{}. {}: **{}**",
            index + 1,
            description,
            if *pass { "PASS" } else { "FAIL" }
        )?;
    }
    writeln!(writer, "\n閾値: signal `>{SIGNAL_TOLERANCE:e}`; trace/Hermiticity `<={TRACE_TOLERANCE:e}`/`<={HERMITICITY_TOLERANCE:e}`; minimum eigenvalue `>=-{POSITIVITY_TOLERANCE:e}`; top level `<{TOP_LEVEL_LIMIT}`; ledger `|r| <= {LEDGER_ABSOLUTE_TOLERANCE:e} + {LEDGER_RELATIVE_TOLERANCE:e}*scale`; convergence `abs<={CONVERGENCE_ABSOLUTE_TOLERANCE:e}` または `rel<={CONVERGENCE_RELATIVE_TOLERANCE:e}`。\n")?;
    writeln!(writer, "## 未確認\n\n同一時刻・同一load energyの公平比較、energy matching、パラメータ探索、連続駆動、仕事抽出、古典比較、量子優位は未確認。")?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams::default();
    println!("running A coarse dt={COARSE_DT}");
    let a_coarse = run_coherent_drive(
        &params,
        quantum_work_network::coherent_drive::CoherentDriveConfig::milestone_5b(0.0, COARSE_DT),
    )?;
    println!("running A fine dt={MAIN_DT}");
    let a_fine = run_coherent_drive(
        &params,
        quantum_work_network::coherent_drive::CoherentDriveConfig::milestone_5b(0.0, MAIN_DT),
    )?;
    println!("running B coarse dt={COARSE_DT}");
    let b_coarse = run_coherent_drive(
        &params,
        quantum_work_network::coherent_drive::CoherentDriveConfig::milestone_5b(0.5, COARSE_DT),
    )?;
    println!("running B fine dt={MAIN_DT}");
    let b_fine = run_coherent_drive(
        &params,
        quantum_work_network::coherent_drive::CoherentDriveConfig::milestone_5b(0.5, MAIN_DT),
    )?;

    assert_eq!(a_fine.samples.len(), b_fine.samples.len());
    assert!(a_fine
        .samples
        .iter()
        .zip(&b_fine.samples)
        .all(|(a, b)| a.time == b.time));
    write_timeseries(
        "coherent_drive_timeseries.csv",
        &[("A", &a_fine), ("B", &b_fine)],
    )?;
    write_summary(
        "coherent_drive_summary.csv",
        &[
            ("A_coarse", &a_coarse),
            ("A", &a_fine),
            ("B_coarse", &b_coarse),
            ("B", &b_fine),
        ],
    )?;
    let convergence_ok = write_convergence(
        "coherent_drive_convergence.csv",
        &[("A", &a_coarse, &a_fine), ("B", &b_coarse, &b_fine)],
    )?;
    write_report("MILESTONE_5B_REPORT.md", &a_fine, &b_fine, convergence_ok)?;
    println!(
        "wrote coherent_drive_timeseries.csv ({} rows)",
        a_fine.samples.len() + b_fine.samples.len()
    );
    println!(
        "A max ergotropy={:.10e}; B={:.10e}; convergence={}",
        a_fine.summary.maximum_load_ergotropy.0,
        b_fine.summary.maximum_load_ergotropy.0,
        convergence_ok
    );
    Ok(())
}
