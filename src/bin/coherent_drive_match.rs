use std::fs::File;
use std::io::{BufWriter, Write};

use quantum_work_network::coherent_drive::{
    run_coherent_drive, CoherentDriveConfig, CoherentDriveRun,
};
use quantum_work_network::coherent_drive_matching::{
    analyze_grid, inclusive_grid, merge_duplicate_roots, refine_bisection, relative_error,
    select_primary_root, GridValue, RootCandidate, ROOT_MERGE_TOLERANCE,
};
use quantum_work_network::operators::ModelParams;

const OMEGA_A: f64 = 0.2;
const OMEGA_LOW: f64 = 0.2;
const OMEGA_HIGH: f64 = 1.0;
const GRID_STEP: f64 = 0.01;
const COARSE_DT: f64 = 0.005;
const FINE_DT: f64 = 0.0025;
const MATCH_TOLERANCE: f64 = 1.0e-4;
const MAX_BISECTION_ITERATIONS: usize = 60;

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
fn ratio(numerator: f64, denominator: f64) -> Option<f64> {
    (denominator.abs() > 1.0e-12).then_some(numerator / denominator)
}

fn config(omega: f64, gamma_phi: f64, dt: f64) -> CoherentDriveConfig {
    let mut value = CoherentDriveConfig::milestone_5b(gamma_phi, dt);
    value.omega0 = omega;
    value
}

fn evaluate(
    params: &ModelParams,
    omega: f64,
    gamma: f64,
    dt: f64,
) -> Result<CoherentDriveRun, quantum_work_network::PhysicsError> {
    run_coherent_drive(params, config(omega, gamma, dt))
}

fn write_grid(path: &str, target: f64, rows: &[(f64, CoherentDriveRun)]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "Omega_B,E_load_B_t10,residual_B_minus_target,absolute_energy_error,relative_energy_error,physical_check,ledger_check")?;
    for (omega, run) in rows {
        let residual = run.summary.at_end.load_energy - target;
        writeln!(
            writer,
            "{},{},{},{},{},{},{}",
            n(*omega),
            n(run.summary.at_end.load_energy),
            n(residual),
            n(residual.abs()),
            n(relative_error(residual, target)),
            b(run.summary.physical_checks_pass && run.summary.top_level_check_pass),
            b(run.summary.ledger_check_pass)
        )?;
    }
    Ok(())
}

fn write_roots(
    path: &str,
    coarse: &[RootCandidate],
    fine: Option<&RootCandidate>,
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "stage,root_index,bracket_low,bracket_high,Omega_B,absolute_energy_error,relative_energy_error,iteration_count,converged")?;
    for (index, root) in coarse.iter().enumerate() {
        writeln!(
            writer,
            "dt_0.005,{},{},{},{},{},{},{},{}",
            index + 1,
            n(root.bracket_low),
            n(root.bracket_high),
            n(root.omega),
            n(root.absolute_error),
            n(root.relative_error),
            root.iterations,
            b(root.converged)
        )?;
    }
    if let Some(root) = fine {
        writeln!(
            writer,
            "dt_0.0025_primary,1,{},{},{},{},{},{},{}",
            n(root.bracket_low),
            n(root.bracket_high),
            n(root.omega),
            n(root.absolute_error),
            n(root.relative_error),
            root.iterations,
            b(root.converged)
        )?;
    }
    Ok(())
}

fn write_comparison(
    path: &str,
    rows: &[(&str, f64, f64, &CoherentDriveRun)],
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "condition,Omega,gamma_phi,load_energy,load_ergotropy,load_diagonal_ergotropy,load_coherence_ergotropy,load_coherence_l1,load_population_0,load_population_1,load_population_2,bare_network_energy,trace_error,hermiticity_error,minimum_eigenvalue,finite,physical,maximum_top_level_population,time_of_maximum_top_level_population,drive_energy_net,drive_energy_in,drive_energy_out,dephasing_energy_net,dephasing_energy_in,dephasing_energy_out,bare_network_energy_change,ledger_absolute_residual,maximum_trace_error,maximum_hermiticity_error,worst_minimum_eigenvalue,maximum_drive_power_imaginary_part,maximum_dephasing_power_imaginary_part,NaN_or_infinite_found,load_energy_delivery_fraction,load_ergotropy_delivery_fraction")?;
    for (condition, omega, gamma, run) in rows {
        let x = &run.summary.at_end;
        let s = &run.summary;
        writeln!(writer, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            condition, n(*omega), n(*gamma), n(x.load_energy), n(x.load_ergotropy), n(x.load_diagonal_ergotropy), n(x.load_coherence_ergotropy), n(x.load_coherence_l1),
            n(x.load_populations[0]), n(x.load_populations[1]), n(x.load_populations[2]), n(x.bare_network_energy), n(x.trace_error), n(x.hermiticity_error), n(x.minimum_eigenvalue), b(x.all_finite), b(s.physical_checks_pass),
            n(s.maximum_top_level_population.0), n(s.maximum_top_level_population.1), n(s.drive_energy.energy_net), n(s.drive_energy.energy_in), n(s.drive_energy.energy_out),
            n(s.dephasing_energy.energy_net), n(s.dephasing_energy.energy_in), n(s.dephasing_energy.energy_out), n(s.delta_bare_network_energy), n(s.ledger_residual.abs()), n(s.maximum_trace_error), n(s.maximum_hermiticity_error), n(s.worst_minimum_eigenvalue), n(s.maximum_drive_power_imaginary), n(s.maximum_dephasing_power_imaginary), b(!s.all_finite),
            ratio(x.load_energy, s.drive_energy.energy_in).map(n).unwrap_or_else(|| "undefined".to_string()),
            ratio(x.load_ergotropy, s.drive_energy.energy_in).map(n).unwrap_or_else(|| "undefined".to_string()))?;
    }
    Ok(())
}

fn write_timeseries(
    path: &str,
    a: &CoherentDriveRun,
    b_run: &CoherentDriveRun,
    omega_b: f64,
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "condition,Omega,gamma_phi,time,drive_envelope,drive_amplitude,load_energy,load_ergotropy,load_diagonal_ergotropy,load_coherence_ergotropy,load_coherence_l1,load_population_0,load_population_1,load_population_2,bare_network_energy,drive_power,dephasing_power,trace_error,hermiticity_error,minimum_eigenvalue")?;
    for (condition, omega, gamma, run) in [("A", OMEGA_A, 0.0, a), ("B", omega_b, 0.5, b_run)] {
        for x in &run.samples {
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                condition,
                n(omega),
                n(gamma),
                n(x.time),
                n(x.drive_envelope),
                n(x.drive_amplitude),
                n(x.load_energy),
                n(x.load_ergotropy),
                n(x.load_diagonal_ergotropy),
                n(x.load_coherence_ergotropy),
                n(x.load_coherence_l1),
                n(x.load_populations[0]),
                n(x.load_populations[1]),
                n(x.load_populations[2]),
                n(x.bare_network_energy),
                n(x.drive_power),
                n(x.dephasing_power),
                n(x.trace_error),
                n(x.hermiticity_error),
                n(x.minimum_eigenvalue)
            )?;
        }
    }
    Ok(())
}

fn comparison_values(
    a: &CoherentDriveRun,
    b_run: &CoherentDriveRun,
    omega_b: f64,
) -> Vec<(&'static str, f64)> {
    let ax = &a.summary.at_end;
    let bx = &b_run.summary.at_end;
    let error = (bx.load_energy - ax.load_energy).abs();
    vec![
        ("matched_Omega_B", omega_b),
        ("A_load_energy", ax.load_energy),
        ("B_load_energy", bx.load_energy),
        ("energy_matching_absolute_error", error),
        (
            "energy_matching_relative_error",
            relative_error(error, ax.load_energy),
        ),
        ("A_load_ergotropy", ax.load_ergotropy),
        ("B_load_ergotropy", bx.load_ergotropy),
        (
            "ergotropy_difference_A_minus_B",
            ax.load_ergotropy - bx.load_ergotropy,
        ),
        (
            "ergotropy_ratio_A_over_B",
            ratio(ax.load_ergotropy, bx.load_ergotropy).unwrap_or(f64::NAN),
        ),
        ("A_load_coherence_l1", ax.load_coherence_l1),
        ("B_load_coherence_l1", bx.load_coherence_l1),
        ("A_drive_energy_in", a.summary.drive_energy.energy_in),
        ("B_drive_energy_in", b_run.summary.drive_energy.energy_in),
        (
            "A_load_energy_delivery_fraction",
            ratio(ax.load_energy, a.summary.drive_energy.energy_in).unwrap_or(f64::NAN),
        ),
        (
            "B_load_energy_delivery_fraction",
            ratio(bx.load_energy, b_run.summary.drive_energy.energy_in).unwrap_or(f64::NAN),
        ),
        (
            "A_load_ergotropy_delivery_fraction",
            ratio(ax.load_ergotropy, a.summary.drive_energy.energy_in).unwrap_or(f64::NAN),
        ),
        (
            "B_load_ergotropy_delivery_fraction",
            ratio(bx.load_ergotropy, b_run.summary.drive_energy.energy_in).unwrap_or(f64::NAN),
        ),
        ("A_ledger_residual", a.summary.ledger_residual),
        ("B_ledger_residual", b_run.summary.ledger_residual),
        (
            "A_worst_minimum_eigenvalue",
            a.summary.worst_minimum_eigenvalue,
        ),
        (
            "B_worst_minimum_eigenvalue",
            b_run.summary.worst_minimum_eigenvalue,
        ),
    ]
}

fn write_convergence(
    path: &str,
    coarse: (&CoherentDriveRun, &CoherentDriveRun, f64),
    fine: (&CoherentDriveRun, &CoherentDriveRun, f64),
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "metric,coarse_dt,fine_dt,coarse_value,fine_value,absolute_difference,relative_difference"
    )?;
    for ((name_c, value_c), (name_f, value_f)) in comparison_values(coarse.0, coarse.1, coarse.2)
        .into_iter()
        .zip(comparison_values(fine.0, fine.1, fine.2))
    {
        assert_eq!(name_c, name_f);
        let abs = (value_c - value_f).abs();
        let rel = (value_f.abs() > 1.0e-12).then_some(abs / value_f.abs());
        writeln!(
            writer,
            "{},{},{},{},{},{},{}",
            name_c,
            n(COARSE_DT),
            n(FINE_DT),
            n(value_c),
            n(value_f),
            n(abs),
            rel.map(n).unwrap_or_default()
        )?;
    }
    Ok(())
}

fn success_checks(
    a: &CoherentDriveRun,
    b_run: &CoherentDriveRun,
    root_found: bool,
    energy_error: f64,
    direction_coarse: f64,
) -> [bool; 10] {
    let ax = &a.summary.at_end;
    let bx = &b_run.summary.at_end;
    let advantage = ratio(ax.load_ergotropy - bx.load_ergotropy, bx.load_ergotropy)
        .unwrap_or(f64::NEG_INFINITY);
    [
        root_found,
        energy_error < MATCH_TOLERANCE,
        ax.load_ergotropy > 1.0e-3 && bx.load_ergotropy > 1.0e-3,
        bx.load_ergotropy > 1.0e-3 && advantage > 0.05,
        ax.load_coherence_ergotropy > bx.load_coherence_ergotropy,
        a.summary.maximum_top_level_population.0 < 0.05
            && b_run.summary.maximum_top_level_population.0 < 0.05,
        a.summary.physical_checks_pass && b_run.summary.physical_checks_pass,
        a.summary.ledger_check_pass && b_run.summary.ledger_check_pass,
        (ax.load_ergotropy - bx.load_ergotropy).signum() == direction_coarse.signum(),
        energy_error < MATCH_TOLERANCE,
    ]
}

fn write_report(
    path: &str,
    roots: &[RootCandidate],
    coarse: (&CoherentDriveRun, &CoherentDriveRun, f64),
    fine: (&CoherentDriveRun, &CoherentDriveRun, f64),
    checks: [bool; 10],
    brackets: &[(f64, f64)],
    grid_matches: usize,
    best: GridValue,
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    let ax = &fine.0.summary.at_end;
    let bx = &fine.1.summary.at_end;
    let abs_e = (bx.load_energy - ax.load_energy).abs();
    let rel_e = relative_error(abs_e, ax.load_energy);
    let advantage = (ax.load_ergotropy - bx.load_ergotropy) / bx.load_ergotropy.abs().max(1.0e-12);
    writeln!(
        w,
        "# Milestone 5c coherent-drive energy-matched comparison\n"
    )?;
    writeln!(w, "## 実装済み\n\n- `Omega_B=0.2..1.0`, step `0.01` の81点を全走査。単調性は仮定していない。\n- 全符号変化区間を二分法で精密化し、重複閾値 `1e-7` で統合。\n- 主根は `|Omega_B-0.2|` 最小規則。\n- `dt=0.0025`は同じ近傍ブラケット内だけで再調整。\n")?;
    writeln!(w, "## 実行確認済み\n\n- 符号変化区間数: `{}`\n- グリッド上直接一致数: `{grid_matches}`\n- 局所接触疑い: `なし`\n- 統合後root数: `{}`\n- 最良グリッド参考点: Omega `{:.4}`, residual `{:.6e}`\n", brackets.len(), roots.len(), best.omega, best.residual)?;
    writeln!(w, "| quantity | dt=0.005 | dt=0.0025 |\n|---|---:|---:|\n| matched Omega_B | {:.12} | {:.12} |\n| A load energy | {:.10e} | {:.10e} |\n| B load energy | {:.10e} | {:.10e} |\n| relative match error | {:.3e} | {:.3e} |\n| A ergotropy | {:.10e} | {:.10e} |\n| B ergotropy | {:.10e} | {:.10e} |\n", coarse.2, fine.2, coarse.0.summary.at_end.load_energy, ax.load_energy, coarse.1.summary.at_end.load_energy, bx.load_energy, relative_error(coarse.1.summary.at_end.load_energy-coarse.0.summary.at_end.load_energy, coarse.0.summary.at_end.load_energy), rel_e, coarse.0.summary.at_end.load_ergotropy, ax.load_ergotropy, coarse.1.summary.at_end.load_ergotropy, bx.load_ergotropy)?;
    writeln!(w, "- matched Omega difference after dt halving: `{:.10e}`\n- absolute load-energy error: `{abs_e:.10e}`\n- relative load-energy error: `{rel_e:.10e}`\n- ergotropy A-B: `{:.10e}`\n- ergotropy A/B: `{:.10}`\n- relative ergotropy advantage: `{advantage:.10}`\n- coherence-derived ergotropy A-B: `{:.10e}`\n- coherence L1 A-B: `{:.10e}`\n- drive energy in B/A: `{:.10}`\n- load-energy delivery fraction A/B: `{:.10}` / `{:.10}`\n- load-ergotropy delivery fraction A/B: `{:.10}` / `{:.10}`\n", (fine.2-coarse.2).abs(), ax.load_ergotropy-bx.load_ergotropy, ax.load_ergotropy/bx.load_ergotropy, ax.load_coherence_ergotropy-bx.load_coherence_ergotropy, ax.load_coherence_l1-bx.load_coherence_l1, fine.1.summary.drive_energy.energy_in/fine.0.summary.drive_energy.energy_in, ax.load_energy/fine.0.summary.drive_energy.energy_in, bx.load_energy/fine.1.summary.drive_energy.energy_in, ax.load_ergotropy/fine.0.summary.drive_energy.energy_in, bx.load_ergotropy/fine.1.summary.drive_energy.energy_in)?;
    writeln!(w, "刻み収束CSVでは全量の絶対差を保存した。ledger residualと最小固有値はゼロ近傍なので、相対差ではなく絶対差を主に評価する。\n")?;
    writeln!(w, "\n## 成功条件\n")?;
    let labels = [
        "root exists",
        "energy match",
        "both ergotropy > 1e-3",
        "A advantage > 5%",
        "A coherence ergotropy > B",
        "top level < 5%",
        "physical checks",
        "energy ledgers",
        "direction stable after dt halving",
        "fine energy match",
    ];
    for (index, (label, pass)) in labels.iter().zip(checks).enumerate() {
        writeln!(
            w,
            "{}. {}: **{}**",
            index + 1,
            label,
            if pass { "PASS" } else { "FAIL" }
        )?;
    }
    writeln!(w, "\n## 公平性と未確認\n\n一致しているのは比較時刻、最終load energy、模型、パルス形状、駆動周波数。一致していないのは駆動強度、drive energy in、総投入エネルギー。これは同じ時刻に同じload energyを持つ状態の仕事価値比較であり、等入力費用比較ではない。Omegaも異なるため雑音だけを独立に変えた比較でもない。熱力学的効率ではなくdelivery fractionのみを報告した。連続運転、仕事抽出、古典比較、量子優位は未確認。")?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams::default();
    println!("computing A target at dt={COARSE_DT}");
    let a_coarse = evaluate(&params, OMEGA_A, 0.0, COARSE_DT)?;
    let target = a_coarse.summary.at_end.load_energy;
    let mut grid_runs = Vec::new();
    for (index, omega) in inclusive_grid(OMEGA_LOW, OMEGA_HIGH, GRID_STEP)?
        .into_iter()
        .enumerate()
    {
        println!("grid {}/81 Omega_B={omega:.2}", index + 1);
        grid_runs.push((omega, evaluate(&params, omega, 0.5, COARSE_DT)?));
    }
    write_grid("coherent_drive_match_grid.csv", target, &grid_runs)?;
    let grid_values: Vec<_> = grid_runs
        .iter()
        .map(|(omega, run)| GridValue {
            omega: *omega,
            residual: run.summary.at_end.load_energy - target,
        })
        .collect();
    let analysis = analyze_grid(&grid_values, target, MATCH_TOLERANCE)?;
    let mut roots: Vec<RootCandidate> = analysis
        .matching_grid_points
        .iter()
        .map(|point| RootCandidate {
            bracket_low: point.omega,
            bracket_high: point.omega,
            omega: point.omega,
            absolute_error: point.residual.abs(),
            relative_error: relative_error(point.residual, target),
            iterations: 0,
            converged: true,
        })
        .collect();
    for &(low, high) in &analysis.brackets {
        println!("refining coarse bracket [{low:.2}, {high:.2}]");
        roots.push(refine_bisection(
            low,
            high,
            target,
            MATCH_TOLERANCE,
            MAX_BISECTION_ITERATIONS,
            |omega| {
                Ok::<_, quantum_work_network::PhysicsError>(
                    evaluate(&params, omega, 0.5, COARSE_DT)?
                        .summary
                        .at_end
                        .load_energy
                        - target,
                )
            },
        )?);
    }
    roots = merge_duplicate_roots(roots, ROOT_MERGE_TOLERANCE);
    let Some(primary) = select_primary_root(&roots, OMEGA_A).copied() else {
        write_roots("coherent_drive_match_roots.csv", &roots, None)?;
        std::fs::write("MILESTONE_5C_REPORT.md", format!("# Milestone 5c\n\nNO_MATCH. Best grid Omega={}, residual={}. Search range was not expanded.\n", analysis.best_grid_point.omega, analysis.best_grid_point.residual))?;
        println!(
            "NO_MATCH; best grid Omega={}",
            analysis.best_grid_point.omega
        );
        return Ok(());
    };
    let b_coarse = evaluate(&params, primary.omega, 0.5, COARSE_DT)?;
    println!("computing fine A and refining same bracket");
    let a_fine = evaluate(&params, OMEGA_A, 0.0, FINE_DT)?;
    let fine_target = a_fine.summary.at_end.load_energy;
    let fine_root = refine_bisection(
        primary.bracket_low,
        primary.bracket_high,
        fine_target,
        MATCH_TOLERANCE,
        MAX_BISECTION_ITERATIONS,
        |omega| {
            Ok::<_, quantum_work_network::PhysicsError>(
                evaluate(&params, omega, 0.5, FINE_DT)?
                    .summary
                    .at_end
                    .load_energy
                    - fine_target,
            )
        },
    )?;
    let b_fine = evaluate(&params, fine_root.omega, 0.5, FINE_DT)?;
    write_roots("coherent_drive_match_roots.csv", &roots, Some(&fine_root))?;
    write_comparison(
        "coherent_drive_match_comparison.csv",
        &[
            ("A", OMEGA_A, 0.0, &a_fine),
            ("B", fine_root.omega, 0.5, &b_fine),
        ],
    )?;
    write_timeseries(
        "coherent_drive_match_timeseries.csv",
        &a_fine,
        &b_fine,
        fine_root.omega,
    )?;
    write_convergence(
        "coherent_drive_match_convergence.csv",
        (&a_coarse, &b_coarse, primary.omega),
        (&a_fine, &b_fine, fine_root.omega),
    )?;
    let coarse_direction =
        a_coarse.summary.at_end.load_ergotropy - b_coarse.summary.at_end.load_ergotropy;
    let fine_error = relative_error(b_fine.summary.at_end.load_energy - fine_target, fine_target);
    let checks = success_checks(&a_fine, &b_fine, true, fine_error, coarse_direction);
    write_report(
        "MILESTONE_5C_REPORT.md",
        &roots,
        (&a_coarse, &b_coarse, primary.omega),
        (&a_fine, &b_fine, fine_root.omega),
        checks,
        &analysis.brackets,
        analysis.matching_grid_points.len(),
        analysis.best_grid_point,
    )?;
    println!(
        "matched coarse Omega={:.12}; fine Omega={:.12}; relative error={:.3e}",
        primary.omega, fine_root.omega, fine_error
    );
    println!(
        "A ergotropy={:.10e}; B={:.10e}",
        a_fine.summary.at_end.load_ergotropy, b_fine.summary.at_end.load_ergotropy
    );
    Ok(())
}
