use std::fs::File;
use std::io::{BufWriter, Write};

use quantum_work_network::diagnostics::integrate_diagnostic_powers;
use quantum_work_network::ergotropy::ergotropy;
use quantum_work_network::operators::ModelParams;
use quantum_work_network::protocol::{
    run_coherent_input_protocol_with_population, LoadTimePoint, PhysicalChecks, ProtocolResult,
};
use quantum_work_network::{ComplexMatrix, C64};

const P_A: f64 = 0.1;
const TIMES: [f64; 4] = [3.0, 5.0, 7.9, 10.0];
const DEPHASING_RATES: [f64; 4] = [0.1, 0.2, 0.5, 1.0];
const RELATIVE_TOLERANCE: f64 = 1.0e-4;

#[derive(Clone)]
struct Metrics {
    load_energy: f64,
    total_ergotropy: f64,
    diagonal_ergotropy: f64,
    coherence_ergotropy: f64,
    coherence_l1: f64,
    populations: Vec<f64>,
    checks: PhysicalChecks,
}

struct OutputRow {
    t_eval: f64,
    gamma_phi: f64,
    comparison_type: &'static str,
    status: &'static str,
    p_b: Option<f64>,
    a: Metrics,
    b: Option<Metrics>,
}

fn relative_difference(actual: f64, target: f64) -> f64 {
    (actual - target).abs() / target.abs().max(1.0e-14)
}

fn load_density_matrix(sample: &LoadTimePoint) -> ComplexMatrix {
    let dim = sample.load_level_populations.len();
    let mut rho = ComplexMatrix::zeros(dim, dim);
    for (level, &population) in sample.load_level_populations.iter().enumerate() {
        rho[(level, level)] = C64::new(population, 0.0);
    }
    let mut index = 0;
    for row in 0..dim {
        for col in (row + 1)..dim {
            let value = sample.load_off_diagonal_values[index];
            rho[(row, col)] = value;
            rho[(col, row)] = value.conj();
            index += 1;
        }
    }
    rho
}

fn local_load_hamiltonian(params: &ModelParams) -> ComplexMatrix {
    let mut h = ComplexMatrix::zeros(params.load_dim, params.load_dim);
    for level in 0..params.load_dim {
        h[(level, level)] = C64::new(level as f64 * params.omega_load, 0.0);
    }
    h
}

fn diagonal_ergotropy(
    params: &ModelParams,
    sample: &LoadTimePoint,
) -> Result<f64, Box<dyn std::error::Error>> {
    let rho = load_density_matrix(sample);
    let mut diagonal = ComplexMatrix::zeros(params.load_dim, params.load_dim);
    for level in 0..params.load_dim {
        diagonal[(level, level)] = rho[(level, level)];
    }
    Ok(ergotropy(&diagonal, &local_load_hamiltonian(params), 1.0e-9)?.ergotropy)
}

fn checks_through(
    result: &ProtocolResult,
    index: usize,
) -> Result<PhysicalChecks, Box<dyn std::error::Error>> {
    let diagnostics = &result.diagnostics[..=index];
    let powers = integrate_diagnostic_powers(diagnostics)?;
    let energy_change = diagnostics[index].total_energy - diagnostics[0].total_energy;
    let trace_ok = diagnostics
        .iter()
        .all(|sample| sample.trace_error <= 1.0e-8);
    let hermiticity_ok = diagnostics
        .iter()
        .all(|sample| sample.hermiticity_error <= 1.0e-8);
    let positivity_ok = diagnostics
        .iter()
        .all(|sample| sample.minimum_eigenvalue >= -1.0e-8);
    let scale = energy_change
        .abs()
        .max(powers.source_energy_net.abs())
        .max(powers.dephasing_energy_net.abs())
        .max(1.0);
    let residual = energy_change - powers.source_energy_net - powers.dephasing_energy_net;
    let energy_balance_ok = residual.abs() <= 5.0e-4 * scale;
    let maximum_top_population = result.load_time_series[..=index]
        .iter()
        .map(|sample| {
            *sample
                .load_level_populations
                .last()
                .expect("load level exists")
        })
        .fold(0.0_f64, f64::max);
    Ok(PhysicalChecks {
        trace_ok,
        hermiticity_ok,
        positivity_ok,
        energy_balance_ok,
        top_level_ok: maximum_top_population < 0.05,
    })
}

fn metrics_at(
    params: &ModelParams,
    result: &ProtocolResult,
    index: usize,
) -> Result<Metrics, Box<dyn std::error::Error>> {
    let sample = &result.load_time_series[index];
    let diagonal_work = diagonal_ergotropy(params, sample)?;
    Ok(Metrics {
        load_energy: sample.load_energy,
        total_ergotropy: sample.load_ergotropy,
        diagonal_ergotropy: diagonal_work,
        coherence_ergotropy: sample.load_ergotropy - diagonal_work,
        coherence_l1: sample.load_off_diagonal_l1,
        populations: sample.load_level_populations.clone(),
        checks: checks_through(result, index)?,
    })
}

fn final_metrics(
    params: &ModelParams,
    result: &ProtocolResult,
) -> Result<Metrics, Box<dyn std::error::Error>> {
    metrics_at(params, result, result.load_time_series.len() - 1)
}

fn bool_text(value: bool) -> &'static str {
    if value {
        "true"
    } else {
        "false"
    }
}

fn number(value: f64) -> String {
    format!("{value:.16e}")
}

fn optional_number(value: Option<f64>) -> String {
    value.map(number).unwrap_or_default()
}

fn write_results(path: &str, rows: &[OutputRow]) -> Result<(), Box<dyn std::error::Error>> {
    let headers = [
        "t_eval",
        "gamma_phi",
        "comparison_type",
        "status",
        "p_A",
        "p_B",
        "initial_energy_ratio_B_over_A",
        "A_load_energy",
        "B_load_energy",
        "energy_relative_difference",
        "A_total_ergotropy",
        "B_total_ergotropy",
        "A_diagonal_ergotropy",
        "B_diagonal_ergotropy",
        "A_coherence_derived_ergotropy",
        "B_coherence_derived_ergotropy",
        "A_coherence_l1",
        "B_coherence_l1",
        "A_level_0",
        "A_level_1",
        "A_level_2",
        "B_level_0",
        "B_level_1",
        "B_level_2",
        "ergotropy_ratio_A_over_B",
        "A_output_ergotropy_over_initial_energy",
        "B_output_ergotropy_over_initial_energy",
        "A_trace_ok",
        "A_hermiticity_ok",
        "A_positivity_ok",
        "A_energy_balance_ok",
        "A_top_level_ok",
        "A_all_physical_checks",
        "B_trace_ok",
        "B_hermiticity_ok",
        "B_positivity_ok",
        "B_energy_balance_ok",
        "B_top_level_ok",
        "B_all_physical_checks",
        "all_physical_checks",
    ];
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "{}", headers.join(","))?;
    for row in rows {
        let b = row.b.as_ref();
        let energy_difference =
            b.map(|value| relative_difference(value.load_energy, row.a.load_energy));
        let ratio = b.and_then(|value| {
            (value.total_ergotropy > 0.0).then_some(row.a.total_ergotropy / value.total_ergotropy)
        });
        let all_checks = b.map(|value| row.a.checks.all_pass() && value.checks.all_pass());
        let values = vec![
            number(row.t_eval),
            number(row.gamma_phi),
            row.comparison_type.to_string(),
            row.status.to_string(),
            number(P_A),
            optional_number(row.p_b),
            optional_number(row.p_b.map(|p| p / P_A)),
            number(row.a.load_energy),
            optional_number(b.map(|value| value.load_energy)),
            optional_number(energy_difference),
            number(row.a.total_ergotropy),
            optional_number(b.map(|value| value.total_ergotropy)),
            number(row.a.diagonal_ergotropy),
            optional_number(b.map(|value| value.diagonal_ergotropy)),
            number(row.a.coherence_ergotropy),
            optional_number(b.map(|value| value.coherence_ergotropy)),
            number(row.a.coherence_l1),
            optional_number(b.map(|value| value.coherence_l1)),
            number(row.a.populations[0]),
            number(row.a.populations[1]),
            number(row.a.populations[2]),
            optional_number(b.map(|value| value.populations[0])),
            optional_number(b.map(|value| value.populations[1])),
            optional_number(b.map(|value| value.populations[2])),
            optional_number(ratio),
            number(row.a.total_ergotropy / P_A),
            optional_number(b.zip(row.p_b).map(|(value, p)| value.total_ergotropy / p)),
            bool_text(row.a.checks.trace_ok).to_string(),
            bool_text(row.a.checks.hermiticity_ok).to_string(),
            bool_text(row.a.checks.positivity_ok).to_string(),
            bool_text(row.a.checks.energy_balance_ok).to_string(),
            bool_text(row.a.checks.top_level_ok).to_string(),
            bool_text(row.a.checks.all_pass()).to_string(),
            b.map(|value| bool_text(value.checks.trace_ok).to_string())
                .unwrap_or_default(),
            b.map(|value| bool_text(value.checks.hermiticity_ok).to_string())
                .unwrap_or_default(),
            b.map(|value| bool_text(value.checks.positivity_ok).to_string())
                .unwrap_or_default(),
            b.map(|value| bool_text(value.checks.energy_balance_ok).to_string())
                .unwrap_or_default(),
            b.map(|value| bool_text(value.checks.top_level_ok).to_string())
                .unwrap_or_default(),
            b.map(|value| bool_text(value.checks.all_pass()).to_string())
                .unwrap_or_default(),
            all_checks
                .map(|value| bool_text(value).to_string())
                .unwrap_or_default(),
        ];
        writeln!(writer, "{}", values.join(","))?;
    }
    Ok(())
}

fn write_summary(
    path: &str,
    rows: &[OutputRow],
    maximum_linearity_error: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let matched: Vec<_> = rows
        .iter()
        .filter(|row| row.comparison_type == "MATCHED" && row.status == "MATCH")
        .collect();
    let no_match_count = rows
        .iter()
        .filter(|row| row.comparison_type == "MATCHED" && row.status == "NO_MATCH")
        .count();
    let a_larger_count = matched
        .iter()
        .filter(|row| row.a.total_ergotropy > row.b.as_ref().unwrap().total_ergotropy)
        .count();
    let reversed_count = matched
        .iter()
        .filter(|row| row.a.total_ergotropy < row.b.as_ref().unwrap().total_ergotropy)
        .count();
    let mut ratios: Vec<f64> = matched
        .iter()
        .filter_map(|row| {
            let b = row.b.as_ref().unwrap();
            (b.total_ergotropy > 0.0).then_some(row.a.total_ergotropy / b.total_ergotropy)
        })
        .collect();
    ratios.sort_by(f64::total_cmp);
    let median = if ratios.is_empty() {
        f64::NAN
    } else if ratios.len() % 2 == 0 {
        0.5 * (ratios[ratios.len() / 2 - 1] + ratios[ratios.len() / 2])
    } else {
        ratios[ratios.len() / 2]
    };
    let matched_physical_failures = matched
        .iter()
        .filter(|row| !row.a.checks.all_pass() || !row.b.as_ref().unwrap().checks.all_pass())
        .count();
    let equal_physical_failures = rows
        .iter()
        .filter(|row| row.comparison_type == "EQUAL_INPUT")
        .filter(|row| !row.a.checks.all_pass() || !row.b.as_ref().unwrap().checks.all_pass())
        .count();
    let max_match_error = matched
        .iter()
        .map(|row| relative_difference(row.b.as_ref().unwrap().load_energy, row.a.load_energy))
        .fold(0.0_f64, f64::max);
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "metric,value")?;
    writeln!(writer, "match_count,{}", matched.len())?;
    writeln!(writer, "A_ergotropy_larger_count,{a_larger_count}")?;
    writeln!(writer, "reversed_count,{reversed_count}")?;
    writeln!(writer, "no_match_count,{no_match_count}")?;
    writeln!(
        writer,
        "minimum_A_over_B_ratio,{:.16e}",
        ratios.first().copied().unwrap_or(f64::NAN)
    )?;
    writeln!(writer, "median_A_over_B_ratio,{median:.16e}")?;
    writeln!(
        writer,
        "maximum_A_over_B_ratio,{:.16e}",
        ratios.last().copied().unwrap_or(f64::NAN)
    )?;
    writeln!(
        writer,
        "matched_physical_check_failure_count,{matched_physical_failures}"
    )?;
    writeln!(
        writer,
        "equal_input_physical_check_failure_count,{equal_physical_failures}"
    )?;
    writeln!(
        writer,
        "maximum_match_relative_error,{max_match_error:.16e}"
    )?;
    writeln!(
        writer,
        "maximum_verified_load_energy_linearity_error,{maximum_linearity_error:.16e}"
    )?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams {
        hopping_j: 1.0,
        coupling_g: 0.25,
        ..ModelParams::default()
    };
    let a_full = run_coherent_input_protocol_with_population(&params, P_A, 0.0, 10.0, 100)?;
    let mut rows = Vec::with_capacity(32);
    let mut maximum_linearity_error = 0.0_f64;

    for gamma in DEPHASING_RATES {
        println!("gamma_phi={gamma:.1}: preparing equal-input and p_B=1 reference runs");
        let equal_full =
            run_coherent_input_protocol_with_population(&params, P_A, gamma, 10.0, 100)?;
        let unit_full =
            run_coherent_input_protocol_with_population(&params, 1.0, gamma, 10.0, 100)?;
        for t_eval in TIMES {
            let index = (t_eval * 10.0).round() as usize;
            let a = metrics_at(&params, &a_full, index)?;
            let equal_b = metrics_at(&params, &equal_full, index)?;
            let unit_b = metrics_at(&params, &unit_full, index)?;
            let expected_equal_energy = P_A * unit_b.load_energy;
            maximum_linearity_error = maximum_linearity_error.max(
                (equal_b.load_energy - expected_equal_energy).abs()
                    / expected_equal_energy.abs().max(1.0e-14),
            );
            rows.push(OutputRow {
                t_eval,
                gamma_phi: gamma,
                comparison_type: "EQUAL_INPUT",
                status: "EQUAL_INPUT",
                p_b: Some(P_A),
                a: a.clone(),
                b: Some(equal_b),
            });

            if unit_b.load_energy <= 1.0e-14 {
                rows.push(OutputRow {
                    t_eval,
                    gamma_phi: gamma,
                    comparison_type: "MATCHED",
                    status: "NO_MATCH",
                    p_b: None,
                    a,
                    b: None,
                });
                continue;
            }
            let p_b = a.load_energy / unit_b.load_energy;
            if !(0.0..=1.0).contains(&p_b) {
                rows.push(OutputRow {
                    t_eval,
                    gamma_phi: gamma,
                    comparison_type: "MATCHED",
                    status: "NO_MATCH",
                    p_b: None,
                    a,
                    b: None,
                });
                continue;
            }
            let steps = index;
            let matched_result =
                run_coherent_input_protocol_with_population(&params, p_b, gamma, t_eval, steps)?;
            let matched_b = final_metrics(&params, &matched_result)?;
            let match_error = relative_difference(matched_b.load_energy, a.load_energy);
            let status = if match_error < RELATIVE_TOLERANCE {
                "MATCH"
            } else {
                "NO_MATCH"
            };
            println!(
                "  t={t_eval:.1}: p_B={p_b:.9}, relative error={match_error:.3e}, status={status}"
            );
            rows.push(OutputRow {
                t_eval,
                gamma_phi: gamma,
                comparison_type: "MATCHED",
                status,
                p_b: (status == "MATCH").then_some(p_b),
                a,
                b: (status == "MATCH").then_some(matched_b),
            });
        }
    }

    rows.sort_by(|left, right| {
        left.t_eval
            .total_cmp(&right.t_eval)
            .then(left.gamma_phi.total_cmp(&right.gamma_phi))
            .then(left.comparison_type.cmp(right.comparison_type))
    });
    write_results("coherent_robustness_results.csv", &rows)?;
    write_summary(
        "coherent_robustness_summary.csv",
        &rows,
        maximum_linearity_error,
    )?;
    Ok(())
}
