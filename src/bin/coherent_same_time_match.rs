use std::fs::File;
use std::io::{BufWriter, Write};

use quantum_work_network::ergotropy::ergotropy;
use quantum_work_network::operators::ModelParams;
use quantum_work_network::protocol::{
    run_coherent_input_protocol_with_population, LoadTimePoint, ProtocolResult,
};
use quantum_work_network::{ComplexMatrix, C64};

const P_A: f64 = 0.2;
const GAMMA_A: f64 = 0.0;
const GAMMA_B: f64 = 0.5;
const T_EVAL: f64 = 7.9;
const TIME_STEPS: usize = 79;
const RELATIVE_TOLERANCE: f64 = 1.0e-4;
const GRID_STEPS: usize = 10;
const MAX_BISECTION_ITERATIONS: usize = 40;

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
    result: &ProtocolResult,
) -> Result<f64, Box<dyn std::error::Error>> {
    let sample = result.load_time_series.last().expect("final sample exists");
    let rho = load_density_matrix(sample);
    let mut diagonal = ComplexMatrix::zeros(params.load_dim, params.load_dim);
    for level in 0..params.load_dim {
        diagonal[(level, level)] = rho[(level, level)];
    }
    Ok(ergotropy(&diagonal, &local_load_hamiltonian(params), 1.0e-9)?.ergotropy)
}

fn write_grid_csv(
    path: &str,
    target: f64,
    grid: &[(f64, ProtocolResult)],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "p_B,B_load_energy,B_minus_A_energy,relative_energy_difference"
    )?;
    for (p, result) in grid {
        writeln!(
            writer,
            "{p:.16e},{:.16e},{:.16e},{:.16e}",
            result.final_load_energy,
            result.final_load_energy - target,
            relative_difference(result.final_load_energy, target),
        )?;
    }
    Ok(())
}

fn write_result_csv(
    path: &str,
    params: &ModelParams,
    p_b: f64,
    a: &ProtocolResult,
    b: &ProtocolResult,
) -> Result<(), Box<dyn std::error::Error>> {
    let a_sample = a.load_time_series.last().expect("final A sample exists");
    let b_sample = b.load_time_series.last().expect("final B sample exists");
    let a_diag = diagonal_ergotropy(params, a)?;
    let b_diag = diagonal_ergotropy(params, b)?;
    let energy_relative_difference = relative_difference(b.final_load_energy, a.final_load_energy);
    let ergotropy_difference = a.final_load_ergotropy - b.final_load_ergotropy;
    let ergotropy_ratio = a.final_load_ergotropy / b.final_load_ergotropy;
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "condition,p,initial_energy,t_eval,load_energy,load_ergotropy,diagonal_ergotropy,coherence_derived_ergotropy,load_coherence_l1,load_level_population_0,load_level_population_1,load_level_population_2,trace_ok,hermiticity_ok,positivity_ok,energy_balance_ok,top_level_ok,all_physical_checks,energy_match_relative_difference,ergotropy_difference_A_minus_B,ergotropy_ratio_A_over_B,load_ergotropy_over_initial_energy")?;
    for (condition, p, result, sample, diagonal_work) in [
        ("A_no_noise", P_A, a, a_sample, a_diag),
        ("B_dephasing", p_b, b, b_sample, b_diag),
    ] {
        let initial_energy = p * params.omega_chain;
        writeln!(
            writer,
            "{condition},{p:.16e},{initial_energy:.16e},{T_EVAL:.16e},{:.16e},{:.16e},{diagonal_work:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{},{},{},{},{},{},{energy_relative_difference:.16e},{ergotropy_difference:.16e},{ergotropy_ratio:.16e},{:.16e}",
            result.final_load_energy,
            result.final_load_ergotropy,
            result.final_load_ergotropy - diagonal_work,
            sample.load_off_diagonal_l1,
            sample.load_level_populations[0],
            sample.load_level_populations[1],
            sample.load_level_populations[2],
            result.checks.trace_ok,
            result.checks.hermiticity_ok,
            result.checks.positivity_ok,
            result.checks.energy_balance_ok,
            result.checks.top_level_ok,
            result.checks.all_pass(),
            result.final_load_ergotropy / initial_energy,
        )?;
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams {
        hopping_j: 1.0,
        coupling_g: 0.25,
        ..ModelParams::default()
    };
    let a = run_coherent_input_protocol_with_population(&params, P_A, GAMMA_A, T_EVAL, TIME_STEPS)?;
    let target = a.final_load_energy;
    println!("A target load energy = {target:.12e}");

    let mut grid = Vec::with_capacity(GRID_STEPS + 1);
    for index in 0..=GRID_STEPS {
        let p = index as f64 / GRID_STEPS as f64;
        let result =
            run_coherent_input_protocol_with_population(&params, p, GAMMA_B, T_EVAL, TIME_STEPS)?;
        println!(
            "grid p_B={p:.1}: load energy={:.12e}",
            result.final_load_energy
        );
        grid.push((p, result));
    }
    write_grid_csv("coherent_same_time_match_grid.csv", target, &grid)?;

    let mut bracket = None;
    for window in grid.windows(2) {
        let f0 = window[0].1.final_load_energy - target;
        let f1 = window[1].1.final_load_energy - target;
        if f0 == 0.0 || f0.signum() != f1.signum() {
            bracket = Some((window[0].0, window[1].0));
            break;
        }
    }
    let Some((mut low, mut high)) = bracket else {
        println!("root status = no root on p_B in [0, 1]");
        return Ok(());
    };
    println!("bracket = [{low:.12e}, {high:.12e}]");

    let mut best: Option<(f64, ProtocolResult)> = None;
    for iteration in 1..=MAX_BISECTION_ITERATIONS {
        let mid = 0.5 * (low + high);
        let result =
            run_coherent_input_protocol_with_population(&params, mid, GAMMA_B, T_EVAL, TIME_STEPS)?;
        let difference = relative_difference(result.final_load_energy, target);
        println!("bisection {iteration}: p_B={mid:.12e}, relative diff={difference:.12e}");
        if best.as_ref().is_none_or(|(_, current)| {
            difference < relative_difference(current.final_load_energy, target)
        }) {
            best = Some((mid, result.clone()));
        }
        if difference < RELATIVE_TOLERANCE {
            best = Some((mid, result));
            break;
        }
        if result.final_load_energy < target {
            low = mid;
        } else {
            high = mid;
        }
    }

    let (p_b, b) = best.expect("bisection evaluates at least one point");
    let match_error = relative_difference(b.final_load_energy, target);
    if match_error >= RELATIVE_TOLERANCE {
        println!(
            "root status = tolerance not reached; best relative difference={match_error:.12e}"
        );
        return Ok(());
    }
    write_result_csv("coherent_same_time_match.csv", &params, p_b, &a, &b)?;
    println!("matched p_B = {p_b:.12e}");
    println!("relative energy difference = {match_error:.12e}");
    println!(
        "A/B ergotropy = {:.12e} / {:.12e}",
        a.final_load_ergotropy, b.final_load_ergotropy
    );
    println!(
        "A/B physical checks = {}/{}",
        a.checks.all_pass(),
        b.checks.all_pass()
    );
    Ok(())
}
