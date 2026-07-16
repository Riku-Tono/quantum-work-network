use quantum_work_network::experiment::{
    summarize, write_csv, write_time_series_csv, write_time_summary_csv,
};
use quantum_work_network::matching::{match_noisy_input, MatchingConfig};
use quantum_work_network::operators::ModelParams;
use quantum_work_network::protocol::ProtocolConfig;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams {
        hopping_j: 1.0,
        coupling_g: 0.25,
        ..ModelParams::default()
    };
    let a = ProtocolConfig {
        input_strength: 0.2,
        dephasing_strength: 0.0,
        end_time: 10.0,
        time_steps: 100,
    };
    let pair = match_noisy_input(&params, a, 0.5, MatchingConfig::default())?;
    write_csv("first_experiment.csv", &pair)?;
    write_time_series_csv("first_experiment_timeseries.csv", &pair)?;
    write_time_summary_csv("first_experiment_time_summary.csv", &pair)?;
    let summary = summarize(&pair);

    println!(
        "matched B input strength = {:.12e}",
        pair.matched.config.input_strength
    );
    println!(
        "A load energy            = {:.12e}",
        pair.reference.final_load_energy
    );
    println!(
        "B load energy            = {:.12e}",
        pair.matched.final_load_energy
    );
    println!(
        "relative energy diff      = {:.12e}",
        summary.energy_relative_difference
    );
    println!(
        "A load ergotropy         = {:.12e}",
        pair.reference.final_load_ergotropy
    );
    println!(
        "B load ergotropy         = {:.12e}",
        pair.matched.final_load_ergotropy
    );
    println!(
        "relative ergotropy diff   = {:.12e}",
        summary.ergotropy_relative_difference
    );
    println!(
        "A source energy net      = {:.12e}",
        pair.reference.source_energy_net
    );
    println!(
        "B source energy net      = {:.12e}",
        pair.matched.source_energy_net
    );
    println!(
        "B dephasing energy net   = {:.12e}",
        pair.matched.dephasing_energy_net
    );
    println!(
        "A max load ergotropy     = {:.12e} at t={:.12e}",
        pair.reference.maximum_load_ergotropy, pair.reference.maximum_load_ergotropy_time
    );
    println!(
        "B max load ergotropy     = {:.12e} at t={:.12e}",
        pair.matched.maximum_load_ergotropy, pair.matched.maximum_load_ergotropy_time
    );
    println!(
        "A max top population     = {:.12e} at t={:.12e}",
        pair.reference.top_level_population, pair.reference.top_level_population_time
    );
    println!(
        "B max top population     = {:.12e} at t={:.12e}",
        pair.matched.top_level_population, pair.matched.top_level_population_time
    );
    println!(
        "A final top population   = {:.12e}",
        pair.reference.final_top_level_population
    );
    println!(
        "B final top population   = {:.12e}",
        pair.matched.final_top_level_population
    );
    println!(
        "physical checks A/B      = {}/{}",
        pair.reference.checks.all_pass(),
        pair.matched.checks.all_pass()
    );
    println!("physical checks A        = {:?}", pair.reference.checks);
    println!("physical checks B        = {:?}", pair.matched.checks);
    for (label, result) in [("A", &pair.reference), ("B", &pair.matched)] {
        let max_trace_error = result
            .diagnostics
            .iter()
            .map(|d| d.trace_error)
            .fold(0.0_f64, f64::max);
        let max_hermiticity_error = result
            .diagnostics
            .iter()
            .map(|d| d.hermiticity_error)
            .fold(0.0_f64, f64::max);
        let minimum_eigenvalue = result
            .diagnostics
            .iter()
            .map(|d| d.minimum_eigenvalue)
            .fold(f64::INFINITY, f64::min);
        println!("{label} max trace error       = {max_trace_error:.12e}");
        println!("{label} max Hermiticity error = {max_hermiticity_error:.12e}");
        println!("{label} minimum eigenvalue    = {minimum_eigenvalue:.12e}");
    }
    println!("success                  = {}", summary.success);
    Ok(())
}
