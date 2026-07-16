use quantum_work_network::experiment::{write_named_time_series_csv, write_named_time_summary_csv};
use quantum_work_network::operators::ModelParams;
use quantum_work_network::protocol::run_coherent_input_protocol;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let params = ModelParams {
        hopping_j: 1.0,
        coupling_g: 0.25,
        ..ModelParams::default()
    };
    let a = run_coherent_input_protocol(&params, 0.0, 10.0, 100)?;
    let b = run_coherent_input_protocol(&params, 0.5, 10.0, 100)?;
    let series = [("A_no_noise", &a), ("B_dephasing", &b)];
    write_named_time_series_csv("coherent_input_timeseries.csv", &series)?;
    write_named_time_summary_csv("coherent_input_summary.csv", &series)?;

    for (condition, result) in series {
        println!("{condition}");
        println!(
            "  maximum load ergotropy    = {:.12e} at t={:.12e}",
            result.maximum_load_ergotropy, result.maximum_load_ergotropy_time
        );
        println!(
            "  maximum off-diagonal abs  = {:.12e} at t={:.12e}",
            result.maximum_load_off_diagonal, result.maximum_load_off_diagonal_time
        );
        println!("  physical checks           = {:?}", result.checks);
    }
    Ok(())
}
