//! CSV output and success criteria for one matched A/B experiment.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::error::PhysicsError;
use crate::matching::MatchedPair;
use crate::protocol::ProtocolResult;

#[derive(Debug, Clone, Copy)]
pub struct ExperimentSummary {
    pub energy_relative_difference: f64,
    pub ergotropy_relative_difference: f64,
    pub ergotropy_absolute_difference: f64,
    pub usable_amount_nonzero: bool,
    pub success: bool,
}

pub fn summarize(pair: &MatchedPair) -> ExperimentSummary {
    let a = &pair.reference;
    let b = &pair.matched;
    let ergotropy_absolute_difference = b.final_load_ergotropy - a.final_load_ergotropy;
    let ergotropy_relative_difference = ergotropy_absolute_difference.abs()
        / a.final_load_ergotropy
            .abs()
            .max(b.final_load_ergotropy.abs())
            .max(1.0e-14);
    let usable_amount_nonzero = a.final_load_ergotropy.max(b.final_load_ergotropy) > 1.0e-6;
    let success = pair.relative_energy_difference < 1.0e-4
        && ergotropy_relative_difference >= 0.05
        && usable_amount_nonzero
        && a.top_level_population < 0.05
        && b.top_level_population < 0.05
        && a.checks.all_pass()
        && b.checks.all_pass();
    ExperimentSummary {
        energy_relative_difference: pair.relative_energy_difference,
        ergotropy_relative_difference,
        ergotropy_absolute_difference,
        usable_amount_nonzero,
        success,
    }
}

pub fn write_csv(path: impl AsRef<Path>, pair: &MatchedPair) -> Result<(), PhysicsError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "quantity,A_no_noise,B_dephasing,difference_B_minus_A"
    )?;
    let rows = [
        (
            "input_strength",
            pair.reference.config.input_strength,
            pair.matched.config.input_strength,
        ),
        (
            "load_energy",
            pair.reference.final_load_energy,
            pair.matched.final_load_energy,
        ),
        (
            "load_ergotropy",
            pair.reference.final_load_ergotropy,
            pair.matched.final_load_ergotropy,
        ),
        (
            "load_passive_energy",
            pair.reference.final_load_passive_energy,
            pair.matched.final_load_passive_energy,
        ),
        (
            "source_energy_net",
            pair.reference.source_energy_net,
            pair.matched.source_energy_net,
        ),
        (
            "source_energy_in",
            pair.reference.source_energy_in,
            pair.matched.source_energy_in,
        ),
        (
            "source_energy_out",
            pair.reference.source_energy_out,
            pair.matched.source_energy_out,
        ),
        (
            "dephasing_energy_net",
            pair.reference.dephasing_energy_net,
            pair.matched.dephasing_energy_net,
        ),
        (
            "dephasing_energy_in",
            pair.reference.dephasing_energy_in,
            pair.matched.dephasing_energy_in,
        ),
        (
            "dephasing_energy_out",
            pair.reference.dephasing_energy_out,
            pair.matched.dephasing_energy_out,
        ),
        (
            "maximum_top_level_population",
            pair.reference.top_level_population,
            pair.matched.top_level_population,
        ),
        (
            "final_top_level_population",
            pair.reference.final_top_level_population,
            pair.matched.final_top_level_population,
        ),
        (
            "energy_balance_residual",
            pair.reference.energy_balance_residual,
            pair.matched.energy_balance_residual,
        ),
    ];
    for (name, a, b) in rows {
        writeln!(writer, "{name},{a:.16e},{b:.16e},{:.16e}", b - a)?;
    }
    let summary = summarize(pair);
    writeln!(
        writer,
        "relative_load_energy_difference,,,{:.16e}",
        summary.energy_relative_difference
    )?;
    writeln!(
        writer,
        "relative_ergotropy_difference,,,{:.16e}",
        summary.ergotropy_relative_difference
    )?;
    writeln!(
        writer,
        "all_physical_checks_A,,,{}",
        pair.reference.checks.all_pass()
    )?;
    writeln!(
        writer,
        "all_physical_checks_B,,,{}",
        pair.matched.checks.all_pass()
    )?;
    writeln!(writer, "success,,,{}", summary.success)?;
    Ok(())
}

fn write_protocol_time_series(
    writer: &mut impl Write,
    condition: &str,
    result: &ProtocolResult,
) -> Result<(), PhysicsError> {
    for sample in &result.load_time_series {
        write!(
            writer,
            "{condition},{:.16e},{:.16e},{:.16e}",
            sample.time, sample.load_energy, sample.load_ergotropy
        )?;
        for population in &sample.load_level_populations {
            write!(writer, ",{population:.16e}")?;
        }
        for (value, magnitude) in sample
            .load_off_diagonal_values
            .iter()
            .zip(&sample.load_off_diagonal_magnitudes)
        {
            write!(
                writer,
                ",{:.16e},{:.16e},{magnitude:.16e}",
                value.re, value.im
            )?;
        }
        writeln!(writer, ",{:.16e}", sample.load_off_diagonal_l1)?;
    }
    Ok(())
}

pub fn write_time_series_csv(
    path: impl AsRef<Path>,
    pair: &MatchedPair,
) -> Result<(), PhysicsError> {
    write_named_time_series_csv(
        path,
        &[
            ("A_no_noise", &pair.reference),
            ("B_dephasing", &pair.matched),
        ],
    )
}

pub fn write_named_time_series_csv(
    path: impl AsRef<Path>,
    series: &[(&str, &ProtocolResult)],
) -> Result<(), PhysicsError> {
    if series.is_empty() {
        return Err(PhysicsError::InvalidParameter(
            "time-series CSV requires at least one protocol result".to_string(),
        ));
    }
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let load_dim = series[0].1.load_time_series[0].load_level_populations.len();
    write!(writer, "condition,time,load_energy,load_ergotropy")?;
    for level in 0..load_dim {
        write!(writer, ",load_level_population_{level}")?;
    }
    for row in 0..load_dim {
        for col in (row + 1)..load_dim {
            write!(
                writer,
                ",load_rho_re_{row}_{col},load_rho_im_{row}_{col},load_rho_abs_{row}_{col}"
            )?;
        }
    }
    writeln!(writer, ",load_off_diagonal_l1")?;
    for &(condition, result) in series {
        write_protocol_time_series(&mut writer, condition, result)?;
    }
    Ok(())
}

pub fn write_time_summary_csv(
    path: impl AsRef<Path>,
    pair: &MatchedPair,
) -> Result<(), PhysicsError> {
    write_named_time_summary_csv(
        path,
        &[
            ("A_no_noise", &pair.reference),
            ("B_dephasing", &pair.matched),
        ],
    )
}

pub fn write_named_time_summary_csv(
    path: impl AsRef<Path>,
    series: &[(&str, &ProtocolResult)],
) -> Result<(), PhysicsError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "condition,maximum_load_ergotropy,maximum_load_ergotropy_time,maximum_load_off_diagonal,maximum_load_off_diagonal_time,maximum_top_level_population,maximum_top_level_population_time,final_top_level_population")?;
    for &(condition, result) in series {
        writeln!(
            writer,
            "{condition},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e},{:.16e}",
            result.maximum_load_ergotropy,
            result.maximum_load_ergotropy_time,
            result.maximum_load_off_diagonal,
            result.maximum_load_off_diagonal_time,
            result.top_level_population,
            result.top_level_population_time,
            result.final_top_level_population,
        )?;
    }
    Ok(())
}
