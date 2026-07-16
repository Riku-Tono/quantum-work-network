//! Match the noisy protocol's final load energy to a reference protocol.

use crate::error::PhysicsError;
use crate::operators::ModelParams;
use crate::protocol::{run_protocol, ProtocolConfig, ProtocolResult};

#[derive(Debug, Clone, Copy)]
pub struct MatchingConfig {
    pub relative_tolerance: f64,
    pub max_iterations: usize,
    pub initial_upper_input: f64,
    pub max_input_strength: f64,
}

impl Default for MatchingConfig {
    fn default() -> Self {
        Self {
            relative_tolerance: 1.0e-4,
            max_iterations: 40,
            initial_upper_input: 0.2,
            max_input_strength: 100.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct MatchedPair {
    pub reference: ProtocolResult,
    pub matched: ProtocolResult,
    pub relative_energy_difference: f64,
    pub iterations: usize,
}

fn relative_difference(actual: f64, target: f64) -> f64 {
    (actual - target).abs() / target.abs().max(1.0e-14)
}

pub fn match_noisy_input(
    params: &ModelParams,
    reference_config: ProtocolConfig,
    noisy_dephasing_strength: f64,
    matching: MatchingConfig,
) -> Result<MatchedPair, PhysicsError> {
    if matching.relative_tolerance <= 0.0 || !matching.relative_tolerance.is_finite() {
        return Err(PhysicsError::InvalidParameter(
            "relative_tolerance must be finite and positive".to_string(),
        ));
    }
    let reference = run_protocol(params, reference_config)?;
    let target = reference.final_load_energy;
    if target <= 1.0e-12 {
        return Err(PhysicsError::MatchingFailure(format!(
            "reference load energy is too small to match reliably: {target:e}"
        )));
    }

    let make_config = |input_strength| ProtocolConfig {
        input_strength,
        dephasing_strength: noisy_dephasing_strength,
        end_time: reference_config.end_time,
        time_steps: reference_config.time_steps,
    };

    let mut low = 0.0;
    let mut low_result = run_protocol(params, make_config(low))?;
    let mut high = matching
        .initial_upper_input
        .max(reference_config.input_strength)
        .max(1.0e-12);
    let mut high_result = run_protocol(params, make_config(high))?;
    let mut iterations = 2;

    while high_result.final_load_energy < target && high < matching.max_input_strength {
        low = high;
        low_result = high_result;
        high = (2.0 * high).min(matching.max_input_strength);
        high_result = run_protocol(params, make_config(high))?;
        iterations += 1;
    }
    if high_result.final_load_energy < target {
        return Err(PhysicsError::MatchingFailure(format!(
            "could not bracket target energy {target:e} below input strength {}",
            matching.max_input_strength
        )));
    }

    let mut best = if relative_difference(low_result.final_load_energy, target)
        < relative_difference(high_result.final_load_energy, target)
    {
        low_result
    } else {
        high_result
    };

    for _ in 0..matching.max_iterations {
        let mid = 0.5 * (low + high);
        let result = run_protocol(params, make_config(mid))?;
        iterations += 1;
        if relative_difference(result.final_load_energy, target)
            < relative_difference(best.final_load_energy, target)
        {
            best = result.clone();
        }
        if relative_difference(result.final_load_energy, target) < matching.relative_tolerance {
            return Ok(MatchedPair {
                reference,
                relative_energy_difference: relative_difference(result.final_load_energy, target),
                matched: result,
                iterations,
            });
        }
        if result.final_load_energy < target {
            low = mid;
        } else {
            high = mid;
        }
    }

    let difference = relative_difference(best.final_load_energy, target);
    if difference >= matching.relative_tolerance {
        return Err(PhysicsError::MatchingFailure(format!(
            "best relative energy difference {difference:e} exceeds tolerance {:e}",
            matching.relative_tolerance
        )));
    }
    Ok(MatchedPair {
        reference,
        matched: best,
        relative_energy_difference: difference,
        iterations,
    })
}
