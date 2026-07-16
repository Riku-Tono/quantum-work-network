//! One complete charging protocol run and its physical checks.

use crate::diagnostics::{
    diagnose_states, integrate_diagnostic_powers, Diagnostics, IntegratedPowers,
};
use crate::error::PhysicsError;
use crate::liouvillian::build_liouvillian;
use crate::matrix::{ComplexMatrix, C64};
use crate::operators::{build_operators, ModelParams};
use crate::partial_trace::partial_trace;
use crate::propagator::{DenseExponentialPropagator, QuantumState};

#[derive(Debug, Clone, Copy)]
pub struct ProtocolConfig {
    /// Lindblad injection rate multiplying D[sigma_1_plus].
    pub input_strength: f64,
    /// Pure-dephasing rate on each of the three chain sites.
    pub dephasing_strength: f64,
    pub end_time: f64,
    /// Number of uniform intervals used for signed-power integration.
    pub time_steps: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct PhysicalChecks {
    pub trace_ok: bool,
    pub hermiticity_ok: bool,
    pub positivity_ok: bool,
    pub energy_balance_ok: bool,
    pub top_level_ok: bool,
}

impl PhysicalChecks {
    pub fn all_pass(&self) -> bool {
        self.trace_ok
            && self.hermiticity_ok
            && self.positivity_ok
            && self.energy_balance_ok
            && self.top_level_ok
    }
}

#[derive(Debug, Clone)]
pub struct LoadTimePoint {
    pub time: f64,
    pub load_energy: f64,
    pub load_ergotropy: f64,
    pub load_level_populations: Vec<f64>,
    /// Independent upper-triangular entries, ordered (0,1), (0,2), ... .
    pub load_off_diagonal_values: Vec<C64>,
    /// Magnitudes of the independent upper-triangular entries, ordered
    /// (0,1), (0,2), ..., (load_dim - 2, load_dim - 1).
    pub load_off_diagonal_magnitudes: Vec<f64>,
    /// Sum of magnitudes over every off-diagonal entry (both triangles).
    pub load_off_diagonal_l1: f64,
}

#[derive(Debug, Clone)]
pub struct ProtocolResult {
    pub config: ProtocolConfig,
    pub final_load_energy: f64,
    pub final_load_ergotropy: f64,
    pub final_load_passive_energy: f64,
    pub source_energy_net: f64,
    pub source_energy_in: f64,
    pub source_energy_out: f64,
    pub dephasing_energy_net: f64,
    pub dephasing_energy_in: f64,
    pub dephasing_energy_out: f64,
    /// Maximum top-level population over the complete sampled time grid.
    pub top_level_population: f64,
    pub top_level_population_time: f64,
    pub final_top_level_population: f64,
    pub maximum_load_ergotropy: f64,
    pub maximum_load_ergotropy_time: f64,
    pub maximum_load_off_diagonal: f64,
    pub maximum_load_off_diagonal_time: f64,
    pub total_energy_change: f64,
    pub energy_balance_residual: f64,
    pub checks: PhysicalChecks,
    pub diagnostics: Vec<Diagnostics>,
    pub load_time_series: Vec<LoadTimePoint>,
}

fn validate_config(config: &ProtocolConfig) -> Result<(), PhysicsError> {
    for (name, value) in [
        ("input_strength", config.input_strength),
        ("dephasing_strength", config.dephasing_strength),
        ("end_time", config.end_time),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(PhysicsError::InvalidParameter(format!(
                "{name} must be finite and nonnegative, got {value}"
            )));
        }
    }
    if config.time_steps == 0 {
        return Err(PhysicsError::InvalidParameter(
            "time_steps must be at least one".to_string(),
        ));
    }
    Ok(())
}

fn vacuum_state(dim: usize) -> ComplexMatrix {
    let mut rho = ComplexMatrix::zeros(dim, dim);
    rho[(0, 0)] = C64::new(1.0, 0.0);
    rho
}

fn coherent_q1_state(
    load_dim: usize,
    excited_population: f64,
) -> Result<ComplexMatrix, PhysicsError> {
    if !excited_population.is_finite() || !(0.0..=1.0).contains(&excited_population) {
        return Err(PhysicsError::InvalidParameter(format!(
            "coherent q1 excited population must be within [0, 1], got {excited_population}"
        )));
    }
    let dim = 2 * 2 * 2 * load_dim;
    let q1_excited_index = 2 * 2 * load_dim;
    let mut rho = ComplexMatrix::zeros(dim, dim);
    let ground_population = 1.0 - excited_population;
    let coherence = (ground_population * excited_population).sqrt();
    rho[(0, 0)] = C64::new(ground_population, 0.0);
    rho[(q1_excited_index, q1_excited_index)] = C64::new(excited_population, 0.0);
    rho[(0, q1_excited_index)] = C64::new(coherence, 0.0);
    rho[(q1_excited_index, 0)] = C64::new(coherence, 0.0);
    Ok(rho)
}

fn make_load_time_series(
    states: &[QuantumState],
    diagnostics: &[Diagnostics],
    dims: &[usize],
    load_dim: usize,
) -> Result<Vec<LoadTimePoint>, PhysicsError> {
    if states.len() != diagnostics.len() {
        return Err(PhysicsError::DimensionMismatch(format!(
            "{} propagated states but {} diagnostic samples",
            states.len(),
            diagnostics.len()
        )));
    }
    states
        .iter()
        .zip(diagnostics)
        .map(|(state, diagnostic)| {
            let rho_load = partial_trace(&state.rho, dims, &[3])?;
            let load_level_populations = (0..load_dim)
                .map(|level| rho_load[(level, level)].re)
                .collect();
            let mut load_off_diagonal_values = Vec::new();
            for row in 0..load_dim {
                for col in (row + 1)..load_dim {
                    load_off_diagonal_values.push(rho_load[(row, col)]);
                }
            }
            let load_off_diagonal_magnitudes = load_off_diagonal_values
                .iter()
                .map(|value| value.norm())
                .collect::<Vec<_>>();
            let load_off_diagonal_l1 = 2.0 * load_off_diagonal_magnitudes.iter().sum::<f64>();
            Ok(LoadTimePoint {
                time: state.time,
                load_energy: diagnostic.load_energy,
                load_ergotropy: diagnostic.load_ergotropy,
                load_level_populations,
                load_off_diagonal_values,
                load_off_diagonal_magnitudes,
                load_off_diagonal_l1,
            })
        })
        .collect()
}

fn maximum_value_and_time(
    samples: &[LoadTimePoint],
    value: impl Fn(&LoadTimePoint) -> f64,
) -> (f64, f64) {
    let first = samples
        .first()
        .expect("uniform propagation returns at least the initial state");
    let mut maximum = value(first);
    let mut time = first.time;
    for sample in &samples[1..] {
        let candidate = value(sample);
        if candidate > maximum {
            maximum = candidate;
            time = sample.time;
        }
    }
    (maximum, time)
}

fn run_protocol_with_initial_state(
    params: &ModelParams,
    config: ProtocolConfig,
    rho0: ComplexMatrix,
    include_lindblad_source: bool,
) -> Result<ProtocolResult, PhysicsError> {
    validate_config(&config)?;
    let operators = build_operators(params)?;
    let dim: usize = operators.dims.iter().product();
    let dephasing: Vec<_> = operators
        .sigma_z_sites
        .iter()
        .map(|sigma_z| sigma_z * C64::new((config.dephasing_strength / 2.0).sqrt(), 0.0))
        .collect();
    let source = include_lindblad_source
        .then(|| &operators.sigma_1_plus * C64::new(config.input_strength.sqrt(), 0.0));
    let source_count = if source.is_some() { 1 } else { 0 };
    let mut all_collapses = Vec::with_capacity(dephasing.len() + source_count);
    if let Some(source_collapse) = &source {
        all_collapses.push(source_collapse.clone());
    }
    all_collapses.extend(dephasing.iter().cloned());
    let liouvillian = build_liouvillian(&operators.h_total, &all_collapses)?;
    let propagator = DenseExponentialPropagator::new(liouvillian, dim)?;
    let states = propagator.propagate_uniform(&rho0, config.end_time, config.time_steps)?;
    let diagnostics = diagnose_states(&states, &operators, params, source.as_ref(), &dephasing)?;
    let load_time_series =
        make_load_time_series(&states, &diagnostics, &operators.dims, params.load_dim)?;
    let powers = integrate_diagnostic_powers(&diagnostics)?;
    let first = diagnostics
        .first()
        .expect("uniform propagation returns initial state");
    let last = diagnostics
        .last()
        .expect("uniform propagation returns final state");
    let final_top_level_population = *load_time_series
        .last()
        .expect("uniform propagation returns final state")
        .load_level_populations
        .last()
        .expect("load dimension is at least two");
    let (top_level_population, top_level_population_time) =
        maximum_value_and_time(&load_time_series, |sample| {
            *sample
                .load_level_populations
                .last()
                .expect("load dimension is at least two")
        });
    let (maximum_load_ergotropy, maximum_load_ergotropy_time) =
        maximum_value_and_time(&load_time_series, |sample| sample.load_ergotropy);
    let (maximum_load_off_diagonal, maximum_load_off_diagonal_time) =
        maximum_value_and_time(&load_time_series, |sample| {
            sample
                .load_off_diagonal_magnitudes
                .iter()
                .copied()
                .fold(0.0_f64, f64::max)
        });
    let total_energy_change = last.total_energy - first.total_energy;
    let energy_balance_residual =
        total_energy_change - powers.source_energy_net - powers.dephasing_energy_net;
    let checks = physical_checks(
        &diagnostics,
        powers,
        total_energy_change,
        top_level_population,
    );

    Ok(ProtocolResult {
        config,
        final_load_energy: last.load_energy,
        final_load_ergotropy: last.load_ergotropy,
        final_load_passive_energy: last.load_passive_energy,
        source_energy_net: powers.source_energy_net,
        source_energy_in: powers.source_energy_in,
        source_energy_out: powers.source_energy_out,
        dephasing_energy_net: powers.dephasing_energy_net,
        dephasing_energy_in: powers.dephasing_energy_in,
        dephasing_energy_out: powers.dephasing_energy_out,
        top_level_population,
        top_level_population_time,
        final_top_level_population,
        maximum_load_ergotropy,
        maximum_load_ergotropy_time,
        maximum_load_off_diagonal,
        maximum_load_off_diagonal_time,
        total_energy_change,
        energy_balance_residual,
        checks,
        diagnostics,
        load_time_series,
    })
}

fn physical_checks(
    diagnostics: &[Diagnostics],
    powers: IntegratedPowers,
    energy_change: f64,
    top_population: f64,
) -> PhysicalChecks {
    let trace_ok = diagnostics.iter().all(|d| d.trace_error <= 1.0e-8);
    let hermiticity_ok = diagnostics.iter().all(|d| d.hermiticity_error <= 1.0e-8);
    let positivity_ok = diagnostics.iter().all(|d| d.minimum_eigenvalue >= -1.0e-8);
    let scale = energy_change
        .abs()
        .max(powers.source_energy_net.abs())
        .max(powers.dephasing_energy_net.abs())
        .max(1.0);
    let residual = energy_change - powers.source_energy_net - powers.dephasing_energy_net;
    let energy_balance_ok = residual.abs() <= 5.0e-4 * scale;
    let top_level_ok = top_population < 0.05;
    PhysicalChecks {
        trace_ok,
        hermiticity_ok,
        positivity_ok,
        energy_balance_ok,
        top_level_ok,
    }
}

pub fn run_protocol(
    params: &ModelParams,
    config: ProtocolConfig,
) -> Result<ProtocolResult, PhysicsError> {
    let dim = 2 * 2 * 2 * params.load_dim;
    let rho0 = vacuum_state(dim);
    run_protocol_with_initial_state(params, config, rho0, true)
}

/// Run the coherent-input sanity check without a Lindblad source.
pub fn run_coherent_input_protocol(
    params: &ModelParams,
    dephasing_strength: f64,
    end_time: f64,
    time_steps: usize,
) -> Result<ProtocolResult, PhysicsError> {
    run_coherent_input_protocol_with_population(
        params,
        0.5,
        dephasing_strength,
        end_time,
        time_steps,
    )
}

/// Run a source-free coherent-input protocol with
/// q1 = sqrt(1-p)|0> + sqrt(p)|1>.
pub fn run_coherent_input_protocol_with_population(
    params: &ModelParams,
    excited_population: f64,
    dephasing_strength: f64,
    end_time: f64,
    time_steps: usize,
) -> Result<ProtocolResult, PhysicsError> {
    let config = ProtocolConfig {
        input_strength: 0.0,
        dephasing_strength,
        end_time,
        time_steps,
    };
    let rho0 = coherent_q1_state(params.load_dim, excited_population)?;
    run_protocol_with_initial_state(params, config, rho0, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(time: f64, ergotropy: f64, top: f64) -> LoadTimePoint {
        LoadTimePoint {
            time,
            load_energy: 0.0,
            load_ergotropy: ergotropy,
            load_level_populations: vec![1.0 - top, 0.0, top],
            load_off_diagonal_values: vec![C64::new(0.0, 0.0); 3],
            load_off_diagonal_magnitudes: vec![0.0; 3],
            load_off_diagonal_l1: 0.0,
        }
    }

    #[test]
    fn extrema_use_all_times_and_keep_the_first_maximum() {
        let samples = vec![
            point(0.0, 0.0, 0.0),
            point(0.1, 0.2, 0.08),
            point(0.2, 0.2, 0.03),
        ];
        assert_eq!(
            maximum_value_and_time(&samples, |s| s.load_ergotropy),
            (0.2, 0.1)
        );
        assert_eq!(
            maximum_value_and_time(&samples, |s| *s.load_level_populations.last().unwrap()),
            (0.08, 0.1)
        );
    }

    #[test]
    fn coherent_q1_state_has_the_requested_superposition() {
        let rho = coherent_q1_state(3, 0.5).unwrap();
        assert_eq!(rho.shape(), (24, 24));
        for &(row, col) in &[(0, 0), (0, 12), (12, 0), (12, 12)] {
            assert_eq!(rho[(row, col)], C64::new(0.5, 0.0));
        }
        assert_eq!(rho.trace(), C64::new(1.0, 0.0));
        assert_eq!(rho.iter().filter(|value| value.norm() > 0.0).count(), 4);
    }

    #[test]
    fn coherent_q1_state_uses_the_requested_population() {
        let rho = coherent_q1_state(3, 0.2).unwrap();
        assert_eq!(rho[(0, 0)], C64::new(0.8, 0.0));
        assert_eq!(rho[(12, 12)], C64::new(0.2, 0.0));
        assert!((rho[(0, 12)].re - 0.4).abs() < 1.0e-15);
        assert!(coherent_q1_state(3, -0.1).is_err());
        assert!(coherent_q1_state(3, 1.1).is_err());
    }
}
