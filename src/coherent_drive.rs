//! Fixed-parameter coherent-drive sanity check for the three-site network.

use nalgebra::linalg::Schur;

use crate::diagnostics::{integrate_signed_power, lindblad_action, SignedEnergyIntegral};
use crate::ergotropy::ergotropy;
use crate::error::PhysicsError;
use crate::matrix::{commutator, expectation, hermiticity_error, ComplexMatrix, C64};
use crate::operators::{build_operators_for_chain, ModelParams, Operators};
use crate::partial_trace::partial_trace;
use crate::propagator::QuantumState;
use crate::time_dependent::{SaveSchedule, TimeDependentRk4};

pub const TRACE_TOLERANCE: f64 = 1.0e-8;
pub const HERMITICITY_TOLERANCE: f64 = 1.0e-8;
pub const POSITIVITY_TOLERANCE: f64 = 1.0e-8;
pub const LEDGER_ABSOLUTE_TOLERANCE: f64 = 5.0e-5;
pub const LEDGER_RELATIVE_TOLERANCE: f64 = 5.0e-4;
pub const SIGNAL_TOLERANCE: f64 = 1.0e-8;
pub const CONVERGENCE_RELATIVE_TOLERANCE: f64 = 5.0e-3;
pub const CONVERGENCE_ABSOLUTE_TOLERANCE: f64 = 1.0e-7;
pub const TOP_LEVEL_LIMIT: f64 = 0.05;

#[derive(Debug, Clone, Copy)]
pub struct CoherentDriveConfig {
    pub omega0: f64,
    pub omega_drive: f64,
    pub tau: f64,
    pub t_end: f64,
    pub dt: f64,
    pub save_interval: f64,
    pub gamma_phi: f64,
}

impl CoherentDriveConfig {
    pub fn milestone_5b(gamma_phi: f64, dt: f64) -> Self {
        Self {
            omega0: 0.2,
            omega_drive: 1.0,
            tau: 3.2,
            t_end: 10.0,
            dt,
            save_interval: 0.01,
            gamma_phi,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CoherentDriveSample {
    pub time: f64,
    pub drive_envelope: f64,
    pub drive_amplitude: f64,
    pub load_energy: f64,
    pub load_ergotropy: f64,
    pub load_diagonal_ergotropy: f64,
    pub load_coherence_ergotropy: f64,
    pub load_coherence_l1: f64,
    pub load_populations: [f64; 3],
    pub bare_network_energy: f64,
    pub drive_power: f64,
    pub drive_power_imaginary: f64,
    pub dephasing_power: f64,
    pub dephasing_power_imaginary: f64,
    pub trace_error: f64,
    pub hermiticity_error: f64,
    pub minimum_eigenvalue: f64,
    pub all_finite: bool,
}

#[derive(Debug, Clone)]
pub struct CoherentDriveSummary {
    pub maximum_load_energy: (f64, f64),
    pub maximum_load_ergotropy: (f64, f64),
    pub maximum_load_coherence_ergotropy: (f64, f64),
    pub maximum_load_coherence_l1: (f64, f64),
    pub maximum_top_level_population: (f64, f64),
    pub at_tau: CoherentDriveSample,
    pub at_end: CoherentDriveSample,
    pub drive_energy: SignedEnergyIntegral,
    pub dephasing_energy: SignedEnergyIntegral,
    pub delta_bare_network_energy: f64,
    pub ledger_residual: f64,
    pub maximum_trace_error: f64,
    pub maximum_hermiticity_error: f64,
    pub worst_minimum_eigenvalue: f64,
    pub maximum_drive_power_imaginary: f64,
    pub maximum_dephasing_power_imaginary: f64,
    pub all_finite: bool,
    pub physical_checks_pass: bool,
    pub ledger_check_pass: bool,
    pub top_level_check_pass: bool,
}

#[derive(Debug, Clone)]
pub struct CoherentDriveRun {
    pub config: CoherentDriveConfig,
    pub samples: Vec<CoherentDriveSample>,
    pub summary: CoherentDriveSummary,
    /// Final full-system density matrix, retained for deterministic downstream diagnostics.
    pub final_state: ComplexMatrix,
    /// Saved full-system states on the same grid as `samples`.
    pub states: Vec<QuantumState>,
}

pub fn drive_envelope(time: f64, tau: f64) -> f64 {
    if time < 0.0 || time > tau {
        0.0
    } else {
        (std::f64::consts::PI * time / tau).sin().powi(2)
    }
}

pub fn drive_hamiltonian(
    time: f64,
    config: &CoherentDriveConfig,
    sigma_plus: &ComplexMatrix,
) -> ComplexMatrix {
    let amplitude = config.omega0 * drive_envelope(time, config.tau);
    let phase_minus = C64::from_polar(1.0, -config.omega_drive * time);
    let phase_plus = phase_minus.conj();
    (sigma_plus * phase_minus + sigma_plus.adjoint() * phase_plus) * C64::new(amplitude, 0.0)
}

pub fn run_coherent_drive(
    params: &ModelParams,
    config: CoherentDriveConfig,
) -> Result<CoherentDriveRun, PhysicsError> {
    run_coherent_drive_with_noise_sites(params, config, &[0, 1, 2])
}

/// Run the coherent-drive model with local dephasing on selected zero-based sites.
///
/// The existing [`run_coherent_drive`] behavior is preserved by selecting all
/// three sites. An empty selection is valid, and zero dephasing creates no
/// collapse operators.
pub fn run_coherent_drive_with_noise_sites(
    params: &ModelParams,
    config: CoherentDriveConfig,
    noise_sites: &[usize],
) -> Result<CoherentDriveRun, PhysicsError> {
    validate_config(&config)?;
    validate_noise_sites(noise_sites)?;
    let mut site_gammas = [0.0; 3];
    for &site in noise_sites {
        site_gammas[site] = config.gamma_phi;
    }
    run_coherent_drive_with_site_gammas(params, config, site_gammas)
}

/// Run the coherent-drive model with an independently specified dephasing rate
/// for each zero-based site. A zero rate creates no collapse operator.
///
/// The scalar `config.gamma_phi` remains part of the historical configuration
/// and is validated for backward compatibility; the dynamics in this entry
/// point use `site_gammas`.
pub fn run_coherent_drive_with_site_gammas(
    params: &ModelParams,
    config: CoherentDriveConfig,
    site_gammas: [f64; 3],
) -> Result<CoherentDriveRun, PhysicsError> {
    run_coherent_drive_general(params, config, 3, &site_gammas)
}

/// Generalized chain-length entry point used by Milestone 8a. It preserves the
/// same drive, load, RK4, diagnostics, and per-site dephasing convention.
pub fn run_coherent_drive_for_chain(
    params: &ModelParams,
    config: CoherentDriveConfig,
    chain_length: usize,
    gamma_phi_per_site: f64,
) -> Result<CoherentDriveRun, PhysicsError> {
    if !gamma_phi_per_site.is_finite() || gamma_phi_per_site < 0.0 {
        return Err(PhysicsError::InvalidParameter(format!(
            "gamma_phi_per_site must be finite and nonnegative, got {gamma_phi_per_site}"
        )));
    }
    let site_gammas = vec![gamma_phi_per_site; chain_length];
    run_coherent_drive_general(params, config, chain_length, &site_gammas)
}

fn run_coherent_drive_general(
    params: &ModelParams,
    config: CoherentDriveConfig,
    chain_length: usize,
    site_gammas: &[f64],
) -> Result<CoherentDriveRun, PhysicsError> {
    validate_config(&config)?;
    validate_site_gamma_slice(site_gammas, chain_length)?;
    let operators = build_operators_for_chain(params, chain_length)?;
    let dim: usize = operators.dims.iter().product();
    let mut rho0 = ComplexMatrix::zeros(dim, dim);
    rho0[(0, 0)] = C64::new(1.0, 0.0);
    let dephasing = dephasing_operators_by_site(&operators, site_gammas);
    let h0 = operators.h_total.clone();
    let sigma_plus = operators.sigma_1_plus.clone();
    let solver = TimeDependentRk4::new(config.dt)?;
    let states = solver.propagate(
        &rho0,
        0.0,
        config.t_end,
        |time| &h0 + drive_hamiltonian(time, &config, &sigma_plus),
        |_| dephasing.clone(),
        SaveSchedule::Interval(config.save_interval),
    )?;

    let load_hamiltonian = local_load_hamiltonian(params);
    let mut samples = Vec::with_capacity(states.len());
    for state in &states {
        let load_index = operators.dims.len() - 1;
        let rho_load = partial_trace(&state.rho, &operators.dims, &[load_index])?;
        let load_result = ergotropy(&rho_load, &load_hamiltonian, 1.0e-9)?;
        let mut diagonal = ComplexMatrix::zeros(params.load_dim, params.load_dim);
        for level in 0..params.load_dim {
            diagonal[(level, level)] = rho_load[(level, level)];
        }
        let diagonal_result = ergotropy(&diagonal, &load_hamiltonian, 1.0e-9)?;
        let coherence_l1 = rho_load
            .iter()
            .enumerate()
            .filter(|(index, _)| index / params.load_dim != index % params.load_dim)
            .map(|(_, value)| value.norm())
            .sum();
        let drive = drive_hamiltonian(state.time, &config, &operators.sigma_1_plus);
        let drive_power_complex =
            expectation(&state.rho, &commutator(&drive, &operators.h_total)) * C64::new(0.0, 1.0);
        let mut dephasing_action = ComplexMatrix::zeros(dim, dim);
        for collapse in &dephasing {
            dephasing_action += lindblad_action(collapse, &state.rho)?;
        }
        let dephasing_power_complex = expectation(&dephasing_action, &operators.h_total);
        let trace_error = (state.rho.trace() - C64::new(1.0, 0.0)).norm();
        let state_hermiticity_error = hermiticity_error(&state.rho);
        let minimum_eigenvalue = minimum_eigenvalue(&state.rho);
        let all_finite = state.rho.iter().all(finite_complex)
            && [
                load_result.energy,
                load_result.ergotropy,
                diagonal_result.ergotropy,
                coherence_l1,
                drive_power_complex.re,
                drive_power_complex.im,
                dephasing_power_complex.re,
                dephasing_power_complex.im,
                trace_error,
                state_hermiticity_error,
                minimum_eigenvalue,
            ]
            .iter()
            .all(|value| value.is_finite());
        samples.push(CoherentDriveSample {
            time: state.time,
            drive_envelope: drive_envelope(state.time, config.tau),
            drive_amplitude: config.omega0 * drive_envelope(state.time, config.tau),
            load_energy: load_result.energy,
            load_ergotropy: load_result.ergotropy,
            load_diagonal_ergotropy: diagonal_result.ergotropy,
            load_coherence_ergotropy: load_result.ergotropy - diagonal_result.ergotropy,
            load_coherence_l1: coherence_l1,
            load_populations: [
                rho_load[(0, 0)].re,
                rho_load[(1, 1)].re,
                rho_load[(2, 2)].re,
            ],
            bare_network_energy: expectation(&state.rho, &operators.h_total).re,
            drive_power: drive_power_complex.re,
            drive_power_imaginary: drive_power_complex.im,
            dephasing_power: dephasing_power_complex.re,
            dephasing_power_imaginary: dephasing_power_complex.im,
            trace_error,
            hermiticity_error: state_hermiticity_error,
            minimum_eigenvalue,
            all_finite,
        });
    }
    let summary = summarize(&samples, &config)?;
    let final_state = states
        .last()
        .expect("time-dependent propagation returns the initial and final states")
        .rho
        .clone();
    Ok(CoherentDriveRun {
        config,
        samples,
        summary,
        final_state,
        states,
    })
}

fn summarize(
    samples: &[CoherentDriveSample],
    config: &CoherentDriveConfig,
) -> Result<CoherentDriveSummary, PhysicsError> {
    let at_tau = sample_at(samples, config.tau)?.clone();
    let at_end = sample_at(samples, config.t_end)?.clone();
    let drive_samples: Vec<_> = samples.iter().map(|s| (s.time, s.drive_power)).collect();
    let dephasing_samples: Vec<_> = samples
        .iter()
        .map(|s| (s.time, s.dephasing_power))
        .collect();
    let drive_energy = integrate_signed_power(&drive_samples)?;
    let dephasing_energy = integrate_signed_power(&dephasing_samples)?;
    let delta = at_end.bare_network_energy - samples[0].bare_network_energy;
    let ledger_residual = delta - drive_energy.energy_net - dephasing_energy.energy_net;
    let scale = delta
        .abs()
        .max(drive_energy.energy_net.abs())
        .max(dephasing_energy.energy_net.abs());
    let ledger_check_pass =
        ledger_residual.abs() <= LEDGER_ABSOLUTE_TOLERANCE + LEDGER_RELATIVE_TOLERANCE * scale;
    let maximum_trace_error = samples.iter().map(|s| s.trace_error).fold(0.0, f64::max);
    let maximum_hermiticity_error = samples
        .iter()
        .map(|s| s.hermiticity_error)
        .fold(0.0, f64::max);
    let worst_minimum_eigenvalue = samples
        .iter()
        .map(|s| s.minimum_eigenvalue)
        .fold(f64::INFINITY, f64::min);
    let all_finite = samples.iter().all(|s| s.all_finite);
    let maximum_top_level_population = maximum(samples, |s| s.load_populations[2]);
    let top_level_check_pass = maximum_top_level_population.0 < TOP_LEVEL_LIMIT;
    let physical_checks_pass = maximum_trace_error <= TRACE_TOLERANCE
        && maximum_hermiticity_error <= HERMITICITY_TOLERANCE
        && worst_minimum_eigenvalue >= -POSITIVITY_TOLERANCE
        && all_finite;
    Ok(CoherentDriveSummary {
        maximum_load_energy: maximum(samples, |s| s.load_energy),
        maximum_load_ergotropy: maximum(samples, |s| s.load_ergotropy),
        maximum_load_coherence_ergotropy: maximum(samples, |s| s.load_coherence_ergotropy),
        maximum_load_coherence_l1: maximum(samples, |s| s.load_coherence_l1),
        maximum_top_level_population,
        at_tau,
        at_end,
        drive_energy,
        dephasing_energy,
        delta_bare_network_energy: delta,
        ledger_residual,
        maximum_trace_error,
        maximum_hermiticity_error,
        worst_minimum_eigenvalue,
        maximum_drive_power_imaginary: samples
            .iter()
            .map(|s| s.drive_power_imaginary.abs())
            .fold(0.0, f64::max),
        maximum_dephasing_power_imaginary: samples
            .iter()
            .map(|s| s.dephasing_power_imaginary.abs())
            .fold(0.0, f64::max),
        all_finite,
        physical_checks_pass,
        ledger_check_pass,
        top_level_check_pass,
    })
}

pub fn sample_at(
    samples: &[CoherentDriveSample],
    time: f64,
) -> Result<&CoherentDriveSample, PhysicsError> {
    samples
        .iter()
        .find(|sample| (sample.time - time).abs() <= 1.0e-12)
        .ok_or_else(|| PhysicsError::InvalidTime(format!("saved time {time} not found")))
}

fn maximum(
    samples: &[CoherentDriveSample],
    value: impl Fn(&CoherentDriveSample) -> f64,
) -> (f64, f64) {
    let mut best = (value(&samples[0]), samples[0].time);
    for sample in &samples[1..] {
        let candidate = value(sample);
        if candidate > best.0 {
            best = (candidate, sample.time);
        }
    }
    best
}

fn validate_noise_sites(noise_sites: &[usize]) -> Result<(), PhysicsError> {
    let mut seen = [false; 3];
    for &site in noise_sites {
        if site >= 3 {
            return Err(PhysicsError::InvalidSubsystem(format!(
                "noise site {site} outside 0..3"
            )));
        }
        if seen[site] {
            return Err(PhysicsError::InvalidSubsystem(format!(
                "duplicate noise site {site}"
            )));
        }
        seen[site] = true;
    }
    Ok(())
}

#[cfg(test)]
fn dephasing_operators(
    operators: &Operators,
    gamma_phi: f64,
    noise_sites: &[usize],
) -> Vec<ComplexMatrix> {
    if gamma_phi == 0.0 {
        return Vec::new();
    }
    noise_sites
        .iter()
        .map(|&site| &operators.sigma_z_sites[site] * C64::new((gamma_phi / 2.0).sqrt(), 0.0))
        .collect()
}

fn validate_site_gamma_slice(site_gammas: &[f64], chain_length: usize) -> Result<(), PhysicsError> {
    if chain_length == 0 || site_gammas.len() != chain_length {
        return Err(PhysicsError::DimensionMismatch(format!(
            "expected {chain_length} site dephasing rates, got {}",
            site_gammas.len()
        )));
    }
    for (site, gamma) in site_gammas.iter().enumerate() {
        if !gamma.is_finite() || *gamma < 0.0 {
            return Err(PhysicsError::InvalidParameter(format!(
                "site {site} dephasing rate must be finite and nonnegative, got {gamma}"
            )));
        }
    }
    Ok(())
}

fn dephasing_operators_by_site(operators: &Operators, site_gammas: &[f64]) -> Vec<ComplexMatrix> {
    site_gammas
        .iter()
        .enumerate()
        .filter(|(_, gamma)| **gamma > 0.0)
        .map(|(site, gamma)| &operators.sigma_z_sites[site] * C64::new((gamma / 2.0).sqrt(), 0.0))
        .collect()
}

fn local_load_hamiltonian(params: &ModelParams) -> ComplexMatrix {
    ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
        params.load_dim,
        (0..params.load_dim).map(|level| C64::new(level as f64 * params.omega_load, 0.0)),
    ))
}

fn minimum_eigenvalue(matrix: &ComplexMatrix) -> f64 {
    let (_, schur) = Schur::new(matrix.clone()).unpack();
    (0..schur.nrows())
        .map(|index| schur[(index, index)].re)
        .fold(f64::INFINITY, f64::min)
}

fn finite_complex(value: &C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

fn validate_config(config: &CoherentDriveConfig) -> Result<(), PhysicsError> {
    for (name, value) in [
        ("omega0", config.omega0),
        ("omega_drive", config.omega_drive),
        ("tau", config.tau),
        ("t_end", config.t_end),
        ("dt", config.dt),
        ("save_interval", config.save_interval),
        ("gamma_phi", config.gamma_phi),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(PhysicsError::InvalidParameter(format!(
                "{name} must be finite and nonnegative, got {value}"
            )));
        }
    }
    if config.tau <= 0.0
        || config.t_end < config.tau
        || config.dt <= 0.0
        || config.save_interval <= 0.0
    {
        return Err(PhysicsError::InvalidParameter(
            "require tau, dt, save_interval > 0 and t_end >= tau".to_string(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> ModelParams {
        ModelParams::default()
    }

    #[test]
    fn drive_is_hermitian_at_multiple_times() {
        let operators = build_operators_for_chain(&params(), 3).unwrap();
        let config = CoherentDriveConfig::milestone_5b(0.0, 0.01);
        for time in [0.0, 0.3, 1.6, 3.2, 4.0] {
            let drive = drive_hamiltonian(time, &config, &operators.sigma_1_plus);
            assert!(hermiticity_error(&drive) < 1.0e-12);
        }
    }

    #[test]
    fn envelope_has_required_boundaries() {
        let tau = 3.2;
        assert_eq!(drive_envelope(0.0, tau), 0.0);
        assert!(drive_envelope(tau, tau).abs() < 1.0e-28);
        assert_eq!(drive_envelope(-0.1, tau), 0.0);
        assert_eq!(drive_envelope(tau + 0.1, tau), 0.0);
    }

    #[test]
    fn zero_drive_leaves_vacuum_unchanged() {
        let mut config = CoherentDriveConfig::milestone_5b(0.0, 0.01);
        config.omega0 = 0.0;
        config.tau = 0.02;
        config.t_end = 0.03;
        config.save_interval = 0.01;
        let run = run_coherent_drive(&params(), config).unwrap();
        assert!(run.summary.at_end.bare_network_energy.abs() < 1.0e-14);
        assert!(run.summary.at_end.load_energy.abs() < 1.0e-14);
    }

    #[test]
    fn zero_dephasing_has_zero_dephasing_power() {
        let mut config = CoherentDriveConfig::milestone_5b(0.0, 0.01);
        config.tau = 0.02;
        config.t_end = 0.03;
        config.save_interval = 0.01;
        let run = run_coherent_drive(&params(), config).unwrap();
        assert!(run.samples.iter().all(|s| s.dephasing_power == 0.0));
    }

    #[test]
    fn all_states_are_finite_and_ab_save_times_match() {
        let mut a = CoherentDriveConfig::milestone_5b(0.0, 0.01);
        a.tau = 0.02;
        a.t_end = 0.03;
        a.save_interval = 0.01;
        let mut b = a;
        b.gamma_phi = 0.5;
        let a_run = run_coherent_drive(&params(), a).unwrap();
        let b_run = run_coherent_drive(&params(), b).unwrap();
        assert!(a_run.samples.iter().all(|s| s.all_finite));
        assert!(b_run.samples.iter().all(|s| s.all_finite));
        assert_eq!(a_run.samples.len(), b_run.samples.len());
        for (left, right) in a_run.samples.iter().zip(&b_run.samples) {
            assert_eq!(left.time, right.time);
        }
    }

    #[test]
    fn local_dephasing_selects_exact_requested_site() {
        let operators = build_operators_for_chain(&params(), 3).unwrap();
        let collapses = dephasing_operators(&operators, 0.5, &[1]);
        assert_eq!(collapses.len(), 1);
        let expected = &operators.sigma_z_sites[1] * C64::new(0.25_f64.sqrt(), 0.0);
        assert_eq!(collapses[0], expected);
        assert!(dephasing_operators(&operators, 0.0, &[1]).is_empty());
    }

    #[test]
    fn local_dephasing_rejects_bad_site_lists() {
        let config = CoherentDriveConfig::milestone_5b(0.5, 0.01);
        assert!(run_coherent_drive_with_noise_sites(&params(), config, &[3]).is_err());
        assert!(run_coherent_drive_with_noise_sites(&params(), config, &[1, 1]).is_err());
        let negative = CoherentDriveConfig::milestone_5b(-0.5, 0.01);
        assert!(run_coherent_drive_with_noise_sites(&params(), negative, &[1]).is_err());
    }

    #[test]
    fn heterogeneous_dephasing_maps_rates_and_omits_zero() {
        let operators = build_operators_for_chain(&params(), 3).unwrap();
        let collapses = dephasing_operators_by_site(&operators, &[0.4, 0.5, 0.0]);
        assert_eq!(collapses.len(), 2);
        assert_eq!(
            collapses[0],
            &operators.sigma_z_sites[0] * C64::new((0.4_f64 / 2.0).sqrt(), 0.0)
        );
        assert_eq!(
            collapses[1],
            &operators.sigma_z_sites[1] * C64::new((0.5_f64 / 2.0).sqrt(), 0.0)
        );
    }

    #[test]
    fn heterogeneous_dephasing_rejects_negative_or_nonfinite_rates() {
        let config = CoherentDriveConfig::milestone_5b(0.5, 0.01);
        assert!(run_coherent_drive_with_site_gammas(&params(), config, [0.5, -0.1, 0.5]).is_err());
        assert!(
            run_coherent_drive_with_site_gammas(&params(), config, [0.5, f64::NAN, 0.5]).is_err()
        );
    }

    #[test]
    fn heterogeneous_all_equal_matches_existing_common_gamma_api() {
        let mut config = CoherentDriveConfig::milestone_5b(0.5, 0.01);
        config.tau = 0.02;
        config.t_end = 0.03;
        config.save_interval = 0.01;
        let existing = run_coherent_drive(&params(), config).unwrap();
        let heterogeneous =
            run_coherent_drive_with_site_gammas(&params(), config, [0.5; 3]).unwrap();
        assert_eq!(existing.samples.len(), heterogeneous.samples.len());
        for (left, right) in existing.states.iter().zip(&heterogeneous.states) {
            assert_eq!(left.time, right.time);
            assert_eq!(left.rho, right.rho);
        }
    }

    #[test]
    fn generalized_three_site_api_matches_existing_api_exactly() {
        let mut config = CoherentDriveConfig::milestone_5b(0.5, 0.01);
        config.tau = 0.02;
        config.t_end = 0.03;
        config.save_interval = 0.01;
        let existing = run_coherent_drive(&params(), config).unwrap();
        let generalized = run_coherent_drive_for_chain(&params(), config, 3, 0.5).unwrap();
        for (left, right) in existing.states.iter().zip(&generalized.states) {
            assert_eq!(left.time, right.time);
            assert_eq!(left.rho, right.rho);
        }
    }
}
