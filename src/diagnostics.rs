//! Static diagnostics for propagated quantum states.
//!
//! This module reads physical quantities from states without modifying,
//! renormalizing, or Hermitian-symmetrizing the supplied density matrix.

use nalgebra::linalg::Schur;

use crate::ergotropy::ergotropy;
use crate::error::PhysicsError;
use crate::matrix::{commutator, expectation, hermiticity_error, ComplexMatrix, C64};
use crate::operators::{ModelParams, Operators};
use crate::partial_trace::partial_trace;
use crate::propagator::QuantumState;

const DIAGNOSTIC_TOLERANCE: f64 = 1.0e-9;

#[derive(Debug, Clone)]
pub struct Diagnostics {
    pub time: f64,
    pub load_energy: f64,
    pub load_ergotropy: f64,
    pub load_passive_energy: f64,
    pub total_energy: f64,
    pub chain_energy: f64,
    pub interaction_energy: f64,
    pub load_energy_current: f64,
    pub source_power: f64,
    pub dephasing_power: f64,
    pub trace_error: f64,
    pub hermiticity_error: f64,
    pub minimum_eigenvalue: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SignedEnergyIntegral {
    pub energy_net: f64,
    pub energy_in: f64,
    pub energy_out: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntegratedPowers {
    pub source_energy_net: f64,
    pub source_energy_in: f64,
    pub source_energy_out: f64,
    pub dephasing_energy_net: f64,
    pub dephasing_energy_in: f64,
    pub dephasing_energy_out: f64,
}

fn validate_square_shape(
    matrix: &ComplexMatrix,
    dim: usize,
    name: &str,
) -> Result<(), PhysicsError> {
    if matrix.shape() != (dim, dim) {
        return Err(PhysicsError::DimensionMismatch(format!(
            "{name} is {:?}, expected ({dim}, {dim})",
            matrix.shape()
        )));
    }
    Ok(())
}

fn validate_finite(matrix: &ComplexMatrix, name: &str) -> Result<(), PhysicsError> {
    if matrix
        .iter()
        .any(|z| !z.re.is_finite() || !z.im.is_finite())
    {
        return Err(PhysicsError::EigenFailure(format!(
            "{name} contains a non-finite value"
        )));
    }
    Ok(())
}

fn minimum_hermitian_eigenvalue(
    matrix: &ComplexMatrix,
    tolerance: f64,
) -> Result<f64, PhysicsError> {
    let error = hermiticity_error(matrix);
    if error > tolerance {
        return Err(PhysicsError::NonHermitian { error, tolerance });
    }

    let schur = Schur::new(matrix.clone());
    let (_q, t) = schur.unpack();
    let mut minimum = f64::INFINITY;
    let mut off_diagonal_norm_sq = 0.0;
    for row in 0..t.nrows() {
        if t[(row, row)].im.abs() > 10.0 * tolerance {
            return Err(PhysicsError::EigenFailure(format!(
                "Hermitian eigenvalue has imaginary part {}",
                t[(row, row)].im
            )));
        }
        minimum = minimum.min(t[(row, row)].re);
        for col in 0..t.ncols() {
            if row != col {
                off_diagonal_norm_sq += t[(row, col)].norm_sqr();
            }
        }
    }
    if off_diagonal_norm_sq.sqrt() > 100.0 * tolerance {
        return Err(PhysicsError::EigenFailure(format!(
            "Schur form is not diagonal enough: {}",
            off_diagonal_norm_sq.sqrt()
        )));
    }
    Ok(minimum)
}

/// Apply one Lindblad dissipator directly to a density matrix.
///
/// The collapse operator must already contain its rate coefficient.
pub fn lindblad_action(
    collapse: &ComplexMatrix,
    rho: &ComplexMatrix,
) -> Result<ComplexMatrix, PhysicsError> {
    if rho.nrows() != rho.ncols() {
        return Err(PhysicsError::DimensionMismatch(
            "density matrix must be square".to_string(),
        ));
    }
    validate_square_shape(collapse, rho.nrows(), "collapse operator")?;
    let l_dag_l = collapse.adjoint() * collapse;
    Ok(
        collapse * rho * collapse.adjoint()
            - (&l_dag_l * rho + rho * &l_dag_l) * C64::new(0.5, 0.0),
    )
}

fn dissipator_power(
    rho: &ComplexMatrix,
    hamiltonian: &ComplexMatrix,
    collapse: &ComplexMatrix,
) -> Result<f64, PhysicsError> {
    Ok(expectation(&lindblad_action(collapse, rho)?, hamiltonian).re)
}

fn local_load_hamiltonian(params: &ModelParams, operators: &Operators) -> ComplexMatrix {
    operators.b_load_local.adjoint() * &operators.b_load_local * C64::new(params.omega_load, 0.0)
}

/// Diagnose one state without modifying it.
pub fn diagnose_state(
    state: &QuantumState,
    operators: &Operators,
    params: &ModelParams,
    source_collapse: Option<&ComplexMatrix>,
    dephasing_collapses: &[ComplexMatrix],
) -> Result<Diagnostics, PhysicsError> {
    let full_dim: usize = operators.dims.iter().product();
    validate_square_shape(&state.rho, full_dim, "density matrix")?;
    validate_square_shape(&operators.h_total, full_dim, "total Hamiltonian")?;
    validate_finite(&state.rho, "density matrix")?;
    if !state.time.is_finite() {
        return Err(PhysicsError::InvalidTime(format!(
            "diagnostic time is not finite: {}",
            state.time
        )));
    }

    let trace = state.rho.trace();
    let trace_error = (trace - C64::new(1.0, 0.0)).norm();
    let state_hermiticity_error = hermiticity_error(&state.rho);
    let minimum_eigenvalue = minimum_hermitian_eigenvalue(&state.rho, DIAGNOSTIC_TOLERANCE)?;
    if trace_error > DIAGNOSTIC_TOLERANCE {
        return Err(PhysicsError::InvalidTrace { trace: trace.re });
    }
    if minimum_eigenvalue < -DIAGNOSTIC_TOLERANCE {
        return Err(PhysicsError::NonPositiveState {
            minimum: minimum_eigenvalue,
        });
    }

    let rho_load = partial_trace(&state.rho, &operators.dims, &[3])?;
    let h_load_local = local_load_hamiltonian(params, operators);
    validate_square_shape(&h_load_local, params.load_dim, "local load Hamiltonian")?;
    let load_result = ergotropy(&rho_load, &h_load_local, DIAGNOSTIC_TOLERANCE)?;

    let total_energy = expectation(&state.rho, &operators.h_total).re;
    let chain_energy = expectation(&state.rho, &operators.h_chain).re;
    let interaction_energy = expectation(&state.rho, &operators.h_interaction).re;

    // Positive sign means energy enters the load:
    // d<H_load>/dt = i Tr(rho [H_interaction, H_load]).
    let load_commutator = commutator(&operators.h_interaction, &operators.h_load);
    let load_energy_current = (expectation(&state.rho, &load_commutator) * C64::new(0.0, 1.0)).re;

    let source_power = match source_collapse {
        Some(collapse) => dissipator_power(&state.rho, &operators.h_total, collapse)?,
        None => 0.0,
    };
    let mut dephasing_power = 0.0;
    for collapse in dephasing_collapses {
        dephasing_power += dissipator_power(&state.rho, &operators.h_total, collapse)?;
    }

    let diagnostics = Diagnostics {
        time: state.time,
        load_energy: load_result.energy,
        load_ergotropy: load_result.ergotropy,
        load_passive_energy: load_result.passive_energy,
        total_energy,
        chain_energy,
        interaction_energy,
        load_energy_current,
        source_power,
        dephasing_power,
        trace_error,
        hermiticity_error: state_hermiticity_error,
        minimum_eigenvalue,
    };

    let scalar_values = [
        diagnostics.load_energy,
        diagnostics.load_ergotropy,
        diagnostics.load_passive_energy,
        diagnostics.total_energy,
        diagnostics.chain_energy,
        diagnostics.interaction_energy,
        diagnostics.load_energy_current,
        diagnostics.source_power,
        diagnostics.dephasing_power,
        diagnostics.trace_error,
        diagnostics.hermiticity_error,
        diagnostics.minimum_eigenvalue,
    ];
    if scalar_values.iter().any(|value| !value.is_finite()) {
        return Err(PhysicsError::EigenFailure(
            "diagnostics contain a non-finite value".to_string(),
        ));
    }

    Ok(diagnostics)
}

/// Diagnose several states in the supplied order.
pub fn diagnose_states(
    states: &[QuantumState],
    operators: &Operators,
    params: &ModelParams,
    source_collapse: Option<&ComplexMatrix>,
    dephasing_collapses: &[ComplexMatrix],
) -> Result<Vec<Diagnostics>, PhysicsError> {
    states
        .iter()
        .map(|state| {
            diagnose_state(
                state,
                operators,
                params,
                source_collapse,
                dephasing_collapses,
            )
        })
        .collect()
}

/// Integrate signed power samples with the trapezoidal rule.
///
/// `energy_in` integrates the positive part and `energy_out` integrates the
/// magnitude of the negative part. When adjacent samples have opposite signs,
/// linear interpolation is assumed and the interval is split at the zero
/// crossing. Time values must be finite and strictly increasing.
pub fn integrate_signed_power(
    samples: &[(f64, f64)],
) -> Result<SignedEnergyIntegral, PhysicsError> {
    if samples.len() < 2 {
        return Err(PhysicsError::InvalidTime(
            "at least two power samples are required".to_string(),
        ));
    }
    let mut result = SignedEnergyIntegral {
        energy_net: 0.0,
        energy_in: 0.0,
        energy_out: 0.0,
    };
    for window in samples.windows(2) {
        let (t0, p0) = window[0];
        let (t1, p1) = window[1];
        if !t0.is_finite() || !t1.is_finite() || !p0.is_finite() || !p1.is_finite() {
            return Err(PhysicsError::InvalidTime(
                "power samples must be finite".to_string(),
            ));
        }
        let dt = t1 - t0;
        if dt <= 0.0 {
            return Err(PhysicsError::InvalidTime(
                "power sample times must be strictly increasing".to_string(),
            ));
        }
        result.energy_net += 0.5 * (p0 + p1) * dt;

        if p0 >= 0.0 && p1 >= 0.0 {
            result.energy_in += 0.5 * (p0 + p1) * dt;
        } else if p0 <= 0.0 && p1 <= 0.0 {
            result.energy_out += -0.5 * (p0 + p1) * dt;
        } else {
            // Assume linear interpolation within the interval and split the
            // positive and negative triangular areas at the zero crossing.
            let alpha = p0.abs() / (p0.abs() + p1.abs());
            let dt_left = alpha * dt;
            let dt_right = dt - dt_left;

            if p0 > 0.0 {
                result.energy_in += 0.5 * p0 * dt_left;
                result.energy_out += 0.5 * (-p1) * dt_right;
            } else {
                result.energy_out += 0.5 * (-p0) * dt_left;
                result.energy_in += 0.5 * p1 * dt_right;
            }
        }
    }
    Ok(result)
}

/// Integrate source and dephasing power columns from a diagnostic time series.
pub fn integrate_diagnostic_powers(
    diagnostics: &[Diagnostics],
) -> Result<IntegratedPowers, PhysicsError> {
    let source_samples: Vec<(f64, f64)> = diagnostics
        .iter()
        .map(|sample| (sample.time, sample.source_power))
        .collect();
    let dephasing_samples: Vec<(f64, f64)> = diagnostics
        .iter()
        .map(|sample| (sample.time, sample.dephasing_power))
        .collect();
    let source = integrate_signed_power(&source_samples)?;
    let dephasing = integrate_signed_power(&dephasing_samples)?;
    Ok(IntegratedPowers {
        source_energy_net: source.energy_net,
        source_energy_in: source.energy_in,
        source_energy_out: source.energy_out,
        dephasing_energy_net: dephasing.energy_net,
        dephasing_energy_in: dephasing.energy_in,
        dephasing_energy_out: dephasing.energy_out,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{eye, kron};
    use crate::operators::build_operators;

    fn basis_density(dim: usize, index: usize) -> ComplexMatrix {
        let mut rho = ComplexMatrix::zeros(dim, dim);
        rho[(index, index)] = C64::new(1.0, 0.0);
        rho
    }

    fn full_state_from_load(load_rho: &ComplexMatrix) -> ComplexMatrix {
        let q0 = basis_density(2, 0);
        kron(&kron(&kron(&q0, &q0), &q0), load_rho)
    }

    fn state(rho: ComplexMatrix) -> QuantumState {
        QuantumState { rho, time: 0.0 }
    }

    #[test]
    fn vacuum_has_zero_load_energy_and_ergotropy() {
        let params = ModelParams::default();
        let ops = build_operators(&params).unwrap();
        let report =
            diagnose_state(&state(basis_density(24, 0)), &ops, &params, None, &[]).unwrap();
        assert!(report.load_energy.abs() < 1.0e-12);
        assert!(report.load_ergotropy.abs() < 1.0e-12);
        assert!(report.load_passive_energy.abs() < 1.0e-12);
    }

    #[test]
    fn pure_excited_load_has_positive_ergotropy() {
        let params = ModelParams::default();
        let ops = build_operators(&params).unwrap();
        let load = basis_density(3, 1);
        let report = diagnose_state(
            &state(full_state_from_load(&load)),
            &ops,
            &params,
            None,
            &[],
        )
        .unwrap();
        assert!((report.load_energy - 1.0).abs() < 1.0e-12);
        assert!(report.load_ergotropy > 0.9);
    }

    #[test]
    fn passive_load_mixture_has_zero_ergotropy() {
        let params = ModelParams::default();
        let ops = build_operators(&params).unwrap();
        let load = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
            C64::new(0.6, 0.0),
            C64::new(0.3, 0.0),
            C64::new(0.1, 0.0),
        ]));
        let report = diagnose_state(
            &state(full_state_from_load(&load)),
            &ops,
            &params,
            None,
            &[],
        )
        .unwrap();
        assert!(report.load_ergotropy.abs() < 1.0e-10);
    }

    #[test]
    fn injection_power_into_vacuum_is_positive() {
        let params = ModelParams::default();
        let ops = build_operators(&params).unwrap();
        let source = &ops.sigma_1_plus * C64::new(0.1_f64.sqrt(), 0.0);
        let report = diagnose_state(
            &state(basis_density(24, 0)),
            &ops,
            &params,
            Some(&source),
            &[],
        )
        .unwrap();
        assert!(report.source_power > 0.0);
    }

    #[test]
    fn absent_dephasing_has_zero_power() {
        let params = ModelParams::default();
        let ops = build_operators(&params).unwrap();
        let report =
            diagnose_state(&state(basis_density(24, 0)), &ops, &params, None, &[]).unwrap();
        assert!(report.dephasing_power.abs() < 1.0e-14);
    }

    #[test]
    fn energy_decomposition_is_consistent() {
        let params = ModelParams::default();
        let ops = build_operators(&params).unwrap();
        let inv_sqrt_two = 1.0 / 2.0_f64.sqrt();
        let mut ket = nalgebra::DVector::from_element(24, C64::new(0.0, 0.0));
        ket[3] = C64::new(inv_sqrt_two, 0.0); // |000,0> and |001,0> coherence below.
        ket[6] = C64::new(0.0, inv_sqrt_two);
        let rho = &ket * ket.adjoint();
        let report = diagnose_state(&state(rho), &ops, &params, None, &[]).unwrap();
        assert!(
            (report.total_energy
                - report.chain_energy
                - report.load_energy
                - report.interaction_energy)
                .abs()
                < 1.0e-10
        );
    }

    #[test]
    fn load_current_sign_matches_short_time_energy_change() {
        let params = ModelParams::default();
        let ops = build_operators(&params).unwrap();
        // Superpose |q3=1, load=0> and |q3=0, load=1> with a phase giving nonzero current.
        let mut ket = nalgebra::DVector::from_element(24, C64::new(0.0, 0.0));
        let index_q3_excited = 1 * 3; // |0,0,1,0>
        let index_load_excited = 1; // |0,0,0,1>
        let inv_sqrt_two = 1.0 / 2.0_f64.sqrt();
        ket[index_q3_excited] = C64::new(inv_sqrt_two, 0.0);
        ket[index_load_excited] = C64::new(0.0, inv_sqrt_two);
        let rho0 = &ket * ket.adjoint();
        let report = diagnose_state(&state(rho0.clone()), &ops, &params, None, &[]).unwrap();

        let dt = 1.0e-6;
        let u = (&ops.h_total * C64::new(0.0, -dt)).exp();
        let rho1 = &u * rho0 * u.adjoint();
        let e0 = report.load_energy;
        let e1 = diagnose_state(
            &QuantumState {
                rho: rho1,
                time: dt,
            },
            &ops,
            &params,
            None,
            &[],
        )
        .unwrap()
        .load_energy;
        let numerical_derivative = (e1 - e0) / dt;
        assert!(report.load_energy_current.abs() > 1.0e-6);
        assert_eq!(
            report.load_energy_current.is_sign_positive(),
            numerical_derivative.is_sign_positive()
        );
        assert!((report.load_energy_current - numerical_derivative).abs() < 1.0e-5);
    }

    #[test]
    fn trapezoidal_integral_handles_constant_and_sign_changing_power() {
        let constant = integrate_signed_power(&[(0.0, 2.0), (1.0, 2.0), (3.0, 2.0)]).unwrap();
        assert!((constant.energy_net - 6.0).abs() < 1.0e-12);
        assert!((constant.energy_in - 6.0).abs() < 1.0e-12);
        assert!(constant.energy_out.abs() < 1.0e-12);

        let single_crossing = integrate_signed_power(&[(0.0, 1.0), (1.0, -1.0)]).unwrap();
        assert!(single_crossing.energy_net.abs() < 1.0e-12);
        assert!((single_crossing.energy_in - 0.25).abs() < 1.0e-12);
        assert!((single_crossing.energy_out - 0.25).abs() < 1.0e-12);

        let changing = integrate_signed_power(&[(0.0, 1.0), (1.0, -1.0), (2.0, 1.0)]).unwrap();
        assert!(changing.energy_net.abs() < 1.0e-12);
        assert!((changing.energy_in - 0.5).abs() < 1.0e-12);
        assert!((changing.energy_out - 0.5).abs() < 1.0e-12);
    }

    #[test]
    fn all_diagnostics_are_finite_and_reduced_trace_is_one() {
        let params = ModelParams::default();
        let ops = build_operators(&params).unwrap();
        let rho = eye(24) * C64::new(1.0 / 24.0, 0.0);
        let report = diagnose_state(&state(rho.clone()), &ops, &params, None, &[]).unwrap();
        let load = partial_trace(&rho, &ops.dims, &[3]).unwrap();
        assert!((load.trace().re - 1.0).abs() < 1.0e-12);
        let values = [
            report.load_energy,
            report.load_ergotropy,
            report.load_passive_energy,
            report.total_energy,
            report.chain_energy,
            report.interaction_energy,
            report.load_energy_current,
            report.source_power,
            report.dephasing_power,
            report.trace_error,
            report.hermiticity_error,
            report.minimum_eigenvalue,
        ];
        assert!(values.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn diagnostic_power_integral_keeps_source_and_dephasing_separate() {
        let samples = vec![
            Diagnostics {
                time: 0.0,
                load_energy: 0.0,
                load_ergotropy: 0.0,
                load_passive_energy: 0.0,
                total_energy: 0.0,
                chain_energy: 0.0,
                interaction_energy: 0.0,
                load_energy_current: 0.0,
                source_power: 2.0,
                dephasing_power: -1.0,
                trace_error: 0.0,
                hermiticity_error: 0.0,
                minimum_eigenvalue: 0.0,
            },
            Diagnostics {
                time: 2.0,
                source_power: 2.0,
                dephasing_power: 1.0,
                ..Diagnostics {
                    time: 0.0,
                    load_energy: 0.0,
                    load_ergotropy: 0.0,
                    load_passive_energy: 0.0,
                    total_energy: 0.0,
                    chain_energy: 0.0,
                    interaction_energy: 0.0,
                    load_energy_current: 0.0,
                    source_power: 0.0,
                    dephasing_power: 0.0,
                    trace_error: 0.0,
                    hermiticity_error: 0.0,
                    minimum_eigenvalue: 0.0,
                }
            },
        ];
        let integrated = integrate_diagnostic_powers(&samples).unwrap();
        assert!((integrated.source_energy_net - 4.0).abs() < 1.0e-12);
        assert!((integrated.source_energy_in - 4.0).abs() < 1.0e-12);
        assert!(integrated.source_energy_out.abs() < 1.0e-12);
        assert!(integrated.dephasing_energy_net.abs() < 1.0e-12);
        assert!((integrated.dephasing_energy_in - 0.5).abs() < 1.0e-12);
        assert!((integrated.dephasing_energy_out - 0.5).abs() < 1.0e-12);
    }
}
