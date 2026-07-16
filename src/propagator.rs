//! Accuracy-first propagation by a dense matrix exponential.

use crate::error::PhysicsError;
use crate::liouvillian::{devectorize_column_major, vectorize_column_major};
use crate::matrix::{ComplexMatrix, C64};

#[derive(Debug, Clone)]
pub struct QuantumState {
    pub rho: ComplexMatrix,
    pub time: f64,
}

#[derive(Debug, Clone)]
pub struct DenseExponentialPropagator {
    liouvillian: ComplexMatrix,
    hilbert_dim: usize,
}

impl DenseExponentialPropagator {
    pub fn new(liouvillian: ComplexMatrix, hilbert_dim: usize) -> Result<Self, PhysicsError> {
        let expected = hilbert_dim * hilbert_dim;
        if liouvillian.shape() != (expected, expected) {
            return Err(PhysicsError::DimensionMismatch(format!(
                "Liouvillian is {:?}, expected ({expected}, {expected})",
                liouvillian.shape()
            )));
        }
        Ok(Self {
            liouvillian,
            hilbert_dim,
        })
    }

    pub fn propagate_times(
        &self,
        rho0: &ComplexMatrix,
        times: &[f64],
    ) -> Result<Vec<QuantumState>, PhysicsError> {
        if rho0.shape() != (self.hilbert_dim, self.hilbert_dim) {
            return Err(PhysicsError::DimensionMismatch(format!(
                "initial state is {:?}, expected ({dim}, {dim})",
                rho0.shape(),
                dim = self.hilbert_dim
            )));
        }
        if times.iter().any(|time| !time.is_finite() || *time < 0.0) {
            return Err(PhysicsError::InvalidTime(
                "all requested times must be finite and nonnegative".to_string(),
            ));
        }

        let rho0_vector = vectorize_column_major(rho0);
        let mut states = Vec::with_capacity(times.len());
        for &time in times {
            let propagator = (&self.liouvillian * C64::new(time, 0.0)).exp();
            let evolved = propagator * &rho0_vector;
            let rho = devectorize_column_major(&evolved, self.hilbert_dim, self.hilbert_dim)?;
            states.push(QuantumState { rho, time });
        }
        Ok(states)
    }
    /// Propagate on a uniform grid using one matrix exponential.
    ///
    /// This is numerically equivalent to repeatedly applying `exp(L * dt)`
    /// and is intended for diagnostic time series where recomputing a dense
    /// exponential at every requested time would be needlessly expensive.
    pub fn propagate_uniform(
        &self,
        rho0: &ComplexMatrix,
        end_time: f64,
        steps: usize,
    ) -> Result<Vec<QuantumState>, PhysicsError> {
        if rho0.shape() != (self.hilbert_dim, self.hilbert_dim) {
            return Err(PhysicsError::DimensionMismatch(format!(
                "initial state is {:?}, expected ({dim}, {dim})",
                rho0.shape(),
                dim = self.hilbert_dim
            )));
        }
        if !end_time.is_finite() || end_time < 0.0 {
            return Err(PhysicsError::InvalidTime(format!(
                "end time must be finite and nonnegative: {end_time}"
            )));
        }
        if steps == 0 {
            return Err(PhysicsError::InvalidTime(
                "uniform propagation requires at least one step".to_string(),
            ));
        }

        let dt = end_time / steps as f64;
        let step_propagator = (&self.liouvillian * C64::new(dt, 0.0)).exp();
        let mut vector = vectorize_column_major(rho0);
        let mut states = Vec::with_capacity(steps + 1);
        states.push(QuantumState {
            rho: rho0.clone(),
            time: 0.0,
        });
        for index in 1..=steps {
            vector = &step_propagator * vector;
            let rho = devectorize_column_major(&vector, self.hilbert_dim, self.hilbert_dim)?;
            states.push(QuantumState {
                rho,
                time: index as f64 * dt,
            });
        }
        Ok(states)
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::linalg::Schur;

    use super::*;
    use crate::liouvillian::build_liouvillian;
    use crate::matrix::{frobenius_norm, hermiticity_error, C64};

    fn qubit_ground() -> ComplexMatrix {
        ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(1.0, 0.0),
                C64::new(0.0, 0.0),
                C64::new(0.0, 0.0),
                C64::new(0.0, 0.0),
            ],
        )
    }

    fn qubit_excited() -> ComplexMatrix {
        ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(0.0, 0.0),
                C64::new(0.0, 0.0),
                C64::new(0.0, 0.0),
                C64::new(1.0, 0.0),
            ],
        )
    }

    fn minimum_hermitian_eigenvalue(matrix: &ComplexMatrix) -> f64 {
        let schur = Schur::new(matrix.clone());
        let (_, t) = schur.unpack();
        (0..t.nrows())
            .map(|i| t[(i, i)].re)
            .fold(f64::INFINITY, f64::min)
    }

    #[test]
    fn t_zero_returns_initial_state() {
        let h = ComplexMatrix::zeros(2, 2);
        let l = build_liouvillian(&h, &[]).unwrap();
        let propagator = DenseExponentialPropagator::new(l, 2).unwrap();
        let rho0 = qubit_excited();
        let states = propagator.propagate_times(&rho0, &[0.0]).unwrap();
        assert!(frobenius_norm(&(&states[0].rho - &rho0)) < 1.0e-13);
    }

    #[test]
    fn zero_hamiltonian_and_no_collapse_leave_state_unchanged() {
        let h = ComplexMatrix::zeros(2, 2);
        let l = build_liouvillian(&h, &[]).unwrap();
        let propagator = DenseExponentialPropagator::new(l, 2).unwrap();
        let rho0 = ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(0.7, 0.0),
                C64::new(0.1, -0.2),
                C64::new(0.1, 0.2),
                C64::new(0.3, 0.0),
            ],
        );
        for state in propagator.propagate_times(&rho0, &[0.0, 0.2, 2.0]).unwrap() {
            assert!(frobenius_norm(&(state.rho - &rho0)) < 1.0e-12);
        }
    }

    #[test]
    fn closed_evolution_matches_unitary_formula() {
        let h = ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(0.0, 0.0),
                C64::new(1.0, 0.0),
                C64::new(1.0, 0.0),
                C64::new(0.0, 0.0),
            ],
        );
        let t = 0.37;
        let l = build_liouvillian(&h, &[]).unwrap();
        let propagator = DenseExponentialPropagator::new(l, 2).unwrap();
        let rho0 = qubit_ground();
        let actual = propagator
            .propagate_times(&rho0, &[t])
            .unwrap()
            .remove(0)
            .rho;
        let u = (&h * C64::new(0.0, -t)).exp();
        let expected = &u * rho0 * u.adjoint();
        assert!(frobenius_norm(&(actual - expected)) < 1.0e-10);
    }

    #[test]
    fn amplitude_damping_matches_analytic_solution() {
        let gamma: f64 = 0.8;
        let sigma_minus = ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(0.0, 0.0),
                C64::new(1.0, 0.0),
                C64::new(0.0, 0.0),
                C64::new(0.0, 0.0),
            ],
        );
        let collapse = sigma_minus * C64::new(gamma.sqrt(), 0.0);
        let h = ComplexMatrix::zeros(2, 2);
        let l = build_liouvillian(&h, &[collapse]).unwrap();
        let propagator = DenseExponentialPropagator::new(l, 2).unwrap();
        let t: f64 = 1.3;
        let rho = propagator
            .propagate_times(&qubit_excited(), &[t])
            .unwrap()
            .remove(0)
            .rho;
        let excited = (-gamma * t).exp();
        assert!((rho[(1, 1)].re - excited).abs() < 1.0e-10);
        assert!((rho[(0, 0)].re - (1.0 - excited)).abs() < 1.0e-10);
    }

    #[test]
    fn pure_dephasing_preserves_populations_and_damps_coherence() {
        let gamma: f64 = 0.6;
        let sigma_z = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
            C64::new(1.0, 0.0),
            C64::new(-1.0, 0.0),
        ]));
        let collapse = sigma_z * C64::new((gamma / 2.0).sqrt(), 0.0);
        let h = ComplexMatrix::zeros(2, 2);
        let l = build_liouvillian(&h, &[collapse]).unwrap();
        let propagator = DenseExponentialPropagator::new(l, 2).unwrap();
        let rho0 = ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(0.5, 0.0),
                C64::new(0.5, 0.0),
                C64::new(0.5, 0.0),
                C64::new(0.5, 0.0),
            ],
        );
        let t: f64 = 1.7;
        let rho = propagator
            .propagate_times(&rho0, &[t])
            .unwrap()
            .remove(0)
            .rho;
        assert!((rho[(0, 0)].re - 0.5).abs() < 1.0e-11);
        assert!((rho[(1, 1)].re - 0.5).abs() < 1.0e-11);
        assert!((rho[(0, 1)].re - 0.5 * (-gamma * t).exp()).abs() < 1.0e-10);
    }

    #[test]
    fn physicality_is_preserved_for_damped_evolution() {
        let gamma: f64 = 0.4;
        let sigma_minus = ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(0.0, 0.0),
                C64::new(1.0, 0.0),
                C64::new(0.0, 0.0),
                C64::new(0.0, 0.0),
            ],
        );
        let collapse = sigma_minus * C64::new(gamma.sqrt(), 0.0);
        let h = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
            C64::new(0.0, 0.0),
            C64::new(1.0, 0.0),
        ]));
        let l = build_liouvillian(&h, &[collapse]).unwrap();
        let propagator = DenseExponentialPropagator::new(l, 2).unwrap();
        for state in propagator
            .propagate_times(&qubit_excited(), &[0.0, 0.1, 1.0, 5.0])
            .unwrap()
        {
            assert!((state.rho.trace().re - 1.0).abs() < 1.0e-10);
            assert!(hermiticity_error(&state.rho) < 1.0e-10);
            assert!(minimum_hermitian_eigenvalue(&state.rho) > -1.0e-10);
        }
    }
}
