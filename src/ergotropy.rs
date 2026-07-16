use nalgebra::linalg::Schur;

use crate::error::PhysicsError;
use crate::matrix::{hermiticity_error, trace_real, ComplexMatrix, C64};

#[derive(Debug, Clone)]
pub struct ErgotropyResult {
    pub energy: f64,
    pub passive_energy: f64,
    pub ergotropy: f64,
    pub state_eigenvalues: Vec<f64>,
    pub hamiltonian_eigenvalues: Vec<f64>,
}

fn hermitian_eigenvalues(matrix: &ComplexMatrix, tolerance: f64) -> Result<Vec<f64>, PhysicsError> {
    let error = hermiticity_error(matrix);
    if error > tolerance {
        return Err(PhysicsError::NonHermitian { error, tolerance });
    }

    // For a Hermitian matrix, complex Schur form is diagonal up to numerical noise.
    let schur = Schur::new(matrix.clone());
    let (_q, t) = schur.unpack();

    let mut values = Vec::with_capacity(t.nrows());
    let mut off_diagonal_norm_sq = 0.0;
    for i in 0..t.nrows() {
        if t[(i, i)].im.abs() > 10.0 * tolerance {
            return Err(PhysicsError::EigenFailure(format!(
                "Hermitian eigenvalue has imaginary part {}",
                t[(i, i)].im
            )));
        }
        values.push(t[(i, i)].re);
        for j in 0..t.ncols() {
            if i != j {
                off_diagonal_norm_sq += t[(i, j)].norm_sqr();
            }
        }
    }

    if off_diagonal_norm_sq.sqrt() > 100.0 * tolerance {
        return Err(PhysicsError::EigenFailure(format!(
            "Schur form of Hermitian matrix is not diagonal enough: {}",
            off_diagonal_norm_sq.sqrt()
        )));
    }

    Ok(values)
}

pub fn ergotropy(
    rho: &ComplexMatrix,
    hamiltonian: &ComplexMatrix,
    tolerance: f64,
) -> Result<ErgotropyResult, PhysicsError> {
    if rho.nrows() != rho.ncols() || hamiltonian.nrows() != hamiltonian.ncols() {
        return Err(PhysicsError::DimensionMismatch(
            "rho and Hamiltonian must be square".to_string(),
        ));
    }
    if rho.shape() != hamiltonian.shape() {
        return Err(PhysicsError::DimensionMismatch(format!(
            "rho is {:?}, Hamiltonian is {:?}",
            rho.shape(),
            hamiltonian.shape()
        )));
    }
    if !(tolerance.is_finite() && tolerance > 0.0) {
        return Err(PhysicsError::DimensionMismatch(
            "tolerance must be finite and positive".to_string(),
        ));
    }

    let rho_h = (rho + rho.adjoint()) * C64::new(0.5, 0.0);
    let h_h = (hamiltonian + hamiltonian.adjoint()) * C64::new(0.5, 0.0);

    let trace = trace_real(&rho_h);
    if !trace.is_finite() || trace.abs() <= tolerance {
        return Err(PhysicsError::InvalidTrace { trace });
    }

    let mut state_eigenvalues = hermitian_eigenvalues(&rho_h, tolerance)?;
    let minimum = state_eigenvalues
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    if minimum < -tolerance {
        return Err(PhysicsError::NonPositiveState { minimum });
    }

    for value in &mut state_eigenvalues {
        if *value < 0.0 {
            *value = 0.0;
        }
    }
    let eigenvalue_sum: f64 = state_eigenvalues.iter().sum();
    if !eigenvalue_sum.is_finite() || eigenvalue_sum <= tolerance {
        return Err(PhysicsError::InvalidTrace {
            trace: eigenvalue_sum,
        });
    }
    for value in &mut state_eigenvalues {
        *value /= eigenvalue_sum;
    }
    state_eigenvalues.sort_by(|a, b| b.total_cmp(a));

    let mut hamiltonian_eigenvalues = hermitian_eigenvalues(&h_h, tolerance)?;
    hamiltonian_eigenvalues.sort_by(|a, b| a.total_cmp(b));

    let normalized_rho = rho_h * C64::new(1.0 / trace, 0.0);
    let energy = (&normalized_rho * &h_h).trace().re;
    let passive_energy: f64 = state_eigenvalues
        .iter()
        .zip(hamiltonian_eigenvalues.iter())
        .map(|(population, energy)| population * energy)
        .sum();
    let raw_ergotropy = energy - passive_energy;
    let ergotropy = if raw_ergotropy.abs() <= 10.0 * tolerance {
        0.0
    } else if raw_ergotropy < 0.0 {
        return Err(PhysicsError::EigenFailure(format!(
            "computed negative ergotropy {raw_ergotropy} beyond tolerance"
        )));
    } else {
        raw_ergotropy
    };

    Ok(ErgotropyResult {
        energy,
        passive_energy,
        ergotropy,
        state_eigenvalues,
        hamiltonian_eigenvalues,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::frobenius_norm;

    fn diag(values: &[f64]) -> ComplexMatrix {
        ComplexMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
            values.len(),
            values.iter().map(|&x| C64::new(x, 0.0)),
        ))
    }

    #[test]
    fn ground_state_has_zero_ergotropy() {
        let rho = diag(&[1.0, 0.0, 0.0]);
        let h = diag(&[0.0, 1.0, 2.0]);
        let result = ergotropy(&rho, &h, 1.0e-12).unwrap();
        assert!(result.ergotropy.abs() < 1.0e-12);
    }

    #[test]
    fn passive_mixture_has_zero_ergotropy() {
        let rho = diag(&[0.6, 0.3, 0.1]);
        let h = diag(&[0.0, 1.0, 2.0]);
        let result = ergotropy(&rho, &h, 1.0e-12).unwrap();
        assert!(result.ergotropy.abs() < 1.0e-12);
    }

    #[test]
    fn inverted_population_is_active() {
        let rho = diag(&[0.1, 0.3, 0.6]);
        let h = diag(&[0.0, 1.0, 2.0]);
        let result = ergotropy(&rho, &h, 1.0e-12).unwrap();
        assert!(result.ergotropy > 0.0);
    }

    #[test]
    fn coherent_pure_state_is_active() {
        let norm = 1.0 / 2.0_f64.sqrt();
        let ket = nalgebra::DVector::from_vec(vec![
            C64::new(norm, 0.0),
            C64::new(norm, 0.0),
            C64::new(0.0, 0.0),
        ]);
        let rho = &ket * ket.adjoint();
        let h = diag(&[0.0, 1.0, 2.0]);
        let result = ergotropy(&rho, &h, 1.0e-12).unwrap();
        assert!((result.ergotropy - 0.5).abs() < 1.0e-10);
    }

    #[test]
    fn rotated_hamiltonian_is_handled_without_diagonal_assumption() {
        let inv_sqrt_two = 1.0 / 2.0_f64.sqrt();
        let u = ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(inv_sqrt_two, 0.0),
                C64::new(-inv_sqrt_two, 0.0),
                C64::new(inv_sqrt_two, 0.0),
                C64::new(inv_sqrt_two, 0.0),
            ],
        );
        let h_diag = diag(&[0.0, 1.0]);
        let h = &u * h_diag * u.adjoint();
        let rho = diag(&[1.0, 0.0]);
        let result = ergotropy(&rho, &h, 1.0e-12).unwrap();
        assert!((result.ergotropy - 0.5).abs() < 1.0e-10);
    }

    #[test]
    fn tiny_negative_eigenvalue_is_clipped_and_renormalized() {
        let rho = diag(&[0.7, 0.3 + 1.0e-13, -1.0e-13]);
        let h = diag(&[0.0, 1.0, 2.0]);
        let result = ergotropy(&rho, &h, 1.0e-10).unwrap();
        assert!(result.state_eigenvalues.iter().all(|&x| x >= 0.0));
        assert!((result.state_eigenvalues.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn clearly_nonpositive_state_is_rejected() {
        let rho = diag(&[0.8, 0.3, -0.1]);
        let h = diag(&[0.0, 1.0, 2.0]);
        let error = ergotropy(&rho, &h, 1.0e-10).unwrap_err();
        assert!(matches!(error, PhysicsError::NonPositiveState { .. }));
    }

    #[test]
    fn passive_energy_uses_only_spectra() {
        let rho = diag(&[0.1, 0.3, 0.6]);
        let h = diag(&[0.0, 1.0, 2.0]);
        let result = ergotropy(&rho, &h, 1.0e-12).unwrap();
        assert!((result.passive_energy - 0.5).abs() < 1.0e-10);
        assert!(frobenius_norm(&(h.clone() - h.adjoint())) < 1.0e-12);
    }
}
