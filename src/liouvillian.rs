//! Liouvillian construction with explicit column-major vectorization.
//!
//! Convention:
//! `vec(rho)` stacks the columns of `rho`, so that
//! `vec(A rho B) = (B^T \otimes A) vec(rho)`.

use nalgebra::DVector;

use crate::error::PhysicsError;
use crate::matrix::{eye, kron, ComplexMatrix, C64};

pub type ComplexVector = DVector<C64>;

pub fn vectorize_column_major(matrix: &ComplexMatrix) -> ComplexVector {
    ComplexVector::from_column_slice(matrix.as_slice())
}

pub fn devectorize_column_major(
    vector: &ComplexVector,
    rows: usize,
    cols: usize,
) -> Result<ComplexMatrix, PhysicsError> {
    if vector.len() != rows * cols {
        return Err(PhysicsError::DimensionMismatch(format!(
            "vector has length {}, expected {} for a {rows}x{cols} matrix",
            vector.len(),
            rows * cols
        )));
    }
    Ok(ComplexMatrix::from_column_slice(
        rows,
        cols,
        vector.as_slice(),
    ))
}

pub fn hamiltonian_superoperator(
    hamiltonian: &ComplexMatrix,
) -> Result<ComplexMatrix, PhysicsError> {
    if hamiltonian.nrows() != hamiltonian.ncols() {
        return Err(PhysicsError::DimensionMismatch(
            "Hamiltonian must be square".to_string(),
        ));
    }
    let dim = hamiltonian.nrows();
    let identity = eye(dim);
    let h_transpose = hamiltonian.transpose();
    Ok((kron(&identity, hamiltonian) - kron(&h_transpose, &identity)) * C64::new(0.0, -1.0))
}

pub fn lindblad_superoperator(collapse: &ComplexMatrix) -> Result<ComplexMatrix, PhysicsError> {
    if collapse.nrows() != collapse.ncols() {
        return Err(PhysicsError::DimensionMismatch(
            "collapse operator must be square".to_string(),
        ));
    }
    let dim = collapse.nrows();
    let identity = eye(dim);
    let l_dag_l = collapse.adjoint() * collapse;
    let jump = kron(&collapse.conjugate(), collapse);
    let left = kron(&identity, &l_dag_l);
    let right = kron(&l_dag_l.transpose(), &identity);
    Ok(jump - (left + right) * C64::new(0.5, 0.0))
}

/// Build the time-independent Lindblad generator.
///
/// Each collapse operator must already include its rate coefficient, e.g.
/// `sqrt(gamma) * L`.
pub fn build_liouvillian(
    hamiltonian: &ComplexMatrix,
    collapse_operators: &[ComplexMatrix],
) -> Result<ComplexMatrix, PhysicsError> {
    if hamiltonian.nrows() != hamiltonian.ncols() {
        return Err(PhysicsError::DimensionMismatch(
            "Hamiltonian must be square".to_string(),
        ));
    }
    let dim = hamiltonian.nrows();
    let mut generator = hamiltonian_superoperator(hamiltonian)?;
    for collapse in collapse_operators {
        if collapse.shape() != (dim, dim) {
            return Err(PhysicsError::DimensionMismatch(format!(
                "collapse operator is {:?}, expected ({dim}, {dim})",
                collapse.shape()
            )));
        }
        generator += lindblad_superoperator(collapse)?;
    }
    Ok(generator)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn column_major_round_trip_is_exact() {
        let matrix = ComplexMatrix::from_row_slice(
            2,
            3,
            &[
                C64::new(1.0, 1.0),
                C64::new(2.0, 2.0),
                C64::new(3.0, 3.0),
                C64::new(4.0, 4.0),
                C64::new(5.0, 5.0),
                C64::new(6.0, 6.0),
            ],
        );
        let vector = vectorize_column_major(&matrix);
        assert_eq!(vector[0], matrix[(0, 0)]);
        assert_eq!(vector[1], matrix[(1, 0)]);
        assert_eq!(vector[2], matrix[(0, 1)]);
        let restored = devectorize_column_major(&vector, 2, 3).unwrap();
        assert_eq!(restored, matrix);
    }

    #[test]
    fn twenty_four_dimensional_state_vectorizes_to_576() {
        let rho = ComplexMatrix::identity(24, 24) * C64::new(1.0 / 24.0, 0.0);
        assert_eq!(vectorize_column_major(&rho).len(), 576);
    }

    #[test]
    fn zero_generator_is_zero() {
        let h = ComplexMatrix::zeros(2, 2);
        let l = build_liouvillian(&h, &[]).unwrap();
        assert!(l.iter().all(|z| z.norm() == 0.0));
    }
}
