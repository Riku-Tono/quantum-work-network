//! Ideal instantaneous local-unitary extraction of load ergotropy.

use nalgebra::{linalg::Schur, DVector};

use crate::error::PhysicsError;
use crate::matrix::{frobenius_norm, hermiticity_error, ComplexMatrix, C64};
use crate::partial_trace::partial_trace;

pub const DEGENERACY_TOLERANCE: f64 = 1.0e-12;
pub const ENTROPY_EIGENVALUE_TOLERANCE: f64 = 1.0e-10;

#[derive(Debug, Clone)]
pub struct PassiveExtraction {
    pub unitary: ComplexMatrix,
    pub passive_state: ComplexMatrix,
    pub transformed_state: ComplexMatrix,
    pub state_eigenvalues_descending: Vec<f64>,
    pub hamiltonian_energies_ascending: Vec<f64>,
    pub energy_before: f64,
    pub passive_energy: f64,
    pub extracted_work: f64,
    pub unitary_error: f64,
    pub passive_mapping_error: f64,
}

pub fn build_passive_extraction(
    rho: &ComplexMatrix,
    hamiltonian: &ComplexMatrix,
) -> Result<PassiveExtraction, PhysicsError> {
    validate_hermitian_pair(rho, hamiltonian)?;
    let dim = rho.nrows();
    let (state_values, state_vectors) = sorted_eigensystem(rho, false)?;
    let (energies, energy_vectors) = sorted_eigensystem(hamiltonian, true)?;
    let unitary = &energy_vectors * state_vectors.adjoint();
    let diagonal = ComplexMatrix::from_diagonal(&DVector::from_iterator(
        dim,
        state_values.iter().map(|&value| C64::new(value, 0.0)),
    ));
    let passive_state = &energy_vectors * diagonal * energy_vectors.adjoint();
    let transformed_state = &unitary * rho * unitary.adjoint();
    let identity = ComplexMatrix::identity(dim, dim);
    let unitary_error = frobenius_norm(&(unitary.adjoint() * &unitary - identity));
    let passive_mapping_error = frobenius_norm(&(&transformed_state - &passive_state));
    let energy_before = (rho * hamiltonian).trace().re;
    let passive_energy = (&passive_state * hamiltonian).trace().re;
    Ok(PassiveExtraction {
        unitary,
        passive_state,
        transformed_state,
        state_eigenvalues_descending: state_values,
        hamiltonian_energies_ascending: energies,
        energy_before,
        passive_energy,
        extracted_work: energy_before - passive_energy,
        unitary_error,
        passive_mapping_error,
    })
}

pub fn von_neumann_entropy(rho: &ComplexMatrix) -> Result<f64, PhysicsError> {
    if rho.nrows() != rho.ncols() || !rho.iter().all(finite) {
        return Err(PhysicsError::InvalidParameter(
            "entropy state must be finite and square".to_string(),
        ));
    }
    if hermiticity_error(rho) > 1.0e-9 {
        return Err(PhysicsError::NonHermitian {
            error: hermiticity_error(rho),
            tolerance: 1.0e-9,
        });
    }
    let (mut values, _) = sorted_eigensystem(rho, true)?;
    if values
        .iter()
        .any(|&value| value < -ENTROPY_EIGENVALUE_TOLERANCE)
    {
        return Err(PhysicsError::NonPositiveState {
            minimum: values.iter().copied().fold(f64::INFINITY, f64::min),
        });
    }
    for value in &mut values {
        if *value < 0.0 {
            *value = 0.0;
        }
    }
    let sum: f64 = values.iter().sum();
    if sum <= 0.0 || !sum.is_finite() {
        return Err(PhysicsError::InvalidTrace { trace: sum });
    }
    Ok(values
        .into_iter()
        .map(|value| value / sum)
        .filter(|&value| value > 0.0)
        .map(|value| -value * value.ln())
        .sum())
}

pub fn mutual_information(
    rho: &ComplexMatrix,
    dims: &[usize],
    chain_subsystems: &[usize],
    load_subsystems: &[usize],
) -> Result<f64, PhysicsError> {
    let chain = partial_trace(rho, dims, chain_subsystems)?;
    let load = partial_trace(rho, dims, load_subsystems)?;
    Ok(von_neumann_entropy(&chain)? + von_neumann_entropy(&load)? - von_neumann_entropy(rho)?)
}

pub fn sorted_hermitian_eigenvalues(matrix: &ComplexMatrix) -> Result<Vec<f64>, PhysicsError> {
    Ok(sorted_eigensystem(matrix, true)?.0)
}

fn validate_hermitian_pair(
    rho: &ComplexMatrix,
    hamiltonian: &ComplexMatrix,
) -> Result<(), PhysicsError> {
    if rho.nrows() == 0 || rho.nrows() != rho.ncols() || rho.shape() != hamiltonian.shape() {
        return Err(PhysicsError::DimensionMismatch(
            "rho and Hamiltonian must be same-size nonempty square matrices".to_string(),
        ));
    }
    if !rho.iter().all(finite) || !hamiltonian.iter().all(finite) {
        return Err(PhysicsError::InvalidParameter(
            "rho and Hamiltonian must contain only finite values".to_string(),
        ));
    }
    for matrix in [rho, hamiltonian] {
        let error = hermiticity_error(matrix);
        if error > 1.0e-9 {
            return Err(PhysicsError::NonHermitian {
                error,
                tolerance: 1.0e-9,
            });
        }
    }
    Ok(())
}

fn sorted_eigensystem(
    matrix: &ComplexMatrix,
    ascending: bool,
) -> Result<(Vec<f64>, ComplexMatrix), PhysicsError> {
    if matrix.nrows() != matrix.ncols() || !matrix.iter().all(finite) {
        return Err(PhysicsError::InvalidParameter(
            "eigensystem matrix must be finite and square".to_string(),
        ));
    }
    let dim = matrix.nrows();
    let (vectors, triangular) = Schur::new(matrix.clone()).unpack();
    let mut entries: Vec<(usize, f64, DVector<C64>)> = (0..dim)
        .map(|index| {
            let mut vector = vectors.column(index).into_owned();
            fix_vector_phase(&mut vector);
            (index, triangular[(index, index)].re, vector)
        })
        .collect();
    entries.sort_by(|left, right| {
        let difference = left.1 - right.1;
        if difference.abs() <= DEGENERACY_TOLERANCE {
            left.0.cmp(&right.0)
        } else if ascending {
            left.1.total_cmp(&right.1)
        } else {
            right.1.total_cmp(&left.1)
        }
    });
    let values = entries.iter().map(|entry| entry.1).collect();
    let columns: Vec<_> = entries.into_iter().map(|entry| entry.2).collect();
    Ok((values, ComplexMatrix::from_columns(&columns)))
}

fn fix_vector_phase(vector: &mut DVector<C64>) {
    let pivot = vector
        .iter()
        .enumerate()
        .max_by(|left, right| left.1.norm().total_cmp(&right.1.norm()))
        .map(|(index, _)| index)
        .unwrap_or(0);
    let value = vector[pivot];
    if value.norm() > 0.0 {
        *vector *= value.conj() / value.norm();
    }
}

fn finite(value: &C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{eye, kron};

    fn diag(values: &[f64]) -> ComplexMatrix {
        ComplexMatrix::from_diagonal(&DVector::from_iterator(
            values.len(),
            values.iter().map(|&value| C64::new(value, 0.0)),
        ))
    }

    #[test]
    fn known_two_level_state_is_made_passive() {
        let result = build_passive_extraction(&diag(&[0.2, 0.8]), &diag(&[0.0, 1.0])).unwrap();
        assert!((result.extracted_work - 0.6).abs() < 1.0e-12);
        assert!(frobenius_norm(&(result.transformed_state - diag(&[0.8, 0.2]))) < 1.0e-12);
    }

    #[test]
    fn passive_state_has_zero_extracted_work() {
        let result = build_passive_extraction(&diag(&[0.8, 0.2]), &diag(&[0.0, 1.0])).unwrap();
        assert!(result.extracted_work.abs() < 1.0e-12);
    }

    #[test]
    fn pure_superposition_matches_analytic_work() {
        let rho = ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(0.5, 0.0),
                C64::new(0.5, 0.0),
                C64::new(0.5, 0.0),
                C64::new(0.5, 0.0),
            ],
        );
        let result = build_passive_extraction(&rho, &diag(&[0.0, 1.0])).unwrap();
        assert!((result.extracted_work - 0.5).abs() < 1.0e-12);
        assert!(result.unitary_error < 1.0e-12);
    }

    #[test]
    fn degenerate_state_is_finite_and_reproducible() {
        let rho = diag(&[0.5, 0.5, 0.0]);
        let h = diag(&[0.0, 1.0, 2.0]);
        let first = build_passive_extraction(&rho, &h).unwrap();
        let second = build_passive_extraction(&rho, &h).unwrap();
        assert!(first.unitary.iter().all(finite));
        assert!(frobenius_norm(&(first.unitary - second.unitary)) < 1.0e-14);
    }

    #[test]
    fn tensor_order_has_load_as_fastest_index() {
        let swap = ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(0.0, 0.0),
                C64::new(1.0, 0.0),
                C64::new(1.0, 0.0),
                C64::new(0.0, 0.0),
            ],
        );
        let full = kron(&eye(8), &swap);
        let mut ket = DVector::from_element(16, C64::new(0.0, 0.0));
        ket[3 * 2] = C64::new(1.0, 0.0);
        let mapped = full * ket;
        assert!((mapped[3 * 2 + 1] - C64::new(1.0, 0.0)).norm() < 1.0e-12);
    }

    #[test]
    fn local_unitary_preserves_chain_and_mutual_information() {
        let norm = 1.0 / 2.0_f64.sqrt();
        let ket = DVector::from_vec(vec![
            C64::new(norm, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(norm, 0.0),
        ]);
        let rho = &ket * ket.adjoint();
        let phase = ComplexMatrix::from_diagonal(&DVector::from_vec(vec![
            C64::new(1.0, 0.0),
            C64::new(0.0, 1.0),
        ]));
        let full = kron(&eye(2), &phase);
        let after = &full * &rho * full.adjoint();
        let chain_before = partial_trace(&rho, &[2, 2], &[0]).unwrap();
        let chain_after = partial_trace(&after, &[2, 2], &[0]).unwrap();
        assert!(frobenius_norm(&(chain_before - chain_after)) < 1.0e-12);
        let before_mi = mutual_information(&rho, &[2, 2], &[0], &[1]).unwrap();
        let after_mi = mutual_information(&after, &[2, 2], &[0], &[1]).unwrap();
        assert!((before_mi - after_mi).abs() < 1.0e-12);
    }

    #[test]
    fn switch_work_sign_convention_is_direct_energy_difference() {
        let rho = diag(&[1.0, 0.0]);
        let h_on = diag(&[1.0, 0.0]);
        let h_off = diag(&[0.0, 0.0]);
        let work = (&rho * (&h_off - &h_on)).trace().re;
        assert_eq!(work, -1.0);
    }

    #[test]
    fn nonfinite_input_is_rejected() {
        let mut rho = diag(&[1.0, 0.0]);
        rho[(0, 0)] = C64::new(f64::NAN, 0.0);
        assert!(build_passive_extraction(&rho, &diag(&[0.0, 1.0])).is_err());
    }
}
