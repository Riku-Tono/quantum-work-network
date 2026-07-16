use std::collections::HashSet;

use crate::error::PhysicsError;
use crate::matrix::{ComplexMatrix, C64};

fn decode_index(mut index: usize, dims: &[usize]) -> Vec<usize> {
    let mut digits = vec![0; dims.len()];
    for position in (0..dims.len()).rev() {
        digits[position] = index % dims[position];
        index /= dims[position];
    }
    digits
}

fn encode_selected(digits: &[usize], dims: &[usize], keep: &[usize]) -> usize {
    let mut index = 0usize;
    for &subsystem in keep {
        index = index * dims[subsystem] + digits[subsystem];
    }
    index
}

pub fn partial_trace(
    rho: &ComplexMatrix,
    dims: &[usize],
    keep: &[usize],
) -> Result<ComplexMatrix, PhysicsError> {
    if dims.is_empty() {
        return Err(PhysicsError::InvalidSubsystem(
            "dims must not be empty".to_string(),
        ));
    }
    if keep.is_empty() {
        return Err(PhysicsError::InvalidSubsystem(
            "keep must contain at least one subsystem".to_string(),
        ));
    }

    let total_dim: usize = dims.iter().product();
    if rho.nrows() != total_dim || rho.ncols() != total_dim {
        return Err(PhysicsError::DimensionMismatch(format!(
            "rho is {}x{}, expected {}x{} from dims {:?}",
            rho.nrows(),
            rho.ncols(),
            total_dim,
            total_dim,
            dims
        )));
    }

    let mut seen = HashSet::new();
    for &subsystem in keep {
        if subsystem >= dims.len() {
            return Err(PhysicsError::InvalidSubsystem(format!(
                "kept subsystem {subsystem} outside 0..{}",
                dims.len()
            )));
        }
        if !seen.insert(subsystem) {
            return Err(PhysicsError::InvalidSubsystem(format!(
                "duplicate kept subsystem {subsystem}"
            )));
        }
    }

    let kept_dim: usize = keep.iter().map(|&i| dims[i]).product();
    let traced: Vec<usize> = (0..dims.len()).filter(|i| !seen.contains(i)).collect();
    let mut reduced = ComplexMatrix::zeros(kept_dim, kept_dim);

    let decoded: Vec<Vec<usize>> = (0..total_dim)
        .map(|index| decode_index(index, dims))
        .collect();

    for row in 0..total_dim {
        for col in 0..total_dim {
            if traced
                .iter()
                .all(|&subsystem| decoded[row][subsystem] == decoded[col][subsystem])
            {
                let reduced_row = encode_selected(&decoded[row], dims, keep);
                let reduced_col = encode_selected(&decoded[col], dims, keep);
                reduced[(reduced_row, reduced_col)] += rho[(row, col)];
            }
        }
    }

    Ok(reduced)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{frobenius_norm, hermiticity_error, kron, trace_real};

    fn ket_density(vector: &[C64]) -> ComplexMatrix {
        let ket = nalgebra::DVector::from_column_slice(vector);
        &ket * ket.adjoint()
    }

    #[test]
    fn product_state_reduces_to_requested_factor() {
        let rho_a = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
            C64::new(0.25, 0.0),
            C64::new(0.75, 0.0),
        ]));
        let rho_b = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
            C64::new(0.6, 0.0),
            C64::new(0.3, 0.0),
            C64::new(0.1, 0.0),
        ]));
        let rho = kron(&rho_a, &rho_b);

        let reduced_b = partial_trace(&rho, &[2, 3], &[1]).unwrap();
        assert!(frobenius_norm(&(reduced_b - rho_b)) < 1.0e-12);
    }

    #[test]
    fn bell_state_reduces_to_maximally_mixed_state() {
        let norm = 1.0 / 2.0_f64.sqrt();
        let rho = ket_density(&[
            C64::new(norm, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(norm, 0.0),
        ]);
        let reduced = partial_trace(&rho, &[2, 2], &[0]).unwrap();
        let expected = ComplexMatrix::identity(2, 2) * C64::new(0.5, 0.0);
        assert!(frobenius_norm(&(reduced - expected)) < 1.0e-12);
    }

    #[test]
    fn trace_and_hermiticity_are_preserved() {
        let norm = 1.0 / 2.0_f64.sqrt();
        let rho = ket_density(&[
            C64::new(norm, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, norm),
            C64::new(0.0, 0.0),
        ]);
        let reduced = partial_trace(&rho, &[2, 2], &[1]).unwrap();
        assert!((trace_real(&reduced) - 1.0).abs() < 1.0e-12);
        assert!(hermiticity_error(&reduced) < 1.0e-12);
    }

    #[test]
    fn keep_order_is_preserved_exactly() {
        // |q, l> = |1, 2>. Keeping [1,0] must produce |2,1> in dimensions [3,2].
        let mut ket = vec![C64::new(0.0, 0.0); 6];
        ket[1 * 3 + 2] = C64::new(1.0, 0.0);
        let rho = ket_density(&ket);
        let reordered = partial_trace(&rho, &[2, 3], &[1, 0]).unwrap();

        let expected_index = 2 * 2 + 1;
        assert!((reordered[(expected_index, expected_index)].re - 1.0).abs() < 1.0e-12);
        assert!((trace_real(&reordered) - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn model_subsystem_shapes_match_contract() {
        let rho = ComplexMatrix::identity(24, 24) * C64::new(1.0 / 24.0, 0.0);
        let load = partial_trace(&rho, &[2, 2, 2, 3], &[3]).unwrap();
        let site3_load = partial_trace(&rho, &[2, 2, 2, 3], &[2, 3]).unwrap();
        assert_eq!(load.shape(), (3, 3));
        assert_eq!(site3_load.shape(), (6, 6));
    }
}
