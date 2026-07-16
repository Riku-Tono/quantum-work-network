use crate::error::PhysicsError;
use crate::matrix::{commutator, eye, kron, ComplexMatrix, C64};

#[derive(Debug, Clone)]
pub struct ModelParams {
    pub omega_chain: f64,
    pub omega_load: f64,
    pub hopping_j: f64,
    pub coupling_g: f64,
    pub load_dim: usize,
}

impl Default for ModelParams {
    fn default() -> Self {
        Self {
            omega_chain: 1.0,
            omega_load: 1.0,
            hopping_j: 1.0,
            coupling_g: 0.25,
            load_dim: 3,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Operators {
    pub dims: Vec<usize>,
    pub h_total: ComplexMatrix,
    pub h_chain: ComplexMatrix,
    pub h_load: ComplexMatrix,
    pub h_site_3: ComplexMatrix,
    pub h_interaction: ComplexMatrix,
    pub sigma_1_plus: ComplexMatrix,
    pub sigma_z_sites: Vec<ComplexMatrix>,
    pub number_sites: Vec<ComplexMatrix>,
    pub number_chain: ComplexMatrix,
    pub number_load: ComplexMatrix,
    pub number_total: ComplexMatrix,
    pub b_load_local: ComplexMatrix,
    pub identity_full: ComplexMatrix,
}

fn qubit_identity() -> ComplexMatrix {
    eye(2)
}

fn sigma_plus() -> ComplexMatrix {
    // |1><0|, with |0> empty and |1> excited.
    ComplexMatrix::from_row_slice(
        2,
        2,
        &[
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
        ],
    )
}

fn sigma_minus() -> ComplexMatrix {
    sigma_plus().adjoint()
}

fn sigma_z() -> ComplexMatrix {
    ComplexMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
        C64::new(1.0, 0.0),
        C64::new(-1.0, 0.0),
    ]))
}

fn number_qubit() -> ComplexMatrix {
    ComplexMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
        C64::new(0.0, 0.0),
        C64::new(1.0, 0.0),
    ]))
}

pub fn truncated_annihilation(dim: usize) -> Result<ComplexMatrix, PhysicsError> {
    if dim < 2 {
        return Err(PhysicsError::DimensionMismatch(
            "load dimension must be at least 2".to_string(),
        ));
    }
    let mut b = ComplexMatrix::zeros(dim, dim);
    for n in 1..dim {
        b[(n - 1, n)] = C64::new((n as f64).sqrt(), 0.0);
    }
    Ok(b)
}

fn embed_local(
    local: &ComplexMatrix,
    subsystem: usize,
    dims: &[usize],
) -> Result<ComplexMatrix, PhysicsError> {
    if subsystem >= dims.len() {
        return Err(PhysicsError::InvalidSubsystem(format!(
            "subsystem {subsystem} outside 0..{}",
            dims.len()
        )));
    }
    if local.nrows() != dims[subsystem] || local.ncols() != dims[subsystem] {
        return Err(PhysicsError::DimensionMismatch(format!(
            "local operator is {}x{}, expected {}x{}",
            local.nrows(),
            local.ncols(),
            dims[subsystem],
            dims[subsystem]
        )));
    }

    let mut out = ComplexMatrix::from_element(1, 1, C64::new(1.0, 0.0));
    for (index, &dim) in dims.iter().enumerate() {
        let factor = if index == subsystem {
            local.clone()
        } else {
            eye(dim)
        };
        out = kron(&out, &factor);
    }
    Ok(out)
}

pub fn build_operators(params: &ModelParams) -> Result<Operators, PhysicsError> {
    build_operators_for_chain(params, 3)
}

/// Build the same nearest-neighbor chain + terminal load model for an arbitrary
/// positive number of two-level chain sites. The historical `build_operators`
/// entry point remains an exact N=3 wrapper.
pub fn build_operators_for_chain(
    params: &ModelParams,
    chain_length: usize,
) -> Result<Operators, PhysicsError> {
    if params.load_dim < 2 {
        return Err(PhysicsError::DimensionMismatch(
            "load_dim must be >= 2".to_string(),
        ));
    }
    if chain_length == 0 {
        return Err(PhysicsError::DimensionMismatch(
            "chain_length must be >= 1".to_string(),
        ));
    }

    // Tensor order: |q1, ..., qN, load>, with load index varying fastest.
    let mut dims = vec![2; chain_length];
    dims.push(params.load_dim);
    let full_dim: usize = dims.iter().product();
    let identity_full = eye(full_dim);

    let sp = sigma_plus();
    let sm = sigma_minus();
    let sz = sigma_z();
    let n = number_qubit();

    let mut sigma_plus_sites = Vec::with_capacity(chain_length);
    let mut sigma_minus_sites = Vec::with_capacity(chain_length);
    let mut sigma_z_sites = Vec::with_capacity(chain_length);
    let mut number_sites = Vec::with_capacity(chain_length);

    for site in 0..chain_length {
        sigma_plus_sites.push(embed_local(&sp, site, &dims)?);
        sigma_minus_sites.push(embed_local(&sm, site, &dims)?);
        sigma_z_sites.push(embed_local(&sz, site, &dims)?);
        number_sites.push(embed_local(&n, site, &dims)?);
    }

    let b_local = truncated_annihilation(params.load_dim)?;
    let b = embed_local(&b_local, chain_length, &dims)?;
    let b_dag = b.adjoint();
    let number_load = &b_dag * &b;

    let number_chain = number_sites
        .iter()
        .fold(ComplexMatrix::zeros(full_dim, full_dim), |acc, op| acc + op);
    let number_total = &number_chain + &number_load;

    let h_onsite = &number_chain * C64::new(params.omega_chain, 0.0);
    let mut h_hopping = ComplexMatrix::zeros(full_dim, full_dim);
    for site in 0..chain_length - 1 {
        h_hopping += (&sigma_plus_sites[site] * &sigma_minus_sites[site + 1]
            + &sigma_minus_sites[site] * &sigma_plus_sites[site + 1])
            * C64::new(params.hopping_j, 0.0);
    }
    let h_chain = h_onsite + h_hopping;
    let h_load = &number_load * C64::new(params.omega_load, 0.0);
    let final_site = chain_length - 1;
    let h_interaction = (&sigma_plus_sites[final_site] * &b
        + &sigma_minus_sites[final_site] * &b_dag)
        * C64::new(params.coupling_g, 0.0);
    let h_total = &h_chain + &h_load + &h_interaction;
    // Historical field name retained for compatibility; for generalized chains
    // it contains the local Hamiltonian of the final chain site.
    let h_site_3 = &number_sites[final_site] * C64::new(params.omega_chain, 0.0);

    Ok(Operators {
        dims,
        h_total,
        h_chain,
        h_load,
        h_site_3,
        h_interaction,
        sigma_1_plus: sigma_plus_sites[0].clone(),
        sigma_z_sites,
        number_sites,
        number_chain,
        number_load,
        number_total,
        b_load_local: b_local,
        identity_full,
    })
}

pub fn excitation_conservation_error(operators: &Operators) -> f64 {
    crate::matrix::frobenius_norm(&commutator(&operators.number_total, &operators.h_total))
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;
    use crate::matrix::{frobenius_norm, hermiticity_error};

    #[test]
    fn operators_are_hermitian_and_dimensions_match() {
        let ops = build_operators(&ModelParams::default()).unwrap();
        assert_eq!(ops.h_total.shape(), (24, 24));
        assert!(hermiticity_error(&ops.h_total) < 1.0e-12);
        assert!(hermiticity_error(&ops.h_chain) < 1.0e-12);
        assert!(hermiticity_error(&ops.h_load) < 1.0e-12);
        assert!(hermiticity_error(&ops.h_interaction) < 1.0e-12);
    }

    #[test]
    fn five_site_operators_have_requested_terminal_mapping() {
        let ops = build_operators_for_chain(&ModelParams::default(), 5).unwrap();
        assert_eq!(ops.dims, vec![2, 2, 2, 2, 2, 3]);
        assert_eq!(ops.h_total.shape(), (96, 96));
        assert_eq!(ops.number_sites.len(), 5);
        assert_eq!(ops.sigma_z_sites.len(), 5);
        assert!(hermiticity_error(&ops.h_total) < 1.0e-12);
        assert!(excitation_conservation_error(&ops) < 1.0e-12);
    }

    #[test]
    fn onsite_number_operators_commute() {
        let ops = build_operators(&ModelParams::default()).unwrap();
        let h0 = &ops.number_chain * C64::new(1.0, 0.0);
        for number in &ops.number_sites {
            assert!(frobenius_norm(&commutator(number, &h0)) < 1.0e-12);
        }
    }

    #[test]
    fn closed_hamiltonian_conserves_total_excitation() {
        let ops = build_operators(&ModelParams::default()).unwrap();
        assert!(excitation_conservation_error(&ops) < 1.0e-12);
    }

    #[test]
    fn truncated_bosonic_commutator_matches_specification() {
        let b = truncated_annihilation(3).unwrap();
        let comm = &b * b.adjoint() - b.adjoint() * &b;
        let mut expected = eye(3);
        expected[(2, 2)] -= C64::new(3.0, 0.0);
        assert!(frobenius_norm(&(comm - expected)) < 1.0e-12);
        assert_abs_diff_eq!(b[(0, 1)].re, 1.0, epsilon = 1.0e-12);
        assert_abs_diff_eq!(b[(1, 2)].re, 2.0_f64.sqrt(), epsilon = 1.0e-12);
    }
}
