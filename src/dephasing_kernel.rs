//! Exact elementwise kernel for computational-basis local sigma-z dephasing.

use crate::error::PhysicsError;
use crate::matrix::{ComplexMatrix, C64};

#[derive(Debug, Clone)]
pub struct DiagonalDephasingKernel {
    chain_length: usize,
    load_dim: usize,
    dim: usize,
    // Column-major, matching nalgebra::DMatrix storage.
    rates: Vec<f64>,
}

impl DiagonalDephasingKernel {
    pub fn new(
        chain_length: usize,
        load_dim: usize,
        site_gammas: &[f64],
    ) -> Result<Self, PhysicsError> {
        if chain_length == 0 || load_dim == 0 || site_gammas.len() != chain_length {
            return Err(PhysicsError::DimensionMismatch(format!(
                "expected {chain_length} site rates and positive dimensions, got {} rates and load_dim={load_dim}",
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
        let chain_dim = 1usize.checked_shl(chain_length as u32).ok_or_else(|| {
            PhysicsError::DimensionMismatch("chain dimension overflow".to_string())
        })?;
        let dim = chain_dim.checked_mul(load_dim).ok_or_else(|| {
            PhysicsError::DimensionMismatch("full dimension overflow".to_string())
        })?;
        let element_count = dim.checked_mul(dim).ok_or_else(|| {
            PhysicsError::DimensionMismatch("kernel element count overflow".to_string())
        })?;
        let mut rates = vec![0.0; element_count];
        for col in 0..dim {
            let col_chain = col / load_dim;
            for row in 0..dim {
                let row_chain = row / load_dim;
                let differing = row_chain ^ col_chain;
                let mut rate = 0.0;
                for (site, gamma) in site_gammas.iter().enumerate() {
                    // Tensor order is |q1,...,qN,load>; q1 is the most-significant chain bit.
                    let shift = chain_length - 1 - site;
                    if differing & (1usize << shift) != 0 {
                        rate += gamma;
                    }
                }
                rates[row + col * dim] = rate;
            }
        }
        Ok(Self {
            chain_length,
            load_dim,
            dim,
            rates,
        })
    }

    pub fn dimension(&self) -> usize {
        self.dim
    }

    pub fn chain_length(&self) -> usize {
        self.chain_length
    }

    pub fn load_dim(&self) -> usize {
        self.load_dim
    }

    pub fn rate(&self, row: usize, col: usize) -> Result<f64, PhysicsError> {
        if row >= self.dim || col >= self.dim {
            return Err(PhysicsError::DimensionMismatch(format!(
                "kernel index ({row},{col}) outside {}x{}",
                self.dim, self.dim
            )));
        }
        Ok(self.rates[row + col * self.dim])
    }

    pub fn apply(&self, rho: &ComplexMatrix) -> Result<ComplexMatrix, PhysicsError> {
        self.validate_matrix(rho, "density matrix")?;
        let mut out = ComplexMatrix::zeros(self.dim, self.dim);
        for ((target, value), rate) in out
            .as_mut_slice()
            .iter_mut()
            .zip(rho.as_slice())
            .zip(&self.rates)
        {
            *target = *value * C64::new(-*rate, 0.0);
        }
        Ok(out)
    }

    pub fn add_to(
        &self,
        rho: &ComplexMatrix,
        derivative: &mut ComplexMatrix,
    ) -> Result<(), PhysicsError> {
        self.validate_matrix(rho, "density matrix")?;
        self.validate_matrix(derivative, "derivative")?;
        for ((target, value), rate) in derivative
            .as_mut_slice()
            .iter_mut()
            .zip(rho.as_slice())
            .zip(&self.rates)
        {
            *target -= *value * C64::new(*rate, 0.0);
        }
        Ok(())
    }

    /// Returns the instantaneous dephasing-kernel-weighted coherence exposure.
    ///
    /// This is `sum_ab Gamma[a,b] * |rho[a,b]|^2`. It is a trajectory
    /// diagnostic, not an energy, work, heat, or entropy-production rate.
    pub fn weighted_coherence_exposure_rate(
        &self,
        rho: &ComplexMatrix,
    ) -> Result<f64, PhysicsError> {
        self.validate_matrix(rho, "density matrix")?;
        if rho
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
        {
            return Err(PhysicsError::InvalidParameter(
                "density matrix contains non-finite values".to_string(),
            ));
        }
        let exposure = self
            .rates
            .iter()
            .zip(rho.as_slice())
            .map(|(rate, value)| rate * value.norm_sqr())
            .sum::<f64>();
        if !exposure.is_finite() {
            return Err(PhysicsError::InvalidParameter(
                "weighted coherence exposure is non-finite".to_string(),
            ));
        }
        Ok(exposure)
    }

    pub fn estimated_bytes(&self) -> usize {
        self.rates.len() * std::mem::size_of::<f64>()
    }

    fn validate_matrix(&self, matrix: &ComplexMatrix, name: &str) -> Result<(), PhysicsError> {
        if matrix.shape() != (self.dim, self.dim) {
            return Err(PhysicsError::DimensionMismatch(format!(
                "{name} is {:?}, expected ({},{})",
                matrix.shape(),
                self.dim,
                self.dim
            )));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::lindblad_action;
    use crate::matrix::{frobenius_norm, hermiticity_error};
    use crate::operators::{build_operators_for_chain, ModelParams};

    fn fixed_hermitian(dim: usize) -> ComplexMatrix {
        let mut out = ComplexMatrix::zeros(dim, dim);
        for col in 0..dim {
            for row in 0..=col {
                let value = C64::new(
                    ((row + 2 * col + 1) as f64).sin() / dim as f64,
                    if row == col {
                        0.0
                    } else {
                        ((3 * row + col + 2) as f64).cos() / dim as f64
                    },
                );
                out[(row, col)] = value;
                out[(col, row)] = value.conj();
            }
        }
        out
    }

    fn dense(
        ops: &crate::operators::Operators,
        gammas: &[f64],
        rho: &ComplexMatrix,
    ) -> ComplexMatrix {
        let mut out = ComplexMatrix::zeros(rho.nrows(), rho.ncols());
        for (z, gamma) in ops.sigma_z_sites.iter().zip(gammas) {
            if *gamma > 0.0 {
                let collapse = z * C64::new((*gamma / 2.0).sqrt(), 0.0);
                out += lindblad_action(&collapse, rho).unwrap();
            }
        }
        out
    }

    #[test]
    fn hamming_rates_exclude_load_and_follow_tensor_order() {
        let k = DiagonalDephasingKernel::new(3, 3, &[0.2, 0.3, 0.5]).unwrap();
        let index = |chain: usize, load: usize| chain * 3 + load;
        assert_eq!(k.rate(index(0b000, 0), index(0b000, 2)).unwrap(), 0.0);
        assert_eq!(k.rate(index(0b000, 0), index(0b100, 0)).unwrap(), 0.2);
        assert_eq!(k.rate(index(0b000, 0), index(0b010, 0)).unwrap(), 0.3);
        assert_eq!(k.rate(index(0b000, 0), index(0b001, 0)).unwrap(), 0.5);
        assert_eq!(k.rate(index(0b101, 1), index(0b010, 2)).unwrap(), 1.0);
    }

    #[test]
    fn basis_mapping_matches_embedded_sigma_z_diagonals() {
        let ops = build_operators_for_chain(&ModelParams::default(), 3).unwrap();
        let k = DiagonalDephasingKernel::new(3, 3, &[0.2, 0.3, 0.5]).unwrap();
        for a in 0..24 {
            for b in 0..24 {
                let expected: f64 = ops
                    .sigma_z_sites
                    .iter()
                    .zip([0.2, 0.3, 0.5])
                    .map(|(z, gamma)| {
                        if z[(a, a)].re != z[(b, b)].re {
                            gamma
                        } else {
                            0.0
                        }
                    })
                    .sum();
                assert!((k.rate(a, b).unwrap() - expected).abs() < 1.0e-15);
            }
        }
    }

    #[test]
    fn kernel_matches_dense_for_n1_and_n3_configurations() {
        for (n, gammas) in [
            (1usize, vec![0.4]),
            (3, vec![0.2, 0.0, 0.7]),
            (3, vec![0.5; 3]),
        ] {
            let ops = build_operators_for_chain(&ModelParams::default(), n).unwrap();
            let rho = fixed_hermitian(ops.h_total.nrows());
            let kernel = DiagonalDephasingKernel::new(n, 3, &gammas).unwrap();
            let difference = kernel.apply(&rho).unwrap() - dense(&ops, &gammas, &rho);
            assert!(difference.iter().map(|z| z.norm()).fold(0.0, f64::max) <= 2.0e-15);
            assert!(frobenius_norm(&difference) <= 2.0e-14);
        }
    }

    #[test]
    fn structure_zero_gamma_trace_diagonal_and_hermiticity() {
        let rho = fixed_hermitian(24);
        let zero = DiagonalDephasingKernel::new(3, 3, &[0.0; 3])
            .unwrap()
            .apply(&rho)
            .unwrap();
        assert!(zero.iter().all(|z| *z == C64::new(0.0, 0.0)));
        let out = DiagonalDephasingKernel::new(3, 3, &[0.2, 0.3, 0.5])
            .unwrap()
            .apply(&rho)
            .unwrap();
        assert!((0..24).all(|i| out[(i, i)] == C64::new(0.0, 0.0)));
        assert_eq!(out.trace(), C64::new(0.0, 0.0));
        assert!(hermiticity_error(&out) < 1.0e-14);
    }

    #[test]
    fn rejects_bad_rates_and_shapes() {
        assert!(DiagonalDephasingKernel::new(3, 3, &[0.1, -0.2, 0.3]).is_err());
        assert!(DiagonalDephasingKernel::new(3, 3, &[0.1, f64::NAN, 0.3]).is_err());
        assert!(DiagonalDephasingKernel::new(3, 3, &[0.1]).is_err());
        let k = DiagonalDephasingKernel::new(1, 3, &[0.1]).unwrap();
        assert!(k.apply(&ComplexMatrix::zeros(2, 2)).is_err());
    }

    #[test]
    fn weighted_exposure_is_zero_for_diagonal_state() {
        let k = DiagonalDephasingKernel::new(1, 2, &[0.7]).unwrap();
        let rho = ComplexMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
            C64::new(0.6, 0.0),
            C64::new(0.4, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
        ]));
        assert_eq!(k.weighted_coherence_exposure_rate(&rho).unwrap(), 0.0);
    }

    #[test]
    fn weighted_exposure_counts_both_hermitian_off_diagonals() {
        let k = DiagonalDephasingKernel::new(1, 1, &[0.7]).unwrap();
        let mut rho = ComplexMatrix::zeros(2, 2);
        rho[(0, 1)] = C64::new(0.3, 0.4);
        rho[(1, 0)] = rho[(0, 1)].conj();
        let expected = k.rate(0, 1).unwrap() * rho[(0, 1)].norm_sqr()
            + k.rate(1, 0).unwrap() * rho[(1, 0)].norm_sqr();
        assert!((k.weighted_coherence_exposure_rate(&rho).unwrap() - expected).abs() < 1.0e-15);
    }

    #[test]
    fn weighted_exposure_is_zero_for_zero_gamma() {
        let k = DiagonalDephasingKernel::new(1, 1, &[0.0]).unwrap();
        let rho = fixed_hermitian(2);
        assert_eq!(k.weighted_coherence_exposure_rate(&rho).unwrap(), 0.0);
    }

    #[test]
    fn weighted_exposure_is_nonnegative_for_finite_hermitian_state() {
        let k = DiagonalDephasingKernel::new(2, 1, &[0.2, 0.5]).unwrap();
        let rho = fixed_hermitian(4);
        assert!(k.weighted_coherence_exposure_rate(&rho).unwrap() >= -1.0e-14);
    }

    #[test]
    fn weighted_exposure_rejects_dimension_mismatch_and_nonfinite_input() {
        let k = DiagonalDephasingKernel::new(1, 1, &[0.5]).unwrap();
        assert!(k
            .weighted_coherence_exposure_rate(&ComplexMatrix::zeros(3, 3))
            .is_err());
        let mut rho = ComplexMatrix::zeros(2, 2);
        rho[(0, 1)] = C64::new(f64::NAN, 0.0);
        assert!(k.weighted_coherence_exposure_rate(&rho).is_err());
    }

    #[test]
    fn weighted_exposure_scales_linearly_with_all_gammas() {
        let rho = fixed_hermitian(4);
        let base = DiagonalDephasingKernel::new(2, 1, &[0.2, 0.5])
            .unwrap()
            .weighted_coherence_exposure_rate(&rho)
            .unwrap();
        let scaled = DiagonalDephasingKernel::new(2, 1, &[0.6, 1.5])
            .unwrap()
            .weighted_coherence_exposure_rate(&rho)
            .unwrap();
        assert!((scaled - 3.0 * base).abs() <= 1.0e-14);
    }
}
