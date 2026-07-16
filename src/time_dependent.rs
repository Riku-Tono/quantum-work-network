//! Fixed-step RK4 propagation for time-dependent Lindblad equations.
//!
//! Unlike [`crate::propagator::DenseExponentialPropagator`], this module
//! updates density matrices directly. It therefore does not use or alter the
//! crate's column-major Liouville-space vectorization convention.

use crate::error::PhysicsError;
use crate::matrix::{ComplexMatrix, C64};
use crate::propagator::QuantumState;

const TIME_EPSILON_FACTOR: f64 = 64.0;

/// Controls which states are retained by the propagator.
#[derive(Debug, Clone, PartialEq)]
pub enum SaveSchedule {
    /// Retain the initial state and every accepted RK4 step.
    EveryStep,
    /// Retain the initial state, points separated by this interval, and the
    /// exact final state. Short RK4 steps are used to land on save points.
    Interval(f64),
    /// Retain the requested in-range times, plus `t0` and `t_end` if absent.
    /// Input times may be unsorted but must be finite and within the interval.
    Times(Vec<f64>),
}

/// Accuracy-first fixed-maximum-step RK4 propagator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimeDependentRk4 {
    /// Maximum step size. A shorter final or save-point step is allowed.
    pub dt: f64,
}

impl TimeDependentRk4 {
    pub fn new(dt: f64) -> Result<Self, PhysicsError> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(PhysicsError::InvalidTime(format!(
                "RK4 dt must be finite and positive, got {dt}"
            )));
        }
        Ok(Self { dt })
    }

    /// Propagate a density matrix under a time-dependent Lindblad equation.
    ///
    /// Collapse operators must include their rate coefficients. No
    /// renormalization, Hermitian symmetrization, eigenvalue clipping, or other
    /// physicality correction is performed.
    pub fn propagate<H, C>(
        &self,
        rho0: &ComplexMatrix,
        t0: f64,
        t_end: f64,
        hamiltonian: H,
        collapse_operators: C,
        save_schedule: SaveSchedule,
    ) -> Result<Vec<QuantumState>, PhysicsError>
    where
        H: Fn(f64) -> ComplexMatrix,
        C: Fn(f64) -> Vec<ComplexMatrix>,
    {
        validate_initial_problem(rho0, t0, t_end)?;
        let save_times = build_save_times(t0, t_end, &save_schedule)?;
        let mut rho = rho0.clone();
        let mut time = t0;
        let mut states = Vec::new();
        states.push(QuantumState {
            rho: rho.clone(),
            time,
        });

        match save_times {
            None => {
                while time < t_end {
                    let step = (t_end - time).min(self.dt);
                    rho = rk4_step(&rho, time, step, &hamiltonian, &collapse_operators)?;
                    time += step;
                    if nearly_equal(time, t_end) {
                        time = t_end;
                    }
                    states.push(QuantumState {
                        rho: rho.clone(),
                        time,
                    });
                }
            }
            Some(targets) => {
                for &target in targets.iter().skip(1) {
                    while time < target && !nearly_equal(time, target) {
                        let step = (target - time).min(self.dt);
                        rho = rk4_step(&rho, time, step, &hamiltonian, &collapse_operators)?;
                        time += step;
                    }
                    time = target;
                    states.push(QuantumState {
                        rho: rho.clone(),
                        time,
                    });
                }
            }
        }
        Ok(states)
    }
}

/// Evaluate the time-dependent Lindblad right-hand side directly on `rho`.
pub fn lindblad_rhs(
    rho: &ComplexMatrix,
    hamiltonian: &ComplexMatrix,
    collapse_operators: &[ComplexMatrix],
) -> Result<ComplexMatrix, PhysicsError> {
    if rho.nrows() == 0 || rho.nrows() != rho.ncols() {
        return Err(PhysicsError::DimensionMismatch(
            "density matrix must be nonempty and square".to_string(),
        ));
    }
    let dim = rho.nrows();
    validate_operator(hamiltonian, dim, "Hamiltonian")?;

    let mut derivative = (hamiltonian * rho - rho * hamiltonian) * C64::new(0.0, -1.0);
    for collapse in collapse_operators {
        validate_operator(collapse, dim, "collapse operator")?;
        let collapse_dagger = collapse.adjoint();
        let l_dagger_l = &collapse_dagger * collapse;
        derivative += collapse * rho * &collapse_dagger
            - (&l_dagger_l * rho + rho * &l_dagger_l) * C64::new(0.5, 0.0);
    }
    Ok(derivative)
}

fn rk4_step<H, C>(
    rho: &ComplexMatrix,
    time: f64,
    dt: f64,
    hamiltonian: &H,
    collapse_operators: &C,
) -> Result<ComplexMatrix, PhysicsError>
where
    H: Fn(f64) -> ComplexMatrix,
    C: Fn(f64) -> Vec<ComplexMatrix>,
{
    let half_dt = C64::new(0.5 * dt, 0.0);
    let full_dt = C64::new(dt, 0.0);

    let h1 = hamiltonian(time);
    let c1 = collapse_operators(time);
    let k1 = lindblad_rhs(rho, &h1, &c1)?;

    let rho2 = rho + &k1 * half_dt;
    let h2 = hamiltonian(time + 0.5 * dt);
    let c2 = collapse_operators(time + 0.5 * dt);
    let k2 = lindblad_rhs(&rho2, &h2, &c2)?;

    let rho3 = rho + &k2 * half_dt;
    let h3 = hamiltonian(time + 0.5 * dt);
    let c3 = collapse_operators(time + 0.5 * dt);
    let k3 = lindblad_rhs(&rho3, &h3, &c3)?;

    let rho4 = rho + &k3 * full_dt;
    let h4 = hamiltonian(time + dt);
    let c4 = collapse_operators(time + dt);
    let k4 = lindblad_rhs(&rho4, &h4, &c4)?;

    Ok(rho
        + (k1 + k2 * C64::new(2.0, 0.0) + k3 * C64::new(2.0, 0.0) + k4) * C64::new(dt / 6.0, 0.0))
}

fn validate_initial_problem(rho0: &ComplexMatrix, t0: f64, t_end: f64) -> Result<(), PhysicsError> {
    if rho0.nrows() == 0 || rho0.nrows() != rho0.ncols() {
        return Err(PhysicsError::DimensionMismatch(
            "initial density matrix must be nonempty and square".to_string(),
        ));
    }
    if !t0.is_finite() || !t_end.is_finite() || t_end < t0 {
        return Err(PhysicsError::InvalidTime(format!(
            "require finite t0 and t_end with t_end >= t0, got {t0}..{t_end}"
        )));
    }
    if rho0
        .iter()
        .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return Err(PhysicsError::InvalidParameter(
            "initial density matrix contains a non-finite value".to_string(),
        ));
    }
    Ok(())
}

fn validate_operator(operator: &ComplexMatrix, dim: usize, name: &str) -> Result<(), PhysicsError> {
    if operator.shape() != (dim, dim) {
        return Err(PhysicsError::DimensionMismatch(format!(
            "{name} is {:?}, expected ({dim}, {dim})",
            operator.shape()
        )));
    }
    if operator
        .iter()
        .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return Err(PhysicsError::InvalidParameter(format!(
            "{name} contains a non-finite value"
        )));
    }
    Ok(())
}

fn build_save_times(
    t0: f64,
    t_end: f64,
    schedule: &SaveSchedule,
) -> Result<Option<Vec<f64>>, PhysicsError> {
    if matches!(schedule, SaveSchedule::EveryStep) {
        return Ok(None);
    }
    let mut times = match schedule {
        SaveSchedule::Interval(interval) => {
            if !interval.is_finite() || *interval <= 0.0 {
                return Err(PhysicsError::InvalidTime(format!(
                    "save interval must be finite and positive, got {interval}"
                )));
            }
            let mut values = vec![t0];
            let mut index = 1usize;
            loop {
                let candidate = t0 + index as f64 * interval;
                if candidate >= t_end || nearly_equal(candidate, t_end) {
                    break;
                }
                values.push(candidate);
                index += 1;
            }
            values.push(t_end);
            values
        }
        SaveSchedule::Times(requested) => {
            if requested.iter().any(|time| {
                !time.is_finite()
                    || *time < t0 && !nearly_equal(*time, t0)
                    || *time > t_end && !nearly_equal(*time, t_end)
            }) {
                return Err(PhysicsError::InvalidTime(format!(
                    "save times must be finite and inside [{t0}, {t_end}]"
                )));
            }
            let mut values = requested.clone();
            values.push(t0);
            values.push(t_end);
            values.sort_by(f64::total_cmp);
            values.dedup_by(|a, b| nearly_equal(*a, *b));
            if let Some(first) = values.first_mut() {
                *first = t0;
            }
            if let Some(last) = values.last_mut() {
                *last = t_end;
            }
            values
        }
        SaveSchedule::EveryStep => unreachable!(),
    };
    times.dedup_by(|a, b| nearly_equal(*a, *b));
    Ok(Some(times))
}

fn nearly_equal(a: f64, b: f64) -> bool {
    let scale = a.abs().max(b.abs()).max(1.0);
    (a - b).abs() <= TIME_EPSILON_FACTOR * f64::EPSILON * scale
}

#[cfg(test)]
mod tests {
    use nalgebra::linalg::Schur;

    use super::*;
    use crate::liouvillian::build_liouvillian;
    use crate::matrix::{frobenius_norm, hermiticity_error};
    use crate::operators::{build_operators, ModelParams};
    use crate::propagator::DenseExponentialPropagator;

    fn ground() -> ComplexMatrix {
        ComplexMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
        ]))
    }

    fn sigma_x() -> ComplexMatrix {
        ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(0.0, 0.0),
                C64::new(1.0, 0.0),
                C64::new(1.0, 0.0),
                C64::new(0.0, 0.0),
            ],
        )
    }

    fn amplitude_damping(gamma: f64) -> ComplexMatrix {
        ComplexMatrix::from_row_slice(
            2,
            2,
            &[
                C64::new(0.0, 0.0),
                C64::new(gamma.sqrt(), 0.0),
                C64::new(0.0, 0.0),
                C64::new(0.0, 0.0),
            ],
        )
    }

    fn minimum_eigenvalue(matrix: &ComplexMatrix) -> f64 {
        let (_, schur) = Schur::new(matrix.clone()).unpack();
        (0..schur.nrows())
            .map(|index| schur[(index, index)].re)
            .fold(f64::INFINITY, f64::min)
    }

    fn constant_reference_error(dt: f64) -> f64 {
        let h = sigma_x() * C64::new(0.7, 0.0);
        let collapses = vec![amplitude_damping(0.3)];
        let liouvillian = build_liouvillian(&h, &collapses).unwrap();
        let reference = DenseExponentialPropagator::new(liouvillian, 2)
            .unwrap()
            .propagate_times(&ground(), &[0.7])
            .unwrap()
            .pop()
            .unwrap()
            .rho;
        let actual = TimeDependentRk4::new(dt)
            .unwrap()
            .propagate(
                &ground(),
                0.0,
                0.7,
                |_| h.clone(),
                |_| collapses.clone(),
                SaveSchedule::Times(vec![0.7]),
            )
            .unwrap()
            .pop()
            .unwrap()
            .rho;
        frobenius_norm(&(actual - reference))
    }

    #[test]
    fn constant_generator_matches_dense_exponential_and_converges() {
        let coarse = constant_reference_error(0.08);
        let medium = constant_reference_error(0.04);
        let fine = constant_reference_error(0.02);
        assert!(medium < coarse, "{medium:e} !< {coarse:e}");
        assert!(fine < medium, "{fine:e} !< {medium:e}");
        assert!(fine < 1.0e-7, "fine-grid error was {fine:e}");
    }

    #[test]
    fn zero_time_dependent_term_matches_existing_network_hamiltonian() {
        let params = ModelParams {
            load_dim: 2,
            ..ModelParams::default()
        };
        let operators = build_operators(&params).unwrap();
        let dim = operators.h_total.nrows();
        let mut rho0 = ComplexMatrix::zeros(dim, dim);
        rho0[(8, 8)] = C64::new(1.0, 0.0);
        let liouvillian = build_liouvillian(&operators.h_total, &[]).unwrap();
        let reference = DenseExponentialPropagator::new(liouvillian, dim)
            .unwrap()
            .propagate_times(&rho0, &[0.05])
            .unwrap()
            .pop()
            .unwrap()
            .rho;
        let actual = TimeDependentRk4::new(0.0025)
            .unwrap()
            .propagate(
                &rho0,
                0.0,
                0.05,
                |_| operators.h_total.clone(),
                |_| Vec::new(),
                SaveSchedule::Times(vec![0.05]),
            )
            .unwrap()
            .pop()
            .unwrap()
            .rho;
        assert!(frobenius_norm(&(actual - reference)) < 1.0e-8);
    }

    #[test]
    fn trace_hermiticity_positivity_and_finiteness_hold_at_saved_times() {
        let h = sigma_x() * C64::new(0.4, 0.0);
        let collapse = amplitude_damping(0.2);
        let states = TimeDependentRk4::new(0.005)
            .unwrap()
            .propagate(
                &ground(),
                0.0,
                1.0,
                |_| h.clone(),
                |_| vec![collapse.clone()],
                SaveSchedule::Interval(0.1),
            )
            .unwrap();
        for state in states {
            assert!((state.rho.trace() - C64::new(1.0, 0.0)).norm() < 1.0e-12);
            assert!(hermiticity_error(&state.rho) < 1.0e-12);
            assert!(minimum_eigenvalue(&state.rho) >= -1.0e-11);
            assert!(state
                .rho
                .iter()
                .all(|value| value.re.is_finite() && value.im.is_finite()));
        }
    }

    #[test]
    fn positivity_defect_does_not_worsen_when_step_is_halved() {
        let h = sigma_x() * C64::new(1.1, 0.0);
        let worst = |dt: f64| {
            TimeDependentRk4::new(dt)
                .unwrap()
                .propagate(
                    &ground(),
                    0.0,
                    1.0,
                    |_| h.clone(),
                    |_| Vec::new(),
                    SaveSchedule::EveryStep,
                )
                .unwrap()
                .iter()
                .map(|state| minimum_eigenvalue(&state.rho))
                .fold(f64::INFINITY, f64::min)
        };
        let coarse = worst(0.2);
        let fine = worst(0.1);
        assert!(coarse >= -1.0e-4, "coarse minimum eigenvalue {coarse:e}");
        assert!(fine >= coarse - 1.0e-12, "fine={fine:e}, coarse={coarse:e}");
    }

    #[test]
    fn t_zero_is_returned_without_modification() {
        let rho0 = ground();
        let states = TimeDependentRk4::new(0.1)
            .unwrap()
            .propagate(
                &rho0,
                2.0,
                2.0,
                |_| sigma_x(),
                |_| Vec::new(),
                SaveSchedule::EveryStep,
            )
            .unwrap();
        assert_eq!(states.len(), 1);
        assert_eq!(states[0].time, 2.0);
        assert_eq!(states[0].rho, rho0);
    }

    #[test]
    fn final_short_step_lands_exactly_on_t_end_without_overshoot() {
        let states = TimeDependentRk4::new(0.3)
            .unwrap()
            .propagate(
                &ground(),
                0.0,
                1.0,
                |_| ComplexMatrix::zeros(2, 2),
                |_| Vec::new(),
                SaveSchedule::EveryStep,
            )
            .unwrap();
        assert_eq!(states.last().unwrap().time, 1.0);
        assert!(states.iter().all(|state| state.time <= 1.0));
        assert_eq!(states.len(), 5);
    }

    #[test]
    fn requested_save_times_are_reached_exactly() {
        let states = TimeDependentRk4::new(0.2)
            .unwrap()
            .propagate(
                &ground(),
                0.0,
                1.0,
                |_| ComplexMatrix::zeros(2, 2),
                |_| Vec::new(),
                SaveSchedule::Times(vec![0.73, 0.11]),
            )
            .unwrap();
        let times: Vec<f64> = states.iter().map(|state| state.time).collect();
        assert_eq!(times, vec![0.0, 0.11, 0.73, 1.0]);
    }

    #[test]
    fn invalid_dimensions_are_reported() {
        let result = TimeDependentRk4::new(0.1).unwrap().propagate(
            &ground(),
            0.0,
            0.1,
            |_| ComplexMatrix::zeros(3, 3),
            |_| Vec::new(),
            SaveSchedule::EveryStep,
        );
        assert!(matches!(result, Err(PhysicsError::DimensionMismatch(_))));
    }
}
