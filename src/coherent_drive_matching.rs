//! Bounded grid and root utilities for Milestone 5c.

use crate::error::PhysicsError;

pub const ROOT_MERGE_TOLERANCE: f64 = 1.0e-7;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridValue {
    pub omega: f64,
    pub residual: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RootCandidate {
    pub bracket_low: f64,
    pub bracket_high: f64,
    pub omega: f64,
    pub absolute_error: f64,
    pub relative_error: f64,
    pub iterations: usize,
    pub converged: bool,
}

#[derive(Debug, Clone)]
pub struct GridAnalysis {
    pub brackets: Vec<(f64, f64)>,
    pub matching_grid_points: Vec<GridValue>,
    pub best_grid_point: GridValue,
}

pub fn inclusive_grid(low: f64, high: f64, step: f64) -> Result<Vec<f64>, PhysicsError> {
    if !low.is_finite() || !high.is_finite() || !step.is_finite() || step <= 0.0 || high < low {
        return Err(PhysicsError::InvalidParameter(
            "grid requires finite low <= high and positive step".to_string(),
        ));
    }
    let count = ((high - low) / step).round() as usize;
    let mut values: Vec<f64> = (0..=count).map(|index| low + index as f64 * step).collect();
    if values
        .last()
        .is_none_or(|value| (*value - high).abs() > 1.0e-12)
    {
        values.push(high);
    } else if let Some(last) = values.last_mut() {
        *last = high;
    }
    Ok(values)
}

pub fn relative_error(residual: f64, target: f64) -> f64 {
    residual.abs() / target.abs().max(1.0e-12)
}

pub fn analyze_grid(
    values: &[GridValue],
    target: f64,
    tolerance: f64,
) -> Result<GridAnalysis, PhysicsError> {
    if values.is_empty() {
        return Err(PhysicsError::MatchingFailure("grid is empty".to_string()));
    }
    if values.windows(2).any(|pair| pair[0].omega >= pair[1].omega) {
        return Err(PhysicsError::MatchingFailure(
            "grid values must be strictly increasing".to_string(),
        ));
    }
    let mut brackets = Vec::new();
    for pair in values.windows(2) {
        if pair[0].residual * pair[1].residual < 0.0 {
            brackets.push((pair[0].omega, pair[1].omega));
        }
    }
    let matching_grid_points = values
        .iter()
        .copied()
        .filter(|value| relative_error(value.residual, target) < tolerance)
        .collect();
    let best_grid_point = values
        .iter()
        .copied()
        .min_by(|left, right| left.residual.abs().total_cmp(&right.residual.abs()))
        .unwrap();
    Ok(GridAnalysis {
        brackets,
        matching_grid_points,
        best_grid_point,
    })
}

pub fn refine_bisection<E>(
    low: f64,
    high: f64,
    target: f64,
    tolerance: f64,
    max_iterations: usize,
    mut residual: impl FnMut(f64) -> Result<f64, E>,
) -> Result<RootCandidate, E> {
    let mut left = low;
    let mut right = high;
    let mut f_left = residual(left)?;
    let f_right = residual(right)?;
    if relative_error(f_left, target) < tolerance {
        return Ok(candidate(left, right, left, f_left, target, 0, true));
    }
    if relative_error(f_right, target) < tolerance {
        return Ok(candidate(left, right, right, f_right, target, 0, true));
    }
    if f_left * f_right > 0.0 {
        let midpoint = 0.5 * (left + right);
        let value = residual(midpoint)?;
        return Ok(candidate(left, right, midpoint, value, target, 0, false));
    }
    let mut midpoint = 0.5 * (left + right);
    let mut value = residual(midpoint)?;
    for iteration in 1..=max_iterations {
        midpoint = 0.5 * (left + right);
        value = residual(midpoint)?;
        if relative_error(value, target) < tolerance {
            return Ok(candidate(
                low, high, midpoint, value, target, iteration, true,
            ));
        }
        if f_left * value <= 0.0 {
            right = midpoint;
        } else {
            left = midpoint;
            f_left = value;
        }
    }
    Ok(candidate(
        low,
        high,
        midpoint,
        value,
        target,
        max_iterations,
        false,
    ))
}

fn candidate(
    low: f64,
    high: f64,
    omega: f64,
    residual: f64,
    target: f64,
    iterations: usize,
    converged: bool,
) -> RootCandidate {
    RootCandidate {
        bracket_low: low,
        bracket_high: high,
        omega,
        absolute_error: residual.abs(),
        relative_error: relative_error(residual, target),
        iterations,
        converged,
    }
}

pub fn merge_duplicate_roots(mut roots: Vec<RootCandidate>, tolerance: f64) -> Vec<RootCandidate> {
    roots.sort_by(|left, right| left.omega.total_cmp(&right.omega));
    let mut merged: Vec<RootCandidate> = Vec::new();
    for root in roots {
        if let Some(previous) = merged.last_mut() {
            if (root.omega - previous.omega).abs() < tolerance {
                if root.absolute_error < previous.absolute_error {
                    *previous = root;
                }
                continue;
            }
        }
        merged.push(root);
    }
    merged
}

pub fn select_primary_root(roots: &[RootCandidate], omega_a: f64) -> Option<&RootCandidate> {
    roots
        .iter()
        .filter(|root| root.converged)
        .min_by(|left, right| {
            let left_distance = (left.omega - omega_a).abs();
            let right_distance = (right.omega - omega_a).abs();
            if (left_distance - right_distance).abs() <= 1.0e-14 {
                left.omega.total_cmp(&right.omega)
            } else {
                left_distance.total_cmp(&right_distance)
            }
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coherent_drive::{run_coherent_drive, CoherentDriveConfig};
    use crate::operators::ModelParams;

    #[test]
    fn grid_includes_endpoints_and_is_unique_and_ascending() {
        let grid = inclusive_grid(0.2, 1.0, 0.01).unwrap();
        assert_eq!(grid.first().copied(), Some(0.2));
        assert_eq!(grid.last().copied(), Some(1.0));
        assert_eq!(grid.len(), 81);
        assert!(grid.windows(2).all(|pair| pair[0] < pair[1]));
    }

    #[test]
    fn all_sign_change_intervals_are_found() {
        let values: Vec<_> = [0.0, 0.25, 0.5, 0.75, 1.0]
            .into_iter()
            .map(|omega| GridValue {
                omega,
                residual: (omega - 0.2) * (omega - 0.6) * (omega - 0.9),
            })
            .collect();
        let analysis = analyze_grid(&values, 1.0, 1.0e-4).unwrap();
        assert_eq!(
            analysis.brackets,
            vec![(0.0, 0.25), (0.5, 0.75), (0.75, 1.0)]
        );
    }

    #[test]
    fn multiple_roots_are_refined() {
        let roots: Vec<_> = [(0.0, 0.5), (0.5, 1.0)]
            .into_iter()
            .map(|(low, high)| {
                refine_bisection(low, high, 1.0, 1.0e-10, 80, |x| {
                    Ok::<_, ()>((x - 0.3) * (x - 0.7))
                })
                .unwrap()
            })
            .collect();
        assert_eq!(roots.len(), 2);
        assert!((roots[0].omega - 0.3).abs() < 1.0e-8);
        assert!((roots[1].omega - 0.7).abs() < 1.0e-8);
    }

    #[test]
    fn duplicate_roots_are_merged() {
        let make = |omega, error| RootCandidate {
            bracket_low: 0.0,
            bracket_high: 1.0,
            omega,
            absolute_error: error,
            relative_error: error,
            iterations: 1,
            converged: true,
        };
        let merged = merge_duplicate_roots(
            vec![make(0.4, 1.0e-5), make(0.4 + 1.0e-8, 1.0e-7)],
            ROOT_MERGE_TOLERANCE,
        );
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].absolute_error, 1.0e-7);
    }

    #[test]
    fn primary_root_uses_distance_then_smaller_omega() {
        let make = |omega| RootCandidate {
            bracket_low: omega,
            bracket_high: omega,
            omega,
            absolute_error: 0.0,
            relative_error: 0.0,
            iterations: 0,
            converged: true,
        };
        let roots = vec![make(0.1), make(0.3), make(0.8)];
        assert_eq!(select_primary_root(&roots, 0.2).unwrap().omega, 0.1);
    }

    #[test]
    fn no_match_keeps_best_grid_reference() {
        let values = vec![
            GridValue {
                omega: 0.2,
                residual: 2.0,
            },
            GridValue {
                omega: 1.0,
                residual: 1.0,
            },
        ];
        let analysis = analyze_grid(&values, 1.0, 1.0e-4).unwrap();
        assert!(analysis.brackets.is_empty());
        assert!(analysis.matching_grid_points.is_empty());
        assert_eq!(analysis.best_grid_point.omega, 1.0);
    }

    #[test]
    fn relative_error_protects_zero_target() {
        assert_eq!(relative_error(1.0e-12, 0.0), 1.0);
    }

    #[test]
    fn matched_style_runs_have_equal_times_and_physical_checks() {
        let mut a = CoherentDriveConfig::milestone_5b(0.0, 0.01);
        a.tau = 0.02;
        a.t_end = 0.03;
        a.save_interval = 0.01;
        let mut b = a;
        b.gamma_phi = 0.5;
        let left = run_coherent_drive(&ModelParams::default(), a).unwrap();
        let right = run_coherent_drive(&ModelParams::default(), b).unwrap();
        assert!(left.summary.physical_checks_pass);
        assert!(right.summary.physical_checks_pass);
        assert_eq!(left.samples.len(), right.samples.len());
        assert!(left
            .samples
            .iter()
            .zip(&right.samples)
            .all(|(x, y)| x.time == y.time));
    }
}
