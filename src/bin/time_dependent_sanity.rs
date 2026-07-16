use nalgebra::linalg::Schur;
use quantum_work_network::liouvillian::build_liouvillian;
use quantum_work_network::matrix::{frobenius_norm, hermiticity_error, ComplexMatrix, C64};
use quantum_work_network::propagator::DenseExponentialPropagator;
use quantum_work_network::time_dependent::{SaveSchedule, TimeDependentRk4};

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

fn damping(gamma: f64) -> ComplexMatrix {
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
    let collapses = vec![damping(0.3)];
    let reference = DenseExponentialPropagator::new(
        build_liouvillian(&h, &collapses).expect("valid constant Liouvillian"),
        2,
    )
    .expect("valid dense propagator")
    .propagate_times(&ground(), &[0.7])
    .expect("dense propagation succeeds")
    .pop()
    .expect("one reference state")
    .rho;
    let actual = TimeDependentRk4::new(dt)
        .expect("positive dt")
        .propagate(
            &ground(),
            0.0,
            0.7,
            |_| h.clone(),
            |_| collapses.clone(),
            SaveSchedule::Times(vec![0.7]),
        )
        .expect("RK4 propagation succeeds")
        .pop()
        .expect("one final state")
        .rho;
    frobenius_norm(&(actual - reference))
}

fn worst_closed_system_minimum_eigenvalue(dt: f64) -> f64 {
    let h = sigma_x() * C64::new(1.1, 0.0);
    TimeDependentRk4::new(dt)
        .expect("positive dt")
        .propagate(
            &ground(),
            0.0,
            1.0,
            |_| h.clone(),
            |_| Vec::new(),
            SaveSchedule::EveryStep,
        )
        .expect("closed-system propagation succeeds")
        .iter()
        .map(|state| minimum_eigenvalue(&state.rho))
        .fold(f64::INFINITY, f64::min)
}

fn main() {
    let pulse_duration = 2.0;
    let omega0 = 0.7;
    let dt = 0.002;
    let x = sigma_x();
    let states = TimeDependentRk4::new(dt)
        .expect("positive dt")
        .propagate(
            &ground(),
            0.0,
            pulse_duration,
            |time| {
                let omega = if (0.0..=pulse_duration).contains(&time) {
                    omega0 * (std::f64::consts::PI * time / pulse_duration).sin().powi(2)
                } else {
                    0.0
                };
                &x * C64::new(omega, 0.0)
            },
            |_| Vec::new(),
            SaveSchedule::Interval(0.01),
        )
        .expect("pulse propagation succeeds");

    let mut maximum_excited = 0.0_f64;
    let mut maximum_coherence = 0.0_f64;
    let mut maximum_trace_error = 0.0_f64;
    let mut maximum_hermiticity_error = 0.0_f64;
    let mut worst_minimum_eigenvalue = f64::INFINITY;
    let mut all_finite = true;
    for state in &states {
        maximum_excited = maximum_excited.max(state.rho[(1, 1)].re);
        maximum_coherence = maximum_coherence.max(state.rho[(0, 1)].norm());
        maximum_trace_error =
            maximum_trace_error.max((state.rho.trace() - C64::new(1.0, 0.0)).norm());
        maximum_hermiticity_error = maximum_hermiticity_error.max(hermiticity_error(&state.rho));
        worst_minimum_eigenvalue = worst_minimum_eigenvalue.min(minimum_eigenvalue(&state.rho));
        all_finite &= state
            .rho
            .iter()
            .all(|value| value.re.is_finite() && value.im.is_finite());
    }

    let zero_states = TimeDependentRk4::new(dt)
        .expect("positive dt")
        .propagate(
            &ground(),
            0.0,
            pulse_duration,
            |_| ComplexMatrix::zeros(2, 2),
            |_| Vec::new(),
            SaveSchedule::Times(vec![pulse_duration]),
        )
        .expect("zero-drive propagation succeeds");
    let zero_drive_change = frobenius_norm(&(&zero_states.last().unwrap().rho - ground()));

    let error_dt = constant_reference_error(0.08);
    let error_half = constant_reference_error(0.04);
    let error_quarter = constant_reference_error(0.02);
    let positivity_dt = worst_closed_system_minimum_eigenvalue(0.2);
    let positivity_half = worst_closed_system_minimum_eigenvalue(0.1);
    let positivity_quarter = worst_closed_system_minimum_eigenvalue(0.05);

    println!("pulse_dt,{dt:.16e}");
    println!(
        "pulse_final_excited_probability,{:.16e}",
        states.last().unwrap().rho[(1, 1)].re
    );
    println!("pulse_maximum_excited_probability,{maximum_excited:.16e}");
    println!("pulse_maximum_coherence,{maximum_coherence:.16e}");
    println!("zero_drive_frobenius_change,{zero_drive_change:.16e}");
    println!("constant_reference_dt,8.0000000000000002e-2");
    println!("constant_reference_error_dt,{error_dt:.16e}");
    println!("constant_reference_error_dt_over_2,{error_half:.16e}");
    println!("constant_reference_error_dt_over_4,{error_quarter:.16e}");
    println!(
        "constant_reference_error_ratio_dt_to_half,{:.16e}",
        error_dt / error_half
    );
    println!(
        "constant_reference_error_ratio_half_to_quarter,{:.16e}",
        error_half / error_quarter
    );
    println!("maximum_trace_error,{maximum_trace_error:.16e}");
    println!("maximum_hermiticity_error,{maximum_hermiticity_error:.16e}");
    println!("worst_minimum_eigenvalue,{worst_minimum_eigenvalue:.16e}");
    println!("positivity_probe_min_eigenvalue_dt,{positivity_dt:.16e}");
    println!("positivity_probe_min_eigenvalue_dt_over_2,{positivity_half:.16e}");
    println!("positivity_probe_min_eigenvalue_dt_over_4,{positivity_quarter:.16e}");
    println!("all_elements_finite,{all_finite}");
}
