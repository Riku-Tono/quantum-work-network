use nalgebra::linalg::Schur;
use quantum_work_network::liouvillian::build_liouvillian;
use quantum_work_network::matrix::{frobenius_norm, hermiticity_error};
use quantum_work_network::operators::{build_operators, ModelParams};
use quantum_work_network::propagator::DenseExponentialPropagator;
use quantum_work_network::{ComplexMatrix, C64};

fn vacuum_density_matrix(dim: usize) -> ComplexMatrix {
    let mut rho = ComplexMatrix::zeros(dim, dim);
    rho[(0, 0)] = C64::new(1.0, 0.0);
    rho
}

fn minimum_hermitian_eigenvalue(matrix: &ComplexMatrix) -> f64 {
    // The propagated state should be Hermitian up to roundoff. Schur is used
    // consistently with the existing Milestone 2 physicality tests.
    let schur = Schur::new(matrix.clone());
    let (_, triangular) = schur.unpack();
    (0..triangular.nrows())
        .map(|i| triangular[(i, i)].re)
        .fold(f64::INFINITY, f64::min)
}

#[test]
#[ignore = "computes a dense 576x576 matrix exponential; run explicitly in release mode"]
fn full_24d_short_time_smoke_test() {
    const HILBERT_DIM: usize = 24;
    const LIOUVILLE_DIM: usize = HILBERT_DIM * HILBERT_DIM;
    const GAMMA_IN: f64 = 0.1;
    const SHORT_TIME: f64 = 0.001;
    const TRACE_TOLERANCE: f64 = 1.0e-9;
    const HERMITICITY_TOLERANCE: f64 = 1.0e-9;
    const POSITIVITY_TOLERANCE: f64 = 1.0e-8;
    const CHANGE_FLOOR: f64 = 1.0e-12;

    // Tensor order is |q1, q2, q3, load>. Basis index zero is
    // |0, 0, 0, 0>, because the load index varies fastest.
    let params = ModelParams::default();
    let operators = build_operators(&params).expect("build 24-dimensional operators");
    assert_eq!(operators.h_total.shape(), (HILBERT_DIM, HILBERT_DIM));

    let rho0 = vacuum_density_matrix(HILBERT_DIM);
    let injection = &operators.sigma_1_plus * C64::new(GAMMA_IN.sqrt(), 0.0);

    let liouvillian = build_liouvillian(&operators.h_total, &[injection])
        .expect("build the 576-dimensional Liouvillian");
    assert_eq!(liouvillian.shape(), (LIOUVILLE_DIM, LIOUVILLE_DIM));

    let propagator = DenseExponentialPropagator::new(liouvillian, HILBERT_DIM)
        .expect("construct dense exponential propagator");
    let states = propagator
        .propagate_times(&rho0, &[0.0, SHORT_TIME])
        .expect("propagate the full model at two short times");

    assert_eq!(states.len(), 2);
    assert!(frobenius_norm(&(&states[0].rho - &rho0)) < 1.0e-12);

    let rho_t = &states[1].rho;
    assert_eq!(rho_t.shape(), (HILBERT_DIM, HILBERT_DIM));

    let trace = rho_t.trace();
    assert!(
        (trace.re - 1.0).abs() < TRACE_TOLERANCE,
        "trace real part was {}",
        trace.re
    );
    assert!(
        trace.im.abs() < TRACE_TOLERANCE,
        "trace imaginary part was {}",
        trace.im
    );

    let hermiticity = hermiticity_error(rho_t);
    assert!(
        hermiticity < HERMITICITY_TOLERANCE,
        "Hermiticity error was {hermiticity:e}"
    );

    let minimum_eigenvalue = minimum_hermitian_eigenvalue(rho_t);
    assert!(
        minimum_eigenvalue > -POSITIVITY_TOLERANCE,
        "minimum eigenvalue was {minimum_eigenvalue:e}"
    );

    assert!(
        rho_t
            .iter()
            .all(|value| value.re.is_finite() && value.im.is_finite()),
        "propagated density matrix contained a non-finite value"
    );

    let change = frobenius_norm(&(rho_t - &rho0));
    assert!(
        change > CHANGE_FLOOR,
        "state did not measurably change; Frobenius distance was {change:e}"
    );
}
