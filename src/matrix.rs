use nalgebra::DMatrix;
use num_complex::Complex64;

pub type C64 = Complex64;
pub type ComplexMatrix = DMatrix<C64>;

pub fn eye(dim: usize) -> ComplexMatrix {
    ComplexMatrix::identity(dim, dim)
}

pub fn dagger(matrix: &ComplexMatrix) -> ComplexMatrix {
    matrix.adjoint()
}

pub fn kron(a: &ComplexMatrix, b: &ComplexMatrix) -> ComplexMatrix {
    let rows = a.nrows() * b.nrows();
    let cols = a.ncols() * b.ncols();
    let mut out = ComplexMatrix::zeros(rows, cols);

    for i in 0..a.nrows() {
        for j in 0..a.ncols() {
            let coefficient = a[(i, j)];
            for k in 0..b.nrows() {
                for l in 0..b.ncols() {
                    out[(i * b.nrows() + k, j * b.ncols() + l)] = coefficient * b[(k, l)];
                }
            }
        }
    }

    out
}

pub fn frobenius_norm(matrix: &ComplexMatrix) -> f64 {
    matrix.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt()
}

pub fn hermiticity_error(matrix: &ComplexMatrix) -> f64 {
    frobenius_norm(&(matrix - matrix.adjoint()))
}

pub fn trace_real(matrix: &ComplexMatrix) -> f64 {
    matrix.trace().re
}

pub fn commutator(a: &ComplexMatrix, b: &ComplexMatrix) -> ComplexMatrix {
    a * b - b * a
}

pub fn expectation(rho: &ComplexMatrix, observable: &ComplexMatrix) -> C64 {
    (rho * observable).trace()
}
