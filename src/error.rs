use thiserror::Error;

#[derive(Debug, Error)]
pub enum PhysicsError {
    #[error("dimension mismatch: {0}")]
    DimensionMismatch(String),

    #[error("invalid subsystem specification: {0}")]
    InvalidSubsystem(String),

    #[error("matrix is not Hermitian within tolerance: error={error:e}, tolerance={tolerance:e}")]
    NonHermitian { error: f64, tolerance: f64 },

    #[error("density matrix trace is invalid: trace={trace:e}")]
    InvalidTrace { trace: f64 },

    #[error("density matrix is not positive semidefinite: minimum eigenvalue={minimum:e}")]
    NonPositiveState { minimum: f64 },

    #[error("eigendecomposition failed physical validation: {0}")]
    EigenFailure(String),

    #[error("invalid propagation time: {0}")]
    InvalidTime(String),

    #[error("invalid protocol parameter: {0}")]
    InvalidParameter(String),

    #[error("matching failed: {0}")]
    MatchingFailure(String),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}
