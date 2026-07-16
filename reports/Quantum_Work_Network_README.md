# Quantum Work Network

Rust reference implementation for the minimal quantum-work-network model.

## Milestone 2

Milestone 1 remains intact and Milestone 2 adds two time-evolution modules:

- `liouvillian`: explicit column-major density-matrix vectorization and construction of Hamiltonian and Lindblad superoperators.
- `propagator`: accuracy-first propagation with a dense matrix exponential at each requested time.

Diagnostics, protocol matching, parameter sweeps, and plotting are deliberately not included yet.

## Basis and vectorization conventions

The physical tensor order used by `operators` is

```text
|q1, q2, q3, load>
```

with the load index varying fastest. Qubit state `|0>` is empty and `|1>` is excited.

Density matrices use **column-major vectorization**:

```text
vec(rho) = [rho(0,0), rho(1,0), ..., rho(0,1), rho(1,1), ...]^T
```

Therefore

```text
vec(A rho B) = (B^T tensor A) vec(rho).
```

A `24 x 24` density matrix becomes a vector of length `576`, and the Liouvillian is `576 x 576`.

## Liouvillian convention

The generator is

```text
L = -i (I tensor H - H^T tensor I)
    + sum_k [L_k* tensor L_k
             - 1/2 I tensor (L_k^dagger L_k)
             - 1/2 (L_k^dagger L_k)^T tensor I].
```

Collapse operators are passed **with their coefficients already included**. Examples:

```text
sqrt(gamma) sigma_minus
sqrt(gamma_phi / 2) sigma_z
```

## Propagator

`DenseExponentialPropagator` computes

```text
vec(rho(t)) = exp(L t) vec(rho(0))
```

independently for every requested nonnegative time. This is a dense, accuracy-first reference implementation. It is intentionally not the eventual parameter-sweep engine; a Krylov `expm_multiply` implementation belongs in a later milestone.

## Tests

Run:

```bash
cargo test
```

Milestone 2 tests cover:

- column-major vectorization and inverse conversion;
- `24 x 24 -> 576` vector dimension;
- `t = 0` identity;
- trace, Hermiticity, and positivity preservation;
- equality with closed-system unitary evolution;
- invariant state for `H = 0` with no collapse operators;
- analytic amplitude damping;
- analytic pure dephasing;
- every Milestone 1 operator, partial-trace, and ergotropy test.

## Numerical scope

The dense exponential is suitable as a correctness baseline. A dense `576 x 576` exponential is computationally expensive, especially at many times. That inefficiency is deliberate at this stage: first make the physics trustworthy, then teach it to run without consuming the afternoon.
