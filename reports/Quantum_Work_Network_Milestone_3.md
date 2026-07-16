# Quantum Work Network

Rust reference implementation for the minimal quantum-work-network model.

## Milestone 3

Milestones 1, 2, and 2.1 remain intact. Milestone 3 adds state diagnostics and signed power accounting without adding protocol matching, efficiency claims, parameter sweeps, or plotting:

- `liouvillian`: explicit column-major density-matrix vectorization and construction of Hamiltonian and Lindblad superoperators.
- `propagator`: accuracy-first propagation with a dense matrix exponential at each requested time.
- `tests/full_24d_short_time.rs`: an opt-in integration test that builds the complete 24-dimensional model, constructs its `576 x 576` Liouvillian, and propagates it from `t = 0` to `t = 0.001`.
- `diagnostics`: reduced-load energy and ergotropy, energy decomposition, load current, source/dephasing power, physicality metrics, and signed trapezoidal power integration.

Protocol matching, efficiency claims, parameter sweeps, and plotting are deliberately not included yet.

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

Run the 26 regular tests:

```bash
cargo test
```

Run the deliberately ignored full 24-dimensional smoke test in release mode:

```bash
cargo test --release full_24d_short_time_smoke_test -- --ignored --nocapture
```

The ignored test computes a dense `576 x 576` matrix exponential and can be substantially slower than the regular suite. It checks the complete `ModelParams::default()` operator model with an injection collapse operator `sqrt(0.1) * sigma_1_plus` at `t = 0` and `t = 0.001`. At the later time it verifies dimensions, trace, Hermiticity, positivity within numerical tolerance, finite entries, and a nonzero change from the vacuum state.

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

Milestone 2.1 adds the ignored integration test described above. It uses only public crate APIs, so it also checks that the modules compose correctly from a downstream user’s perspective.

## Numerical scope

The dense exponential is suitable as a correctness baseline. A dense `576 x 576` exponential is computationally expensive, especially at many times. That inefficiency is deliberate at this stage: first make the physics trustworthy, then teach it to run without consuming the afternoon.
