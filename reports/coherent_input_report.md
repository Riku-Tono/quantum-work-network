# Coherent-input sanity check

## Scope

This is a single sanity check, not a parameter search. No energy matching or Hamiltonian/model change was performed.

Initial state:

- q1 = (|0> + |1>) / sqrt(2)
- q2 = |0>
- q3 = |0>
- load = |0>
- no Lindblad injection collapse operator

Conditions:

- A: gamma_phi = 0.0
- B: gamma_phi = 0.5
- J = 1.0, g = 0.25
- t = 0.0 through 10.0 at intervals of 0.1
- 101 sampled times per condition

## Sampled extrema

| Condition | Maximum load ergotropy | Time | Maximum off-diagonal magnitude | Component | Time |
|---|---:|---:|---:|---|---:|
| A: no noise | 0.24479125936715285 | 7.9 | 0.39237850954350106 | rho_01 | 7.9 |
| B: dephasing | 0.0093142147516385304 | 3.8 | 0.091611812049261859 | rho_01 | 3.7 |

At the maximum off-diagonal sample:

- A, t=7.9: rho_01 = 0.39196311483247470 + 0.018050245494275006 i
- B, t=3.7: rho_01 = -0.048539248957815126 + 0.077695980704026918 i

These maxima are maxima on the sampled 0.1-spaced time grid.

## Verdict

The model does produce nonzero load ergotropy from the phase-aligned coherent initial input.

Pure dephasing does not reduce the result exactly to zero in this run, but the sampled maximum falls from about 0.24479 in A to about 0.009314 in B.

The load never reaches level 2 in either condition, as expected for this source-free initial state containing at most one excitation. All physical checks pass for both A and B.

## CSV contents

The time-series CSV contains `load_energy`, `load_ergotropy`, all three load-level populations, and the real part, imaginary part, and magnitude of each independent off-diagonal density-matrix component (`rho_01`, `rho_02`, `rho_12`).

## Verification

- `cargo fmt --all -- --check`: passed
- `cargo test --release --offline`: 38 passed, 0 failed, 1 explicitly ignored dense smoke test
- `cargo run --release --offline --bin coherent_input_sanity`: completed
