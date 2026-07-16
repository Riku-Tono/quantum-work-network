# Milestone 4: first A/B experiment

## Conditions

- J = 1.0
- g = 0.25
- A input strength = 0.2
- A dephasing strength = 0.0
- B dephasing strength = 0.5
- End time = 10
- Time steps = 100
- Matched B input strength = 0.76484375

## Results

| Quantity | A: no noise | B: dephasing | B - A |
|---|---:|---:|---:|
| Load energy | 0.4539387586036373 | 0.4539114669191970 | -0.0000272916844403 |
| Load ergotropy | 0 | 0 | 0 |
| Source energy, net | 1.5218145006801025 | 2.7557440397959345 | 1.2339295391158320 |
| Dephasing energy, net | 0 | 1.8616679943449102e-17 | 1.8616679943449102e-17 |
| Top-level population | 5.0373974182078764% | 8.6239344045699398% | 3.5865369863620634 percentage points |

Relative load-energy difference: 0.0060121952406573816%.

## Success criteria

| Criterion | Result | Verdict |
|---|---:|---|
| Relative load-energy difference below 0.01% | 0.0060121952406573816% | PASS |
| Relative ergotropy difference at least 5% | 0% | FAIL |
| Ergotropy not approximately zero | A = 0, B = 0 | FAIL |
| Top-level population below 5% | A = 5.037397%, B = 8.623934% | FAIL |
| Trace check | A/B true/true | PASS |
| Hermiticity check | A/B true/true | PASS |
| Positivity check | A/B true/true | PASS |
| Energy-balance check | A/B true/true | PASS |

Overall success: **false**.

The aggregate physical check is false for both protocols only because the top-level population check failed. The trace, Hermiticity, positivity, and energy-balance checks all passed.

## Numerical diagnostics

| Diagnostic | A | B |
|---|---:|---:|
| Maximum trace error | 7.105427357601e-15 | 6.883382752676e-15 |
| Maximum Hermiticity error | 2.709486393190e-15 | 5.062806078492e-16 |
| Minimum eigenvalue | -1.946559121972e-20 | -2.375851389861e-20 |
| Energy-balance residual | -3.4075232014707524e-5 | -4.7840089765973287e-4 |

## Verification

- `cargo test --release --offline`: 36 passed, 0 failed, 1 ignored.
- `cargo run --release --offline --bin first_experiment`: completed successfully.
- The ignored test is the explicitly ignored dense 24-dimensional short-time smoke test.
