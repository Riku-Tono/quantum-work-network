# Milestone 4: current A/B time-series result

## Scope

No parameter search and no model change were performed. The existing matched pair was rerun on the existing uniform grid:

- A: input strength 0.2, dephasing strength 0.0
- B: matched input strength 0.76484375, dephasing strength 0.5
- J = 1.0, g = 0.25
- time = 0.0 through 10.0
- interval = 0.1
- 101 samples per condition

## Extrema

| Condition | Maximum load ergotropy | Time | Maximum top-level population | Time | Final top-level population |
|---|---:|---:|---:|---:|---:|
| A: no noise | 0 | 0.0 | 0.050373974182078764 | 10.0 | 0.050373974182078764 |
| B: dephasing | 0 | 0.0 | 0.086239344045699398 | 10.0 | 0.086239344045699398 |

The reported ergotropy time is 0.0 because every one of the 101 sampled values is exactly zero and the summary keeps the first occurrence of a tied maximum.

## Time-series checks

- `load_ergotropy` is zero at every sampled time for both A and B.
- `|rho_01|`, `|rho_02|`, `|rho_12|`, and the full off-diagonal L1 sum are zero at every sampled time for both A and B.
- The load populations remain ordered `p0 >= p1 >= p2` at every sampled time, so the diagonal load state remains passive on this grid.
- The maximum absolute population-sum error is 7.11e-15 for A and 6.88e-15 for B.
- The maximum top-level population occurs at the final sample for both conditions. Therefore the corrected all-time check has the same numerical value as the former final-time-only check for this run.

## Corrected truncation verdict

The physical check now compares the 5% threshold against the maximum top-level population over all sampled times, not only the final value.

- A: 5.037397% -> FAIL
- B: 8.623934% -> FAIL

## CSV columns

The time-series CSV uses long format with one row per condition and time. It contains:

- `load_energy`
- `load_ergotropy`
- `load_level_population_0`, `_1`, `_2`
- `load_rho_abs_0_1`, `_0_2`, `_1_2`
- `load_off_diagonal_l1`, defined as the sum of magnitudes over both off-diagonal triangles

## Verification

- `cargo fmt --all -- --check`: passed
- `cargo test --release --offline`: 37 passed, 0 failed, 1 explicitly ignored dense smoke test
- `cargo run --release --offline --bin first_experiment`: completed
