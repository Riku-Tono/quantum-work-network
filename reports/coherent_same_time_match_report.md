# Coherent input: same-time, same-load-energy comparison

## Scope

- Evaluation time: t = 7.9
- No Lindblad injection
- No Hamiltonian, model, or dephasing-definition change
- A: p_A = 0.2, gamma_phi = 0
- B: gamma_phi = 0.5, with p_B matched within [0, 1]
- This is a small root-finding calculation, not a large parameter search

The initial state is `sqrt(1-p)|0> + sqrt(p)|1>` on q1, with q2, q3, and the load in their ground states.

## Root check and matching

The p_B grid used 0.0, 0.1, ..., 1.0. The target was bracketed between p_B=0.6 and p_B=0.7:

- p_B=0.6: load_energy = 0.11616636569349417
- A target: load_energy = 0.12316871580126314
- p_B=0.7: load_energy = 0.13552742664240999

Bisection stopped after 9 iterations:

- matched p_B = 0.6361328125
- absolute load-energy difference = 0.00000665425708631928
- relative load-energy difference = 0.000054025545716139065
- required tolerance = 0.0001

The match passes.

## Matched result at t=7.9

| Quantity | A: no noise | B: dephasing |
|---|---:|---:|
| p | 0.2 | 0.6361328125 |
| Initial energy | 0.2 | 0.6361328125 |
| Load energy | 0.12316871580126314 | 0.12316206154417682 |
| Load ergotropy | 0.11361421674084743 | 0.008089034375357115 |
| Diagonal ergotropy | 0 | 0 |
| Coherence-derived ergotropy | 0.11361421674084743 | 0.008089034375357115 |
| Load coherence L1 | 0.62780561526960177 | 0.15699608343320479 |
| Level 0 population | 0.87683128419873746 | 0.87683793845581570 |
| Level 1 population | 0.12316871580126322 | 0.12316206154417575 |
| Level 2 population | 0 | 0 |
| Trace check | PASS | PASS |
| Hermiticity check | PASS | PASS |
| Positivity check | PASS | PASS |
| Energy-balance check | PASS | PASS |
| Top-level check | PASS | PASS |

## Differences

- Ergotropy difference A - B: 0.10552518236549031
- Ergotropy ratio A / B: 14.045460986909402
- Load ergotropy / initial energy:
  - A: 0.56807108370423709, or 56.8071%
  - B: 0.012715952103723806, or 1.2716%

For both conditions, diagonal ergotropy is zero. Therefore the reported load ergotropy is entirely coherence-derived on this matched sample.

## Interpretation boundary

This comparison uses the same evaluation time and nearly the same load energy. However, the initial input energies are not the same: B starts with 0.6361328125 while A starts with 0.2. B's initial energy is about 3.1807 times A's.

Therefore this is not an equal-input performance comparison. It establishes that, even after increasing B's initial excitation enough to match the load energy at t=7.9, the matched B load state has much less coherence and ergotropy in this single condition.

## Verification

- `cargo fmt --all -- --check`: passed
- `cargo test --release --offline`: 39 passed, 0 failed, 1 explicitly ignored dense smoke test
- Grid bracketing and bisection were performed by `coherent_same_time_match`
