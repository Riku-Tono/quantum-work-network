# Quantum Work Network: Static Foundations

Milestone 1 implements and tests the three static modules that must be correct before any Liouvillian time evolution is trusted:

1. `operators`
2. `partial_trace`
3. `ergotropy`

## Frozen conventions

- Tensor order: `|q1, q2, q3, load>`.
- The rightmost load index varies fastest in flattened basis indices.
- Qubit `|0>` is empty and `|1>` is excited.
- Chain onsite Hamiltonian is `omega * sum_i |1><1|_i`.
- Load basis is `|0>, |1>, ..., |d-1>`.
- The truncated load annihilation operator satisfies `[b,b†] = I - d |d-1><d-1|`.
- Column-major vectorization is reserved for the later Liouvillian milestone.

## Run

```bash
cargo test
```

The current execution environment used to generate this project did not contain a Rust toolchain, so the source and tests were produced but not compiled here. Run the command above in a Rust environment before treating the milestone as certified.
