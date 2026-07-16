# Milestone 4 confirmed coherent-input result

This file freezes the coherent-input results established before Milestone 5a.
The existing result CSV files are the numerical record and are not modified by
the time-dependent-propagator work.

## Initial state and protocol

```text
q1 = sqrt(1-p)|0> + sqrt(p)|1>
q2 = |0>
q3 = |0>
load = |0>
```

- Lindblad injection was not used.
- This was a single prepared-initial-state transport experiment.
- The matched comparison used the same time and the same `load_energy`.
- Evaluation times were `3.0`, `5.0`, `7.9`, and `10.0`.
- B dephasing strengths were `0.1`, `0.2`, `0.5`, and `1.0`.
- All 16 combinations produced an energy match.

## Confirmed results

- A ergotropy was greater than B ergotropy in all 16 matched conditions.
- The matched A/B ergotropy ratio ranged from `1.229` to `49.318`.
- Diagonal ergotropy was zero in all matched conditions.
- The observed ergotropy difference was therefore coherence-derived
  ergotropy in these conditions.
- In the separate equal-input comparison (`p_B = p_A`), A ergotropy was also
  greater than B ergotropy in all 16 conditions.
- There were zero physical-check failures.

## Limits of the result

- This was single-shot initial-state transport.
- It did not create the state from vacuum using an external drive.
- It was not continuous supply.
- Some comparisons, including the energy-matched comparison, did not use equal
  initial input energy.
- No quantum advantage over a classical method was tested.
- No claim about real-world power-transmission performance is made.
