# Coherent-input comparison at equal load energy

## Scope and interpretation

Only the existing `coherent_input_timeseries.csv` was used. No model change, new time evolution, energy matching run, or parameter search was performed.

This is an **existence check at equal load energy on different time points**. It is not a performance comparison at the same time or under the same supply conditions.

- A first rising branch: t=0.0 to 3.9, load_energy=0 to 0.1291409314244132
- B first rising branch: t=0.0 to 10.0, load_energy=0 to 0.10792882657075184
- Common target grid: 0.01 to 0.10 in steps of 0.01

For each target, the complex 3x3 load density matrix was linearly interpolated between the two neighboring samples. Energy, ergotropy, coherence L1, and dephased-state ergotropy were then recalculated from that interpolated density matrix.

## Result

| target load energy | A time | B time | A ergotropy | B ergotropy | A - B | A / B |
|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 1.649284 | 1.930763 | 0.00503604 | 0.00273906 | 0.00229698 | 1.8386 |
| 0.02 | 1.915409 | 2.347795 | 0.01027345 | 0.00490011 | 0.00537334 | 2.0966 |
| 0.03 | 2.105632 | 2.703258 | 0.01567855 | 0.00673901 | 0.00893954 | 2.3265 |
| 0.04 | 2.262573 | 3.071493 | 0.02116594 | 0.00818675 | 0.01297918 | 2.5854 |
| 0.05 | 2.404115 | 3.530580 | 0.02695307 | 0.00914145 | 0.01781162 | 2.9484 |
| 0.06 | 2.535577 | 4.319788 | 0.03275973 | 0.00908165 | 0.02367808 | 3.6072 |
| 0.07 | 2.663351 | 5.502045 | 0.03882197 | 0.00848285 | 0.03033912 | 4.5765 |
| 0.08 | 2.790982 | 6.369773 | 0.04514556 | 0.00857846 | 0.03656710 | 5.2627 |
| 0.09 | 2.923190 | 7.191722 | 0.05152518 | 0.00855935 | 0.04296582 | 6.0198 |
| 0.10 | 3.065451 | 8.352317 | 0.05812139 | 0.00776160 | 0.05035980 | 7.4883 |

At every target energy, A has larger load ergotropy and larger load coherence L1 than B.

The dephased-state ergotropy `W(Delta[rho_load])` is zero for both A and B at all ten targets. Therefore, on this grid:

`W_coh = W(rho_load) - W(Delta[rho_load]) = W(rho_load)`

## Interpolation checks

- Maximum absolute A target-energy error: 1.3877787807814457e-17
- Maximum absolute B target-energy error: 1.3877787807814457e-17
- Maximum absolute actual A/B energy difference: 2.0816681711721685e-17
- Minimum eigenvalue among interpolated density matrices: 0

The same-energy condition is therefore satisfied to floating-point precision, and the interpolated density matrices remain positive semidefinite.
