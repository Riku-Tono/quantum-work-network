# Milestone 5b coherent-drive sanity check

## 実装済み

- 真空初期状態、既存24次元bare network `H0`。
- `Omega0=0.2`, `omega_drive=1`, `tau=3.2`, `t_end=10` の単一 `sin^2` パルス。
- A: `gamma_phi=0`; B: `gamma_phi=0.5`。Bは3サイトすべてへ `sqrt(gamma_phi/2) sigma_z,j`。loadへの直接雑音なし。
- Lindblad励起注入、energy matching、探索、最適化は未使用。

## 実行確認済み

主計算は最大RK4刻み `0.005`、CSV間隔 `0.01`。

### A

- 最大load energy: `5.6859734359e-2` at `t=9.50`
- 最大load ergotropy: `5.5424064639e-2` at `t=9.48`
- 最大coherence-derived ergotropy: `5.5424064639e-2` at `t=9.48`
- 最大coherence L1: `4.7204015323e-1` at `t=9.58`
- 最大最上段占有率: `6.0632991982e-5` at `t=10.00`
- drive energy net/in/out: `7.5315933440e-2` / `7.5315933440e-2` / `0.0000000000e0`
- dephasing energy net/in/out: `0.0000000000e0` / `0.0000000000e0` / `0.0000000000e0`
- ledger absolute residual: `1.4342804722e-11`
- max trace/Hermiticity errors: `1.776e-15` / `0.000e0`
- worst minimum eigenvalue: `-1.049e-10`
- max power imaginary parts drive/dephasing: `8.707e-18` / `0.000e0`
- finite=true, physical=true, ledger=true, top-level=true

### B

- 最大load energy: `1.2596874861e-2` at `t=10.00`
- 最大load ergotropy: `3.0302005933e-3` at `t=5.63`
- 最大coherence-derived ergotropy: `3.0302005933e-3` at `t=5.63`
- 最大coherence L1: `1.1145862109e-1` at `t=5.63`
- 最大最上段占有率: `1.6920767383e-5` at `t=10.00`
- drive energy net/in/out: `5.9618618775e-2` / `5.9618618775e-2` / `0.0000000000e0`
- dephasing energy net/in/out: `1.1327655621e-13` / `1.1547305350e-12` / `1.0414539788e-12`
- ledger absolute residual: `9.9229088294e-12`
- max trace/Hermiticity errors: `1.332e-15` / `0.000e0`
- worst minimum eigenvalue: `-2.465e-18`
- max power imaginary parts drive/dephasing: `3.447e-18` / `4.321e-18`
- finite=true, physical=true, ledger=true, top-level=true

## 同時刻A/B確認

- 最大load ergotropy A/B: `5.5424064639e-2` / `3.0302005933e-3` (比 `18.290560`)
- 最大load coherence L1 A/B: `4.7204015323e-1` / `1.1145862109e-1` (比 `4.235116`)
- t=tau: load energy A `2.2449518483e-3`, B `8.6907555399e-4`; ergotropy A `2.1805953633e-3`, B `5.6588983861e-4`; coherence-derived ergotropy A `2.1805953633e-3`, B `5.6588983861e-4`; coherence L1 A `9.4071838695e-2`, B `4.7795200327e-2`
- t=10: load energy A `5.4450767879e-2`, B `1.2596874861e-2`; ergotropy A `5.2798274953e-2`, B `2.3652476827e-3`; coherence-derived ergotropy A `5.2798274953e-2`, B `2.3652476827e-3`; coherence L1 A `4.6453837011e-1`, B `9.7568729205e-2`
- A最大ergotropy時刻 t=9.48: load energy A `5.6854985096e-2`, B `1.2204746895e-2`; ergotropy A `5.5424064639e-2`, B `2.5123340537e-3`; coherence-derived ergotropy A `5.5424064639e-2`, B `2.5123340537e-3`; coherence L1 A `4.7159598955e-1`, B `1.0062483350e-1`

A/Bそれぞれの最大値は時刻が異なる可能性があり、その比は公平な同時刻比較ではなくsanity check用。今回は同一load energyへのmatchingもしていない。

刻み収束の全16指標はCSVに絶対差と相対差を保存した。ledger residualは分母がほぼゼロなので、相対差ではなく絶対差を主判定に使用した。

## 成功確認

1. Aの非ゼロload coherence: **PASS**
2. Aの非ゼロload ergotropy: **PASS**
3. 同時刻でBのcoherenceまたはergotropy低下: **PASS**
4. A/Bの最上段占有率5%未満: **PASS**
5. A/Bの物理チェック: **PASS**
6. A/Bの固定H0エネルギー台帳: **PASS**
7. dt半減による主要量収束: **PASS**

閾値: signal `>1e-8`; trace/Hermiticity `<=1e-8`/`<=1e-8`; minimum eigenvalue `>=-1e-8`; top level `<0.05`; ledger `|r| <= 5e-5 + 5e-4*scale`; convergence `abs<=1e-7` または `rel<=5e-3`。

## 未確認

同一時刻・同一load energyの公平比較、energy matching、パラメータ探索、連続駆動、仕事抽出、古典比較、量子優位は未確認。
