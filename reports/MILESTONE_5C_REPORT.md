# Milestone 5c coherent-drive energy-matched comparison

## 実装済み

- `Omega_B=0.2..1.0`, step `0.01` の81点を全走査。単調性は仮定していない。
- 全符号変化区間を二分法で精密化し、重複閾値 `1e-7` で統合。
- 主根は `|Omega_B-0.2|` 最小規則。
- `dt=0.0025`は同じ近傍ブラケット内だけで再調整。

## 実行確認済み

- 符号変化区間数: `1`
- グリッド上直接一致数: `0`
- 局所接触疑い: `なし`
- 統合後root数: `1`
- 最良グリッド参考点: Omega `0.4300`, residual `-4.419149e-4`

| quantity | dt=0.005 | dt=0.0025 |
|---|---:|---:|
| matched Omega_B | 0.431953125000 | 0.431953125000 |
| A load energy | 5.4450767879e-2 | 5.4450767878e-2 |
| B load energy | 5.4452946588e-2 | 5.4452946589e-2 |
| relative match error | 4.001e-5 | 4.001e-5 |
| A ergotropy | 5.2798274953e-2 | 5.2798274942e-2 |
| B ergotropy | 8.2846362486e-3 | 8.2846362481e-3 |

- matched Omega difference after dt halving: `0.0000000000e0`
- absolute load-energy error: `2.1787115244e-6`
- relative load-energy error: `4.0012503208e-5`
- ergotropy A-B: `4.4513638694e-2`
- ergotropy A/B: `6.3730347792`
- relative ergotropy advantage: `5.3730347792`
- coherence-derived ergotropy A-B: `4.4513638694e-2`
- coherence L1 A-B: `2.8506417830e-1`
- drive energy in B/A: `3.4278523291`
- load-energy delivery fraction A/B: `0.7229647884` / `0.2109174044`
- load-ergotropy delivery fraction A/B: `0.7010239737` / `0.0320896128`

刻み収束CSVでは全量の絶対差を保存した。ledger residualと最小固有値はゼロ近傍なので、相対差ではなく絶対差を主に評価する。


## 成功条件

1. root exists: **PASS**
2. energy match: **PASS**
3. both ergotropy > 1e-3: **PASS**
4. A advantage > 5%: **PASS**
5. A coherence ergotropy > B: **PASS**
6. top level < 5%: **PASS**
7. physical checks: **PASS**
8. energy ledgers: **PASS**
9. direction stable after dt halving: **PASS**
10. fine energy match: **PASS**

## 公平性と未確認

一致しているのは比較時刻、最終load energy、模型、パルス形状、駆動周波数。一致していないのは駆動強度、drive energy in、総投入エネルギー。これは同じ時刻に同じload energyを持つ状態の仕事価値比較であり、等入力費用比較ではない。Omegaも異なるため雑音だけを独立に変えた比較でもない。熱力学的効率ではなくdelivery fractionのみを報告した。連続運転、仕事抽出、古典比較、量子優位は未確認。
