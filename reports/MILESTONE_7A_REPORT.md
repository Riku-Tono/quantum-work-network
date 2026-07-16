# Milestone 7a: 雑音位置の比較

## 1. 目的

固定した有限3サイト模型で、位相雑音を1サイトだけに置き、最後のloadへの有限時間内の影響を比較した。原因機構の断定は行わない。

## 2. 既存模型から変更した点

位相雑音を入れる0始まりsite集合を指定できる入口を追加した。既存の全3サイト雑音APIは全site `[0,1,2]` を渡すラッパーとして保持した。Milestone 7a専用binで比較・CSV・診断を実装した。

## 3. 変更していない条件

Hamiltonian、3つの二準位site、3準位load、全系真空、J=1、g=0.25、各角周波数=1、tau=3.2、Omega=0.2、t_max=10、pulse、load定義、ergotropy、RK4を変更していない。

## 4. 雑音位置の定義

site 1/2/3は内部index 0/1/2で、それぞれ入口・中央・出口側と呼ぶ。雑音あり条件は `sqrt(gamma_phi/2) sigma_z` を指定した1サイトにだけ置き、gamma_phi=0.5とした。noise_freeはcollapse operatorを作らない。

## 5. Milestone 5cとの違い

Milestone 5cの雑音ありBは3サイトすべてに位相雑音があった。今回は1サイトだけであり、数値を直接比較しない。

## 6. 数値手法

Milestone 5c本計算と同じ基準刻み `dt=0.0025`、保存間隔 `0.01`、固定刻みRK4を使用した。半減確認は `dt=0.00125`。積分は共通保存グリッド上の台形則。比の分母許容値は既存 `SIGNAL_TOLERANCE=1e-8`。分類の `relative_tol=5e-3` は既存の収束相対許容値を採用した。

## 7. 数値品質チェック

4条件のtrace、Hermiticity、positivity、load縮約整合、population bounds、top-level、energy ledger、予期しない非有限値、共通時間グリッド、共通初期状態はすべてPASS。ledger基準は `|r| <= 5e-5 + 5e-4*scale`。詳細は `local_noise_placement_checks.csv`。

## 8. t=10の比較

| condition | noise site | E(t=10) | W(t=10) | usable fraction | W_max | W_time_area |
|---|---|---:|---:|---:|---:|---:|
| noise_free | none | 5.4450767878e-2 | 5.2798274942e-2 | 9.6965161374e-1 | 5.5424064638e-2 | 1.9511404530e-1 |
| noise_entrance | site1 | 2.0554525850e-2 | 9.6153984705e-3 | 4.6779957566e-1 | 1.0762811785e-2 | 5.6246438529e-2 |
| noise_middle | site2 | 4.8178995305e-2 | 4.5346118640e-2 | 9.4120100165e-1 | 4.5346118640e-2 | 1.6829515255e-1 |
| noise_exit | site3 | 2.4577202789e-2 | 1.0192261965e-2 | 4.1470390478e-1 | 1.1318513177e-2 | 5.8513109845e-2 |

## 9. 最大ergotropyの比較

各条件の `W_max`, `t_at_W_max`, `E_at_W_max` はsummary CSVに保存した。最小順位は `site1`。

## 10. W_time_areaとW_time_meanの比較

`W_time_area`は0〜10にergotropyがどの程度存在したかの補助値であり、累積仕事・総仕事・供給された仕事ではない。最小順位は `site1`。

## 11. 雑音なし条件に対する比

| condition | R_E | R_W | R_use | R_Wmax | R_Warea | R_Wmean |
|---|---:|---:|---:|---:|---:|---:|
| noise_entrance | 0.37748826 | 0.18211577 | 0.48244088 | 0.19419023 | 0.28827468 | 0.28827468 |
| noise_middle | 0.88481756 | 0.85885606 | 0.97065893 | 0.81816660 | 0.86254761 | 0.86254761 |
| noise_exit | 0.45136559 | 0.19304157 | 0.42768341 | 0.20421658 | 0.29989184 | 0.29989184 |

## 12. 壊れ方の便宜的分類

- noise_entrance: `transport_and_quality_loss`
- noise_middle: `transport_and_quality_loss`
- noise_exit: `transport_and_quality_loss`

これは今回だけの整理で、一般法則や物理定理ではない。

## 13. 雑音位置の順位

- t=10でload_ergotropyを最も小さくした位置: `site1`
- t=10でusable_fractionを最も小さくした位置: `site3`
- W_maxを最も小さくした位置: `site1`
- W_time_areaを最も小さくした位置: `site1`

## 14. 時間刻み半減の結果

既存基準に従い、絶対差 `<= 1e-7` または相対差 `<= 5e-3` をPASSとした。

| condition | metric | abs diff | rel diff |
|---|---|---:|---:|
| noise_free | E_at_t10 | 8.654188e-14 | 1.589360e-12 |
| noise_free | W_at_t10 | 6.723233e-13 | 1.273381e-11 |
| noise_free | usable_fraction_at_t10 | 1.080624e-11 | 1.114446e-11 |
| noise_free | W_max | 9.919843e-14 | 1.789808e-12 |
| noise_free | W_time_area | 8.311130e-13 | 4.259626e-12 |
| noise_entrance | E_at_t10 | 1.661518e-14 | 8.083466e-13 |
| noise_entrance | W_at_t10 | 1.294798e-13 | 1.346588e-11 |
| noise_entrance | usable_fraction_at_t10 | 5.921152e-12 | 1.265746e-11 |
| noise_entrance | W_max | 2.677893e-14 | 2.488098e-12 |
| noise_entrance | W_time_area | 1.216319e-13 | 2.162481e-12 |

W_at_t10が最小の雑音位置は刻み半減後も `site1` で変わらず、全10比較量がPASSした。

## 15. 直接確認できたこと

固定した有限3サイト模型、gamma_phi=0.5、Omega=0.2、固定単発pulse、0<=t<=10で、1サイトだけに雑音を置いた場合のload energy・ergotropy等の位置別差を直接確認した。t=10のW最小は `site1`、usable fraction最小は `site3` だった。

## 16. 確認できていないこと

他のgamma_phi、Omega、pulse、長いnetwork、長時間、連続運転、放電、抽出操作、現実装置、古典比較、文献上の新規性は確認していない。

## 17. 主張してはいけないこと

雑音位置の普遍則、任意パラメータや長いnetworkへの一般化、量子優位、現実の送電・装置効率、雑音が常にergotropyを減らす一般論、各site雑音が注入・輸送・受け渡しだけを壊したという因果断定はできない。

## 18. 作成ファイル一覧

- `src/bin/local_noise_placement.rs`
- `local_noise_placement_timeseries.csv`
- `local_noise_placement_summary.csv`
- `local_noise_placement_ratios.csv`
- `local_noise_placement_checks.csv`
- `local_noise_placement_convergence.csv`
- `MILESTONE_7A_REPORT.md`

