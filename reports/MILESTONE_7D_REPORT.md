# Milestone 7d: Partial end-protection strength sweep

## 1. 目的

中央siteのgamma=0.5を固定し、両端gammaだけを0.5から0へ下げたときの回復曲線を調べた。

## 2. Milestone 7cから進めた問い

7cの完全両端雑音除去という端点比較から、途中の有限gammaで回復がどのようにつながるかへ進めた。

## 3. 部分保護の定義

site1/site3のdephasing rateを同時に数値的に弱める理想感度試験。現実的保護装置、制御、cost、不完全装置モデルではない。

## 4. sweep条件

実行gamma_end: `0.50, 0.40, 0.30, 0.20, 0.15, 0.10, 0.05, 0.00`。基本6点を先に実行。追加: max normalized second difference=0.354119 (threshold=0.20), max sensitivity ratio=21.420200 (threshold=3.0)

## 5. 変更していない物理条件

3site+3準位load、J=1、g=0.25、周波数=1、真空初期状態、tau=3.2、Omega=0.2、t_max=10、中央gamma=0.5、load無雑音を固定した。

## 6. site別gamma実装

新APIは `[gamma_site1, gamma_site2, gamma_site3]` を受け、負値・非有限値を拒否し、gamma=0のcollapse operatorを除外する。既存共通gamma APIは新APIへ委譲し、短時間密度行列の完全一致をunit testした。

## 7. 数値手法

RK4、基準dt=0.0025、半減dt=0.00125、保存間隔0.01。時間面積は台形則。感度と二階差分は離散点だけから計算し、連続微分・物理的感受率・臨界指数とは呼ばない。

## 8. 数値品質チェック

全152項目PASS。trace、Hermiticity、positivity、population、縮約状態、load top-level、ledger、有限性、W<=E、usable範囲、drive整合、gamma mapping、collapse数、共通grid、初期状態を確認した。

## 9. 端点再現

gamma_end=0.5は7c all_noisy、gamma_end=0は7c protect_both_ends、別計算のnoise_freeは7c noise_freeの8比較量を絶対誤差1e-9以内で再現した。

## 10. t=10の回復曲線

| gamma_end | E | W | normalized W recovery |
|---:|---:|---:|---:|
| 0.50 | 1.2596874861e-2 | 2.3652476826e-3 | 0.00000000 |
| 0.40 | 1.5020610827e-2 | 3.5423888746e-3 | 0.02738756 |
| 0.30 | 1.8447315010e-2 | 5.7168843200e-3 | 0.07797973 |
| 0.20 | 2.3649039203e-2 | 1.0137269800e-2 | 0.18082514 |
| 0.15 | 2.7373125074e-2 | 1.4070100199e-2 | 0.27232702 |
| 0.10 | 3.2283452181e-2 | 2.0131518735e-2 | 0.41335298 |
| 0.05 | 3.8931384397e-2 | 2.9733407828e-2 | 0.63675211 |
| 0.00 | 4.8178995305e-2 | 4.5346118640e-2 | 1.00000000 |

## 11. W_maxの回復曲線

| gamma_end | W_max | t_at_W_max | normalized recovery |
|---:|---:|---:|---:|
| 0.50 | 3.0302005931e-3 | 5.63 | 0.00000000 |
| 0.40 | 4.1847112328e-3 | 5.85 | 0.02728313 |
| 0.30 | 6.3184273762e-3 | 8.51 | 0.07770662 |
| 0.20 | 1.0756371404e-2 | 8.95 | 0.18258308 |
| 0.15 | 1.4600721000e-2 | 9.19 | 0.27343187 |
| 0.10 | 2.0465310304e-2 | 9.46 | 0.41202249 |
| 0.05 | 2.9800196232e-2 | 9.80 | 0.63262235 |
| 0.00 | 4.5346118640e-2 | 10.00 | 1.00000000 |

## 12. E_time_areaとW_time_area

両者は状態量の時間面積で、累積流入エネルギー・累積抽出仕事ではない。gamma低下時の判定はE=`monotonic_nondecreasing`、W=`monotonic_nondecreasing`。

## 13. usable fractionの回復曲線

0.5で0.18776464、0で0.94120100。0.4で完全両端保護回復の6.38%、0.3で16.21%へ到達した。

## 14. 隣接感度

- E_at_t10: 最大区間 `0.05->0.00`, sensitivity=1.84952218e-1
- W_at_t10: 最大区間 `0.05->0.00`, sensitivity=3.12254216e-1
- usable_fraction_at_t10: 最大区間 `0.05->0.00`, sensitivity=3.54924453e0
- E_time_area: 最大区間 `0.05->0.00`, sensitivity=5.60883556e-1
- W_time_area: 最大区間 `0.05->0.00`, sensitivity=8.90681418e-1

## 15. 離散曲率

絶対二階差分最大は `usable_fraction_at_t10` の中心gamma=0.30（隣接区間 0.40->0.30）で 4.46831663e-2。追加判断で使った正規化二階差分最大は `W_at_t10` のgamma=0.10で 0.35411918。相転移や閾値現象とは断定しない。

## 16. 単調性

| metric | status |
|---|---|
| E_at_t10 | monotonic_nondecreasing |
| W_at_t10 | monotonic_nondecreasing |
| usable_fraction_at_t10 | monotonic_nondecreasing |
| W_max | monotonic_nondecreasing |
| E_time_area | monotonic_nondecreasing |
| W_time_area | monotonic_nondecreasing |

## 17. 小さな雑音低減の効果

| metric | recovery at gamma=0.4 | recovery at gamma=0.3 |
|---|---:|---:|
| W_at_t10 | 2.7388% | 7.7980% |
| usable_fraction_at_t10 | 6.3802% | 16.2109% |
| W_time_area | 4.6965% | 12.3360% |

## 18. 回復水準到達点

補間は使わず、離散点で初めて超えた最大gammaを保存した。

| metric | target | gamma | observed recovery |
|---|---:|---:|---:|
| W_at_t10 | 10% | 0.20 | 0.18082514 |
| W_at_t10 | 25% | 0.15 | 0.27232702 |
| W_at_t10 | 50% | 0.05 | 0.63675211 |
| W_at_t10 | 75% | 0.00 | 1.00000000 |
| W_at_t10 | 90% | 0.00 | 1.00000000 |
| usable_fraction_at_t10 | 10% | 0.30 | 0.16210883 |
| usable_fraction_at_t10 | 25% | 0.20 | 0.31972173 |
| usable_fraction_at_t10 | 50% | 0.10 | 0.57844523 |
| usable_fraction_at_t10 | 75% | 0.05 | 0.76446289 |
| usable_fraction_at_t10 | 90% | 0.00 | 1.00000000 |
| W_time_area | 10% | 0.30 | 0.12335983 |
| W_time_area | 25% | 0.20 | 0.25554933 |
| W_time_area | 50% | 0.10 | 0.50168169 |
| W_time_area | 75% | 0.00 | 1.00000000 |
| W_time_area | 90% | 0.00 | 1.00000000 |

## 19. 時間窓別感度

- pulse_interval: W_time_area最大感度区間 0.05->0.00, 9.29827024e-4
- early_post_pulse: W_time_area最大感度区間 0.05->0.00, 3.73815113e-2
- middle_interval: W_time_area最大感度区間 0.05->0.00, 2.49957005e-1
- late_interval: W_time_area最大感度区間 0.05->0.00, 6.02413075e-1

## 20. 時系列感度

各保存時刻でE、W、usable fractionの絶対回復、正規化回復、隣接感度、最大感度gamma区間をtimeseries CSVへ保存した。初期の小分母はNaNのままとした。

## 21. 刻み幅整合性

最終sweep全点とnoise-freeを半減刻みで再計算し、要約値、正規化回復、感度、単調性、回復水準順序、最大感度区間、追加点判断の全行がPASS。

## 22. 直接確認できたこと

固定模型・Omega・中央gamma内で、両端gamma低減に対する回復曲線、離散単調性、有限差分感度、回復水準、時間窓差を確認した。

## 23. 確認できていないこと

現実的実装、cost、異なる中央gamma/Omega/network、長時間・連続運転、因果機構、新規性は未確認。

## 24. 主張してはいけないこと

実装可能な必要保護強度、物理的感受率・臨界指数・相転移、一般最適、量子優位、実用送電性能は主張しない。

## 25. load雑音を扱わない理由

loadは3準位でsiteの二準位sigma_zを流用できず、同じgammaが公平な強度を意味しないため別検証が必要。

## 26. 次段階への判断材料

Wの最大離散感度区間は `0.05->0.00`。追加点判定の最大感度比は 21.42020012。これは候補情報であり、自動的に次Milestoneへ進まない。

## 27. 生成ファイル一覧

- `src/bin/partial_end_protection.rs`
- `partial_end_protection_timeseries.csv`
- `partial_end_protection_summary.csv`
- `partial_end_protection_recovery.csv`
- `partial_end_protection_sensitivity.csv`
- `partial_end_protection_thresholds.csv`
- `partial_end_protection_windows.csv`
- `partial_end_protection_checks.csv`
- `partial_end_protection_convergence.csv`
- `MILESTONE_7D_REPORT.md`

