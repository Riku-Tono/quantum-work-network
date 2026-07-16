# Milestone 7c: Ideal partial protection upper-bound test

## 1. 目的

全3site雑音から指定siteの位相雑音演算子だけを理想的に除去した反実仮想条件で、load状態量の回復上限を比較した。

## 2. Milestone 7aとの問いの違い

7aは1siteだけへ雑音を置く有害配置比較。7cは全site雑音から雑音項を選択的に除く回復比較であり、順位一致を仮定しない。

## 3. 理想保護の定義

指定siteの `sqrt(gamma_phi/2) sigma_z` collapse operatorを完全に除去する。現実的装置、制御pulse、cost、誤り訂正、有限精度保護ではない。

## 4. 比較条件

`all_noisy=[0,1,2]`, `protect_entrance=[1,2]`, `protect_exit=[0,1]`, `protect_both_ends=[1]`, `noise_free=[]`。protect_both_endsでは中央雑音だけが残る。

## 5. 変更していない物理条件

3site+3準位load、Hamiltonian、真空初期状態、J=1、g=0.25、各周波数=1、tau=3.2、Omega=0.2、gamma_phi=0.5、t_max=10、pulse、load、ergotropy、RK4を固定した。

## 6. 数値手法

基準dt=0.0025、半減dt=0.00125、保存間隔0.01。状態量時間面積とpowerは台形則。比のepsilon=1e-8。非加算性の約0判定はabs<=1e-7またはnormalized abs<=5e-3。

## 7. 数値品質チェック

全75項目PASS。collapse数は3/2/2/1/0、site mapping、trace、Hermiticity、positivity、population、load縮約、top-level、ledger、有限性、共通grid/初期状態、W<=E、usable範囲、drive整合、既存baselineを確認した。

## 8. t=10の5条件比較

| condition | E | W | usable fraction | W/Ein |
|---|---:|---:|---:|---:|
| all_noisy | 1.2596874861e-2 | 2.3652476826e-3 | 1.8776464073e-1 | 3.9672970145e-2 |
| protect_entrance | 2.2717839144e-2 | 9.0447521802e-3 | 3.9813435260e-1 | 1.1520713324e-1 |
| protect_exit | 1.9779100573e-2 | 7.2883177108e-3 | 3.6848580066e-1 | 1.2166726969e-1 |
| protect_both_ends | 4.8178995305e-2 | 4.5346118640e-2 | 9.4120100165e-1 | 5.7488202578e-1 |
| noise_free | 5.4450767878e-2 | 5.2798274942e-2 | 9.6965161374e-1 | 7.0102397372e-1 |

## 9. all_noisyからの回復量

| condition | metric | point | absolute | normalized |
|---|---|---|---:|---:|
| protect_entrance | E | t10 | 1.01209643e-2 | 0.24181656 |
| protect_entrance | E | time_area | 3.85130623e-2 | 0.26146431 |
| protect_entrance | W | t10 | 6.67950450e-3 | 0.13244306 |
| protect_entrance | W | W_max | 6.48247650e-3 | 0.12372587 |
| protect_entrance | W | time_area | 3.29629852e-2 | 0.18544293 |
| protect_entrance | usable_fraction | t10 | 2.10369712e-1 | 0.26905386 |
| protect_entrance | usable_fraction | at_W_max | 5.09988235e-2 | 0.08854317 |
| protect_exit | E | t10 | 7.18222571e-3 | 0.17160233 |
| protect_exit | E | time_area | 3.52582321e-2 | 0.23936734 |
| protect_exit | W | t10 | 4.92307003e-3 | 0.09761599 |
| protect_exit | W | W_max | 5.10317428e-3 | 0.09740023 |
| protect_exit | W | time_area | 2.83874420e-2 | 0.15970187 |
| protect_exit | usable_fraction | t10 | 1.80721160e-1 | 0.23113463 |
| protect_exit | usable_fraction | at_W_max | 4.80075869e-2 | 0.08334984 |
| protect_both_ends | E | t10 | 3.55821204e-2 | 0.85015080 |
| protect_both_ends | E | time_area | 1.26806771e-1 | 0.86088832 |
| protect_both_ends | W | t10 | 4.29808710e-2 | 0.85223659 |
| protect_both_ends | W | W_max | 4.23159180e-2 | 0.80765026 |
| protect_both_ends | W | time_area | 1.50933830e-1 | 0.84912246 |
| protect_both_ends | usable_fraction | t10 | 7.53436361e-1 | 0.96361288 |
| protect_both_ends | usable_fraction | at_W_max | 5.42345838e-1 | 0.94161034 |

## 10. noise-freeまでの残留損失

| condition | metric | t=10 residual loss |
|---|---|---:|
| protect_entrance | E | 0.58278202 |
| protect_entrance | W | 0.82869228 |
| protect_entrance | usable_fraction | 0.58940474 |
| protect_exit | E | 0.63675259 |
| protect_exit | W | 0.86195917 |
| protect_exit | usable_fraction | 0.61998124 |
| protect_both_ends | E | 0.11518244 |
| protect_both_ends | W | 0.14114394 |
| protect_both_ends | usable_fraction | 0.02934107 |

符号付きで保存した。負値はnoise-free超過であり0へ丸めていない。

## 11. W_max比較

W_max回復最大: `protect_both_ends`。

## 12. E_time_areaとW_time_area比較

E_time_area回復最大: `protect_both_ends`。W_time_area回復最大: `protect_both_ends`。両者は状態量の時間面積で、累積流入/抽出ではない。

## 13. 時間窓別回復

- pulse_interval: W_time_area回復最大 `protect_both_ends` (3.05833350e-4)
- early_post_pulse: W_time_area回復最大 `protect_both_ends` (9.76660121e-3)
- middle_interval: W_time_area回復最大 `protect_both_ends` (4.92226516e-2)
- late_interval: W_time_area回復最大 `protect_both_ends` (9.16387437e-2)

## 14. 入口保護と出口保護の順位変化

- E: t=1.80: tie->protect_entrance; t=3.89: protect_entrance->tie; t=4.27: tie->protect_exit; t=7.20: protect_exit->tie; t=7.33: tie->protect_entrance
- W: t=4.59: tie->protect_entrance
- usable fraction: t=3.97: protect_exit->tie; t=4.25: tie->protect_entrance

5点持続した変化だけを確定した。

## 15. 両端保護の結果

中央雑音だけが残る。t=10 W=4.5346118640e-2, usable=9.4120100165e-1。noise-freeへのW残留損失=0.14114394、all-noisyからのW回復率=0.85223659。

| t=10 metric | both - entrance | both - exit |
|---|---:|---:|
| E | 2.54611562e-2 | 2.83998947e-2 |
| W | 3.63013665e-2 | 3.80578009e-2 |
| usable fraction | 5.43066649e-1 | 5.72715201e-1 |

中央雑音が無害/保護不要とは言わない。

## 16. 回復の非加算性

| metric | point | synergy | normalized | class |
|---|---|---:|---:|---|
| E | t10 | 1.82789304e-2 | 0.43673191 | positive_nonadditivity |
| W | t10 | 3.13782964e-2 | 0.62217753 | positive_nonadditivity |
| usable_fraction | t10 | 3.62345489e-1 | 0.46342438 | positive_nonadditivity |
| E | time_area | 5.30354765e-2 | 0.36005666 | positive_nonadditivity |
| W | time_area | 8.95834026e-2 | 0.50397767 | positive_nonadditivity |

Synergyは相互作用エネルギーではなく、選択的雑音除去に対する観測量応答の非加算性。

## 17. 刻み幅整合性

5条件×6量、t=10 synergy E/W/use、5種類の保護順位がPASS。順位とsynergy符号は不変。基準/半減刻みのラベルもconvergence CSVへ保存した。

## 18. 直接確認できたこと

固定条件で全3site雑音から入口/出口/両端雑音を理想除去したときのt=10、最大値、時間面積、時間窓、非加算性、中央雑音残存時の残留損失を確認した。

## 19. 確認できていないこと

現実的実装、cost、不完全保護、他gamma/Omega/network、長時間/反復/定常、因果機構、新規性は未確認。

## 20. 主張してはいけないこと

保護の一般則・実用最適、synergyを物理相互作用エネルギーとする表現、中央雑音が無害、量子優位、現実送電性能は主張しない。

## 21. load雑音を今回扱わない理由

loadは3準位でsiteの二準位sigma_zを流用できず、同じgammaが公平な強度を意味しない。演算子、decay rate、規格化、強度対応を別途検証する必要がある。

## 22. 次段階への判断材料

t=10 W回復最大 `protect_both_ends`、usable fraction回復最大 `protect_both_ends`、E回復最大 `protect_both_ends`。これは次の比較候補であり、自動的に次Milestoneへ進まない。

## 23. 生成ファイル一覧

- `src/bin/ideal_partial_protection.rs`
- `ideal_partial_protection_timeseries.csv`
- `ideal_partial_protection_summary.csv`
- `ideal_partial_protection_recovery.csv`
- `ideal_partial_protection_synergy.csv`
- `ideal_partial_protection_windows.csv`
- `ideal_partial_protection_rankings.csv`
- `ideal_partial_protection_checks.csv`
- `ideal_partial_protection_convergence.csv`
- `MILESTONE_7C_REPORT.md`

