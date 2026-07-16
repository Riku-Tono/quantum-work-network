# Milestone 7b: Local-noise damage onset analysis

## 1. 目的

Milestone 7aの固定済み時系列から、局所位相雑音条件とnoise-freeとの差がいつ持続的に現れたかを記述する。因果分解は行わない。

## 2. 使用した既存データ

主要入力は `local_noise_placement_timeseries.csv` のみ。`local_noise_placement_summary.csv` と `local_noise_placement_ratios.csv` はt=10、W_time_area、比の整合性確認だけに使用した。各条件1001点、合計4004行。

## 3. 新規時間発展を行っていないこと

このbinはCSV reader、算術、台形積分、順位判定、CSV/Markdown writerだけを含む。Hamiltonian、Lindblad、RK4、collapse operator生成、量子状態伝播は呼び出していない。

## 4. onsetの定義

基準量がminimum_reference以上で、`reference - noisy`が絶対閾値以上、かつ相対損失が指定閾値以上となる状態が連続5点続いた最初の点をsustained onsetとした。単一点の超過は採用しない。

## 5. 閾値と持続条件

推奨値を変更せず採用した。E/Wの絶対閾値`1e-5`、usable fractionの絶対閾値`1e-3`。minimum referenceはE/W `1e-4`、usable fraction `1e-3`。相対閾値はweak `1%`、medium `5%`、strong `10%`。時間幅0.01の5点持続を約0.05時間として記録した。

## 6. 数値品質チェック

全40項目がPASS。必須列、4条件、点数、時間グリッド、単調性、重複、ラベル、定数、有限値/許容NaN、population、E/W/useの範囲、分母処理、summary/ratios整合を確認した。

## 7. entrance条件の被害開始

| quantity | weak | medium | strong |
|---|---:|---:|---:|
| E | 2.2500000000000000e0 | 2.2500000000000000e0 | 2.2500000000000000e0 |
| W | 2.2500000000000000e0 | 2.2500000000000000e0 | 2.2500000000000000e0 |
| usable_fraction | 9.7999999999999998e-1 | 1.5700000000000001e0 | 2.9300000000000002e0 |

時刻と同時点のpopulationの対応であり、特定過程だけの損傷や因果原因を意味しない。

## 8. middle条件の被害開始

| quantity | weak | medium | strong |
|---|---:|---:|---:|
| E | 2.2500000000000000e0 | 2.2500000000000000e0 | 2.2500000000000000e0 |
| W | 2.2500000000000000e0 | 2.2500000000000000e0 | 2.2500000000000000e0 |
| usable_fraction | 9.7999999999999998e-1 | 1.6200000000000001e0 | NaN |

時刻と同時点のpopulationの対応であり、特定過程だけの損傷や因果原因を意味しない。

## 9. exit条件の被害開始

| quantity | weak | medium | strong |
|---|---:|---:|---:|
| E | 2.2500000000000000e0 | 2.2500000000000000e0 | 2.2500000000000000e0 |
| W | 2.2500000000000000e0 | 2.2500000000000000e0 | 2.2500000000000000e0 |
| usable_fraction | 9.7999999999999998e-1 | 1.2700000000000000e0 | 2.5400000000000000e0 |

時刻と同時点のpopulationの対応であり、特定過程だけの損傷や因果原因を意味しない。

## 10. E、W、usable fractionの開始時刻比較

medium閾値での持続onsetをまとめる。

| condition | E | W | usable fraction |
|---|---:|---:|---:|
| noise_entrance | 2.2500000000000000e0 | 2.2500000000000000e0 | 1.5700000000000001e0 |
| noise_middle | 2.2500000000000000e0 | 2.2500000000000000e0 | 1.6200000000000001e0 |
| noise_exit | 2.2500000000000000e0 | 2.2500000000000000e0 | 1.2700000000000000e0 |

EとWは3条件・3閾値すべてでt=2.25となった。これは同じminimum referenceと絶対閾値を使った診断上の同時開始であり、energy lossとquality lossの物理機構が同時だという因果主張ではない。usable fractionは閾値依存性があり、weakでは3条件ともt=0.98、mediumではexit 1.27、entrance 1.57、middle 1.62、strongではexit 2.54、entrance 2.93、middleは未検出だった。entranceとexitの順序反転はなく、medium/strongではexitが早いが、単一の正確な開始時刻へ圧縮しない。

## 11. 最大損失時刻

| condition | E loss max | W loss max | use loss max |
|---|---:|---:|---:|
| noise_entrance | 9.20 | 10.00 | 0.98 |
| noise_middle | 3.92 | 3.87 | 3.70 |
| noise_exit | 9.75 | 10.00 | 0.98 |

最大値の詳細と同時点のpopulationは `local_noise_damage_extrema.csv` に保存した。

## 12. onset時のsite population

medium W onsetでの雑音site populationを示す。

| condition | onset | noisy-site population | chain内比 | dominant site |
|---|---:|---:|---:|---|
| noise_entrance | 2.2500000000000000e0 | 2.7699954847829177e-2 | 5.8401015613112817e-1 | site1 |
| noise_middle | 2.2500000000000000e0 | 1.6476035301427060e-2 | 2.6687137699981156e-1 | site1 |
| noise_exit | 2.2500000000000000e0 | 6.5886858113160841e-3 | 1.1008522248295555e-1 | site1 |

これは励起分布の同時記録であり、populationが被害の原因だとは示さない。

## 13. 時間窓別比較

各窓の平均と時間面積は `local_noise_damage_windows.csv` に保存した。W_time_areaが最小の条件を窓ごとに示す。

- pulse_interval: `noise_exit` (W_time_area=4.99993445e-4)
- early_post_pulse: `noise_entrance` (W_time_area=9.91809215e-3)
- middle_interval: `noise_entrance` (W_time_area=2.13576986e-2)
- late_interval: `noise_entrance` (W_time_area=2.44692058e-2)

E/W_time_areaは状態量の時間面積であり、累積流入エネルギーや累積抽出仕事ではない。

## 14. 順位の時間変化

- worst Eの確定切替: t=1.78: noise_entrance+noise_middle+noise_exit -> noise_entrance+noise_middle; t=1.85: noise_entrance+noise_middle -> noise_entrance; t=3.90: noise_entrance -> noise_entrance+noise_exit; t=4.31: noise_entrance+noise_exit -> noise_exit; t=6.44: noise_exit -> noise_entrance+noise_exit; t=6.51: noise_entrance+noise_exit -> noise_entrance
- worst Wの確定切替: t=1.82: noise_entrance+noise_middle+noise_exit -> noise_entrance+noise_exit; t=4.75: noise_entrance+noise_exit -> noise_entrance
- worst usable fractionの確定切替: t=4.00: noise_exit -> noise_entrance+noise_exit; t=4.31: noise_entrance+noise_exit -> noise_entrance; t=6.80: noise_entrance -> noise_entrance+noise_exit; t=6.86: noise_entrance+noise_exit -> noise_exit
- entrance/exit Wの入替: t=4.75: tie -> entrance_worse
- entrance/exit usable fractionの入替: t=4.00: exit_worse -> tie; t=4.31: tie -> entrance_worse; t=6.80: entrance_worse -> tie; t=6.86: tie -> exit_worse


Tieを含めて各条件が最悪集合に入った時間割合を示す（tieのため合計は100%を超えうる）。

| condition | E | W | usable fraction |
|---|---:|---:|---:|
| noise_entrance | 0.7872 | 1.0000 | 0.3178 |
| noise_middle | 0.1848 | 0.1818 | 0.0000 |
| noise_exit | 0.4386 | 0.4745 | 0.7243 |
短い反転を除くため、同じ順位状態が5点続いた時だけ切替とした。全時刻の順位は `local_noise_damage_rankings.csv`。

## 15. middle条件の被害が小さいことへの手がかり

| condition | max noisy-site p | p time-area | pulse mean p | post-pulse mean p | p at W onset | p at max W loss |
|---|---:|---:|---:|---:|---:|---:|
| noise_entrance | 2.773840e-2 | 1.277890e-1 | 1.223039e-2 | 1.302873e-2 | 2.769995e-2 | 1.318198e-2 |
| noise_middle | 2.366341e-2 | 1.215187e-1 | 8.502065e-3 | 1.386725e-2 | 1.647604e-2 | 8.197671e-3 |
| noise_exit | 5.416329e-2 | 1.769388e-1 | 7.084881e-3 | 2.268701e-2 | 6.588686e-3 | 1.507433e-2 |

middleの雑音site population時間面積はentranceに近く、exitより小さかった。pulse中平均はentranceより小さくexitより大きい3条件中間だった一方、7aのW損失はmiddleが明らかに軽かった。したがって、少なくとも単純な最大populationまたは時間面積だけではW損失差を説明し切れない。時刻のずれや未計算の量も候補として残る。これは手がかりであり原因説明ではない。

## 16. 直接確認できたこと

固定された7aデータ内で、E/W/usable fractionの差が各閾値で持続的に現れた時刻、その時刻のsite population、時間窓別の状態量、順位切替を直接確認した。

## 17. 確認できていないこと

特定物理過程への因果分解、populationだけによる説明、coherence current、site間current、保護による回復、他パラメータ・長いnetworkへの一般化は確認していない。

## 18. 主張してはいけないこと

entrance/middle/exit雑音が注入/輸送/受け渡しだけを壊すという断定、populationが大きいことを原因とする断定、先に変化した量を原因とする断定、保護効果の予測はできない。

## 19. 次の保護実験を選ぶための判断材料

- entrance protection候補: 7aでt=10のW、W_max、W_time_areaが最小であり、今回の時系列でも入口雑音条件の持続損失が確認された。
- exit protection候補: 7aでt=10のusable fractionが最小であり、今回その開始時刻と順位変化を確認した。
- 同時比較候補: entranceは仕事量、exitは使える割合で異なる悪化指標を持つため、同じ固定条件で並べる価値がある。

これらは候補選定の材料であり、保護による回復を予測するものではない。保護機能は実装していない。

## 20. 生成ファイル一覧

- `src/bin/local_noise_damage_analysis.rs`
- `local_noise_damage_timeseries.csv`
- `local_noise_damage_onsets.csv`
- `local_noise_damage_extrema.csv`
- `local_noise_damage_windows.csv`
- `local_noise_damage_rankings.csv`
- `local_noise_damage_checks.csv`
- `MILESTONE_7B_REPORT.md`

