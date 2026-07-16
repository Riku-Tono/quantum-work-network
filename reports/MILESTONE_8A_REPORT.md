# Milestone 8a: Chain-length reachability and usable-work degradation

## 1. 目的

N=3からN=5へchain長だけを変え、load energy、ergotropy、usable fractionを比較した。

## 2. 今回の問い

N=5の有限時間到達、W最大値、ピーク遅延、共通時刻と個別ピーク、全site雑音損失を評価した。

## 3. N=3からN=5への模型一般化

既存N=3 APIを保存し、演算子構築とcoherent-drive実行にchain_length引数の追加APIを設けた。dim=2^N*3、drive site=0、load coupling site=N-1である。

## 4. 変更していない物理条件

局所周波数1、J=1、g=0.25、3準位load、Omega=0.2、drive周波数1、tau=3.2、真空初期状態、RK4、load無雑音を固定した。

## 5. chain長比較で同時に変わるもの

伝播距離、bond数、Hilbert次元、all-site noisyではcollapse operator数と総雑音寄与が同時に変わる。等総散逸比較ではない。

## 6. N=3回帰確認

noise-freeとall-site noisyのt=10基準値を絶対誤差2e-9以内で再現した。

## 7. N=5 noise-free到達試験

max E>1e-5、max W>1e-6を満たし、loadへの数値的到達を確認した。

## 8. 共通観測時間の決定

N=5 noise-freeのW/Eピークが終端から必要な余白を持つ最小候補を採用した。

## 9. 数値手法

dense complex matrix、Lindblad master equation、time-dependent RK4。基準dt=0.0025、保存0.01。

## 10. 数値品質チェック

全チェックPASS。

## 11. 4条件の共通時刻比較

|condition|E(t10)|W(t10)|usable(t10)|W/Ein(t10)|E(tmax)|W(tmax)|
|---|---:|---:|---:|---:|---:|---:|
|N3_noise_free|5.4450767878e-2|5.2798274942e-2|9.6965161374e-1|7.0102397372e-1|5.4450767878e-2|5.2798274942e-2|
|N3_all_site_noisy|1.2596874861e-2|2.3652476826e-3|1.8776464073e-1|3.9672970145e-2|1.2596874861e-2|2.3652476826e-3|
|N5_noise_free|2.2157813813e-2|2.1085865709e-2|9.5162211793e-1|2.8096700365e-1|2.2157813813e-2|2.1085865709e-2|
|N5_all_site_noisy|6.1088332712e-3|8.0518142592e-4|1.3180608967e-1|1.3531676843e-2|6.1088332712e-3|8.0518142592e-4|

## 12. 各条件のW最大値比較

|condition|W_max|t_at_W_max|E_at_W_max|usable_at_W_max|
|---|---:|---:|---:|---:|
|N3_noise_free|5.5424064638e-2|9.48|5.6854985097e-2|9.7483210211e-1|
|N3_all_site_noisy|3.0302005931e-3|5.63|7.5972454924e-3|3.9885516351e-1|
|N5_noise_free|2.1974463817e-2|6.63|2.3377531195e-2|9.3998222623e-1|
|N5_all_site_noisy|1.0968175558e-3|6.58|3.8915410559e-3|2.8184658470e-1|

## 13. load energy比較

N=5/N=3のE(t10)比はfree `4.0693299059e-1`、noisy `4.8494831763e-1`。

## 14. usable fraction比較

t=10のN=3/N=5はfree `9.6965161374e-1` / `9.5162211793e-1`、noisy `1.8776464073e-1` / `1.3180608967e-1`。

## 15. W/Ein比較

t=10のN=3/N=5はfree `7.0102397372e-1` / `2.8096700365e-1`、noisy `3.9672970145e-2` / `1.3531676843e-2`。制御費用を含む総合効率ではない。

## 16. 到達時刻とピーク遅延

N=5 freeのenergy/W閾値到達時刻は `3.18` / `2.66`。Wピーク遅延N5-N3はfree `-2.85`、noisy `0.95`。閾値は数値診断用である。

## 17. N=3からN=5への有限差

Wmax比N5/N3はfree `3.9647874908e-1`、noisy `3.6196202928e-1`。2点だけなので減衰則や距離scalingを推定しない。

## 18. Nごとの全site雑音損失

W(t10) noisy/freeはN=3 `4.4797821239e-2`、N=5 `3.8185836760e-2`。Wmax noisy/freeはN=3 `5.4673012759e-2`、N=5 `4.9913279567e-2`。後者は異なるピーク時刻の比較である。

## 19. 時間窓別比較

固定窓と延長窓のmean/area/maxを `chain_length_reachability_windows.csv` に保存した。時間面積は状態量の面積であり累積仕事ではない。
## 20. site populationの時間発展

可変長long形式を `chain_length_site_populations.csv` に保存した。
## 21. 刻み幅整合性

基準dtと半減dtの主要量比較は **PASS**。N=3 free、N=5 free、N=5 noisyを半減刻みで再計算し、N=5到達、freeのN=3/N=5 Wmax順位、N=5のfree/noisy順位、freeのpeak-delay符号を確認した。N=3 noisyの半減計算は行っていないため、noisy条件のN=3/N=5 peak-delay符号は半減刻みで独立確認していない。

## 22. 計算負荷

|condition|dt|dim|saved points|seconds|
|---|---:|---:|---:|---:|
|N3_noise_free|0.0025|24|1001|1.115|
|N3_all_site_noisy|0.0025|24|1001|5.304|
|N5_noise_free|0.0025|96|1001|43.849|
|N5_all_site_noisy|0.0025|96|1001|403.443|
|N3_noise_free|0.00125|24|1001|1.127|
|N5_noise_free|0.00125|96|1001|62.823|
|N5_all_site_noisy|0.00125|96|1001|819.015|

## 23. 直接確認できたこと

- N=5 noise-freeで有限時間内にload ergotropyが生成された。
- WmaxのN5/N3比はfreeで0.39648、noisyで0.36196だった。
- t=10と個別Wピークの両方で、N=5のWはN=3より小さく、all-site noisyはfreeより小さい。同じ大小結論だが比率は同一ではない。
- N=5 freeは閾値到達がN=3 freeより遅い一方、W最大時刻は2.85早かった。振動する有限系では到達遅延と最大時刻遅延を同一視できない。
- siteごとのgammaを固定したall-site noise損失をN別に測定した。

## 24. 確認できていないこと

N>5、連続N sweep、位置別弱点、一般scaling、総雑音一致、保護費用は確認していない。N=3 noisyの半減刻みは未実行である。
## 25. 主張してはいけないこと

指数/べき減衰、熱力学極限、距離だけの純粋因果、実機効率、量子優位、新規性、N=5位置別弱点の予測。
## 26. 次段階への判断材料

N=5 freeの到達と数値品質・指定3条件の刻み幅整合性が合格したため、位置別雑音比較は次段階候補にできる。ただし今回は実装していない。
## 27. 生成ファイル一覧

- `src/bin/chain_length_reachability.rs`
- `chain_length_reachability_timeseries.csv`
- `chain_length_site_populations.csv`
- `chain_length_reachability_summary.csv`
- `chain_length_reachability_ratios.csv`
- `chain_length_reachability_arrivals.csv`
- `chain_length_reachability_windows.csv`
- `chain_length_reachability_checks.csv`
- `chain_length_reachability_convergence.csv`
- `chain_length_reachability_performance.csv`
- `MILESTONE_8A_REPORT.md`
