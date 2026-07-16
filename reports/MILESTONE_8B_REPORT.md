# Milestone 8b: N=7 computational feasibility and short-time reachability probe

## 1. 目的

N=7のdense Lindblad RK4が現在環境で計算可能か、短時間probeだけで評価した。

## 2. 今回は本計算ではないこと

`t=10`本計算、半減刻み本計算、位置別雑音、最適化は実行していない。

## 3. N=7模型構成

7個の二準位siteと3準位load。J=1、g=0.25、各周波数1、drive site=0、load coupling site=6、Omega=0.2、tau=3.2、真空初期状態、load無雑音。

## 4. 次元とoperator mapping

Hilbert次元384、density matrix 384x384（147456要素）、bond数6、drive site 0、load coupling site 6。collapse数はfree 0、noisy 7。

## 5. construction-only結果

全construction検査PASS。詳細は `n7_feasibility_construction.csv` と checks CSV。

## 6. one-step結果

両条件でRK4 1 stepが完了し、trace、Hermiticity、positivity、finite、load縮約、ledger検査はPASS。個別時間は checks CSV。

## 7. short-time benchmark

|condition|longest probe|seconds/step|estimated t=10 hours|
|---|---:|---:|---:|
|N7_noise_free|0.100|1.100105|1.222|
|N7_all_site_noisy|0.100|21.317895|23.687|

`t=0.1` noisy の時点で推定が9時間を大きく超え、現行dense法の infeasible 基準が確定した。このため、判断を変えず約1時間超を追加消費する `t=0.5` と、条件付きの `t=1.0` は実行しなかった。

## 8. noise-free実行時間推定

dt=0.0025で約1.222時間。最長probeのpropagation step時間を使用し、startupは除外した。

## 9. all-site noisy実行時間推定

dt=0.0025で約23.687時間。最長probeのpropagation step時間を使用し、startupは除外した。

## 10. 半減刻み実行時間推定

未実行。step数が2倍としてfree約2.445時間、noisy約47.373時間と概算した。

## 11. メモリ推定

Complex64は16 bytes。1行列は2359296 bytes（約2.25 MiB）。保守的peakはfree約83.2 MiB、noisy約108.0 MiB。construction時の実測working setはfree約60.7 MiB、noisy約76.6 MiB、one-step後の実測peakはnoisy約106.6 MiBだった。実行時の利用可能物理メモリは約18.6 GiBで、メモリ不足は見られなかった。

## 12. N=5実測との比較

Hilbert次元比N7/N5=4、density element比=16。seconds/step比はfree 100.35、noisy 211.36。N=5は1001保存点の診断を含み、N=7 probeはpropagation主体なので厳密な同条件比較ではない。collapse数比noisy=7/5=1.4。2点から一般scalingは推定しない。

## 13. 数値品質チェック

全記録検査: **PASS**。

## 14. load初期応答

最長probe終点の値は steps CSV末尾に保存した。t<=1でEやWがほぼ0でも到達不能とは判定しない。今回確認するのは初期兆候だけである。

- N7_noise_free,2.5000000000000001e-3,40,1.0000000000000001e-1,4.4004186599999997e1,1.1001046649999999e0,7.0219357908988491e-36,0.0000000000000000e0,5.3077481825799702e-18,NaN,4.2135088241534798e-9,0.0000000000000000e0,0.0000000000000000e0,-3.1706867892292381e-18,-1.0220358751523007e-10
- N7_all_site_noisy,2.5000000000000001e-3,40,1.0000000000000001e-1,8.5271581509999999e2,2.1317895377500001e1,6.6760129276378057e-36,0.0000000000000000e0,5.1416566515921428e-18,NaN,4.1683294557653592e-9,0.0000000000000000e0,0.0000000000000000e0,-5.6887848173237804e-30,-1.0068939883314184e-10

## 15. 本計算へ進む基準

construction、one-step、short-time数値品質、メモリ、free<=1h、noisy<=6h、probe長で推定が大きく発散しないこと。

## 16. feasibility判定

**infeasible_with_current_dense_method**。これは現在環境・現在dense実装での計算可否分類であり、物理的可能性ではない。`t=10`本計算は自動実行していない。

## 17. 最適化候補

必要なら、一時行列削減、collapse共通部分再利用、sparse表現、正確性を保つ別積分法、Krylov、exponential integrator、並列行列積、BLAS/release設定確認を比較候補にする。今回は実装していない。

## 18. 直接確認できたこと

- N=7の構築とmapping。
- 両条件の1 stepおよび実行済み短時間probeの安定性。
- 現在実装でのstep時間、t=10概算、memory概算。

## 19. 確認できていないこと

N=7のt=10最終値、半減刻み実測、長時間安定性、最終到達、距離依存則、最適化後性能。

## 20. 主張してはいけないこと

N=7が到達不能、N=7の最終W低下、指数/べきscaling、物理的限界、量子優位、新規性。

## 21. 次段階への判断材料

現行dense法のまま`t=10`本計算へ進む候補とはしない。まず一時行列削減、collapse共通部分再利用、BLAS/release設定など、物理模型を変えない最適化候補を別Milestoneで比較する判断材料になった。今回は最適化を実装していない。

## 22. 生成ファイル一覧

- `src/bin/n7_feasibility_probe.rs`
- `n7_feasibility_construction.csv`
- `n7_feasibility_steps.csv`
- `n7_feasibility_benchmarks.csv`
- `n7_feasibility_estimates.csv`
- `n7_feasibility_memory.csv`
- `n7_feasibility_checks.csv`
- `MILESTONE_8B_REPORT.md`
