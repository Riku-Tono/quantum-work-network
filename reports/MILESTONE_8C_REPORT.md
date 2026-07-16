# Milestone 8c: Exact dephasing-kernel optimization and equivalence validation

## 1. 目的

局所sigma_z位相雑音のLindblad項だけを、物理模型を変えずに厳密高速化した。

## 2. Milestone 8bでのボトルネック

N=7 all-site noisyは旧dense collapse実装で1 step約21.317895秒、t=10推定約23.69時間だった。

## 3. 物理模型を変更していないこと

Hamiltonian、drive、RK4、dt、load、初期状態、gamma、ergotropy、Hilbert空間、dense density matrixを変更していない。

## 4. dephasing kernelの数式

各要素へ -Gamma[a,b] rho[a,b] を加え、Gammaはchain site bitが異なるsiteのgamma和として一度だけ構築した。これは同じLindblad dissipatorの厳密な成分表示である。

## 5. basis mapping

既存tensor順序 |q1,...,qN,load> をoperator対角要素と全basis pairで照合した。load levelはGammaへ含めない。

## 6. 実装

旧Dense pathを残し、DiagonalDephasingKernelを独立moduleとして追加した。kernelはf64連続配列で、N=7では1,179,648 bytes。

## 7. unit checks

Hamming rate、load除外、site順序、diagonal zero、gamma zero、trace、Hermiticity、入力検証を追加した。

## 8. dissipator単体一致

N=3/N=5のall-siteとsite-dependent partial条件でdense dissipatorと比較した。詳細はrhs equivalence CSV。

## 9. trajectory一致

N=3 all-site/partialはt=0.1,1,10、N=5 all-siteはt=0.1,1をdense再計算し、N=5 t=10はMilestone 8a既存値と照合した。

## 10. N=3回帰

noise-freeとall-site noisyのt=10主要量、W peak、時間面積をMilestone 8aへ回帰した。

## 11. N=5回帰

all-site noisyのE、W、usable fraction、W peak、時間面積、drive inputをMilestone 8aへ回帰した。

## 12. N=7 benchmark

最長probeについてwarmup 1回、measurement 3回をrelease buildで実行した。詳細はbenchmarks CSV。

## 13. speedup

N=7 noisyのmedianは `0.711565` s/step、旧値に対するspeedupは **29.96x**。

## 14. t=10実行時間再推定

N=7 noisyはmedianから約 **0.791時間**。分類は **feasible_candidate**。今回はt=10本計算を実行していない。

## 15. memory

N=7 Gamma kernelは1,179,648 bytes（1.125 MiB）。dephasing評価用の7個のcollapse行列はkernel hot pathで使用しない。既存operator bundleは診断・reference互換性のため維持した。

## 16. 数値品質チェック

checks CSVは **140 PASS / 0 FAIL**。許容値はRHS max abs 1e-12、trajectory主要量2e-9、trace/Hermiticity/positivity 1e-8、ledger 5e-5。

## 17. 直接確認できたこと

- sigma_z dephasingをelementwise kernelで再現した。
- N=3/N=5で旧dense path・既存値と許容内一致した。
- N=7 noisyの短時間速度とt=10再推定を測定した。

## 18. 確認できていないこと

N=7 t=10の最終物理量、長時間の実測性能、N scaling、他のLindblad operatorへの適用。

## 19. 主張してはいけないこと

物理近似、新現象、一般Lindbladへの普遍適用、N=7最終結果、scaling則、GPU/sparse法より優れるという主張。

## 20. N=7本計算へ進む判断

数値等価性、speedup>=5、推定<=6時間を満たす場合だけ **次段階候補** とする。実行承認ではなく、今回自動実行していない。現在判定: **feasible_candidate**。

## 21. 生成ファイル一覧

- `src/dephasing_kernel.rs`
- `src/bin/dephasing_kernel_benchmark.rs`
- `dephasing_kernel_unit_checks.csv`
- `dephasing_kernel_rhs_equivalence.csv`
- `dephasing_kernel_trajectory_equivalence.csv`
- `dephasing_kernel_regression.csv`
- `dephasing_kernel_benchmarks.csv`
- `dephasing_kernel_estimates.csv`
- `dephasing_kernel_checks.csv`
- `MILESTONE_8C_REPORT.md`

