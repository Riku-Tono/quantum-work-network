# Milestone 9c diagnostic: N=7 t=0.02 minimum-eigenvalue NaN

## 1. 目的

Milestone 9cのN=7 t=0.02 minimum-eigenvalue NaNを、短時間再計算だけで切り分けた。

## 2. Milestone 9cで発生した症状

t=0.02のdirect SymmetricEigen 1点だけがNaNとなり、positivityとfinite_valuesがFAILした。

## 3. 再計算範囲

t=0から0.03まで12 RK4 stepsを3回実行し、t=0,0.01,0.02,0.03を保存した。t=10本計算は再実行していない。

## 4. 変更していない模型

N=7、total gamma=1.5、gamma_site=1.5/7、all-site noise、Omega=0.2、tau=3.2、J=1、g=0.25、dt=0.0025、RK4、exact dephasing kernel、真空初期状態を維持した。

## 5. t=0.01結果

保存済み9c CSVと8主要値は差0で一致した。directおよびHermitianized SymmetricEigenには一部非有限固有値があったが、有限値だけから作られたminimumは -2.808e-22 だったため、9cのminimum列では異常が表面化していなかった。Schurは全固有値有限だった。

## 6. t=0.02結果

rho finite=true、trace=9.9999999999999989e-1+0.0000000000000000e0i、Hermiticity=0.000e0。

## 7. t=0.03結果

保存済み9c CSVと8主要値を照合した。

## 8. density matrix有限性

t=0.02の全要素finite=true、Frobenius norm=9.9999999999999956e-1、max abs element=9.9999999999973599e-1。

## 9. traceとHermiticity

trace error=1.110e-16、Hermiticity error=0.000e0、Hermitian化補正norm=0.000e0。

## 10. direct SymmetricEigen

raw出力には非有限固有値があり、minimum集約結果=-Inf、finite=false。9cのCSV formatterは非有限値を一律 `NaN` と記録するため、元CSVのNaNに対応する現象は3回すべてで再現した。

## 11. Hermitianized SymmetricEigen

minimum=-Inf、finite=false、eigenvalue sum-trace差=NaN。

## 12. independent solver

Complex Schur minimum=-3.2372313550345495e-24、max eigenvalue imag=1.436e-30、sum-trace差=9.992e-16。

## 13. solver間比較

HermitianizedとSchurのminimum差=inf、許容値=1.0e-10。

## 14. 保存済み9c結果との一致

24 scalar比較の最大絶対差=0.000e0、許容値=1.0e-12。9cにはfull density matrixが保存されていないため、density matrix要素のrun-to-saved直接比較はできない。代わりに3回の決定論的rho要約を比較した。

## 15. 3回再現性

t=0.02の決定論的rho要約一致=true。

## 16. minimum eigenvalue判定

独立solver minimum=-3.2372313550345495e-24、基準 >= -1.0e-8。

## 17. 直接確認できたこと

短時間状態の有限性、trace、Hermiticity、3 solver、保存済み主要値、3回再現性を直接確認した。状態側の検査とSchur結果は正常だった一方、rawとHermitianizedの両SymmetricEigenで非有限固有値が再現した。

## 18. 確認できていないこと

t=0.03より後のdensity matrix再構築、別dt、別gamma、別solver crate、N=7 t=10再計算は行っていない。

## 19. Milestone 9c判定の更新可否

更新可否=不可。元レポートは上書きせず、この補足だけを追加した。

## 20. 最終判定

指定判定規則による診断判定は **state_level_numerical_issue**。これは「Hermitianized SymmetricEigenも非有限ならCase B」という停止規則によるラベルであり、rho自体にNaN、trace異常、Hermiticity異常、非決定性が見つかったという意味ではない。Milestone 9c更新後判定は **numerical_issue_stop** のまま。全必須チェック=false。

## 21. 生成ファイル一覧

- `n7_t002_eigen_diagnostic.csv`
- `n7_t002_state_summary.csv`
- `n7_t002_reproducibility.csv`
- `n7_t002_saved_value_comparison.csv`
- `n7_t002_diagnostic_checks.csv`
- `MILESTONE_9C_DIAGNOSTIC.md`
