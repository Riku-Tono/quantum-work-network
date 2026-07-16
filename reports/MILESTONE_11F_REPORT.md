# Milestone 11f: 局所線形補間点の単発評価と等入力matching判定

## 1. 目的

Omega_trial=0.18748395731510084を1本だけ評価した。反復二分法は実行していない。

## 2. 新規軌道

N=7、TOTAL_GAMMA=1.5、dt=0.0025、T=10、4000 steps、1001保存点の1本だけ。

## 3. reference target

11c正式値 target_E_drive_in=5.9618618770136536e-2。reference軌道は再計算していない。

## 4. 旧bracket

[0.1870, 1.8770762191709489e-1]。端点は11e・11dから読み込んだ。

## 5. trial入力会計

E_drive_in=5.9618450901925489e-2、E_drive_out=0.0000000000000000e0、E_drive_net=5.9618450901925489e-2、identity residual=0.000e0。

## 6. matching残差

F_trial=-1.6786821104702865e-7、absolute mismatch=1.6786821104702865e-7、relative mismatch=2.8157011099880636e-6、matching tolerance passed=true。

## 7. matching判定

**matched_input_found_with_fallback_diagnostic**

更新bracket=[NaN, NaN]、次補間候補=NaN、中点=NaN。matching成功時はこれらを生成していない。

## 8. 数値品質

checks=PASS、trace max=1.776e-15、Hermiticity max=0.000e0、worst eigenvalue=-5.061e-18、primary success/failure=1000/1, fallback success/attempt=1/1, solver failure=0、ledger max=5.197e-7。

## 9. matching成功時の等入力比較

matching成功時のみ `equal_input_N3_vs_N7_comparison.csv` にN=3正式値とN=7 trial値を記述的に保存した。

## 10. 入力正規化診断

W(t10)/Ein=3.9022615843454733e-2、Wmax/Ein=5.6781920195615787e-2、W time area/Ein=2.0153131959434209e-1、Eload(t10)/Ein=1.0575057241498595e-1。これらは装置効率・総合効率ではない。

## 11. 直接確認できたこと

この1条件の入力matching、保存物理量、数値品質だけを直接確認した。

## 12. 確認できていないこと

唯一root、広域唯一解、全Omega単調性、matching条件のdt収束、N=5、TOTAL_GAMMA=3.0、N>7。

## 13. 主張してはいけないこと

装置効率、送電効率、総合効率、実機優位、量子優位、長い鎖の一般的優位は主張しない。

## 14. 最終判定

**matched_input_found_with_fallback_diagnostic**

## 15. 次段階

等入力比較結果を確認し、matched条件のdt半減検証が必要か判断する。

自動実行していない。

## 16. 実行記録

fmt PASS、release tests 119 passed / 0 failed / 1 ignored、single trial PASS。実測 2348.114s（propagation 1930.718s、diagnostics 413.804s）。追加Omega、grid、二分法、dt半減は未実行。
