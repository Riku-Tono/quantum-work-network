# Milestone 9c validation: robust positivity diagnostic

## 1. 目的

N=7 fixed-total-noiseをrobust positivity診断でt=10まで再検証した。

## 2. 元9cの停止理由

SymmetricEigenの非有限出力をstate finite失敗と混同しnumerical_issue_stopとなった。

## 3. diagnosticで判明したこと

rhoは有限、traceとHermiticityは正常、Schurはt=0.02でminimum -3.237e-24を返した。

## 4. 判定規則の問題

state finitenessとsolver finitenessを独立に判定する必要があった。

## 5. robust solver policy

Hermitian化SymmetricEigenをprimaryとし、不合格時だけHermitian化Complex Schurへfallbackした。

## 6. state finiteとsolver finiteの分離

primary失敗だけでは物理状態をFAILにせず、両solver失敗時だけsolver_failureとした。

## 7. unit tests

必須10検査をunit testとruntime CSVで確認した。

## 8. N=7再検証条件

N=7、total gamma=1.5、gamma_site=1.5/7、dt=0.0025、4000 RK4 steps、1001保存点。物理模型は不変。

## 9. 実行時間

total 2703.425s、propagation 1893.772s、diagnostics 805.168s。

## 10. primary solver結果

成功999時刻、失敗2時刻。

## 11. fallback結果

attempt 2、success 2、両solver失敗 0。max imag=1.436e-30。

## 12. t=0.01診断

selected=complex_schur_fallback minimum=-2.8079689841970771e-22 fallback=true。

## 13. t=0.02診断

selected=complex_schur_fallback minimum=-3.2372313550345495e-24 fallback=true。

## 14. t=0.03診断

selected=symmetric_eigen minimum=-2.5580541622425445e-24 fallback=false。

## 15. 全1001時刻のpositivity

worst selected minimum=-5.2780482602309927e-18、solver_failure=0。

## 16. traceとHermiticity

max trace=3.997e-15、max Hermiticity=0.000e0。

## 17. energy ledger

max abs ledger=5.892e-7。

## 18. 既存9c trajectoryとの一致

全指定物理量を1001時刻、許容値1e-12で比較。all checks=true。

## 19. t=10結果

E10=7.1549705124731794e-3 W10=2.6169190648232064e-3 usable=3.6574840668611575e-1 W/Ein=3.8664482773736401e-2。

## 20. W最大値

Wmax=3.8080717406769921e-3 at t=7.7000000000000002e0、E at Wmax=7.3641120898461504e-3、usable=5.1711213710715664e-1。

## 21. N=3/N=5/N=7比較

N7/N5=1.0918825139330242e0、N7/N3=1.2567061564545519e0。N7/N5はsimilar within 10 percent descriptive band。

## 22. fixed-per-site比較との違い

今回はtotal gammaを1.5へ固定した比較であり、siteごと0.5固定とは総雑音が異なる。

## 23. 中心結果

判定 **completed_comparison_with_fallback_diagnostic**。fallbackは診断層だけで、物理時間発展を変更していない。

## 24. 直接確認できたこと

N=7 t=10 trajectory、全保存点のselected positivity、fallback成功、既存9c物理量再現を確認した。

## 25. 確認できていないこと

別dt、別gamma、t>10、N>7、別外部solver crateは未確認。

## 26. 主張してはいけないこと

10%帯は統計的有意差でなく、距離だけの因果証明でもない。

## 27. Milestone 9cの正式判定

**completed_comparison_with_fallback_diagnostic**。物理比較結果を正式採用してよい。

## 28. 生成ファイル一覧

- `robust_eigen_diagnostic_unit_checks.csv`
- `n7_fixed_total_validation_timeseries.csv`
- `n7_fixed_total_validation_eigen_diagnostics.csv`
- `n7_fixed_total_validation_summary.csv`
- `n7_fixed_total_validation_trajectory_comparison.csv`
- `n7_fixed_total_validation_checks.csv`
- `n7_fixed_total_validation_performance.csv`
- `fixed_total_noise_final_comparison.csv`
- `MILESTONE_9C_VALIDATION.md`

selected eigenvalue sum-trace maximum=2.442e-15。
