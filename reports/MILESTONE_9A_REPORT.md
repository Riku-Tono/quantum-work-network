# Milestone 9a: N=7 noise-free full reachability run

## 1. 目的

N=7 noise-freeをdt=0.0025でt=10まで本計算し、load energy・ergotropy到達とN=3/N=5との差を確認した。

## 2. Milestone 8aから8cまでの位置づけ

8aのN=3/N=5結果、8bのN=7 feasibility、8cの厳密dephasing kernel検証を前提にした。

## 3. 今回はN=7 noise-freeだけであること

N=7 all-site noisy、細刻み、t>10、追加最適化は実行していない。

## 4. 変更していない物理模型

Hamiltonian、drive、RK4、dt、基底、load、初期真空、観測量を既存実装から変更していない。

## 5. N=7模型構成

7二準位site、3準位load、dim=384、bond=6、drive site=0、load coupling site=6、collapse=0。

## 6. 数値手法

dense complex density matrixのtime-dependent RK4。4000 step、保存間隔0.01、1001点。各保存点で縮約・ergotropy・全系最小固有値を診断した。

## 7. 構成検査

次元、mapping、真空、drive端点、Hamiltonian Hermiticityを時間発展前に検査し全合格した。

## 8. 実行時間とメモリ

construction 3.049s、propagation 1973.170s、diagnostics 976.635s、total 2953.536s。peak working set 99332096 bytes。

## 9. 数値品質チェック

全チェックPASS。max trace=2.109e-15、max Hermiticity=0.000e0、min eigenvalue=-6.578e-12、max ledger=6.566e-7。

## 10. load energy到達

持続閾値1e-4の到達時刻は 4.1299999999999999e0。

## 11. load ergotropy到達

持続閾値1e-5の到達時刻は 3.5500000000000003e0。

## 12. coherence生成

coherence L1最大値は 3.1015172865e-1、時刻 7.78。

## 13. t=10結果

E=1.3960933446e-2、W=1.3085005957e-2、usable=9.3725867314e-1、W/Ein=1.7435814453e-1。

## 14. W最大値と時刻

W_max=2.2435519997e-2、t=7.71、その時のE=2.3850844709e-2。

## 15. 終端挙動とピーク判定

分類 `peak_resolved`。W(10)-W(9.9)=-1.070e-4、E(10)-E(9.9)=-1.252e-4、最終10点W slope=-1.049e-3、E slope=-1.229e-3。

## 16. usable fraction

t=10で 9.3725867314e-1、W最大時で 9.4065934647e-1。

## 17. W/Ein

t=10で 1.7435814453e-1、W最大時で 2.9895405863e-1。制御費用を含む総合効率ではない。

## 18. site population時間発展

1001時刻×7site=7007行をlong形式で保存した。site_indexは0始まり、site_labelは1始まり。

## 19. 時間窓解析

pulse、early post-pulse、middle、lateの4窓をwindows CSVへ保存した。時間面積は状態量の面積であり累積仕事ではない。

## 20. N=3、N=5、N=7比較

WmaxはN3 5.5424064638e-2、N5 2.1974463817e-2、N7 2.2435519997e-2。N7/N5=1.0209814530e0、N7/N3=4.0479744933e-1。

## 21. N=3からN=7への有限差

3点の有限長比較のみであり、指数・べき則・漸近scalingは推定しない。Nとともに距離、次元、bond数、干渉構造が同時に変わる。

## 22. 直接確認できたこと

固定模型でN=7のt=10までのenergy/W到達、W最大値、usable、W/Ein、有限長3点差、実測時間を確認した。

## 23. 確認できていないこと

N=7 noisy、細刻み、t>10、最終到達上限、N>7、scaling、実機性能は未確認。

## 24. 主張してはいけないこと

指数/べき減衰、熱力学極限、N>7外挿、物理的輸送限界、noisy予測、量子優位、新規性、距離だけの純粋因果。

## 25. N=7 noisy本計算へ進む判断

判定 **proceed_candidate**。8c noisy推定は約0.791時間。ただしnoisy本計算は自動実行していない。

## 26. 生成ファイル一覧

`src/bin/n7_noise_free_full.rs` と指定9成果物を新規作成した。

