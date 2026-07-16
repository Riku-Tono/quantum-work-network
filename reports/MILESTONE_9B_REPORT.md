# Milestone 9b: N=7 all-site noisy full run

## 1. 目的

N=7の全7siteへgamma_phi=0.5の位相雑音を入れ、t=10までenergyとergotropyがloadへ届くか確認した。

## 2. 今回の範囲

新規本計算はN=7 all-site noisy、dt=0.0025、t=10の1条件だけである。noise-free再計算、細刻み、時間延長、sweep、次Milestoneは実行していない。

## 3. 参照成果物

Milestone 8aのN=3/N=5、8cのkernel検証、9aのN=7 noise-free成果物を読み取り専用で参照した。

## 4. 変更していない物理模型

Hamiltonian、drive、RK4、dt、基底、load、初期真空、観測量を既存実装から変更していない。

## 5. 雑音模型

全7siteにL_phi,j=sqrt(0.5/2) sigma_z,jを入れ、loadには直接雑音を入れていない。

## 6. dephasing kernel

Milestone 8cでdense collapse pathと等価性確認済みのDiagonalDephasingKernelを使用した。物理近似ではなく同じLindblad項の成分表示である。

## 7. N=7模型構成

7二準位site、3準位load、dim=384、bond=6、drive site=0、load coupling site=6、noisy site count=7。

## 8. 数値手法

dense complex density matrixのtime-dependent RK4。4000 step、保存間隔0.01、1001点。各保存点で縮約、ergotropy、全系最小固有値、dephasing powerを診断した。

## 9. 構成検査

次元、全gamma、kernel寸法・対称性・非負性・mapping・load除外、真空、drive端点、Hamiltonian Hermiticityを本計算前に検査した。

## 10. 実行時間とメモリ

construction 4.160s、propagation 2030.362s、diagnostics 863.969s、total 2899.253s。8c推定 2846.262s。peak working set 99262464 bytes。

## 11. 数値品質チェック

全チェックPASS。max trace=1.443e-15、max Hermiticity=0.000e0、min eigenvalue=-7.444e-18、max ledger=5.157e-7。

## 12. load energy到達

持続閾値1e-4の到達時刻は 4.8500000000000005e0。

## 13. load ergotropy到達

持続閾値1e-5の到達時刻は 4.2199999999999998e0。

## 14. t=10結果

E=3.3040555691e-3、W=3.1157407105e-4、usable=9.4300493600e-2、W/Ein=5.2362601887e-3、coherence L1=3.5398341400e-2。

## 15. W最大値と時刻

W_max=4.0437475514e-4、t=7.65、その時のE=2.0692277651e-3、usable=1.9542302784e-1。

## 16. 終端挙動とピーク判定

分類 `peak_resolved`。W(10)-W(9.9)=-3.130e-6、E(10)-E(9.9)=4.694e-5、最終10点W slope=-3.122e-5、E slope=4.693e-4。

## 17. energy ledger

t=10でdrive net=5.9503168257e-2、dephasing net=1.1166652695e-14。ledgerは両方を含み、最大残差は 5.157e-7。

## 18. site populationと時間窓

1001時刻×7site=7007行を保存し、pulse、early post-pulse、middle、lateの4窓を集計した。

## 19. N=7 noisy/free比較

t=10でE noisy/free=2.3666437362e-1、W noisy/free=2.3811534520e-2。Wmax noisy/free=1.8023863730e-2。Wmaxは異なる時刻同士の最大値比較である。

## 20. N=3、N=5、N=7 noisy比較

WmaxはN3 3.0302005931e-3、N5 1.0968175558e-3、N7 4.0437475514e-4。N7/N5=3.6868005348e-1、N7/N3=1.3344818031e-1。

## 21. 中心問題への有限計算上の答え

N=7 all-site noisyでenergy到達=true、ergotropy到達=true。noise-freeで見えたN7 Wmax>N5 Wmaxという特徴はnoisyでは残らなかった。N7/N5=3.6868005348e-1。

## 22. 解釈上の制限

Nとともに距離と雑音site数が同時に増える。3点の有限長比較から距離だけの因果、指数/べき則、漸近scaling、輸送限界を主張しない。

## 23. 未確認

dt半減、t>10、位置別雑音、保護、gamma/Omega sweep、N>7、実機性能は未確認である。

## 24. 最終判定

判定 **completed_comparison**。どの判定でも時間延長や次Milestoneは自動実行していない。

## 25. 生成ファイル

`src/bin/n7_all_site_noisy_full.rs` と指定11成果物を新規作成した。既存成果物は上書きしていない。
