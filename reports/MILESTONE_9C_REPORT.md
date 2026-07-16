# Milestone 9c: Fixed-total-dephasing comparison across N=3 N=5 N=7

## 1. 目的

全site gammaの単純和を1.5へ固定し、N増加と総雑音増加の交絡を部分的に切り分けた。

## 2. 9b比較の交絡

fixed-per-site gamma=0.5ではtotal gammaがN3=1.5、N5=2.5、N7=3.5と同時に増えていた。

## 3. fixed-total-noise設計

TOTAL_GAMMA=1.5を全siteへ均等配分した。これは総dephasing rateの単純和を固定する記述的比較である。

## 4. 今回の新規計算

N=5 gamma_site=0.3とN=7 gamma_site=1.5/7だけを新規本計算した。N=3は既存gamma_site=0.5結果を参照した。

## 5. 変更していない物理模型

Hamiltonian、drive、RK4、dt=0.0025、t=10、load、初期真空、観測量を変更していない。

## 6. 雑音定義

各chain siteにL_phi,j=sqrt(gamma_j/2) sigma_z,jを適用し、loadへ直接雑音を入れていない。

## 7. gamma配分

N5 gamma_site=2.9999999999999999e-1、N7 gamma_site=2.1428571428571427e-1。両条件のsum gammaは1.5。

## 8. dephasing kernel

Milestone 8cでdense pathとの等価性を確認したDiagonalDephasingKernelを使用した。新しい近似ではない。

## 9. 数値手法

time-dependent dense density-matrix RK4、4000 steps、1001保存点。各保存点で最小固有値とpower ledgerを診断した。

## 10. 構成検査

chain長、次元、bond、drive/load mapping、gamma総和、kernel mapping、load除外、真空、Hermiticityを本計算前に検査した。

## 11. N=5実行結果

E10=8.7957620049e-3、W10=2.5945447799e-3、Wmax=3.4876204098e-3 at t=6.62、peak=peak_resolved。

## 12. N=7実行結果

E10=7.1549705125e-3、W10=2.6169190648e-3、Wmax=3.8080717407e-3 at t=7.70、peak=peak_resolved。

## 13. 数値品質

N5 checks=true、N7 checks=false。N7では `t=0.02` の最小固有値診断1点だけが `NaN` となり、`positivity` と `finite_values` の2検査がFAILした。他の1000保存点の最小固有値は有限で、その最小値は -5.278e-18。max traceは3.997e-15、max Hermiticityは0、max ledgerは5.892e-7だった。仕様に従い品質合格とはせず、issueとして停止した。

## 14. 実行時間

N5 total 50.418s、N7 total 3041.451s。時間差は性能診断であり物理結果ではない。

## 15. fixed-total N=3/5/7比較

WmaxはN3 3.0302005931e-3、N5 3.4876204098e-3、N7 3.8080717407e-3。

## 16. fixed-per-site N=3/5/7比較

既存WmaxはN3 3.0302005931e-3、N5 1.0968175558e-3、N7 4.0437475514e-4。

## 17. N=5での正規化効果

Wmax fixed-total/fixed-per-site=3.1797634814e0、absolute gain=2.3908028540e-3。

## 18. N=7での正規化効果

Wmax fixed-total/fixed-per-site=9.4171846593e0、absolute gain=3.4036969855e-3。

## 19. noise-free基準の残存率

metricごとのnoise-free基準残存率をfixed_total_noise_recovery.csvへ保存した。巨大な比率だけでなく絶対値も併記した。

## 20. Wmax比較

fixed-total N7/N5=1.0918825139e0。分類 `similar within 10 percent descriptive band`。10%帯は統計的有意差ではない。

## 21. t=10比較

N5 E10=8.7957620049e-3 W10=2.5945447799e-3、N7 E10=7.1549705125e-3 W10=2.6169190648e-3。

## 22. arrival time比較

energy arrival N5=3.4100000000000001e0 N7=4.4000000000000004e0、ergotropy arrival N5=2.8599999999999999e0 N7=3.7900000000000000e0。絶対閾値が異なるため単純な到着順とは解釈しない。

## 23. usable fraction比較

t10 usableはN5 2.9497669200e-1、N7 3.6574840669e-1。

## 24. W/Ein比較

t10 W/EinはN5 3.9885039584e-2、N7 3.8664482774e-2。制御費用を含む総合効率ではない。

## 25. 時間面積比較

E area N5 3.7378742342e-2 N7 2.7453340088e-2、W area N5 1.5515466186e-2 N7 1.3522095341e-2。

## 26. 時間窓解析

pulse、early post-pulse、middle、lateの4窓を両条件で同じ定義により保存した。

## 27. 中心判定

判定 **numerical_issue_stop**。原因はN7の `t=0.02` における最小固有値診断1点の `NaN`。fixed-totalでの回復とN5/N7差は参考値として併記するが、品質合格した中心判定とは扱わない。

## 28. 直接確認できたこと

total gamma=1.5のN5/N7本計算、N3参照との有限長比較、fixed-per-site差、noise-free基準残存率を確認した。

## 29. 確認できていないこと

距離だけの純粋因果、他のtotal gamma、dt半減、t>10、N>7、位置別雑音、occupation-weighted exposureは未確認。

## 30. 主張してはいけないこと

指数/べき則、熱力学極限、一般的輸送限界、統計的有意差、total noiseまたはlengthだけの単独因果を主張しない。

## 31. 次段階候補

total gamma sweep、noise exposure integral、site occupation weighted dephasingは候補のみ。自動実行していない。

## 32. 生成ファイル一覧

`src/bin/fixed_total_noise_comparison.rs` と指定11成果物を新規作成した。既存成果物は上書きしていない。
