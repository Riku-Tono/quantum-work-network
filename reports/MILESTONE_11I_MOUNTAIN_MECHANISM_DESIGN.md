# Milestone 11i: 等入力N=7ピークの物理機構診断設計と保存量監査

## 1. 実行範囲

11hの正式成果物5点を読み取り専用で監査した。新しい時間発展、物理軌道、parameter scanは実行していない。

## 2. 既存データ監査

保存済みなのはload energy、load ergotropy、load coherenceの集約量、power・数値健全性診断、N3/N7比較、passive energy、最良変形残差である。site別population、bond/site-load current、site-site coherence、密度行列、相互情報、negativity、mode occupationは保存されていない。

## 3. 現時点で可能な解析

R_W(t)の記述統計と候補周波数抽出、load量の比較は既存時系列から可能。ただし周波数ピークだけを固有modeや機構とは呼べない。

## 4. 群速度・有限鎖mode設計

長鎖基準 epsilon(k)=omega+2J cos(k)、v_g=-2J sin(k)、最大速度2Jをarrival scaleにだけ使う。N=7開放鎖の離散k_m=m*pi/8と、終端loadを含む拡張sectorの静的固有modeを分ける。drive終了t=3.2でsingle-excitation sector probabilityとmode occupationを保存し、sector probabilityが小さい場合はsingle-particle解釈を停止する。

## 5. 反射診断

entrance、middle、exit/loadのforward arrival後に、current sign reversal、backward-moving front、距離と整合する再到着が揃うことを候補条件とする。current反転だけでは反射を証明しない。

## 6. 相関・entanglement診断

相互情報 I(A:B)=S(A)+S(B)-S(AB) は古典・量子相関を含み、entanglementそのものとは呼ばない。exit-site:loadとlast-two-sites:loadでは negativity=(||rho^{T_A}||_1-1)/2 を選択時刻だけ計算する。

## 7. 残差5.5%との対応

R_Wとexit-load current、反射current、相互情報、negativity、mode beating候補を比較する。peak time差、cross-correlation lag、Pearson、Spearmanを保存するが、相関を因果としない。

## 8. 事前固定時刻

t=0、3.2、3.83000000（ergotropy arrival）、6.871718418689（crossing 1）、7.70000000（N7 W最大）、9.378580405004（crossing 2）、10。結果を見て時刻を自動増殖させない。

## 9. 最小の次回計算候補

N=7、TOTAL_GAMMA=1.5、Omega=0.18748395731510084、dt=0.0025、T=10のmatched軌道1本だけ。物理模型は変えず、全1001時刻に7 site populations、6 bond currents、site-load current、6 nearest-neighbor coherences、single-excitation-sector weightを追加する。選択7時刻だけ小さいreduced statesと相関量を計算する。

## 10. データ量監査

N=7の全Hilbert次元は2^7*3=384。384x384 complex128は1時刻2,359,296 bytes（2.25 MiB）、1001時刻で約2.20 GiBのbinary相当で、CSVはさらに大きい。そのため全密度行列時系列は保存しない。site-resolved scalar表は約1001行・約28列で、概算1 MiB未満。6x6と12x12のreduced statesだけを7選択時刻に保存する。

## 11. 現在の候補機構判定

群速度、分散、境界反射、correlation/entanglement front、load局所応答、mode beatingのいずれも未確定。既存load時系列だけでは必要な空間情報と状態情報が不足する。

## 12. Checks

15/15 PASS。既存データ監査、定義分離、事前時刻固定、全密度行列時系列禁止、単一再計算制限、新規時間発展なしを確認した。

## 13. 最終判定

**completed_design_with_targeted_recomputation_required**

## 14. 実行記録

`cargo fmt --all -- --check` PASS、`cargo test --release --offline` 119 passed / 0 failed / 1 ignored、`cargo run --release --offline --bin mountain_mechanism_design_audit` PASS。

## 15. 次段階

Milestone 11iの監査後、
N=7 matched軌道を同一物理条件で1回だけ再計算し、
site-resolved current・population・選択時刻の相関量を追加保存するか判断する。
