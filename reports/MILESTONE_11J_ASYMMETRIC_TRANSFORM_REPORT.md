# Milestone 11j: 等入力W曲線の左右非対称時間変形

## 1. 実行範囲

11hの正式成果物と、11hが使用したN=3/N=7等入力時系列だけを使用した。新しいRK4時間発展、N=7再計算、site-resolved診断はない。

## 2. 入力

`fixed_total_gamma_1_5_xgamma_timeseries.csv`（N=3、TOTAL_GAMMA=1.5）と`input_matching_interpolated_trial_timeseries.csv`（N=7、TOTAL_GAMMA=1.5、Omega=0.18748395731510084）を使用。11hのmodels、residuals、summary、reportも監査した。

## 3. 固定モデル

N3正式ピークt3=5.63を境界とし、t_b=delta+s_rise*t3。上り側u=(t-delta)/s_rise、下り側u=t3+(t-t_b)/s_fall。このため時間写像とWモデルは境界で連続し、境界は追加パラメータではない。自由度はA、delta、s_rise、s_fallの4個。

## 4. 探索

A=[0.5,2.0]解析解、delta=[-1,4]、s_rise/s_fall=[0.4,1.6]。粗gridはdelta 0.05・scale 0.04、細gridは粗最良点の周囲をdelta 0.005・scale 0.005で探索。結果を見た探索範囲変更はない。

## 5. 公平比較

主区間は11hと同じpost-arrival t=3.83～10、同じ618点、線形補間、外挿除外。N3正式ピーク=5.63、N7正式ピーク=7.70、最良写像の対応境界=7.75295000。

## 6. 結果

Baseline: A=1.0234364536、delta=2.49000、s=0.91000、normalized RMSE=0.05484958、BIC-like=-10599.800045。

Asymmetric: A=1.0567451682、delta=2.32000、s_rise=0.96500、s_fall=0.51500、|差|=0.45000、normalized RMSE=0.03898819、BIC-like=-11015.264522。scale差は細探索刻み0.005より大きい。

## 7. 判定

事前RMSE規則: **asymmetric_time_scaling_partially_supported**。

複雑度ペナルティ: **improvement_supported_after_complexity_penalty**。AIC-like/BIC-likeは記述的比較であり、生成モデルや物理機構の証明ではない。

## 8. 残差

残差符号変化=3、最大正=1.847604e-4 at t=7.65、最大負=-2.826439e-4 at t=9.50。詳細とearly/peak/late MAEはCSVに保存した。

## 9. 解釈限界

左右非対称変形が単一scaleよりどこまで記述的に改善するかだけを検査した。反射、群速度、mode beating、entanglement、因果機構は判定していない。残差が残るため11iは将来候補として保存するが、自動実行しない。

## 10. Checks

23/23 PASS。7入力は解析前の正式11h成果物と解析後の作業コピーでSHA-256が一致した。

- `fixed_total_gamma_1_5_xgamma_timeseries.csv`: `aa582f956aa8942bc886c537d02bbd6eadd2c75984cce10190079a81d917f3f5`
- `input_matching_interpolated_trial_timeseries.csv`: `09996dd78fc360a8535acaba0c8521ae80fb820f5222888bd657de7162e83164`
- `input_matching_interpolated_trial_summary.csv`: `91205d98ceb7360bcd46a20ba01fe5c2721c1952122022f3a920ceffd9eaa46e`
- `equal_input_curve_transform_models.csv`: `f35dc80f9d161845b84c051d06e2752047d598ffc781fbf0d3aa01b7161f3cfb`
- `equal_input_curve_transform_residuals.csv`: `c35561822eac83978b51b205ce225ed39b2aed923a72d3efa1413c0277498356`
- `equal_input_transform_decomposition_summary.csv`: `1ef6ba9d2499c71ab5a94c3a41be2adfbf6a0aec918f5fb942f18bd76985a940`
- `MILESTONE_11H_REPORT.md`: `30ed88468e38c7ac8632723c1b7f52bd54809bc186a7e04a82a08a75dd92b850`

実行記録：`cargo fmt --all -- --check` PASS。`cargo test --release --offline` は119 passed / 0 failed / 1 ignored、実測115.328秒。`cargo run --release --offline --bin equal_input_asymmetric_transform_analysis` PASS。

## 11. 最終判定

**asymmetric_time_scaling_partially_supported**

補助判定: **improvement_supported_after_complexity_penalty**

## 12. 停止

N=7再計算、dt半減、追加N/Omega、物理機構診断、Milestone 11kへ進まない。
