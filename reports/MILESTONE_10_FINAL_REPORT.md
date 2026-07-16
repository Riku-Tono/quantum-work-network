# Milestone 10 Final Report

## 1. Milestone 10の目的

既存比較を整理し、固定総雑音TOTAL_GAMMA=1.5と3.0を同じ診断系で比較し、dephasing-kernel-weighted coherence exposure XGammaを導入した。10dでは10a・10b・10cの正式成果物だけを読み、新しい時間発展やXGamma再計算をしていない。

## 2. Milestone 10a

新しい物理計算を行わず既存結果を比較し、fixed-total 1.5のN=3・5に残っていた欠損を明示した。最終判定は `completed_with_explicit_missing_values`。

## 3. Milestone 10b

TOTAL_GAMMA=3.0をN=3・5・7で計算し、XGammaを初導入した。fallback 0、solver failure 0で数値検査を通過し、最終判定は `completed_fixed_total_gamma_3_comparison`。

## 4. Milestone 10c

TOTAL_GAMMA=1.5を同じXGamma診断付きでN=3・5・7について再計算し、10aのfixed-total欠損を正式補完した。N=7は9c正本と7物理量×1001時刻で最大差0。N=7のprimary診断2回不合格に対してfallbackは2/2成功、solver failure 0で、最終判定は `completed_with_fallback_diagnostic`。

## 5. 統合された主要結果

| TOTAL_GAMMA | metric | ranking |
|---:|---|---|
| 1.5 | W_max | N=7 > N=5 > N=3 |
| 1.5 | W(t=10) | N=7 > N=5 > N=3 |
| 1.5 | usable fraction | N=7 > N=5 > N=3 |
| 1.5 | W_time_area | N=3 > N=5 > N=7 |
| 1.5 | ergotropy arrival | N=3 fastest, N=5 middle, N=7 slowest |
| 1.5 | XGamma | N=7 > N=5 > N=3 |
| 3.0 | W_max | N=7 > N=5 > N=3 |
| 3.0 | W(t=10) | N=7 > N=5 > N=3 |
| 3.0 | usable fraction | N=7 > N=5 > N=3 |
| 3.0 | W_time_area | N=3 > N=5 > N=7 |
| 3.0 | ergotropy arrival | N=3 fastest, N=5 middle, N=7 slowest |
| 3.0 | XGamma | N=7 > N=5 > N=3 |

## 6. 最大値と時間面積の違い

固定総雑音1.5と3.0の両方で、長い鎖ほどW_max、W(t=10)、usable fractionは大きかったが、W_time_areaは短い鎖ほど大きかった。これはこの模型・有限条件における記述結果で、一般法則ではない。長短を一括評価せず、**評価指標によって順位が逆転する**と結論する。W_time_areaはergotropy状態量の時間積分であり、累積抽出仕事ではない。E_time_areaもload energy状態量の時間積分であり、累積入力エネルギーではない。

## 7. 到着時刻

両fixed-total条件でNが長いほどergotropy arrivalは遅かった。これを輸送速度の普遍則やballistic/diffusive scalingへ結びつけない。

## 8. TOTAL_GAMMA倍増の有限比較

1.5から3.0への24個の直接比較を示す。W系指標は各Nで大幅に減少し、ergotropy arrivalの比は約1.08〜1.09と相対的に小さい変化だった。XGammaは各Nで減少した。

| N | metric | gamma 1.5 | gamma 3.0 | ratio 3/1.5 | absolute difference |
|---:|---|---:|---:|---:|---:|
| 3 | W_max | 3.030201e-3 | 5.464453e-4 | 1.803330e-1 | 2.483755e-3 |
| 3 | W_at_t10 | 2.365248e-3 | 3.963319e-4 | 1.675646e-1 | 1.968916e-3 |
| 3 | W_time_area | 1.736132e-2 | 3.180859e-3 | 1.832152e-1 | 1.418046e-2 |
| 3 | E_at_t10 | 1.259687e-2 | 6.777269e-3 | 5.380119e-1 | 5.819606e-3 |
| 3 | E_time_area | 5.551910e-2 | 2.697790e-2 | 4.859211e-1 | 2.854120e-2 |
| 3 | usable_fraction_at_t10 | 1.877646e-1 | 5.847959e-2 | 3.114515e-1 | 1.292851e-1 |
| 3 | ergotropy_arrival_time | 1.950000e0 | 2.130000e0 | 1.092308e0 | 1.800000e-1 |
| 3 | XGamma | 5.459024e-2 | 4.652299e-2 | 8.522217e-1 | 8.067254e-3 |
| 5 | W_max | 3.487620e-3 | 6.283973e-4 | 1.801794e-1 | 2.859223e-3 |
| 5 | W_at_t10 | 2.594545e-3 | 4.644843e-4 | 1.790234e-1 | 2.130060e-3 |
| 5 | W_time_area | 1.551547e-2 | 2.978707e-3 | 1.919831e-1 | 1.253676e-2 |
| 5 | E_at_t10 | 8.795762e-3 | 5.177189e-3 | 5.886004e-1 | 3.618573e-3 |
| 5 | E_time_area | 3.737874e-2 | 1.797228e-2 | 4.808154e-1 | 1.940647e-2 |
| 5 | usable_fraction_at_t10 | 2.949767e-1 | 8.971746e-2 | 3.041510e-1 | 2.052592e-1 |
| 5 | ergotropy_arrival_time | 2.860000e0 | 3.110000e0 | 1.087413e0 | 2.500000e-1 |
| 5 | XGamma | 5.912513e-2 | 5.434015e-2 | 9.190704e-1 | 4.784974e-3 |
| 7 | W_max | 3.808072e-3 | 6.991746e-4 | 1.836033e-1 | 3.108897e-3 |
| 7 | W_at_t10 | 2.616919e-3 | 5.274431e-4 | 2.015512e-1 | 2.089476e-3 |
| 7 | W_time_area | 1.352210e-2 | 2.657171e-3 | 1.965059e-1 | 1.086492e-2 |
| 7 | E_at_t10 | 7.154971e-3 | 3.987531e-3 | 5.573092e-1 | 3.167439e-3 |
| 7 | E_time_area | 2.745334e-2 | 1.174549e-2 | 4.278346e-1 | 1.570785e-2 |
| 7 | usable_fraction_at_t10 | 3.657484e-1 | 1.322731e-1 | 3.616505e-1 | 2.334753e-1 |
| 7 | ergotropy_arrival_time | 3.790000e0 | 4.100000e0 | 1.081794e0 | 3.100000e-1 |
| 7 | XGamma | 6.018345e-2 | 5.831420e-2 | 9.689409e-1 | 1.869243e-3 |

これはこの模型・初期条件・N=3・5・7・t<=10の2点比較である。gamma倍増の普遍倍率、一般的な関数形、非線形応答は主張しない。

## 9. XGamma

```text
x_gamma(t) = sum_ab Gamma[a,b] |rho[a,b](t)|^2
XGamma(T) = integral_0^T x_gamma(t) dt
```

XGammaはkernelが重み付けしたcoherence exposureという診断量であり、失われた仕事、散逸エネルギー、dephasing power、熱、entropy production、効率、損傷量ではない。強いdephasingがcoherence自体を早く抑え、kernelが重み付けする対象が減った可能性は候補説明になり得るが、今回の計算では因果機構として確認していない。

## 10. 数値品質

10bはfallback 0。10cのN=7はfallback 2/2成功。全fixed-total条件でsolver failure 0。N=7・TOTAL_GAMMA=1.5軌道は9c正本と最大差0で、trace、Hermiticity、positivity、ledger検査はPASSした。

## 11. 直接確認できたこと

この模型、vacuum初期状態、N=3・5・7、t<=10、指定されたnoise-free・fixed-per-site・fixed-total 1.5・3.0条件について、正式成果物に保存された最大値、最終値、時間面積、到着時刻、usable fractionを横断整理した。XGammaの直接比較は同一定義で計算されたfixed-total 1.5と3.0に限る。

## 12. 確認できていないこと

中間gamma、連続gamma sweep、TOTAL_GAMMA依存の関数形、臨界値、XGammaの因果機構、XGamma一致比較、dt半減、t>10、N>7、等入力費用比較、連続運転、抽出サイクル、scaling law、量子優位、実機性能は確認していない。

## 13. 主張してはいけないこと

長い鎖が一般に優れる、短い鎖が一般に優れる、XGammaがW損失を引き起こした、XGammaが損失量そのもの、2つのgamma点から関数形を決定、指数則・べき則・相転移、N>7への外挿、量子優位は主張しない。

## 14. Milestone 10最終判定

10a・10b・10cの正式成果物間に矛盾はなく、closure checksは全件PASSした。noise-freeとfixed-per-siteのE系およびXGammaは正式10a統合表に存在しないため、0や推測値で埋めず `not_available` とした。

最終判定: **completed_with_explicit_missing_values**

実行・検証記録：

```text
cargo fmt --all -- --check
cargo test --release --offline
cargo run --release --offline --bin milestone_10_closure
```

closure binはstd-onlyのCSV／report統合処理で、Hamiltonian、Liouvillian、RK4、dephasing kernel、時間発展モジュールを呼んでいない。Cargo testは107 passed、0 failed、1 ignored。

## 15. 次段階判断

Milestone 10はここで完了する。
追加gamma点、N>7、XGamma matching、等入力費用比較の
どれを次に行うかは、研究目的を再確認してから決定する。
