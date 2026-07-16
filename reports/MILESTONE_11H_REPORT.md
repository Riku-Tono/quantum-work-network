# Milestone 11h: 等入力W曲線の制約付き時間変形適合とenergy-ergotropy分解

## 1. 目的

保存済みN=3・N=7等入力軌道だけを使い、単純な振幅・時間移動・時間拡縮適合とE-W-passive分解を行った。新規時間発展はない。

## 2. 入力matching確認

相対入力差=2.8157011099880636e-6で1e-4以内。

## 3. 使用した正式時系列

両条件1001点、t=0〜10、同一0.01 grid。

## 4. 変形モデル

Wmodel(t)=A W3((t-delta)/s)。Model 0〜3だけを決定論的粗→細探索し、外挿は除外した。

## 5. 各モデルの適合結果

full最良（BIC-like）=model_1_amplitude、post-arrival最良=model_3_amplitude_shift_scale。

## 6. 複雑度ペナルティ付き比較

AIC-like/BIC-likeは生成モデル推論ではなく記述的複雑度比較である。

## 7. 最良モデル

post: A=1.02343645、delta=2.49000000、s=0.91000000、normalized RMSE=0.054850、max residual=4.673381e-4。

## 8. 残差構造

符号変化=3、最大正残差=2.897797e-4 at t=7.78、最大負残差=-4.673381e-4 at t=10.00、正/負面積=4.504434e-4/4.083505e-4。

## 9. 単純変形で説明できる範囲

分類: **partial_transform_with_structured_residual**。閾値は記述規則で物理法則ではない。

## 10. E-W-passive energy分解

t10 P3=1.023163e-2、P7=3.978217e-3、比=0.388816。

## 11. usable fraction差のexact分解

t10 DeltaU=1.812415e-1 = W差項 -6.150946e-3 + energy分母項 1.873925e-1。N=3基準の一順序で唯一の寄与分解ではない。

## 12. 一次分解

t10 W成分=-3.078524e-3、E成分=9.378919e-2、非線形残差=9.053088e-2。有限差の一次近似で因果分解ではない。

## 13. 時間窓別の差

4窓のE/W/P面積とmean UをCSVに保存し、総面積との一致を確認した。

## 14. 直接確認できたこと

単純変形の記述精度、残差時間構造、passive energy差、usable差の代数分解だけを直接確認した。

## 15. 確認できていないこと

物理機構、因果寄与、matched dt半減、N=5/TOTAL_GAMMA=3.0等入力、追加Omega、N>7。

## 16. 主張してはいけないこと

浄化、選別、時間変形による機構証明、独立した性能向上、分解項の因果寄与。

## 17. 最終判定

**completed_equal_input_transform_decomposition**

## 18. 次段階

単純変形で説明できない残差構造とenergy分解を確認後、matched N=7条件のdt半減検証を行う価値を判断する。自動実行していない。

## 19. 実行記録

`cargo fmt --all -- --check` PASS、`cargo test --release --offline` 119 passed / 0 failed / 1 ignored、`cargo run --release --offline --bin equal_input_transform_decomposition` PASS。動的時間伸縮・高自由度warping・新規時間発展なし。
