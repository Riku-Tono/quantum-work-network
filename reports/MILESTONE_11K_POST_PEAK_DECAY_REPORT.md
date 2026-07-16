# Milestone 11k: 等入力Wピーク後の減衰率・尾部残差解析

## 1. 範囲

11j正式成果物とN=3/N=7等入力保存時系列のみを使用。新規RK4・軌道再計算なし。主区間t=7.70～10.00、補助区間7.70～9.00と9.00～10.00。11j境界t_b=7.75295。

## 2. 瞬間減衰率

k_Q=-d ln(Q)/dtを中央差分、端点のみ片側差分で計算。floor=1e-10。N7実測W、11jモデルW、N3 Wと、N7のE/P/coherenceをCSV化した。

## 3. 低自由度モデル

Model A単一指数、Model B連続二段階指数（switch 8.5/9.0/9.5のみ）、Model C指数×線形だけを評価。AIC-like/BIC-likeは記述的比較で物理過程の証明ではない。

## 4. 主結果

BIC-like最良=model_B_two_stage_exponential、normalized RMSE=0.00809806、BIC-like=-4836.732603。Model A nRMSE=0.01014188、Model B nRMSE=0.00809806（switch=8.50, lambda1=0.15000, lambda2=0.20000）、Model C nRMSE=0.00913268。

二段階指数は単一指数よりBIC-likeが98.530改善し、|lambda2-lambda1|=0.05だった。このため二段階減衰は記述的には支持される。ただし、二つの物理過程を意味しない。

## 5. 11j残差の尾部

post-peak absolute area=4.075384e-4、9～10 absolute fraction=0.596160、符号変化=1、最大正=1.839769e-4 at 7.70、最大負=-2.826439e-4 at 9.50。分類: **tail_residual_concentrated_late**。

## 6. E/W/passive/coherence

各量をt=7.70値で規格化し、対数減衰率とともに保存した。これは同時変化の記述であり、coherenceがW減衰を引き起こしたとは主張しない。t=8～10のDelta E/W/passive/usableと11j残差も同一表に保存した。usable fraction上昇を独立性能向上とは扱わない。

t=10の規格化値はE=0.9713、W=0.6872、passive=1.2809、coherence=0.8289だった。したがってこの区間では、E全体よりWが速く低下し、passive部分は低下せず増加した。N7-N3 usable差は正のままだが、7.70の0.2444からt=10の0.1812へ縮小しており、終盤残差と同時に独立したusable向上が起きたとは読まない。

## 7. Checks

13/13 PASS。正式入力は実行前後および11j提出版とのSHA-256一致を確認した。

- N7時系列: `09996dd78fc360a8535acaba0c8521ae80fb820f5222888bd657de7162e83164`
- N3時系列: `aa582f956aa8942bc886c537d02bbd6eadd2c75984cce10190079a81d917f3f5`
- 11j残差: `79d4a1c6bd8ee7d4eacae686784f8a087746078c1ca972306bcc973a3178d81a`
- 11jモデル: `c58b50cd045e85afd2ec4b780736111cb1d40c7b3dd2bcb310da3a0036963f15`
- 11jレポート: `605246c7f0a6898b1d14ba27514025257c2c9fb3a45dd3633635b4e2a8cea525`

`cargo fmt --all -- --check` PASS。`cargo test --release --offline`は119 passed / 0 failed / 1 ignored、実測17.529秒。解析bin PASS。

## 8. 判定

**late_tail_structure_remains**

併記：`two_stage_tail_decay_supported`（減衰fit内部の記述的比較）。11j残差の9～10集中判定を優先し、最終判定は上記とした。

## 9. 解釈限界

減衰fitは形状診断であり、境界反射、群速度、二つの物理過程、因果機構を証明しない。

## 10. 次段階

終盤構造が残る場合、11iのsite-resolved診断を同一N=7 matched軌道1本に限定して検討できる。ただし自動実行しない。
