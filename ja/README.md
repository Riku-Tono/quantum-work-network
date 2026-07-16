# Quantum Work Network

## プロジェクト概要

Quantum Work Network は、小さな量子ネットワークを流れるエネルギーが、最終的に「取り出せる仕事（ergotropy）」としてどれだけ受け皿（load）に届くかを、Rustで数値的に計算するための参照実装です。

物理系は、鎖状につながった複数の二準位サイト（qubit）と、その端につながった1個の受け皿（3準位のload）からなります。片端のサイトを外部パルスで駆動してエネルギーを注入し、それが鎖を伝わってloadに蓄えられる様子を、密度行列の時間発展として追跡します。位相雑音（dephasing）を加えた条件と加えない条件を比較し、「同じ状況で、雑音があると局所的に取り出せる仕事がどう変わるか」を、有限サイズ・有限時間の範囲で直接測定することが目的です。

このリポジトリが計算するのは、あくまで「この物理模型・この初期条件・この有限時間・このRust実装で直接観測された数値」です。量子優位性や普遍的なスケーリング則を証明するものではありません。原因の説明（なぜそうなるか）と、観測された結果（何が起きたか）は区別して記述しています。

開発はMilestone 1から9c（およびその検証）まで段階的に進みました。各段階は、前段階の物理模型・CSV・レポートを変更せずに、新しい実装や比較を積み増す方針で作られています。

---

## 模型と数値規約

以下は Milestone 1 で凍結され、その後の全段階で維持されている規約です。

**物理模型**

- サイト構成：N個の二準位サイト（`|0>`＝空、`|1>`＝励起）＋1個の3準位load。
- 既定パラメータ（ℏ=1）：サイト間結合 `J=1`、サイト–load結合 `g=0.25`、各サイト・load・駆動の角周波数=1。
- 鎖のオンサイトHamiltonian：`omega * sum_i |1><1|_i`。
- 駆動site=0（片端）、load coupling site=N-1（反対端）。
- coherent drive：`H_drive(t) = Omega f(t){ exp(-i omega t) sigma_1^+ + exp(+i omega t) sigma_1^- }`、包絡 `f(t)=sin^2(pi t / tau)`（`0<=t<=tau`、それ以外0）、既定 `tau=3.2`、`Omega=0.2`。
- 位相雑音（dephasing）：各chain siteに `L_phi,j = sqrt(gamma_phi/2) sigma_z,j`。loadには直接雑音を入れない。
- 時間発展：Lindbladマスター方程式 `d rho/dt = -i[H(t), rho] + sum_k D[L_k] rho`。

**基底順序とベクトル化**

- 一般のNではテンソル順序は `|q1, q2, ..., qN, load>`。右端のload indexが最も速く動く。
- N=3の場合は `|q1, q2, q3, load>`。
- 密度行列は **column-major vectorization**：`vec(rho) = [rho(0,0), rho(1,0), ..., rho(0,1), ...]^T`。したがって `vec(A rho B) = (B^T ⊗ A) vec(rho)`。
- N=3では `24 x 24` 密度行列 → 長さ576ベクトル、Liouvillianは `576 x 576`。

**Liouvillian規約**

```
L = -i (I ⊗ H - H^T ⊗ I)
    + sum_k [ L_k* ⊗ L_k
            - 1/2 I ⊗ (L_k^dagger L_k)
            - 1/2 (L_k^dagger L_k)^T ⊗ I ]
```

collapse operatorは、係数（例：`sqrt(gamma) sigma_minus`、`sqrt(gamma_phi/2) sigma_z`）を**含めた形**で渡す規約です。

**局所ergotropy（取り出せる仕事）**

```
W(rho_L) = Tr(rho_L H_L) - min_U Tr(U rho_L U^dagger H_L)
```

loadの縮約状態から計算される、そのloadに対する局所的な最大抽出可能仕事です。

---

## Milestone一覧

各段階について「目的」「追加・変更したもの」「実行または検証した内容」「その段階ではまだ確認していないこと」を簡潔に示します。数値の詳細は各MILESTONEレポートおよびCSVを参照してください。

### Milestone 1
- **目的**：時間発展を信頼する前に必要な静的モジュール（operators / partial_trace / ergotropy）を正しく実装する。
- **追加・変更**：基底順序・load演算子・ergotropy計算などの規約を凍結。
- **実行・検証**：各モジュールの単体テスト（このプロジェクトを生成した環境にはRustツールチェーンが無かったため、ソースとテストは生成されたがこの場ではコンパイルされていない旨が当時の資料に明記されている）。
- **未確認**：時間発展（Liouvillian）の正しさはこの段階の対象外。

### Milestone 2
- **目的**：時間発展の2モジュールを追加する。
- **追加・変更**：`liouvillian`（column-majorベクトル化とHamiltonian/Lindblad superoperator構築）、`propagator`（各時刻でのdense行列指数による精度優先伝播）。
- **実行・検証**：26個の通常テスト（ベクトル化、`t=0`恒等、trace/Hermiticity/positivity保存、閉系ユニタリとの一致、振幅減衰・純粋dephasingの解析解一致など）。
- **未確認**：診断量、protocol matching、パラメータ探索、作図は未導入。

### Milestone 2.1
- **目的**：完全な24次元模型が下流ユーザー視点で正しく組み合わさることを確認する。
- **追加・変更**：`tests/full_24d_short_time.rs`（`576 x 576` Liouvillianを構築し `t=0` から `t=0.001` まで伝播する、opt-in統合smoke test）。
- **実行・検証**：`ModelParams::default()` に注入collapse operator `sqrt(0.1) sigma_1_plus` を与え、次元・trace・Hermiticity・positivity・有限性・真空からの非ゼロ変化を確認。
- **未確認**：診断量やmatchingはまだ無い。

### Milestone 3
- **目的**：状態診断と符号付きpower会計を追加する。
- **追加・変更**：`diagnostics`（縮約load energy・ergotropy、エネルギー分解、load current、source/dephasing power、physicality指標、符号付き台形power積分）。
- **実行・検証**：診断量を含む形での既存テスト維持。protocol matching・効率主張・パラメータ探索・作図は意図的に未導入。
- **未確認**：効率の主張、matching、sweepはこの段階では行わない。

### Milestone 3.1
- **目的**：符号付きpower積分の扱いを明文化・修正する（Milestone 3とは別の新しい物理実験ではなく、診断機能を維持したままの整備）。
- **追加・変更**：`integrate_signed_power` で、台形積分区間の両端でpowerの符号が変わる場合に、線形補間したゼロ交差点で区間を分割する処理を明示。例：`[(0.0, 1.0), (1.0, -1.0)]` は `energy_net=0`、`energy_in=0.25`、`energy_out=0.25`。
- **実行・検証**：`energy_in`/`energy_out` がゼロ交差分割で正しく分離されることの確認。
- **未確認**：物理模型そのものは変更していないため、新しい物理的主張は無い。

### Milestone 4
- **目的**：初期コヒーレント状態を用いた単発輸送で、同一時刻・同一load energyに一致させたときのA（雑音なし）とB（位相雑音あり）のergotropyを比較する。
- **追加・変更**：coherent-input実験と、確定結果を凍結する `MILESTONE_4_RESULT.md`。
- **実行・検証**：評価時刻 `3.0/5.0/7.9/10.0`、B dephasing強度 `0.1/0.2/0.5/1.0` の全16条件でenergy一致。全条件でA ergotropy > B ergotropy（一致比 `1.229`〜`49.318`）。diagonal ergotropyは全条件0で、差はcoherence由来。等入力（`p_B=p_A`）比較でも全16条件でA>B。物理チェック失敗0。
- **未確認**：単発輸送であり、真空からの生成・連続供給ではない。等入力費用比較でない条件を含む。古典比較・量子優位は未検証。

### Milestone 5a
- **目的**：時間依存Lindblad方程式のRK4伝播器を実装・検証する。
- **追加・変更**：`src/time_dependent.rs`（固定最大刻みRK4、時刻依存 `H(t)` と collapse operators、保存スケジュール）、公開モジュール化、1量子ビットsanity check bin。
- **実行・検証**：時間不変参照問題でdense指数解に対し刻み半減で単調収束（誤差 `7.80e-7 → 4.99e-8 → 3.19e-9`、誤差比≈15.6、RK4期待16に近い）。`cargo test --release`：47 passed / 0 failed / 1 ignored。
- **未確認**：本番ネットワークでの時間依存駆動、A/B比較、energy matching、連続供給、RK4が任意の大刻み・強駆動でpositivityを保つこと。

### Milestone 5b
- **目的**：真空初期状態から有限パルス駆動でloadのcoherenceと非ゼロergotropyが生成されるかを確認する（sanity check）。
- **追加・変更**：既存24次元 `H0` に単一 `sin^2` パルスを適用。A：`gamma_phi=0`、B：`gamma_phi=0.5`（3サイト全て）。
- **実行・検証**：最大load ergotropy A `5.5424e-2`（t=9.48）、B `3.0302e-3`（t=5.63）。同時刻比較の成功条件7項目PASS。物理チェック・energy ledger・dt半減収束を確認。
- **未確認**：同一時刻・同一load energyでの公平比較、energy matching、探索、連続駆動、仕事抽出、古典比較、量子優位。各最大値は時刻が異なりうるためsanity check用。

### Milestone 5c（中心的な比較結果）
- **目的**：t=10でA/Bのload energyを条件付きで一致させ、そのときのergotropyを比較する。
- **追加・変更**：`Omega_B` を `0.2〜1.0`（81点）走査し、符号変化区間を二分法で精密化。単調性は仮定しない。
- **実行・検証**：唯一の符号変化区間から `Omega_B=0.431953125`。t=10でload energy相対差 `4.001e-5`（`<1e-4`）に一致。ergotropy A `5.2798e-2` / B `8.2846e-3`、**A/B比 = 6.373**。成功条件10項目すべてPASS。dt半減で方向不変。
- **未確認**：一致しているのは比較時刻・最終load energy・模型・パルス形状・駆動周波数のみ。**駆動強度Omegaと総投入エネルギーは一致していない**ため、等入力費用比較でも、雑音だけを独立に変えた因果比較でもない。連続運転・仕事抽出・古典比較・量子優位は未確認。

### Milestone 6a
- **目的**：5cの終状態から、ergotropyの定義に対応する理想局所unitaryで実際に仕事を回収する実装検算（新しい物理発見ではない）。
- **追加・変更**：5c設定を `dt=0.0025` で決定論的に再構成し、`H_drive(10)=0` を検算のうえ、突然切断後に理想的瞬時load-local unitaryを適用。
- **実行・検証**：gross extracted work = load ergotropy（A `5.2798e-2`、B `8.2846e-3`）を回収、抽出後ergotropyは両方0。switch workは数値的に0。18/18項目がA/Bともに合格。gross work A-B比 `6.373`。
- **未確認**：抽出unitaryの制御費用、有限時間切断・抽出、繰り返し放電、連続運転、相関からのglobal extraction、古典比較、量子優位。switch workの定義は理想化で装置費用を完全評価していない。

### Milestone 7a
- **目的**：固定N=3模型で、位相雑音を1サイトだけに置いたときのloadへの影響を位置別に比較する。
- **追加・変更**：雑音を入れるsite集合を指定できる入口を追加（既存の全3サイトAPIはラッパーとして保持）。
- **実行・検証**：`gamma_phi=0.5` で入口/中央/出口の3配置とnoise-freeを比較。t=10でW最小は `site1`（入口）、usable fraction最小は `site3`（出口）。中央雑音は他2配置より損失が小さい。全物理チェックPASS、dt半減でW最小位置は不変。
- **未確認**：他の `gamma_phi`・Omega・パルス・長い鎖・長時間・連続運転・抽出・古典比較。雑音位置の普遍則や因果断定はしない。

### Milestone 7b
- **目的**：7aの固定済み時系列から、雑音条件とnoise-freeの差が「いつ持続的に現れたか」を記述する（新しい時間発展は行わない）。
- **追加・変更**：CSV読み込みと解析だけのbin（Hamiltonian/RK4等は呼ばない）。持続onset・最大損失時刻・onset時population・時間窓比較・順位切替を算出。
- **実行・検証**：全40チェックPASS。EとWは3条件・3閾値でt=2.25の診断上同時開始。usable fractionは閾値依存（mediumでexit 1.27 < entrance 1.57 < middle 1.62）。middleの被害が軽い理由は最大populationや時間面積だけでは説明し切れないという「手がかり」を提示。
- **未確認**：特定物理過程への因果分解、populationだけによる説明、coherence/site間current、保護による回復、一般化。

### Milestone 7c
- **目的**：全3サイト雑音から、指定サイトの位相雑音演算子だけを理想的に除去した反実仮想条件で、load状態量の回復上限を比較する。
- **追加・変更**：`all_noisy=[0,1,2]`、`protect_entrance=[1,2]`、`protect_exit=[0,1]`、`protect_both_ends=[1]`、`noise_free=[]` の5条件。理想保護（collapse operator完全除去）であり、現実装置・制御・cost・誤り訂正ではない。
- **実行・検証**：全75チェックPASS。t=10 W回復・usable fraction回復・E回復いずれも `protect_both_ends`（中央雑音のみ残す）が最大。入口保護と出口保護の回復量は近く、順位は時間で入れ替わる。保護の非加算性（synergy）を positive_nonadditivity と分類（synergyは相互作用エネルギーではなく観測量応答の非加算性）。
- **未確認**：現実的実装・cost・不完全保護・他パラメータ・長時間・因果機構。中央雑音が無害という主張はしない。

### Milestone 7d
- **目的**：中央site gamma=0.5を固定し、両端gammaだけを `0.5→0` へ下げたときの回復曲線を調べる。
- **追加・変更**：site別gamma API（負値・非有限を拒否、gamma=0のcollapseを除外）。sweep点 `0.50/0.40/0.30/0.20/0.15/0.10/0.05/0.00`。
- **実行・検証**：全152チェックPASS。端点は7cを絶対誤差 `1e-9` 以内で再現。W_maxもusable fractionも離散的に単調非減少。最大感度区間はいずれも `0.05→0.00`。曲率・感度は離散点のみから計算し、臨界指数・相転移とは呼ばない。
- **未確認**：実装可能な必要保護強度、物理的感受率・臨界指数・相転移、他パラメータ、長時間・連続運転、因果機構。

### Milestone 8a
- **目的**：N=3からN=5へ鎖長だけを変え、load energy・ergotropy・usable fractionを比較する。
- **追加・変更**：`chain_length` 引数を取る演算子構築・駆動実行API（既存N=3 APIは保存）。`dim=2^N*3`。
- **実行・検証**：N=3回帰を絶対誤差 `2e-9` 以内で再現。N=5 noise-freeでload ergotropy生成を確認。W_max比 N5/N3 は free `0.3965`、noisy `0.3620`。t=10・個別ピークともにN=5 W < N=3 W、noisy < free（大小結論は同じだが比率は異なる）。指定3条件のdt半減整合PASS（N=3 noisyの半減は未実行）。
- **未確認**：N>5、連続Nスイープ、位置別弱点、一般scaling、総雑音一致、保護費用。鎖長を変えると距離・次元・bond数・all-site時の総雑音が同時に変わる（等総散逸比較ではない）。

### Milestone 8b
- **目的**：N=7のdense Lindblad RK4が現在の計算環境で実行可能かを、短時間probeだけで評価する。
- **追加・変更**：構築のみ・1ステップ・短時間probeを行うfeasibility probe bin。
- **実行・検証**：Hilbert次元384、密度行列 `384x384`。構築・1ステップ・短時間probeの数値品質PASS。`t=0.1` noisyのstep時間 `21.32 s/step` から t=10 推定 約23.7時間となり、現行dense法での **infeasible_with_current_dense_method** を確定。
- **未確認**：N=7のt=10最終値、半減刻み実測、長時間安定性、最適化後性能。これは「現在環境・現在実装での可否」であり物理的可能性ではない。

### Milestone 8c
- **目的**：局所sigma_z位相雑音のLindblad項だけを、物理模型を変えずに厳密高速化する。
- **追加・変更**：`DiagonalDephasingKernel`（各要素へ `-Gamma[a,b] rho[a,b]` を加える、同じLindblad dissipatorの厳密な成分表示。物理近似ではない）。旧dense pathも保持。
- **実行・検証**：checks CSV **140 PASS / 0 FAIL**。N=3/N=5でdense pathおよび既存値と許容内一致。N=7 noisyのmedian `0.7116 s/step`、旧値比 **29.96x**。t=10再推定 約0.791時間、分類 **feasible_candidate**。
- **未確認**：N=7 t=10の最終物理量、長時間の実測性能、N scaling、他のLindblad operatorへの適用。

### Milestone 9a
- **目的**：N=7 noise-freeを `dt=0.0025` でt=10まで本計算し、到達とN=3/N=5との差を確認する。
- **追加・変更**：`n7_noise_free_full.rs` と成果物一式。
- **実行・検証**：全チェックPASS（max trace `2.109e-15`、min eigenvalue `-6.578e-12` 等）。t=10で E `1.3961e-2`、W `1.3085e-2`、usable `0.9373`。W_max `2.2436e-2`（t=7.71）。W_max比 N7/N5=`1.0210`、N7/N3=`0.4048`。分類 `peak_resolved`。
- **未確認**：N=7 noisy、細刻み、t>10、最終到達上限、N>7、scaling、実機性能。

### Milestone 9b
- **目的**：N=7の全7サイトに `gamma_phi=0.5` を入れた all-site noisy を t=10 まで本計算する。
- **追加・変更**：`n7_all_site_noisy_full.rs`（8cのkernelを使用）と成果物一式。
- **実行・検証**：全チェックPASS。t=10で E `3.3041e-3`、W `3.1157e-4`。W_max `4.0437e-4`（t=7.65）、分類 `peak_resolved`。fixed-per-site noisyでの W_max は N3 `3.0302e-3` > N5 `1.0968e-3` > N7 `4.0437e-4`（N7/N5=`0.3687`）。**noise-freeで見えた「N7 W_max > N5 W_max」という特徴は、この fixed-per-site noisy 条件では残らなかった。**
- **未確認**：dt半減、t>10、位置別雑音、保護、gamma/Omega sweep、N>7、実機性能。判定 **completed_comparison**。この段階では、Nとともに距離と雑音site数（したがって総雑音）が同時に増えることに注意（等総雑音比較ではない）。

### Milestone 9c
- **目的**：全site gammaの単純和を `TOTAL_GAMMA=1.5` に固定し、「N増加」と「総雑音増加」の交絡を部分的に切り分ける（fixed-total-noise比較）。
- **追加・変更**：N=5（gamma_site=0.3）とN=7（gamma_site=1.5/7）を新規本計算。N=3は既存 gamma_site=0.5 結果を参照。
- **実行・検証（この段階の途中経過）**：N=5はchecks=true。N=7では `t=0.02` の最小固有値診断1点が `NaN` となり `positivity` と `finite_values` の2検査がFAIL。**この段階のレポート `MILESTONE_9C_REPORT.md` は途中経過として `numerical_issue_stop` と判定した**（後述のとおり、これは最終結論ではない）。他の1000保存点の最小固有値は有限で最小 `-5.278e-18`、trace/Hermiticity/ledgerは正常だった。
- **未確認（この段階）**：`NaN` の原因切り分け、および比較結果を正式採用してよいかの判断は、この段階では未完了（診断・検証で解決）。

### Milestone 9c diagnostic
- **目的**：N=7 `t=0.02` の最小固有値 `NaN` を、短時間再計算だけで切り分ける（`t=10` 本計算は再実行しない）。
- **追加・変更**：`t=0`〜`0.03` を12 RK4 stepで3回再計算する診断bin。
- **実行・検証**：`t=0.02` で密度行列 rho は**全要素有限**、trace `≈1`、Hermiticity誤差0、3回再現性一致。原因は**固有値ソルバー側**：raw / Hermitianized SymmetricEigen が非有限固有値を返し、9cのCSV formatterがそれを一律 `NaN` と記録していた。独立solver（Complex Schur）は全固有値有限で minimum `-3.237e-24`。24 scalar比較の最大絶対差0。
- **未確認（この段階）**：`t=0.03` 以降の再構築、別dt/gamma、t=10再計算。**この診断はあくまで補足記録であり、元レポートを上書きしていない。診断ラベルは `state_level_numerical_issue`（停止規則によるラベル）で、rho自体にNaN・trace異常・Hermiticity異常・非決定性があったという意味ではない。**

### Milestone 9c validation（9cの最終正本）
- **目的**：state finiteness と solver finiteness を独立に判定するrobust positivity診断で、N=7 fixed-total-noise を t=10 まで再検証する。
- **追加・変更**：primaryをHermitian化SymmetricEigen、不合格時のみHermitian化Complex Schurへfallbackするsolver policy。物理模型・時間発展は不変（fallbackは診断層のみ）。
- **実行・検証**：条件はN=7、total gamma=1.5、gamma_site=1.5/7、`dt=0.0025`、4000 RK4 steps、1001保存点。下記「現在の主要結果」に検証値を掲載。**最終判定 `completed_comparison_with_fallback_diagnostic`。物理比較結果を正式採用してよい。**
- **未確認**：別dt、別gamma、t>10、N>7、別の外部solver crate。10%帯は統計的有意差ではなく、距離だけの因果証明でもない。

---

## Milestone 1–3：基礎実装

Milestone 1〜3.1は、時間発展を信頼する前の土台づくりです。まず静的モジュール（演算子・部分トレース・ergotropy）を確定し（M1）、次にcolumn-majorベクトル化に基づくLiouvillianとdense行列指数伝播器を追加し（M2）、24次元の完全模型が組み合わさることを統合テストで確認しました（M2.1）。その後、状態診断と符号付きpower会計を導入し（M3）、符号反転区間をゼロ交差で分割するsigned-power積分の扱いを明文化しました（M3.1）。

この段階の伝播器（`DenseExponentialPropagator`）は各時刻で `exp(L t) vec(rho(0))` を独立に計算する精度優先の実装で、`576 x 576` の行列指数は重い一方、正しさの基準として使えます。効率的な計算はより後の段階の課題として意図的に先送りされています。

---

## Milestone 4–6：比較実験と仕事抽出

ここからが本題の比較です。Milestone 4では初期コヒーレント状態の単発輸送で、同一時刻・同一load energyに一致した16条件すべてで雑音なしAが位相雑音ありBを上回りました。Milestone 5aで時間依存RK4伝播器を検証し、5bで真空からの有限パルスが非ゼロのload ergotropyを生むことを確認しました。

中心的な結果はMilestone 5cです。t=10でA/Bのload energyを相対誤差 `4.001e-5` に一致させたとき、ergotropy比は **A/B=6.373** でした。ただしこの一致は「最終load energy」についての条件付き一致であり、駆動強度Omegaと総投入エネルギーはA/Bで異なります。したがってこれは等入力費用比較でも、雑音だけを単独で変えた因果比較でもありません。

Milestone 6aは新しい物理発見ではなく、ergotropyの定義に対応する理想局所unitaryを構成して予測どおりの仕事量が回収できることを示した**実装検算**です。

（この4–6a区間の要約は `Milestone_4-6a_研究結果ノート.pdf` にもまとめられています。）

---

## Milestone 7：雑音位置と部分保護

固定N=3模型で、雑音を「どこに置くか」「どこから外すか」を調べた区間です。7aは雑音を1サイトだけに置く有害配置比較で、t=10のW最小は入口（site1）、usable fraction最小は出口（site3）でした。7bはその時系列から差が持続的に現れる時刻を記述しました（新しい時間発展はしていません）。

7cは逆に、全サイト雑音から特定サイトの雑音を理想的に除去したときの回復上限を比較し、両端保護（中央雑音のみ残す）がW・usable fraction・Eいずれの回復でも最大でした。7dは中央gammaを固定して両端gammaを `0.5→0` に下げ、回復が離散的に単調非減少でつながることを確認しました。いずれも理想的な雑音除去であり、現実の保護装置・cost・不完全保護ではありません。

---

## Milestone 8：鎖長一般化とN=7実行可能性

鎖長を一般化した区間です。8aでN=3からN=5へ鎖長だけを変え、load ergotropyの到達とW_maxの縮小を確認しました。8bはN=7のdense法での実行可能性をprobeし、現行dense実装ではt=10まで約23.7時間かかる（infeasible）と判定しました。8cは物理模型を変えずにsigma_z dephasing項だけを厳密な成分表示（`DiagonalDephasingKernel`）に置き換え、N=3/N=5でdense pathと一致することを確認したうえで約30倍の高速化を得て、N=7 t=10を feasible_candidate に更新しました。

鎖長を変えると、伝播距離・Hilbert次元・bond数・（all-siteでは）総雑音が同時に変わる点に注意が必要です。2〜3点の有限長比較から距離だけの因果やscaling則は主張しません。

---

## Milestone 9：N=7本計算と固定総雑音比較

8cの高速化を使って、N=7の本計算に進んだ区間です。

9aはN=7 noise-freeをt=10まで本計算し、load energy・ergotropyの到達を確認しました。9bはN=7 all-site noisy（各site gamma=0.5、fixed-per-site）を本計算しました。この fixed-per-site 条件では総雑音がN3=1.5、N5=2.5、N7=3.5と鎖長とともに増えるため、noise-freeで見えた「N7 W_max > N5 W_max」という特徴は残りませんでした（N7/N5=0.369）。

9cは、この「N増加」と「総雑音増加」の交絡を部分的に切り分けるため、全site gammaの単純和を `TOTAL_GAMMA=1.5` に固定した比較（fixed-total-noise）を行いました。この9cについては、以下の順序で経緯を整理します。

1. **当初の診断で固有値ソルバーの失敗が見つかった。** 9c本計算のN=7で、`t=0.02` の最小固有値診断1点が `NaN` となり、`positivity` と `finite_values` の2検査がFAILした。このため `MILESTONE_9C_REPORT.md` は途中経過として `numerical_issue_stop` と記録した。

2. **厳密診断を実施した。** `MILESTONE_9C_DIAGNOSTIC.md` で `t=0.02` を短時間再計算したところ、密度行列 rho は全要素有限、trace `≈1`、Hermiticity誤差0、3回再現一致だった。`NaN` は**状態ではなく固有値ソルバー側**の非有限出力を、CSV formatterが一律 `NaN` と記録していたことに由来する。独立のComplex Schurでは全固有値が有限だった。

3. **対称固有値計算に失敗した2点ではSchur分解によるfallbackを使用した。** `MILESTONE_9C_VALIDATION.md` で、Hermitian化SymmetricEigenをprimaryとし、失敗した2時刻（`t=0.01`、`t=0.02`）だけHermitian化Complex Schurへfallbackするrobust診断を導入した。fallbackは診断層のみで、物理的な時間発展は変更していない。

4. **1001時刻すべてでpositivity判定が完了した。** primary成功999、失敗2、fallback成功2/2、両solver失敗（solver_failure）0。worst selected minimum eigenvalue `-5.278e-18`。

5. **既存9c軌道との比較は最大差0だった。** 指定物理量を1001時刻・許容値 `1e-12` で比較し、最大絶対差はすべて0（all checks=true）。

6. **最終状態は `completed_comparison_with_fallback_diagnostic`。**

7. **したがって、固定総位相雑音条件でのN=3、N=5、N=7比較は、正式な比較結果として採用できる。**

以上より、9cの**最終的な正本は `MILESTONE_9C_VALIDATION.md`** です。`MILESTONE_9C_REPORT.md` の `numerical_issue_stop` は途中経過であり、最終結論ではありません。`MILESTONE_9C_DIAGNOSTIC.md` は途中診断の記録です。

---

## 現在の主要結果

いずれも「この模型・この初期条件・この有限時間・このRust実装で直接確認した数値」です。一般法則ではありません。

**Milestone 5c（N=3、条件付きload-energy一致、t=10）**
- ergotropy A/B比：`6.373`（A `5.2798e-2` / B `8.2846e-3`）、load energy相対差 `4.001e-5`。
- 一致しているのは最終load energyのみ。Omegaと総投入エネルギーは不一致。

**Milestone 9a/9b（N=7、fixed-per-site 参考）**
- noise-free W_max：`2.2436e-2`（t=7.71）。
- all-site noisy（gamma=0.5）W_max：`4.0437e-4`（t=7.65）。
- fixed-per-site の W_max：N3 `3.0302e-3` > N5 `1.0968e-3` > N7 `4.0437e-4`。noise-freeの「N7 > N5」はこの条件では残らない。

**Milestone 9c validation（N=7、fixed-total-noise、total gamma=1.5、最終正本）**
- positivity診断完了：**1001 / 1001時刻**、solver_failure=0。
- SymmetricEigen（primary）成功：**999**、失敗：**2**（`t=0.01`、`t=0.02`）。
- Complex Schur fallback成功：**2 / 2**。
- 既存9c軌道との最大差：**0**（1001時刻、許容 `1e-12`）。
- N=7 W_max：**0.0038080717406769921**（t=7.70）。
- fixed-total の W_max：N3 `3.0302e-3`、N5 `3.4876e-3`、N7 `3.8081e-3`。
- N=7 / N=5 W_max比：**1.0918825139**。この比は**記述的な「10%帯」の範囲内**。

なお、この「10%帯」は理論法則ではなく、**今回の有限サイズ比較のための記述的な目安**にすぎません。統計的有意差でも、距離だけの因果証明でもありません。

---

## ビルドとテスト

Rustツールチェーンのある環境で実行してください（`Cargo.toml`：edition 2021、依存 `nalgebra`, `num-complex`, `thiserror`, dev-dependency `approx`）。

```bash
# フォーマット確認
cargo fmt --all -- --check

# 通常テスト（リリース、オフライン）
cargo test --release --offline

# ビルド
cargo build --release --offline
```

意図的にignoredな `576 x 576` 24次元smoke testは、明示的に実行します。

```bash
cargo test --release full_24d_short_time_smoke_test -- --ignored --nocapture
```

---

## 主要な実行コマンド

各Milestoneは専用binで実行します（代表例。詳細な引数・出力は各レポート末尾の「生成ファイル一覧」を参照）。

```bash
# 5a: 時間依存RK4のsanity check
cargo run --release --offline --bin time_dependent_sanity

# 8c: 厳密dephasing kernelのベンチマーク・等価性
cargo run --release --offline --bin dephasing_kernel_benchmark

# 9a: N=7 noise-free 本計算
cargo run --release --offline --bin n7_noise_free_full

# 9b: N=7 all-site noisy 本計算
cargo run --release --offline --bin n7_all_site_noisy_full

# 9c: fixed-total-noise 比較
cargo run --release --offline --bin fixed_total_noise_comparison

# 9c diagnostic: t=0.02固有値NaNの切り分け
cargo run --release --offline --bin n7_t002_eigen_diagnostic

# 9c validation: robust positivity診断による最終再検証
cargo run --release --offline --bin n7_fixed_total_validation
```

N=7本計算は数十分規模の計算時間になります（9a total ≈ 2953s、9b total ≈ 2899s、9c validation total ≈ 2703s）。

---

## 出力ファイル

各段階はレポート（`MILESTONE_*.md`）と複数のCSVを生成します。READMEは入口であり、詳細な表はレポートとCSVを参照してください。代表的なもの：

- **5c**：`coherent_drive_match_{grid,roots,comparison,timeseries,convergence}.csv`
- **6a**：`explicit_load_extraction_{results,checks,mapping}.csv`
- **7a/7b**：`local_noise_placement_*.csv` / `local_noise_damage_*.csv`
- **7c/7d**：`ideal_partial_protection_*.csv` / `partial_end_protection_*.csv`
- **8a/8b/8c**：`chain_length_reachability_*.csv` / `n7_feasibility_*.csv` / `dephasing_kernel_*.csv`
- **9a/9b**：N=7 noise-free / all-site noisy の timeseries・summary・checks 等
- **9c validation（正本）**：
  - `n7_fixed_total_validation_summary.csv`（t=10値、W_max、solver会計、最終判定）
  - `n7_fixed_total_validation_checks.csv`（各チェックのPASS記録）
  - `n7_fixed_total_validation_trajectory_comparison.csv`（既存9c軌道との差、全0）
  - `n7_fixed_total_validation_eigen_diagnostics.csv`、`..._timeseries.csv`、`..._performance.csv`
  - `fixed_total_noise_final_comparison.csv`（N3/N5/N7 W_max比較）
  - `robust_eigen_diagnostic_unit_checks.csv`

**CSVの読み方の要点**：`W_time_area` や `E_time_area` は状態量の時間面積であり、累積流入エネルギーや累積抽出仕事ではありません。`W/Ein` は制御費用を含む総合効率ではありません。`usable_fraction` はload energyに対するergotropyの割合です。

---

## 数値上の注意

- 位相雑音時間発展はdense密度行列RK4で計算します。8cで導入した `DiagonalDephasingKernel` はsigma_z dephasing項の**厳密な成分表示**であり、物理近似ではありません（dense pathと許容内一致を確認済み）。
- RK4は完全正値写像を厳密には保証しないため、最小固有値を明示的に検査しています。観測される負値は丸め誤差の範囲（例：9c validationのworst `-5.278e-18`）で、固有値補正はしていません。
- 固有値の**診断**で非有限値が出た場合でも、それは密度行列そのものの異常を意味しません。9cでは、状態の有限性とソルバーの有限性を分離し、primary失敗時のみSchur fallbackを使う方針で1001時刻すべてのpositivity判定を完了しています（fallbackは診断層のみ）。
- ledger residualは分母がほぼゼロのため、相対差ではなく絶対差を主判定に使う場面があります。
- 鎖長Nを変えると距離・次元・bond数・（all-siteでは）総雑音が同時に変わります。有限個のNの比較から距離だけの因果やscaling則は導けません。

---

## 確認できていないこと

以下は、このプロジェクトが**まだ示していない**ことです。READMEで断定してはいけない事項でもあります。

- 量子優位性・量子超越、または古典方式より高性能であること。
- 普遍的なスケーリング則、指数/べき減衰、熱力学極限、N→∞の挙動。
- 「雑音が性能を改善する一般的機構」や、雑音位置・保護強度の普遍則。
- 現実の送電・蓄電・装置効率、実機での有効性。
- 抽出unitaryの制御費用、有限時間切断・抽出、繰り返しサイクル、連続供給、定常出力、長時間安定性。
- 別のdt・別のgamma・t>10・N>7・別の外部solver crate（9c validationの範囲外）。
- 古典的波動・確率模型との公平な基準比較。
- 文献に対する新規性・優先権（文献調査は未実施）。
- レポートに書かれていない因果説明（例：特定サイトの雑音が「注入/輸送/受け渡しだけ」を壊したという断定）。

9c validationでは、`n7_fixed_total_validation_checks.csv` のruntime検査47項目がすべてPASSした。また、`robust_eigen_diagnostic_unit_checks.csv` に記録されたrobust固有値診断の必須10項目もすべてPASSした。Cargoテスト全体について別途言及される「90/90」という総数は、今回README作成に使用した資料からは確認できないため、本READMEの検証件数には採用していない。

---

## リポジトリ構成

`Cargo.toml`（`quantum_work_network`、edition 2021）と `src/lib.rs` が公開する主なモジュール：

```
src/
  operators.rs            # 演算子（M1）
  partial_trace.rs        # 部分トレース（M1）
  ergotropy.rs            # 局所ergotropy（M1）
  matrix.rs               # ComplexMatrix / C64
  error.rs                # PhysicsError
  liouvillian.rs          # column-majorベクトル化・superoperator（M2）
  propagator.rs           # dense行列指数伝播（M2）
  diagnostics.rs          # 状態診断・signed power会計（M3 / M3.1）
  time_dependent.rs       # 時間依存RK4伝播器（M5a）
  coherent_drive.rs       # coherent駆動（M5b/5c 系）
  coherent_drive_matching.rs
  matching.rs / protocol.rs / experiment.rs
  load_extraction.rs      # 理想局所unitary抽出（M6a）
  dephasing_kernel.rs     # 厳密dephasing kernel（M8c）
  bin/
    time_dependent_sanity.rs        # M5a
    local_noise_placement.rs        # M7a
    local_noise_damage_analysis.rs  # M7b
    ideal_partial_protection.rs     # M7c
    partial_end_protection.rs       # M7d
    chain_length_reachability.rs    # M8a
    n7_feasibility_probe.rs         # M8b
    dephasing_kernel_benchmark.rs   # M8c
    n7_noise_free_full.rs           # M9a
    n7_all_site_noisy_full.rs       # M9b
    fixed_total_noise_comparison.rs # M9c
    n7_t002_eigen_diagnostic.rs   # M9c diagnostic
    n7_fixed_total_validation.rs  # M9c validation（最終正本）
tests/
  full_24d_short_time.rs            # M2.1（ignored smoke test）
MILESTONE_*.md                      # 各段階のレポート
*.csv                               # 各段階の出力
```

上記の `src/lib.rs` のモジュール宣言は実ファイルに基づきます。`bin/` 内の対応付けは各レポートの「生成ファイル一覧」に記載されたファイル名に基づく整理です。個々のCSVの完全な一覧は、対応するMILESTONEレポート末尾を参照してください。

---

## Milestone 10：fixed-total比較とXGamma診断

Milestone 10では、既存結果を整理したうえで、固定総位相雑音 `TOTAL_GAMMA=1.5` と `3.0` をN=3・5・7について同じ診断系で比較しました。

- **Milestone 10a**：新しい物理計算を行わず、既存結果を比較整理し、fixed-total条件の欠損を明示しました。
- **Milestone 10b**：`TOTAL_GAMMA=3.0` をN=3・5・7で計算し、XGammaを初導入しました。
- **Milestone 10c**：`TOTAL_GAMMA=1.5` をXGamma付きでN=3・5・7について再計算し、10aの欠損を補完しました。N=7の物理軌道は9c正本と1001時刻で一致しました。

fixed-total 1.5と3.0の両方で、`W_max`、`W(t=10)`、usable fractionは `N=7 > N=5 > N=3`、`W_time_area` は `N=3 > N=5 > N=7`、ergotropy arrivalはN=3が最速でした。したがって、鎖長の優劣を一括して決めるのではなく、**評価指標によって順位が逆転する**という有限条件の結果です。`W_time_area` はergotropy状態量の時間積分であり、累積抽出仕事ではありません。

XGammaは次のdephasing-kernel-weighted coherence exposureです。

```text
x_gamma(t) = sum_ab Gamma[a,b] |rho[a,b](t)|^2
XGamma(T) = integral_0^T x_gamma(t) dt
```

XGammaは診断量であり、失われた仕事、散逸エネルギー、dephasing power、熱、entropy production、効率、損傷量ではありません。今回の有限比較は、長い鎖の一般的優位性、scaling law、XGammaの因果機構、gamma倍増の普遍倍率、量子優位を示しません。
