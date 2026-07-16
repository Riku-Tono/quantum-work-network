# Milestone 5a implementation report

## 実装済み

- `src/time_dependent.rs`
  - 密度行列を直接更新する固定最大刻みRK4伝播器
  - 時刻依存の `H(t)` と時刻依存collapse operators `L_k(t)`
  - 保存方式：毎ステップ、一定間隔、指定時刻
  - 終了時刻・保存時刻へ到達するための短い最終ステップ
  - 時間依存Lindblad右辺を行列形式で実装
- `src/lib.rs`
  - `time_dependent` モジュールを公開
- `src/bin/time_dependent_sanity.rs`
  - 1量子ビットの滑らかな `sin^2` パルスによる動作確認
  - 時間不変問題でのDense exponentialとの誤差測定
- `MILESTONE_4_RESULT.md`
  - Milestone 4 coherent-input確定結果と限界を固定
- `src/time_dependent.rs` 内の8テスト
  - 時間不変generatorとの一致と収束
  - 時間依存項ゼロで既存ネットワークHamiltonianとの一致
  - トレース、エルミート性、正値性、有限値
  - `t=0`、指定保存時刻、非整数倍終了時刻、次元エラー

既存の `DenseExponentialPropagator`、物理モデル、演算子、診断量、
coherent-input実験コードおよび既存CSVは変更していない。

## 新しいAPI

```rust
let solver = TimeDependentRk4::new(dt)?;
let states = solver.propagate(
    &rho0,
    t0,
    t_end,
    |t| hamiltonian_at(t),
    |t| collapse_operators_at(t),
    SaveSchedule::Interval(save_interval),
)?;
```

返り値は既存の `QuantumState` の列で、各要素に保存時刻と密度行列を含む。
collapse operatorは既存規約と同じく係数を含めて渡す。

標準動作には、再正規化、エルミート化、固有値clip、対角要素補正、純化を
一切含めていない。

## 実行確認済み

### 時間不変参照問題

比較対象：`DenseExponentialPropagator`、終了時刻 `0.7`。

| RK4最大刻み | Frobenius差 |
|---:|---:|
| `dt = 0.08` | `7.8023058757e-7` |
| `dt/2 = 0.04` | `4.9923346793e-8` |
| `dt/4 = 0.02` | `3.1902121214e-9` |

- 3刻み中の最大差：`7.8023058757e-7`
- 誤差比 `error(dt) / error(dt/2)`：`15.6286`
- 誤差比 `error(dt/2) / error(dt/4)`：`15.6489`
- 刻み半減で単調に誤差が減少し、RK4の期待値16に近い。

既存ネットワークと同じ時間不変Hamiltonianに時間依存項ゼロを与えた
別テストも、Dense exponentialとの差 `1e-8` 未満で合格した。

### 1量子ビットsanity check

```text
H(t) = Omega_0 sin^2(pi t / T) sigma_x
Omega_0 = 0.7
T = 2.0
RK4 dt = 0.002
```

- 最終励起確率：`0.4150164285`
- 最大非対角成分：`0.4927248650`
- ゼロ駆動時の状態変化：`0.0`
- 最大トレース誤差：`1.3322676296e-15`
- 最大エルミート性誤差：`0.0`
- 最小固有値の最悪値：`-1.5388686396e-16`
- NaN / 無限値：なし

正値性のテスト許容値は通常診断で `-1e-11`。RK4は完全正値写像を厳密には
保証しないため明示的に検査した。sanity checkの負値は丸め誤差の範囲であり、
固有値補正はしていない。別の刻み幅確認 `0.2, 0.1, 0.05` では、保存時刻中の
最小固有値はいずれも `0.0` で、刻み半減による悪化はなかった。

### コマンド

- `cargo fmt --all -- --check`：合格
- `cargo test --release --offline`：47 passed / 0 failed / 1 ignored
- `cargo build --release --offline`：合格
- `cargo run --release --offline --bin time_dependent_sanity`：合格

ignoredの1件は、既存の明示実行用 `576 x 576` dense exponential smoke test。
既存の未使用importとdead-codeの警告2件は残るが、新規エラー・警告はない。

## 未確認

- 本番量子ネットワークでの時間依存coherent drive
- 本番A/B比較、energy matching、`p_B`探索
- 連続供給、パラメータスイープ、適応刻み
- RK4が任意の大きな刻み・強い駆動で正値性を保つこと
- 現実の送電性能または古典方式に対する量子優位

これらはMilestone 5aの成功条件には含めていない。
