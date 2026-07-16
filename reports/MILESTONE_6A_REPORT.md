# Milestone 6a explicit load work extraction

## 実装済み

Milestone 5cではfull density matrixを保存していなかったため、同一設定と`dt=0.0025`で確定状態を決定論的に再構成した。新しい探索・matching・実験ではない。`H_drive(10)=0`を明示検算し、突然切断後に理想的瞬時load-local unitaryを適用した。

再構成差 A energy/ergotropy `1.015e-13`/`4.463e-13`, B `4.228e-13`/`1.262e-13`。許容値 `1e-9`。

| quantity | A | B |
|---|---:|---:|
| load energy before | 5.4450767878e-2 | 5.4452946589e-2 |
| ergotropy before | 5.2798274942e-2 | 8.2846362481e-3 |
| switch work | 4.6917468130e-15 | -5.1513600007e-16 |
| gross extracted work | 5.2798274942e-2 | 8.2846362481e-3 |
| signed net work | 5.2798274942e-2 | 8.2846362481e-3 |
| conservative net work | 5.2798274942e-2 | 8.2846362481e-3 |
| load energy after | 1.6524929355e-3 | 4.6168310341e-2 |
| ergotropy after | 0.0000000000e0 | 0.0000000000e0 |

- gross work A-B `4.4513638694e-2`, ratio `6.37303478`
- switch work A-B `5.2068828131e-15`
- signed net work A-B `4.4513638694e-2`, ratio `6.37303478`
- conservative net work A-B `4.4513638694e-2`, ratio `6.37303478`
- gross work recovery fraction A/B `0.70102397` / `0.03208961`
- conservative work recovery fraction A/B `0.70102397` / `0.03208961`

全検算: A `18/18 PASS`, B `18/18 PASS`。相互情報量変化 A `9.714e-16`, B `9.992e-16`。

閾値: unitarity/trace/Hermiticity `1e-10`; positivity `-1e-8`; state/eigenvalue and energy agreement `1e-9`; post ergotropy `1e-9`; mutual information `1e-9`; degeneracy `1e-12`; entropy negative-eigenvalue clip `1e-10`。

## 未確認

抽出unitaryの制御費用、有限時間切断・抽出、繰り返し放電、連続運転、正味周期仕事、相関からのglobal extraction、古典比較、量子優位は未確認。switch workのsigned/conservative定義は理想化であり、装置費用を完全評価していない。
