# He 1D PIC/MCC simulator

He ガス用の 1D3V 静電 PIC/MCC コードです。Rust 単体実装なので外部クレートなしでビルドできます。

## 実装内容

- 1D 空間、3D 速度の electrostatic PIC
- Dirichlet 電極境界、左電極は DC + RF 正弦波を設定可能
- 冷たい中性 He ガス近似、密度は `pressure_pa / (kB * gas_temperature_k)`
- `xsec/Cross section.txt` の LXCat 断面積を読み込み
  - 電子: elastic, excitation, ionization
  - He+: Phelps Backscat, Isotropic
- イオン壁衝突時の二次電子放出
- 電子壁衝突時の確率的反射
- Fowler-Nordheim 放出
- `output/diagnostics.csv` と `output/fields.csv` へ出力

## 実行

```powershell
cargo run --release -- config.example.ini
```

短い動作確認だけなら:

```powershell
cargo run --release -- config.smoke.ini
```

独自設定を作る場合:

```powershell
cargo run --release -- --write-config config.ini
cargo run --release -- config.ini
```

## 物理モデル上の注意

このコードは研究用の出発点として、保守しやすさと速度のバランスを優先しています。電子衝突は LXCat の断面積から直接 MCC 選択し、励起はしきい値分のエネルギー損失、電離は残余エネルギーを一次電子と二次電子へランダム分配します。He+ の Backscat は共鳴電荷交換に相当する近似として冷たい中性粒子速度へリセットし、Isotropic はエネルギー保存の等方散乱にしています。

時間刻みはプラズマ振動、セル幅は Debye 長に対して十分小さくしてください。粒子数が増えすぎる場合は `particle_weight`、`initial_particles_per_species`、`max_particles`、FN 放出係数を調整してください。
