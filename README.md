# 🌀 Typhoon Monitor — WN1 路徑 + dH_curl 強度 + 鞍點環 一體化追蹤

**WN1 Shadow** 用 500 hPa 徑向風 wavenumber-1 相位（出流非對稱影子）預測活躍颱風移動方向，GitHub Actions 每 6 小時自動測量。

**Paper / 論文：** [EN](paper_EN/paper_wn1_shadow_en.md) · [ZH](paper_ZH/paper_wn1_shadow_zh.md) · PDF ([EN](paper_EN/paper_wn1_shadow_en.pdf) / [ZH](paper_ZH/paper_wn1_shadow_zh.pdf)) · [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21862294.svg)](https://doi.org/10.5281/zenodo.21862294)

## 🔬 方法

颱風移動主要由 **500 hPa 引導流**控制（76 樣本前瞻驗證：mean 23.9°，隨機 90°；ellipt≤0.4 後 21.1°，<45° 命中 90%）。WN1 = 徑向風傅立葉分解一階模態相位，代表颱風結構非對稱（影子）方向；**相位指向移動方向**時可作前瞻訊號。

**UQ 門檻：** ellipt ≤ 0.4（ellipt = WN2/WN1，細 = 結構單純）

**多維度監測（2026-08-11 起）：** 每 6h 同時收集三類數據 —

| 維度 | 指標 | 物理意義 |
|---|---|---|
| 🧭 路徑 | WN1 相位（500 hPa） | 出流非對稱影子 → 移動方向 |
| 💪 強度 | dH_curl = H_shell − H_core（850 hPa ζ，core=5°/shell=8°） | 負 = 組織化結構；正 = 發散/genesis |
| 🔗 鞍點環 | 眼牆環帶 (1.5°, 3.0°) 鞍點點集 box-count D_fold | D_fold≈0 且 z_1D 顯著 = 鞍點凝聚成環（渦旋特異結構） |
| 🔮 轉向 onset | WN1 相位連續扭轉（amp gate + 6-18h 滯後跟隨） | 結構 re-orientation 起動 → 預測移動 6-18h 內跟隨轉向 |

## 🔮 轉向 onset 測試方法原理（2026-09-01 起）

**研究問題：** 颱風「回馬槍」/「急轉向」點樣提前分辨？主流模型視颱風移動為「環境引導流推著走」，本 repo 嘅核心框架反轉呢個假設——

> **颱風移動由結構慣性主導，大氣（steering）只係干擾源。**（框架轉移，§13）

### 核心發現（三層）

**① 移動跟 WN1（結構），唔跟 steering（大氣）**
76 樣本 70 對 6h 段量度：移動方向 ↔ WN1 相位同期相關 **r=0.695 / dCor 0.521（p=0.002）**；移動方向 ↔ steering 三個 lag 全部 ≈0。即係話颱風「行邊度」由自身結構非對稱（WN1）決定，唔係環境風直接推。

**② 轉向三步次序**（19 個大轉事件 composite）

```
環境先轉 → 結構抵抗 0-6h（順時針更明顯 −20°）→ WN1 6-18h 追 46° → 移動最後跟
```

關鍵：onset 觸發後**即刻 6h 移動反而反方向**（64% 反號，Fisher p=0.018）——呢個「抵抗期」係實錘，意味住用即刻窗判「跟唔跟隨」係方法錯誤，必須用 **6-18h 滯後窗**。

**③ 弱 amp 偽影**（SAOLA 案例）
WN1 phase 嘅可信度由 amplitude 決定。Bootstrap 校準：amp<1.0 時 σ_phase ≈ 32°（phase 幾乎係雜訊）；amp≥1.5 時 σ ≈ 8-9°。SAOLA 三段「+84°/+84°/+122° 爆轉」實錘係弱 amp 偽影——因為同一時段強 amp 嘅 850hPa core 同移動都連續逆時針轉，而弱場 500hPa 順時針爆跳，方向矛盾。**故 onset 偵測必須加 amp gate ≥2.0。**

### 偵測規則（`detect_onset_live`，turn_onset_detector v2）

```
onset = |ΔWN1| > 15° 連續 2 段（同向）+ 兩端 amp ≥ 2.0
命中 = onset 後 6-18h（t+1 或 t+2）移動同方向大轉 > 15°
```

### 驗證結果（amp gate 後）

| 窗 | onset 同號命中 | 對照 | Fisher p |
|---|---|---|---|
| t+0 即刻 | 0%（反號 60%） | 0% | — |
| t+6h | 20% | 8% | 0.387 |
| t+12h | 50% | 14% | 0.128 |
| **合併 6-18h** | **75%（3/4）** | 16% | **0.028 ✅** |

⚠️ 樣本細（onset n=5、合併窗 n=4），p=0.028 係方向性提示，非最終實錘——live 數據累積中，預期逐 cycle 驗證「跟隨中」onset 嘅命中率會否收斂到歷史 75%。

## 📊 Key Numbers（論文核心數字）

| 指標 | 數值 |
|---|---|
| 樣本（6 颱風：HATO/MANGKHUT/MERANTI/HAGIBIS/SAOLA/GONI） | 76 |
| 12h 前瞻平均誤差（500 hPa） | **23.9°**（隨機 90°） |
| UQ 過濾後（ellipt≤0.4） | **21.1°**，<45° 命中 **90%** |
| Lead 退化（6h→48h） | 20.1° → 26.6°（≈0.16°/h，可用窗口 ≥48h） |
| 層級對比 | 850=52.0° / 200=31.8° / **500=23.9°（最優）** |
| 轉向警報（相位跳變） | ρ(Δφ,Δmotion)=0.589 |

## ⚙️ 自動化流程

每 6h（03:40/09:40/15:40/21:40 UTC）GitHub Actions 自動：

1. **cyclocane** → 活躍颱風列表 + JTWC advisory 位置
2. **NOMADS GFS** 0.25° 分析場（500 + 850 hPa 一次過）
3. 對每個颱風計算：
   - **WN1 相位 + ellipt**（500 hPa，路徑）
   - **dH_curl**（850 hPa ζ，強度）
   - **鞍點環 D_fold**（850 hPa ζ Morse 分類 + box-count，結構）
   - **轉向 onset**（WN1 連續扭轉 + amp gate + 6-18h 滯後跟隨驗證）
4. 冪等寫入 `typhoon_history.json`
5. 自動 commit push（5 次 retry）

## 🗂️ 資料結構

```
├── paper_EN/          # 英文論文（md + pdf）
├── paper_ZH/          # 中文論文（md + pdf）
├── scripts/           # 可重現分析 scripts（前瞻驗證/UQ/lead-time/層級對比/track）
├── output/            # 結果 JSON（76 樣本完整記錄）
├── figures/           # 論文圖（CDF/lead-time/層級對比/dφ）
├── typhoon_monitor.py # GitHub Actions 一體化追蹤 script（WN1 + dH_curl + 鞍點環）
├── typhoon_history.json # 實時三維數據時間序列（颱風 × cycle）
├── typhoon_latest.md  # 最新一輪摘要
├── wn1_track.py       # 舊版單維度追蹤（保留參考）
├── wn1_history.json   # 舊版 WN1 相位時間序列（歷史兼容）
├── wn1_latest.md      # 舊版摘要
└── .github/workflows/wn1_track.yml
```

## 📡 數據出處（Data Provenance）

| 數據 | 出處 | 用途 |
|---|---|---|
| ERA5 再分析（0.25°，850/500/200 hPa u/v） | Copernicus CDS（Hersbach et al. 2020） | 76 樣本前瞻驗證 |
| 颱風最佳路徑（6 颱風） | JTWC / HKO 綜合 | 移動方向真值 |
| GFS 作業分析（0.25° 500/850 hPa） | NOMADS（NCEP） | 實時自動追蹤 |
| 活躍颱風列表 | cyclocane（JTWC advisory） | 自動追蹤目標 |

原始 ERA5/GFS 場唔打包（容量），可經 scripts 重拉重現。

## 🚀 手動觸發

```bash
# 指定 cycle
gh workflow run wn1_track.yml -f cycle=20260811,00

# 覆寫颱風位置（測試）
gh workflow run wn1_track.yml -f storms='CHAN-HOM:36.5,142.0;PEILOU:26.0,158.0'
```

## 📋 相關

- 方法詳情：`scripts/wn1_forward_validation.py` / `scripts/wn1_uq_shape.py`（ERA5 76 樣本驗證）
- 轉向 onset 偵測：`turn_onset_detector.py`（v2，amp gate + 6-18h 滯後窗）；`typhoon_monitor.py::detect_onset_live`（live 自動輸出）
- 論文：`paper_EN/` `paper_ZH/`（中英雙語、數據出處全列明）
- Dolphin 自動化參考：`MMDR10/dolphin-watch`（dH_curl 追蹤）
- 鞍點環方法：Ô-HAT 框架 phase15（眼牆鞍點環 D_fold 簽名）

## 📄 License

CC BY 4.0
