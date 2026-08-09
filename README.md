# 🌀 WN1 Shadow — 活躍颱風 500 hPa WN1 相位前瞻自動追蹤

**WN1 Shadow** 用 500 hPa 徑向風 wavenumber-1 相位（出流非對稱影子）預測活躍颱風移動方向，GitHub Actions 每 6 小時自動測量。

## 🔬 方法

颱風移動主要由 **500 hPa 引導流**控制（76 樣本前瞻驗證：mean 23.9°，隨機 90°；ellipt≤0.4 後 21.1°，<45° 命中 91%）。WN1 = 徑向風傅立葉分解一階模態相位，代表颱風結構非對稱（影子）方向；**相位指向移動方向**時可作前瞻訊號。

**UQ 門檻：** ellipt ≤ 0.4（ellipt = WN2/WN1，細 = 結構單純）

## ⚙️ 自動化流程

每 6h（03:40/09:40/15:40/21:40 UTC）GitHub Actions 自動：

1. **cyclocane** → 活躍颱風列表 + JTWC advisory 位置
2. **NOMADS GFS** 0.25° 500 hPa 分析場
3. 計算每個颱風 **WN1 相位 + ellipt**
4. 冪等寫入 `wn1_history.json`
5. 自動 commit push（5 次 retry）

## 📊 數據

- `wn1_history.json` — 完整時間序列（颱風 × cycle）
- `wn1_latest.md` — 最新一輪摘要

## 🚀 手動觸發

```bash
# 指定 cycle
gh workflow run wn1_track.yml -f cycle=20260809,12

# 覆寫颱風位置（測試）
gh workflow run wn1_track.yml -f storms='CHAN-HOM:32.3,150.8;PEILOU:21.2,142.7'
```

## 📋 相關

- 方法詳情：`wn1_forward_validation.py` / `wn1_uq_shape.py`（ERA5 76 樣本驗證）
- Dolphin 自動化參考：`MMDR10/dolphin-watch`（dH_curl 追蹤）
