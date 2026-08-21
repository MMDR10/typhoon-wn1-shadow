# 🌀 Typhoon Monitor 最新追蹤

**GFS 2026-08-20T18:00** — 自動更新 (2026-08-20 21:57 UTC)

| 颱風 | 位置 | WN1 相位 | ellipt | UQ | dH_curl | 強度模式 | 鞍點 n | D_fold |
|------|------|---------|--------|-----|---------|---------|--------|--------|
| EIGHTEEN | 20.1N 108.4E | 273.5° | 4.449 | 🔴 LOW | -2.32e-05 | organized | 179 | 0.003 |
| SAUDEL | 13.7N 149.6E | 305.2° | 0.572 | 🔴 LOW | -3.20e-05 | collapse | 202 | 0.012 |
| TWO-C | 11.7N -142.6E | 335.6° | 1.0 | 🔴 LOW | -1.90e-05 | organized | — | — |
| HURRICANE LALA | 23.8N -172.0E | 47.7° | 0.505 | 🔴 LOW | -2.99e-05 | organized | — | — |

**UQ 機制 v4**（2026-08-21 backtest 修正）：🟢🟢 Very High = ellipt≤0.4 + Amp≥緯度門檻（<20°→10 / 20-30°→7 / ≥30°→5）；🟡 Medium = ellipt≤0.6 + Amp≥7（v4 新增：backtest 證明 0.4-0.6 組 mean 11.1°、100% <45°）；🔴 Low = 其他。Amp≥7 係最強 UQ 指標（Live 95% <45° vs 回溯 94%）。

