# 📌 總結記錄 — WN1 颱風監測 × 環境引導流（2026-08-15）

> **以本文件為準**（MKP 指示 2026-08-15：下次 session 以總結記錄為準，唔好重做已確認嘅測試）
> 完整技術細節：`notes/2026-08-15-wn1-steering-combined-test.md`、`notes/2026-08-15-wn1-user-manual.md`、`notes/2026-08-15-wn1-nonlinear-redo.md`

---

## 一、WN1 係咩（測量本質）

**WN1 = 500 hPa 內環（r=1-5°）徑向風嘅 wavenumber-1 FFT 相位** — 量度颱風內核偶極子嘅方位（結構量，唔係移動量）。

- 純幾何測量：**100% 準**（量到咩就係咩）
- **唔係直接測移動方向**：bias +19.8° 證明（結構北 ≠ 移動北）
- 2D 單層測量（500 hPa 單層 = 垂直平均場嘅良好代理，82% <30° 一致）
- 三層（850/500/200）垂直一致率僅 13% — 唔係統一 3D 結構，500 hPa 係最佳單層

## 二、測得準咩（數據結論）

### 移動方向預測（條件性）

| 方法 | Mean | <45° | 註 |
|------|------|------|-----|
| **WN1 + LOO bias 校正** | **19.4°** | **90.5%** | 🏆 目前最佳 |
| WN1 高品質（ellipt≤0.4 + amp≥10） | ~15-20° | 94% | UQ 門檻內 |
| WN1 全體（500 hPa） | 22.6-23.9° | 83.8% | 74-76 樣本 |
| 純 Steering（r=8-12° 平均風向量） | 54.3° | 56.8% | 有災難長尾 |

### 診斷（新發現）

- **WN1-Steer 分離 >30° = 環境引導失效警號**（Mann-Whitney p=0.0008）：
  - WN1 ≈ Steering（<30°）：兩者都準（~14°），打和
  - WN1 ≠ Steering（>30°）：Steering 完全失效（56-134°），WN1 仍 carry 訊號（18-42°）
- **三種失敗模式有 signature 可提前認出**：
  - GONI 型（急增強 RI）：ellipt 高 + 相位穩定但偏 60-80°
  - SAOLA 型（急轉向/數字 6 路徑）：ellipt 高 + 相位跳動 + turn >40°/12h
  - CHAN-HOM 型（轉向點）：環境流快速旋轉（3.3× 快過移動），WN1 超前 90°+（ellipt UQ 盲點）

## 三、WN1 × Steering 合併測試結論（今日核心）

**「加」冇用，「取代」先有用。**

1. **逐樣本對決：WN1 贏 52/74 (70%)**；WN1 喺 5/6 颱風贏（MANGKHUT 打和 15.5 vs 16.0°）
2. **合併策略反而差過純 WN1**：Conditional (ellipt≤0.4→WN1) = 41.8°、Weighted vector = 36.3°，都差過純 WN1 22.6° — ellipt 唔係好嘅「揀邊個」準則（ellipt>0.4 組 WN1 贏率反而更高 81%，因 SAOLA 全喺嗰組而 Steering 平均錯 122°）
3. **Steering 有極端長尾**（p90=149.9°、11 樣本 >120°）；WN1 全部 <90°（max 89.0°）→ WN1 更準 + 更穩健
4. Oracle 上限 18.8°（每樣本揀啱嗰個）→ 兩者有理論互補空間，但現有 selector 搵唔到

**最終策略：直接用「WN1 + LOO bias 校正」（19.4°/90.5%）就係最佳**，唔使做合併 selector。
**Steering 嘅真正價值 = 診斷工具**：WN1-Steer >30° = 唔好信 steering 為主嘅預報。

## 四、LOO Bias 校正（可以落地）

- 全樣本 signed bias = **17.7°**（circ std 25.6°）
- **LOO 校正：22.6° → 19.4° / 90.5%**；LOO ≈ in-sample（19.1°）→ **唔係 look-ahead 假象，係真改善**
- 42/74 樣本改善，平均 +3.2°
- 修正之前「bias 唔穩定唔落地」結論（seq 28349）：bias 按颱風/時間分組有 spread，但**整體平均 bias 有系統性**，校正穩健改善

## 五、測咗但係負結果 / 已撤回

| 項目 | 狀態 | 結論 |
|------|------|------|
| 地磁場調制 WN1 amplitude | ⛔ 撤回 | r=0.752 係緯度偽相關（corr(F,lat)=0.924，partial r≈0 甚至反號） |
| 鞍點環 → 移動方向 | ❌ 負結果 | 四假說全不顯著，環對稱穩定，同移動無關 |
| 850 hPa WN1 | ❌ 淘汰 | 太近內核被對流污染（52.0°），200 hPa fallback |
| SST 梯度 → 移動 | ❌ 負結果 | 颱風唔係向暖水移動（前方冷 0.28°C，梯度垂直移動方向 122.5°） |
| 大氣波動假說 | ❌ 否定 | WN1 係準靜態偶極子，WN2-4 冇方向資訊 |
| dphi vs dsteer 零相關 | ⚠️ 修正 | 500 hPa 實為弱正相關 dCor≈0.30 (p≈0.04) |

## 六、使用方式（實戰 5 步速查）

1. 攞 500 hPa 風場（GFS 分析場 live / ERA5 歷史）
2. 計算 WN1 相位（r=1-5° 徑向風 FFT）+ ellipt + amp
3. **LOO bias 校正：相位 − 17.7°**
4. 過 UQ 燈號先可信：🟢 ellipt≤0.4 + amp≥緯度門檻（<20°→10 / 20-30°→7 / ≥30°→5）；🔴 ellipt>0.4 唔好用
5. 交叉檢查 WN1-Steer 分離：>30° = 環境引導失效，信 WN1；<30° = 兩者皆可

## 七、檔案索引

| 檔案 | 內容 |
|------|------|
| `notes/2026-08-15-wn1-steering-combined-test.md` | 合併測試 + LOO 校正完整細節 |
| `notes/2026-08-15-wn1-user-manual.md` | 使用說明書（速查表 + 失敗模式） |
| `wn1_steering_combined_test.py` / `.json` | 策略比較 script + 數據 |
| `wn1_steering_bias_loo.py` / `.json` | LOO 校正 + 特徵分析 script + 數據 |
| `wn1_vs_steering_validation_500.json` | 74 樣本原始數據（WN1/Steering/Move） |
| `wn1_uq_shape_v2.json` | ellipt UQ 數據 |

## 八、未解決 / 下一步候選

1. Oracle 18.8° vs WN1 22.6° 中間嘅現實 selector（「Steering 贏」樣本嘅特徵，除 WN1-Steer 角度差外）
2. MANGKHUT 組 Steering 贏 — 拆解「咩情況 Steering 反而好」
3. SAOLA 離群原因（bias +37.9°）
4. 框架 §7 ERRATUM 未標註（dphi-dsteer 修正）
5. 中/東太平洋擴展未決；數據源交叉（ERA5 vs GFS）留待 GFS 歷史場
