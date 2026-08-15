# 2026-08-15 WN1 × Steering 合併增量測試

> 📌 **以總結記錄為準**：`notes/SUMMARY-2026-08-15-wn1-steering.md`（MKP 2026-08-15 指示：下次 session 以總結記錄為準，唔好重做已確認測試）

**任務來源：** MKP 問「環境引導流 加WN1或相位沒用?」→ DR 跑合併測試（74 樣本，500 hPa）

## 直接答案

**有用，但唔係「加」— 係「取代」。** WN1/相位 carry 獨立資訊，單獨已經贏 Steering；而「合併」策略暫時贏唔到純 WN1。

## 數據

- 74 樣本（`wn1_vs_steering_validation_500.json` + `wn1_uq_shape_v2.json` merge）
- WN1 相位：500 hPa 徑向風 r=1-5° W1 FFT
- Steering：500 hPa u,v r=8-12° area-averaged vector direction
- Move：t→t+12h IBTrACS 移動方向

## 結果

### 策略比較

| 策略 | Mean | Median | <45° |
|------|------|--------|------|
| ① 純 Steering | 54.3° | 31.2° | 56.8% |
| ② 純 WN1 | **22.6°** | **17.2°** | **83.8%** |
| ③ Conditional (ellipt≤0.4→WN1) | 41.8° | 22.7° | 68.9% |
| ④ Weighted vector (w=1-ellipt) | 36.3° | 22.0° | 73.0% |
| ⑤ Oracle 上限 | 18.8° | 12.4° | 87.8% |

### 關鍵發現

1. **逐樣本對決：WN1 贏 52/74 (70%)** — 唔係打和，係壓倒性
2. **⚠️ 合併策略反而差過純 WN1！** Conditional 41.8°、Weighted 36.3° 都比 WN1 單獨 22.6° 差
3. **ellipt 唔係好嘅「揀邊個」準則**：
   - ellipt≤0.4（n=43）：WN1 贏 27/43 (63%)
   - ellipt>0.4（n=31）：WN1 贏 25/31 (**81%**) — 反直覺！ellipt 差嗰組 WN1 贏率更高
   - 原因：ellipt>0.4 嗰組 Steering mean 73.2°（災難），WN1 27.5° — 因為 SAOLA 全喺嗰組
4. **Steering 有極端長尾**：p90=149.9°、11 個 >120°、17 個 >90°；WN1 全部 <90°（max 89.0°）
   - WN1 贏 70% 但 mean 贏咁多嘅原因：Steering 輸嗰陣輸得好甘（SAOLA 122°），WN1 輸嗰陣輸得少（max 89°）

### Per-storm

| 颱風 | n | Steer | WN1 | Cond | WVec | Orcl |
|------|---|-------|-----|------|------|------|
| GONI | 10 | 75.9 | **29.1** | 55.8 | 48.0 | 27.7 |
| HAGIBIS | 15 | 34.3 | **14.9** | 27.9 | 22.1 | 13.4 |
| HATO | 5 | 38.4 | **9.9** | 9.9 | 15.3 | 9.9 |
| MANGKHUT | 18 | **15.5** | 16.0 | 17.0 | 14.2 | 9.5 |
| MERANTI | 10 | 31.8 | **15.4** | 20.2 | 21.7 | 11.2 |
| SAOLA | 16 | 122.1 | **41.7** | 97.2 | 82.7 | 36.0 |

- WN1 喺 5/6 颱風贏；MANGKHUT 打和（Steer 15.5 vs WN1 16.0）
- Oracle 每個颱風都 ≤ 兩者單獨 → 兩者有互補空間（18.8° vs 22.6°），但現實規則未搵到

## 詮釋修正（重要）

之前（seq 28363）話 WN1 同 steering「完全獨立（零相關）互補」— 呢個測試精確化咗：

- **WN1 唔係要「加」落 Steering，係要「優先」**：任何用 ellipt 做 selector 嘅合併都輸俾純 WN1
- 真正嘅增量關係：WN1 捕捉到 Steering 嘅大部分訊號 + Steering 完全 miss 嘅內部動力學（MANGKHUT 85.9° 案例）
- **「環境引導流 + WN1」冇用；「環境引導流被 WN1 取代/校正」先有用** — 或者「Steering 做 fallback，WN1 做主」

## 下一步（未做）

1. ~~搵「Oracle 18.8° vs WN1 22.6°」中間嘅現實 selector~~ ✅ 已做（見下節）
2. ~~Steering 做 WN1 嘅 bias 校正~~ ✅ 已做（見下節）
3. 合併誤差函數（兩個方向都錯嗰陣點算）

---

## 追加：LOO bias 校正 + Steering 診斷（同日完成）

### A. LOO Bias 校正 — 有效，可以落地

| 方法 | Mean | Median | <45° |
|------|------|--------|------|
| 純 WN1（無校正） | 22.6° | 17.2° | 83.8% |
| **WN1 + LOO 校正** | **19.4°** | **14.8°** | **90.5%** |
| WN1 + in-sample 校正 | 19.1° | 14.5° | 90.5% |
| 純 Steering（對照） | 54.3° | 31.2° | 56.8% |

- LOO（19.4°）≈ in-sample（19.1°）→ **唔係 look-ahead 假象，係真改善**
- 全樣本 signed bias = 17.7°（circ std 25.6°）
- 42/74 樣本改善，平均 +3.2°
- **補充 seq 28349「bias 唔穩定唔落地」**：bias 按颱風/時間分組有 spread，但整體平均 bias 仍然有系統性，用嚟校正穩健改善

### B. 「Steering 贏」特徵 — WN1-Steer 角度差係唯一顯著 selector

Mann-Whitney p=0.0008（唯一 <0.05 特徵；steering_mag p=0.12、wn1_amp p=0.16、ellipt p=0.16 全部唔顯著）

| WN1-Steer 差 | n | Steering 贏 | WN1 err | Steering err |
|-------------|---|------------|---------|--------------|
| <30° 一致 | 37 | 51% | 14.4° | 14.1° |
| 30-60° | 15 | 7% | 17.7° | 56.5° |
| 60-90° | 9 | 11% | 36.0° | 100.6° |
| >90° 相反 | 13 | 8% | 42.3° | 133.9° |

**物理意義**：
- WN1 ≈ Steering（<30°）→ 兩者都準，打和（MANGKHUT 主場）
- WN1 ≠ Steering（>30°）→ Steering 完全失效（56-134°），WN1 仍 carry 訊號（18-42°）
- **WN1-Steer 分離 = 環境引導流失效警號**（颱風靠內部動力學移動，MANGKHUT 85.9° 案例）

### 最終策略建議

**直接用「WN1 + LOO bias 校正」（19.4°/90.5%）就係最佳** — WN1 喺每個分層都唔差過 Steering（<30° 打和 14.4 vs 14.1，≥30° 大勝），唔使做 selector。
Steering 嘅真正價值係**診斷**：WN1-Steer >30° = 環境引導失效 → 唔好信 steering 為主嘅預報。

- script: `wn1_steering_bias_loo.py`
- data: `wn1_steering_bias_loo.json`

## 檔案

- script: `wn1_steering_combined_test.py` / `wn1_steering_combined_deep.py`
- data: `wn1_steering_combined_test.json`
