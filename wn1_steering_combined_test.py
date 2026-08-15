#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌀 增量價值測試：環境引導流 + WN1/相位 合併 vs 單獨
=====================================================
核心問題：喺已經有環境引導流（steering）嘅情況下，加 WN1/相位有冇增量價值？

策略比較（同一 74 樣本，500 hPa）：
1. 純 Steering：steering_dir vs move_fwd
2. 純 WN1：wn1_phi vs move_fwd
3. Conditional：ellipt ≤ 0.4 → 信 WN1；ellipt > 0.4 → 信 Steering
4. Weighted vector：w_wn1 = 1-ellipt, w_steer = ellipt，向量合成方向
5. Oracle：每樣本揀誤差細嗰個（理論上限，唔係可實現）

指標：mean error、median、<45° 比例；per-storm 拆分
"""
import json
import numpy as np
from pathlib import Path

BASE = Path("/app/working/workspaces/tygtDc/projects/typhoon-dh-curl")


def ang_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def vec_dir(w1, a1, w2, a2):
    """兩個方向按權重合成向量方向（處理 wraparound）"""
    r1 = np.radians(a1)
    r2 = np.radians(a2)
    x = w1 * np.cos(r1) + w2 * np.cos(r2)
    y = w1 * np.sin(r1) + w2 * np.sin(r2)
    return (np.degrees(np.arctan2(y, x)) + 360) % 360


def main():
    steer_rows = json.load(open(BASE / "wn1_vs_steering_validation_500.json"))
    uq_rows = json.load(open(BASE / "wn1_uq_shape_v2.json"))
    uq_map = {(r["storm"], r["iso"]): r for r in uq_rows}

    rows = []
    for r in steer_rows:
        u = uq_map.get((r["storm"], r["iso"]))
        if u is None:
            continue
        rows.append(dict(
            storm=r["storm"], iso=r["iso"],
            wn1_phi=r["wn1_phi"], wn1_amp=r["wn1_amp"],
            steering_dir=r["steering_dir"], steering_mag=r["steering_mag"],
            move_fwd=r["move_fwd"], ellipt=u["ellipt"],
        ))
    n = len(rows)
    print(f"總樣本: {n}\n")

    # ── 各策略誤差 ──
    err_steer = [ang_diff(r["steering_dir"], r["move_fwd"]) for r in rows]
    err_wn1 = [ang_diff(r["wn1_phi"], r["move_fwd"]) for r in rows]

    # Conditional: ellipt ≤ 0.4 → WN1，否則 Steering
    err_cond = []
    for r in rows:
        if r["ellipt"] <= 0.4:
            err_cond.append(ang_diff(r["wn1_phi"], r["move_fwd"]))
        else:
            err_cond.append(ang_diff(r["steering_dir"], r["move_fwd"]))

    # Weighted vector: w_wn1 = max(1-ellipt, 0), w_steer = ellipt
    err_wvec = []
    for r in rows:
        ww = max(1.0 - r["ellipt"], 0.0)
        ws = r["ellipt"]
        comb = vec_dir(ww, r["wn1_phi"], ws, r["steering_dir"])
        err_wvec.append(ang_diff(comb, r["move_fwd"]))

    # Oracle（理論上限）
    err_oracle = [min(e1, e2) for e1, e2 in zip(err_steer, err_wn1)]

    strategies = [
        ("① 純 Steering", err_steer),
        ("② 純 WN1", err_wn1),
        ("③ Conditional (ellipt≤0.4→WN1)", err_cond),
        ("④ Weighted vector (w=1-ellipt)", err_wvec),
        ("⑤ Oracle 上限", err_oracle),
    ]

    print("=" * 78)
    print("📊 策略比較（vs 實際移動方向）")
    print("=" * 78)
    print(f"{'策略':<32}{'n':>4}{'Mean':>9}{'Median':>9}{'<45°':>8}")
    print("-" * 78)
    for name, errs in strategies:
        pct45 = sum(1 for e in errs if e < 45) / n * 100
        print(f"{name:<32}{n:>4}{np.mean(errs):>9.1f}{np.median(errs):>9.1f}{pct45:>7.1f}%")

    # ── 增量價值判斷 ──
    print("\n" + "=" * 78)
    print("🔬 增量價值判斷")
    print("=" * 78)
    m_steer, m_wn1 = np.mean(err_steer), np.mean(err_wn1)
    m_cond, m_wvec, m_oracle = np.mean(err_cond), np.mean(err_wvec), np.mean(err_oracle)
    print(f"  Steering 單獨:     {m_steer:.1f}°")
    print(f"  WN1 單獨:          {m_wn1:.1f}°")
    print(f"  Conditional 合併:  {m_cond:.1f}°  (改善 vs Steering: {m_steer-m_cond:+.1f}°)")
    print(f"  Weighted 合併:     {m_wvec:.1f}°  (改善 vs Steering: {m_steer-m_wvec:+.1f}°)")
    print(f"  Oracle 上限:       {m_oracle:.1f}°  (改善 vs Steering: {m_steer-m_oracle:+.1f}°)")

    # 每樣本邊個贏（err_wn1 < err_steer → WN1 贏）
    wn1_wins = sum(1 for e1, e2 in zip(err_wn1, err_steer) if e1 < e2)
    steer_wins = sum(1 for e1, e2 in zip(err_wn1, err_steer) if e2 < e1)
    print(f"\n  每樣本對決: WN1 贏 {wn1_wins}/{n} ({wn1_wins/n*100:.0f}%), Steering 贏 {steer_wins}/{n} ({steer_wins/n*100:.0f}%)")

    # Conditional 有幾多樣本用 WN1
    n_wn1_used = sum(1 for r in rows if r["ellipt"] <= 0.4)
    print(f"  Conditional 中 ellipt≤0.4 用 WN1 嘅樣本: {n_wn1_used}/{n} ({n_wn1_used/n*100:.0f}%)")

    # ── per-storm 拆分 ──
    print("\n" + "=" * 78)
    print("🌀 Per-storm 拆分")
    print("=" * 78)
    storms = {}
    for i, r in enumerate(rows):
        storms.setdefault(r["storm"], []).append(i)
    print(f"{'颱風':<10}{'n':>4}{'Steer':>9}{'WN1':>9}{'Cond':>9}{'WVec':>9}{'Orcl':>9}")
    print("-" * 60)
    for storm, idxs in sorted(storms.items()):
        e_s = [err_steer[i] for i in idxs]
        e_w = [err_wn1[i] for i in idxs]
        e_c = [err_cond[i] for i in idxs]
        e_v = [err_wvec[i] for i in idxs]
        e_o = [err_oracle[i] for i in idxs]
        print(f"{storm:<10}{len(idxs):>4}{np.mean(e_s):>9.1f}{np.mean(e_w):>9.1f}"
              f"{np.mean(e_c):>9.1f}{np.mean(e_v):>9.1f}{np.mean(e_o):>9.1f}")

    # ── 儲存 ──
    out = BASE / "wn1_steering_combined_test.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(dict(
            n=n,
            strategies={name: dict(mean=float(np.mean(e)), median=float(np.median(e)),
                                   pct_lt45=float(sum(1 for x in e if x < 45) / n * 100))
                        for name, e in strategies},
            per_storm={s: dict(n=len(idxs),
                               steer=float(np.mean([err_steer[i] for i in idxs])),
                               wn1=float(np.mean([err_wn1[i] for i in idxs])),
                               cond=float(np.mean([err_cond[i] for i in idxs])),
                               wvec=float(np.mean([err_wvec[i] for i in idxs])),
                               oracle=float(np.mean([err_oracle[i] for i in idxs])))
                       for s, idxs in sorted(storms.items())},
            per_sample=[dict(storm=rows[i]["storm"], iso=rows[i]["iso"],
                             err_steer=float(err_steer[i]), err_wn1=float(err_wn1[i]),
                             ellipt=rows[i]["ellipt"],
                             wn1_wins=err_wn1[i] < err_steer[i])
                        for i in range(n)],
        ), f, indent=2, ensure_ascii=False)
    print(f"\n💾 saved: {out}")


if __name__ == "__main__":
    main()
