#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌀 深入拆解：點解 mean Steering 大但逐樣本 Steering 贏 70%？
=============================================================
補充分析：
1. 誤差分佈（極端值影響 mean）
2. ellipt 分層：ellipt≤0.4 / >0.4 兩組入面 WN1 vs Steering 邊個贏
3. 有冇「現實可揀」嘅準則（唔係 Oracle）可以接近 Oracle 上限
"""
import json
import numpy as np
from pathlib import Path

BASE = Path("/app/working/workspaces/tygtDc/projects/typhoon-dh-curl")


def ang_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def main():
    rows = json.load(open(BASE / "wn1_steering_combined_test.json"))["per_sample"]

    err_s = np.array([r["err_steer"] for r in rows])
    err_w = np.array([r["err_wn1"] for r in rows])
    ellipt = np.array([r["ellipt"] for r in rows])
    n = len(rows)

    # ── 1. 誤差分佈 ──
    print("=" * 70)
    print("📈 誤差分佈（極端值對 mean 嘅影響）")
    print("=" * 70)
    for name, e in [("Steering", err_s), ("WN1", err_w)]:
        print(f"\n{name}:")
        print(f"  mean={e.mean():.1f}°  median={np.median(e):.1f}°  p90={np.percentile(e,90):.1f}°  max={e.max():.1f}°")
        print(f"  >90° 樣本數: {(e>90).sum()} ({(e>90).mean()*100:.0f}%)")
        print(f"  >120° 樣本數: {(e>120).sum()}")

    # ── 2. ellipt 分層 ──
    print("\n" + "=" * 70)
    print("🔍 ellipt 分層：每組入面邊個贏")
    print("=" * 70)
    for mask_name, mask in [("ellipt ≤ 0.4", ellipt <= 0.4), ("ellipt > 0.4", ellipt > 0.4)]:
        sub_s, sub_w = err_s[mask], err_w[mask]
        wn1_wins = (sub_w < sub_s).sum()
        print(f"\n{mask_name} (n={mask.sum()}):")
        print(f"  Steering mean={sub_s.mean():.1f}°  WN1 mean={sub_w.mean():.1f}°")
        print(f"  WN1 贏: {wn1_wins}/{mask.sum()} ({wn1_wins/mask.sum()*100:.0f}%)")
        print(f"  Steering 贏: {mask.sum()-wn1_wins}/{mask.sum()} ({(mask.sum()-wn1_wins)/mask.sum()*100:.0f}%)")

    # ── 3. 其他候選揀選準則 ──
    print("\n" + "=" * 70)
    print("🎯 候選揀選準則（現實可實現，非 Oracle）")
    print("=" * 70)
    # 全部樣本
    print(f"\n全部 (n={n}):  Steering={err_s.mean():.1f}°  WN1={err_w.mean():.1f}°")

    # 準則 A：ellipt ≤ 0.4 → WN1（已知 41.8°）
    e_cond = np.where(ellipt <= 0.4, err_w, err_s)
    print(f"準則A ellipt≤0.4→WN1: {e_cond.mean():.1f}°")

    # 準則 B：ellipt ≤ 0.4 且 WN1 贏 → 但現實唔知邊個贏，用「WN1 有 UQ 信號」代替：amp≥10 → WN1
    # 冇 amp 喺 per_sample，暫時跳過

    # ── 4. SAOLA 影響 ──
    print("\n" + "=" * 70)
    print("🌀 SAOLA（Steering 災難源頭）影響")
    print("=" * 70)
    saola = np.array([r["storm"] == "SAOLA" for r in rows])
    non_saola = ~saola
    print(f"SAOLA (n={saola.sum()}):  Steering={err_s[saola].mean():.1f}°  WN1={err_w[saola].mean():.1f}°")
    print(f"非SAOLA (n={non_saola.sum()}):  Steering={err_s[non_saola].mean():.1f}°  WN1={err_w[non_saola].mean():.1f}°")
    wn1_wins_non = (err_w[non_saola] < err_s[non_saola]).sum()
    print(f"非SAOLA 逐樣本: WN1 贏 {wn1_wins_non}/{non_saola.sum()} ({wn1_wins_non/non_saola.sum()*100:.0f}%)")


if __name__ == "__main__":
    main()
