#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌀 下一步：LOO bias 校正 + 「Steering 贏」特徵分析
====================================================
背景：合併測試發現 WN1 贏 52/74 (70%)，但合併策略反輸純 WN1。
已知 WN1 bias = +19.83°（WN1 相位 vs 移動方向嘅 signed bias，in-sample）。

本 script：
A. LOO bias 校正 — 每個樣本用「其他 n-1 個樣本嘅 signed bias」校正自己
   → 誠實版（冇 look-ahead），睇校正後 WN1 誤差有冇改善
B. 對照：in-sample bias 校正（之前嘅做法，有 look-ahead，做對比）
C. 「Steering 贏」特徵分析 — 22 個 Steering 贏樣本 vs 52 個 WN1 贏樣本
   → storm / steering_mag / wn1_amp / ellipt / wind / WN1-Steering 角度差
"""
import json
import numpy as np
from pathlib import Path

BASE = Path("/app/working/workspaces/tygtDc/projects/typhoon-dh-curl")


def ang_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def signed_diff(a, b):
    """a - b 嘅 signed circular difference（-180, 180]"""
    d = (a - b) % 360
    if d > 180:
        d -= 360
    return d


def circular_mean(angles_deg):
    """角度嘅 circular mean（向量平均）"""
    r = np.radians(angles_deg)
    return (np.degrees(np.arctan2(np.sin(r).sum(), np.cos(r).sum())) + 360) % 360


def main():
    rows = json.load(open(BASE / "wn1_vs_steering_validation_500.json"))
    uq = json.load(open(BASE / "wn1_uq_shape_v2.json"))
    uq_map = {(r["storm"], r["iso"]): r for r in uq}

    recs = []
    for r in rows:
        u = uq_map.get((r["storm"], r["iso"]))
        if u is None:
            continue
        recs.append(dict(
            storm=r["storm"], iso=r["iso"],
            wn1_phi=r["wn1_phi"], wn1_amp=r["wn1_amp"],
            steering_dir=r["steering_dir"], steering_mag=r["steering_mag"],
            move_fwd=r["move_fwd"], ellipt=u["ellipt"],
        ))
    n = len(recs)

    # ── A. LOO bias 校正 ──
    # signed bias = WN1 相位 - 移動方向（用 signed_diff）
    signed_biases = np.array([signed_diff(r["wn1_phi"], r["move_fwd"]) for r in recs])

    # LOO：每個樣本用其他 n-1 個樣本嘅 circular mean signed bias 校正
    err_wn1_raw = []
    err_wn1_loo = []
    err_wn1_insample = []
    for i, r in enumerate(recs):
        err_raw = ang_diff(r["wn1_phi"], r["move_fwd"])
        err_wn1_raw.append(err_raw)

        # LOO bias（排除自己）
        others = np.delete(signed_biases, i)
        bias_loo = circular_mean(others)
        # 校正：WN1 相位 - bias（如果 WN1 系統性領先 move +bias，減 bias 拉返近）
        corr_loo = (r["wn1_phi"] - bias_loo) % 360
        err_wn1_loo.append(ang_diff(corr_loo, r["move_fwd"]))

        # in-sample bias（全部樣本，有 look-ahead）
        bias_in = circular_mean(signed_biases)
        corr_in = (r["wn1_phi"] - bias_in) % 360
        err_wn1_insample.append(ang_diff(corr_in, r["move_fwd"]))

    print("=" * 70)
    print("🔬 A. WN1 bias 校正（signed bias = WN1 - Move）")
    print("=" * 70)
    print(f"全樣本 signed bias: {circular_mean(signed_biases):.1f}° (n={n})")
    print(f"bias 嘅散佈 (circ std): {np.std(np.radians(signed_biases))*180/np.pi:.1f}°")
    print(f"\n{'方法':<28}{'Mean':>9}{'Median':>9}{'<45°':>8}")
    print("-" * 60)
    for name, errs in [
        ("純 WN1（無校正）", err_wn1_raw),
        ("WN1 + LOO 校正", err_wn1_loo),
        ("WN1 + in-sample 校正", err_wn1_insample),
        ("純 Steering（對照）", [ang_diff(r["steering_dir"], r["move_fwd"]) for r in recs]),
    ]:
        pct = sum(1 for e in errs if e < 45) / n * 100
        print(f"{name:<28}{np.mean(errs):>9.1f}{np.median(errs):>9.1f}{pct:>7.1f}%")

    # 每樣本 LOO 校正改善幾多
    impro = np.array(err_wn1_raw) - np.array(err_wn1_loo)
    print(f"\nLOO 校正 per-sample 改善: mean={impro.mean():+.1f}° (正=改善, 負=變差)")
    print(f"改善嘅樣本: {(impro > 0).sum()}/{n}, 變差: {(impro < 0).sum()}/{n}")

    # ── C. 「Steering 贏」特徵分析 ──
    print("\n" + "=" * 70)
    print("🔍 C. 「Steering 贏」樣本特徵（n=22 vs WN1 贏 n=52）")
    print("=" * 70)
    err_s = np.array([ang_diff(r["steering_dir"], r["move_fwd"]) for r in recs])
    err_w = np.array([ang_diff(r["wn1_phi"], r["move_fwd"]) for r in recs])
    steer_wins = err_s < err_w
    wn1_wins = err_w < err_s

    print(f"\nSteering 贏: {steer_wins.sum()}  WN1 贏: {wn1_wins.sum()}")

    # storm 分佈
    print("\n① Storm 分佈:")
    storms = sorted(set(r["storm"] for r in recs))
    for s in storms:
        idxs = [i for i, r in enumerate(recs) if r["storm"] == s]
        sw = sum(1 for i in idxs if steer_wins[i])
        print(f"   {s}: Steering 贏 {sw}/{len(idxs)} ({sw/len(idxs)*100:.0f}%)")

    # 特徵對比（Mann-Whitney 用 scipy？簡單比較 median/mean）
    print("\n② 特徵對比（Steering 贏 vs WN1 贏）:")
    feats = [
        ("steering_mag", [r["steering_mag"] for r in recs]),
        ("wn1_amp", [r["wn1_amp"] for r in recs]),
        ("ellipt", [r["ellipt"] for r in recs]),
        ("wind", [r.get("wind", np.nan) for r in recs]),
        ("WN1-Steer 角度差", [ang_diff(r["wn1_phi"], r["steering_dir"]) for r in recs]),
    ]
    print(f"   {'特徵':<20}{'Steer贏 mean':>14}{'WN1贏 mean':>14}{'Steer贏 med':>14}{'WN1贏 med':>14}")
    for fname, vals in feats:
        v = np.array(vals)
        print(f"   {fname:<20}{v[steer_wins].mean():>14.2f}{v[wn1_wins].mean():>14.2f}"
              f"{np.median(v[steer_wins]):>14.2f}{np.median(v[wn1_wins]):>14.2f}")

    # 用 scipy 做 Mann-Whitney U（非參數，適合 n 細）
    try:
        from scipy.stats import mannwhitneyu
        print("\n   Mann-Whitney U (p 值):")
        for fname, vals in feats:
            v = np.array(vals)
            u, p = mannwhitneyu(v[steer_wins], v[wn1_wins], alternative="two-sided")
            print(f"   {fname:<20}U={u:>7.1f}  p={p:.4f}  {'⭐' if p<0.05 else ''}")
    except ImportError:
        print("\n   (scipy 唔可用，跳過 Mann-Whitney)")

    # ③ WN1-Steer 角度差分層：兩者差得遠時邊個贏？
    print("\n③ WN1 vs Steering 角度差分層（兩者差得遠 = 邊個啱？）:")
    ws_diff = np.array([ang_diff(r["wn1_phi"], r["steering_dir"]) for r in recs])
    for lo, hi, label in [(0, 30, "<30° 一致"), (30, 60, "30-60°"), (60, 90, "60-90°"), (90, 181, ">90° 相反")]:
        mask = (ws_diff >= lo) & (ws_diff < hi)
        if mask.sum() == 0:
            continue
        sw = steer_wins[mask].sum()
        ww = wn1_wins[mask].sum()
        print(f"   {label:<12} n={mask.sum():>3}  Steering 贏 {sw:>2} ({sw/mask.sum()*100:.0f}%)  "
              f"WN1 贏 {ww:>2} ({ww/mask.sum()*100:.0f}%)  "
              f"Steer err={err_s[mask].mean():.1f}°  WN1 err={err_w[mask].mean():.1f}°")

    # ④ Steering 贏嗰啲樣本嘅共同點（逐個列）
    print("\n④ Steering 贏嘅 22 個樣本（睇有冇 pattern）:")
    print(f"   {'storm':<10}{'wn1_amp':>9}{'steer_mag':>10}{'ellipt':>8}{'wind':>6}"
          f"{'WN1-Steer':>11}{'Steer err':>11}{'WN1 err':>9}")
    for i, r in enumerate(recs):
        if steer_wins[i]:
            print(f"   {r['storm']:<10}{r['wn1_amp']:>9.1f}{r['steering_mag']:>10.2f}"
                  f"{r['ellipt']:>8.2f}{r.get('wind',0):>6.0f}"
                  f"{ws_diff[i]:>11.1f}{err_s[i]:>11.1f}{err_w[i]:>9.1f}")

    # 儲存
    out = BASE / "wn1_steering_bias_loo.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(dict(
            n=n,
            signed_bias_full=float(circular_mean(signed_biases)),
            bias_circ_std=float(np.std(np.radians(signed_biases))*180/np.pi),
            methods=dict(
                raw=dict(mean=float(np.mean(err_wn1_raw)), median=float(np.median(err_wn1_raw)),
                         pct_lt45=float(sum(1 for e in err_wn1_raw if e < 45)/n*100)),
                loo=dict(mean=float(np.mean(err_wn1_loo)), median=float(np.median(err_wn1_loo)),
                         pct_lt45=float(sum(1 for e in err_wn1_loo if e < 45)/n*100)),
                insample=dict(mean=float(np.mean(err_wn1_insample)), median=float(np.median(err_wn1_insample)),
                              pct_lt45=float(sum(1 for e in err_wn1_insample if e < 45)/n*100)),
                steering=dict(mean=float(np.mean(err_s)), median=float(np.median(err_s)),
                              pct_lt45=float(sum(1 for e in err_s if e < 45)/n*100)),
            ),
            loo_improve_mean=float(impro.mean()),
            loo_improve_n=dict(better=int((impro > 0).sum()), worse=int((impro < 0).sum())),
            per_sample=[dict(storm=recs[i]["storm"], iso=recs[i]["iso"],
                             err_steer=float(err_s[i]), err_wn1=float(err_w[i]),
                             wn1_amp=recs[i]["wn1_amp"], steering_mag=recs[i]["steering_mag"],
                             ellipt=recs[i]["ellipt"], wind=recs[i].get("wind"),
                             wn1_steer_diff=float(ws_diff[i]),
                             steer_wins=bool(steer_wins[i]))
                        for i in range(n)],
        ), f, indent=2, ensure_ascii=False)
    print(f"\n💾 saved: {out}")


if __name__ == "__main__":
    main()
