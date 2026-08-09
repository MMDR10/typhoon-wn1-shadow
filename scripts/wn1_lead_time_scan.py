#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌀 Lead-Time 掃描：φ(t) 對 t+6h/+12h/+18h/+24h 移動嘅預測能力退化曲線
======================================================================
問題：「早 12 小時」到底有幾可靠？早 6h 定早 24h 會點？
方法：對每個樣本 t（76 個），用 φ(t)（200hPa WN1 相位）對比 track 上
      t → t+Lh 嘅移動方位角，L = 6, 12, 18, 24。
      同時計隨機 baseline + 同期對照（t-6→t）。
"""
import json
import numpy as np
from datetime import datetime, timedelta
from pathlib import Path

TRACKS = Path("era5_downloader/tracks_6typhoon.json")
ROWS = Path("wn1_forward_validation.json")
OUT = Path("wn1_lead_time.json")


def load_tracks():
    with open(TRACKS) as f:
        return json.load(f)


def bearing(lat1, lon1, lat2, lon2):
    """兩點方位角（度，0-360）"""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(dlon)
    return (np.degrees(np.arctan2(x, y)) + 360) % 360


def ang_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def find_point_at(track_points, target_dt):
    best, best_delta = None, None
    for p in track_points:
        if not p["iso"]:
            continue
        try:
            t = datetime.strptime(p["iso"][:16], "%Y-%m-%d %H:%M")
        except ValueError:
            continue
        delta = abs((t - target_dt).total_seconds())
        if best_delta is None or delta < best_delta:
            best_delta, best = delta, p
    return best


def move_direction(points, t0, lead_h):
    """track 上 t0 → t0+lead_h 嘅移動方位角（用 track 3h 內插）"""
    p0 = find_point_at(points, t0)
    if p0 is None:
        return None
    t1 = t0 + timedelta(hours=lead_h)
    p1 = find_point_at(points, t1)
    if p1 is None:
        return None
    return bearing(float(p0["lat"]), float(p0["lon"]), float(p1["lat"]), float(p1["lon"]))


def main():
    tracks = load_tracks()
    rows = json.load(open(ROWS))

    leads = [6, 12, 18, 24]
    results = []

    # 隨機 baseline（每個樣本固定 random seed 對比）
    rng = np.random.default_rng(42)
    random_baselines = {L: rng.uniform(0, 180, len(rows)) for L in leads}

    for i, r in enumerate(rows):
        storm, iso = r["storm"], r["iso"]
        t0 = datetime.strptime(iso, "%Y-%m-%d %H:%M")
        wn1 = r["wn1"]
        pts = tracks[storm]["points"]

        row = {"storm": storm, "iso": iso, "wn1": wn1,
               "wind": r["wind"], "amp": r["wn1_amp"]}
        # 同期（t-6 → t）：t-6 相位 vs t-6→t 移動 = 診斷
        move_6back = move_direction(pts, t0 - timedelta(hours=6), 6)
        row["diag_6h"] = None if move_6back is None else round(ang_diff(wn1, move_6back), 1)

        for L in leads:
            mv = move_direction(pts, t0, L)
            row[f"move_{L}h"] = None if mv is None else round(mv, 1)
            row[f"err_{L}h"] = None if mv is None else round(ang_diff(wn1, mv), 1)
            row[f"rand_{L}h"] = round(float(random_baselines[L][i]), 1)
        results.append(row)

    # ── 統計輸出 ──
    print("=== Lead-Time 掃描（φ(t) vs 移動 t→t+Lh）===")
    print(f"{'lead':>6} {'n':>4} {'mean':>7} {'median':>7} {'<30':>5} {'<45':>5} {'<60':>5}  | 隨機 mean")
    for L in leads:
        errs = [x[f"err_{L}h"] for x in results if x[f"err_{L}h"] is not None]
        rands = [x[f"rand_{L}h"] for x in results if x[f"err_{L}h"] is not None]
        if errs:
            print(f"  {L:>4}h {len(errs):>4} {np.mean(errs):>7.1f} {np.median(errs):>7.1f} "
                  f"{np.mean([e<30 for e in errs])*100:>5.0f} {np.mean([e<45 for e in errs])*100:>5.0f} "
                  f"{np.mean([e<60 for e in errs])*100:>5.0f}  | {np.mean(rands):.1f}")
    print(f"\n  同期（診斷）t-6→t: ", end="")
    diag = [x["diag_6h"] for x in results if x["diag_6h"] is not None]
    print(f"n={len(diag)} mean={np.mean(diag):.1f} <45={np.mean([e<45 for e in diag])*100:.0f}%")

    # 效能提升（vs 隨機）
    print("\n=== 效能提升 vs 隨機 baseline ===")
    for L in leads:
        errs = [x[f"err_{L}h"] for x in results if x[f"err_{L}h"] is not None]
        rands = [x[f"rand_{L}h"] for x in results if x[f"err_{L}h"] is not None]
        if errs:
            impr = (np.mean(rands) - np.mean(errs)) / np.mean(rands) * 100
            print(f"  lead {L:>2}h: 前瞻 {np.mean(errs):5.1f}° vs 隨機 {np.mean(rands):5.1f}° → 提升 {impr:.0f}%")

    # 分颱風（lead 12 同 24）
    print("\n=== 每颱風（lead 12h vs 24h）===")
    for s in sorted({x["storm"] for x in results}):
        sub = [x for x in results if x["storm"] == s]
        e12 = [x["err_12h"] for x in sub if x["err_12h"] is not None]
        e24 = [x["err_24h"] for x in sub if x["err_24h"] is not None]
        if e12 and e24:
            print(f"  {s}: n={len(e12):>2}  L12={np.mean(e12):5.1f}°  L24={np.mean(e24):5.1f}°")

    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nsaved: {OUT}")


if __name__ == "__main__":
    main()
