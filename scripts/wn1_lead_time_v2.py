#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lead-time 掃描 v2：UQ 收緊 + 轉向段分層 + 移動速度控制"""
import json
import numpy as np
from datetime import datetime, timedelta

TRACKS = "era5_downloader/tracks_6typhoon.json"
ROWS = "wn1_forward_validation.json"

def load_tracks():
    with open(TRACKS) as f:
        return json.load(f)

def bearing(lat1, lon1, lat2, lon2):
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
    p0 = find_point_at(points, t0)
    if p0 is None:
        return None
    p1 = find_point_at(points, t0 + timedelta(hours=lead_h))
    if p1 is None:
        return None
    return bearing(float(p0["lat"]), float(p0["lon"]), float(p1["lat"]), float(p1["lon"]))

def main():
    tracks = load_tracks()
    rows = json.load(open(ROWS))
    uq = json.load(open("wn1_uq_shape_v2.json"))
    # UQ map
    uq_map = {(x["storm"], x["iso"]): x for x in uq}

    leads = [6, 12, 18, 24, 36, 48]
    results = []

    for r in rows:
        storm, iso = r["storm"], r["iso"]
        t0 = datetime.strptime(iso, "%Y-%m-%d %H:%M")
        wn1 = r["wn1"]
        pts = tracks[storm]["points"]
        u = uq_map.get((storm, iso), {})
        ellipt = u.get("ellipt")
        row = {"storm": storm, "iso": iso, "wn1": wn1, "wind": r["wind"],
               "ellipt": ellipt, "turn": r["turn"]}
        # 轉向段定義：|Δmove| t-6→t+6 大 = 轉向
        m_back = move_direction(pts, t0 - timedelta(hours=6), 6)
        m_fwd = move_direction(pts, t0, 6)
        row["turn_rate"] = None if (m_back is None or m_fwd is None) else round(ang_diff(m_back, m_fwd), 1)
        for L in leads:
            mv = move_direction(pts, t0, L)
            row[f"err_{L}h"] = None if mv is None else round(ang_diff(wn1, mv), 1)
        results.append(row)

    def stat(rows, label):
        if not rows:
            print(f'  {label}: n=0')
            return
        e6 = [x["err_6h"] for x in rows if x["err_6h"] is not None]
        e12 = [x["err_12h"] for x in rows if x["err_12h"] is not None]
        e24 = [x["err_24h"] for x in rows if x["err_24h"] is not None]
        e36 = [x["err_36h"] for x in rows if x["err_36h"] is not None]
        e48 = [x["err_48h"] for x in rows if x["err_48h"] is not None]
        print(f'  {label}: n={len(e12):>3} | L6={np.mean(e6):5.1f} L12={np.mean(e12):5.1f} L24={np.mean(e24):5.1f} L36={np.mean(e36):5.1f} L48={np.mean(e48):5.1f}')

    print("=== Lead-time 掃描 v2（含 UQ 收緊 + 轉向分層）===")
    stat(results, '全部 76')

    # UQ 收緊 ellipt<=0.4
    good = [x for x in results if x["ellipt"] is not None and x["ellipt"] <= 0.4]
    bad = [x for x in results if x["ellipt"] is not None and x["ellipt"] > 0.4]
    stat(good, 'ellipt<=0.4 (可信)')
    stat(bad, 'ellipt>0.4 (剔除)')
    print()

    # 轉向分層（turn_rate = 12h 內移動方向變化）
    turning = [x for x in results if x["turn_rate"] is not None and x["turn_rate"] >= 20]
    straight = [x for x in results if x["turn_rate"] is not None and x["turn_rate"] < 20]
    stat(turning, f'轉向段 (|Δmove|>={20}°): n={len(turning)}')
    stat(straight, f'直線段 (<{20}°): n={len(straight)}')
    print()

    # 交叉：可信 + 轉向
    good_turn = [x for x in good if x["turn_rate"] is not None and x["turn_rate"] >= 20]
    good_str = [x for x in good if x["turn_rate"] is not None and x["turn_rate"] < 20]
    stat(good_turn, '可信 + 轉向段')
    stat(good_str, '可信 + 直線段')
    print()

    # 轉向段詳細 lead 退化
    print("=== 轉向段 lead 退化（最重要）===")
    for thr in [15, 25, 35]:
        sub = [x for x in results if x["turn_rate"] is not None and x["turn_rate"] >= thr]
        if len(sub) >= 5:
            stat(sub, f'turn>={thr}°')

    # 保存
    with open("wn1_lead_time_v2.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print("\nsaved: wn1_lead_time_v2.json")

if __name__ == "__main__":
    main()
