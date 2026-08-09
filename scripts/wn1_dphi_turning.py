#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌀 dφ/dt 轉向點偵測分析
========================
假說：轉向點（moving direction 快速變化）嘅 WN1 相位都會快速轉動
      （dφ/dt 大）→ 可以作為新 UQ 指標。

用 6 颱風 76 樣本（500 hPa WN1），計算：
1. 每樣本嘅 WN1 相位 φ(t)（已經有）
2. 相鄰樣本相位差 dφ = |φ(t+12h) - φ(t)|
3. 同時期移動方向變化 dmove = |move(t+12h) - move(t)|
4. 測試：dφ 大 係咪對應 dmove 大？（轉向）
5. 測試：dφ 大 時前瞻誤差係咪大？（UQ 有效性）
"""
import xarray as xr
import numpy as np
import json
from datetime import datetime, timedelta
from pathlib import Path

TRACKS = Path("era5_downloader/tracks_6typhoon.json")
DATA_500 = Path("era5_500")
ROWS = Path("wn1_forward_validation.json")
OUT = Path("wn1_dphi_turning.json")


def shape_at(f, clat, clon, level, rmin=1.0, rmax=5.0):
    ds = xr.open_dataset(f)
    lat = ds["latitude"].values; lon = ds["longitude"].values
    LON, LAT = np.meshgrid(lon, lat)
    dlat = (LAT - clat) * np.pi / 180.0
    dlon = (LON - clon) * np.pi / 180.0
    r_deg = np.degrees(np.arccos(np.clip(
        np.sin(clat*np.pi/180)*np.sin(LAT*np.pi/180) +
        np.cos(clat*np.pi/180)*np.cos(LAT*np.pi/180)*np.cos(dlon), -1, 1)))
    theta = np.degrees(np.arctan2(
        np.sin(dlon)*np.cos(LAT*np.pi/180),
        np.cos(clat*np.pi/180)*np.sin(LAT*np.pi/180) -
        np.sin(clat*np.pi/180)*np.cos(LAT*np.pi/180)*np.cos(dlon)))
    theta = (theta + 360) % 360
    levs = ds["pressure_level"].values
    li = np.argmin(np.abs(levs - level))
    u = ds["u"].values[0, li]; v = ds["v"].values[0, li]
    m = (r_deg >= rmin) & (r_deg <= rmax)
    th = np.radians(theta[m])
    uu = u[m]; vv = v[m]
    te = np.radians(90.0 - theta[m])
    u_r = uu * np.cos(te) + vv * np.sin(te)
    A1 = 2.0*np.mean(u_r*np.cos(th)); B1 = 2.0*np.mean(u_r*np.sin(th))
    wn1_phi = (np.degrees(np.arctan2(B1, A1))+360) % 360
    A2 = 2.0*np.mean(u_r*np.cos(2*th)); B2 = 2.0*np.mean(u_r*np.sin(2*th))
    wn2_r = np.hypot(A2, B2)
    wn1_r = np.hypot(A1, B1)
    ellipt = wn2_r / max(wn1_r, 1e-9)
    ds.close()
    return dict(wn1_phi=wn1_phi, ellipt=ellipt, wn1_amp=wn1_r)


def find_point_at(track_points, target_iso):
    target = datetime.strptime(target_iso, "%Y-%m-%d %H:%M")
    best, best_delta = None, None
    for p in track_points:
        if not p["iso"]:
            continue
        try:
            t = datetime.strptime(p["iso"][:16], "%Y-%m-%d %H:%M")
        except ValueError:
            continue
        delta = abs((t - target).total_seconds())
        if best_delta is None or delta < best_delta:
            best_delta, best = delta, p
    return best


def ang_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def bearing(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(dlon)
    return (np.degrees(np.arctan2(x, y)) + 360) % 360


def move_direction(points, t0, lead_h):
    p0 = find_point_at(points, t0.strftime("%Y-%m-%d %H:%M"))
    if p0 is None:
        return None
    p1 = find_point_at(points, (t0 + timedelta(hours=lead_h)).strftime("%Y-%m-%d %H:%M"))
    if p1 is None:
        return None
    return bearing(float(p0["lat"]), float(p0["lon"]), float(p1["lat"]), float(p1["lon"]))


def main():
    with open(TRACKS) as f:
        tracks = json.load(f)
    with open(ROWS) as f:
        rows = json.load(f)

    # 第一步：計算每個樣本嘅 500 hPa WN1 相位 + 移動方向
    samples = []
    for r in rows:
        storm, iso = r["storm"], r["iso"]
        fname = f"{storm}_{iso[:10]}{iso[11:13]}.nc".replace("-", "")
        fpath = DATA_500 / storm / fname
        if not fpath.exists():
            continue
        c = find_point_at(tracks[storm]["points"], iso)
        if c is None:
            continue
        t0 = datetime.strptime(iso, "%Y-%m-%d %H:%M")
        mv = move_direction(tracks[storm]["points"], t0, 12)
        if mv is None:
            continue
        sh = shape_at(fpath, c["lat"], c["lon"], 500)
        samples.append(dict(storm=storm, t=t0, iso=iso,
                            phi=sh["wn1_phi"], ellipt=sh["ellipt"], move=mv))

    # 第二步：按颱風排序，計算相鄰 12h 樣本嘅 dφ 同 dmove
    results = []
    for storm in sorted({s["storm"] for s in samples}):
        sub = sorted([s for s in samples if s["storm"] == storm], key=lambda x: x["t"])
        for i in range(len(sub) - 1):
            a, b = sub[i], sub[i+1]
            dt_h = (b["t"] - a["t"]).total_seconds() / 3600
            if dt_h != 12:  # 只要 12h 間隔
                continue
            dphi = ang_diff(a["phi"], b["phi"])
            dmove = ang_diff(a["move"], b["move"])
            # b 嘅前瞻誤差（b 時刻 WN1 → b+12h 移動）
            t_b = b["t"]
            mv_fwd = move_direction(tracks[storm]["points"], t_b, 12)
            if mv_fwd is None:
                continue
            err_b = ang_diff(b["phi"], mv_fwd)
            results.append(dict(
                storm=storm, iso_b=b["iso"],
                dphi=dphi, dmove=dmove,
                phi_a=a["phi"], phi_b=b["phi"],
                move_a=a["move"], move_b=b["move"],
                ellipt_b=b["ellipt"], err_b=err_b,
            ))

    print(f'=== dφ/dt 轉向點偵測（12h 間隔樣本對, n={len(results)}）===\n')
    print('--- 基礎統計 ---')
    dphis = [r["dphi"] for r in results]
    dmoves = [r["dmove"] for r in results]
    print(f'dφ (12h 相位變化): mean={np.mean(dphis):5.1f}°  median={np.median(dphis):5.1f}°')
    print(f'dmove (12h 移動方向變化): mean={np.mean(dmoves):5.1f}°  median={np.median(dmoves):5.1f}°')
    print()

    print('--- 測試 1: dφ 大 ⇔ dmove 大（轉向）？ ---')
    corr = np.corrcoef(dphis, dmoves)[0, 1]
    print(f'相關係數: {corr:.3f}')
    # 分組
    hi_dphi = [r for r in results if r["dphi"] > 60]
    lo_dphi = [r for r in results if r["dphi"] <= 60]
    if hi_dphi and lo_dphi:
        print(f'dφ>60° 組 (n={len(hi_dphi)}): dmove mean={np.mean([r["dmove"] for r in hi_dphi]):.1f}°')
        print(f'dφ≤60° 組 (n={len(lo_dphi)}): dmove mean={np.mean([r["dmove"] for r in lo_dphi]):.1f}°')
    print()

    print('--- 測試 2: dφ 大 ⇔ 前瞻誤差大（UQ 有效性）？ ---')
    corr2 = np.corrcoef([r["dphi"] for r in results], [r["err_b"] for r in results])[0, 1]
    print(f'相關係數 (dφ vs err): {corr2:.3f}')
    if hi_dphi and lo_dphi:
        print(f'dφ>60° 組: err_b mean={np.mean([r["err_b"] for r in hi_dphi]):.1f}°  <45={np.mean([r["err_b"]<45 for r in hi_dphi])*100:.0f}%')
        print(f'dφ≤60° 組: err_b mean={np.mean([r["err_b"] for r in lo_dphi]):.1f}°  <45={np.mean([r["err_b"]<45 for r in lo_dphi])*100:.0f}%')
    print()

    print('--- 每颱風 dφ 統計 ---')
    for s in sorted({r["storm"] for r in results}):
        sub = [r for r in results if r["storm"] == s]
        print(f'  {s}: n={len(sub):>2}  dφ mean={np.mean([r["dphi"] for r in sub]):5.1f}°  max={np.max([r["dphi"] for r in sub]):5.1f}°')

    print()
    print('--- 最大 dφ 樣本（轉向候選） ---')
    top = sorted(results, key=lambda r: -r["dphi"])[:6]
    for r in top:
        print(f'  {r["storm"]} {r["iso_b"]}: dφ={r["dphi"]:.0f}°  dmove={r["dmove"]:.0f}°  err_b={r["err_b"]:.0f}°')

    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f'\nsaved: {OUT}')


if __name__ == "__main__":
    main()
