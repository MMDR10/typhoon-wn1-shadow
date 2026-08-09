#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌀 500 hPa WN1 對比（真正引導流層）
=====================================
用 era5_500/ 嘅獨立 500 hPa 檔案，計算 WN1 相位 + ellipt，
同 200 hPa 對比前瞻誤差。
"""
import xarray as xr
import numpy as np
import json
from datetime import datetime, timedelta
from pathlib import Path

TRACKS = Path("era5_downloader/tracks_6typhoon.json")
DATA_500 = Path("era5_500")
DATA_200 = Path("era5_6typhoon")
ROWS = Path("wn1_forward_validation.json")
OUT = Path("wn1_500_vs_200.json")


def shape_at(f, clat, clon, level, rmin=1.0, rmax=5.0):
    """指定層數嘅徑向風傅立葉分解：WN1 相位 + ellipt"""
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
    wn0_r = np.mean(u_r)
    A1 = 2.0*np.mean(u_r*np.cos(th)); B1 = 2.0*np.mean(u_r*np.sin(th))
    wn1_r = np.hypot(A1, B1); wn1_phi = (np.degrees(np.arctan2(B1, A1))+360) % 360
    A2 = 2.0*np.mean(u_r*np.cos(2*th)); B2 = 2.0*np.mean(u_r*np.sin(2*th))
    wn2_r = np.hypot(A2, B2)
    ellipt = wn2_r / max(wn1_r, 1e-9)
    ds.close()
    return dict(wn1_phi=wn1_phi, ellipt=ellipt, wn1_amp=wn1_r, asym=wn1_r/max(abs(wn0_r), 1e-9))


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

    results = []
    for r in rows:
        storm, iso = r["storm"], r["iso"]
        fname = f"{storm}_{iso[:10]}{iso[11:13]}.nc".replace("-", "")
        f500 = DATA_500 / storm / fname
        if not f500.exists():
            continue
        c = find_point_at(tracks[storm]["points"], iso)
        if c is None:
            continue
        t0 = datetime.strptime(iso, "%Y-%m-%d %H:%M")
        mv = move_direction(tracks[storm]["points"], t0, 12)
        if mv is None:
            continue
        sh500 = shape_at(f500, c["lat"], c["lon"], 500)
        # 200 hPa 由 850/200 檔案
        f200 = DATA_200 / storm / fname
        if not f200.exists():
            continue
        sh200 = shape_at(f200, c["lat"], c["lon"], 200)
        results.append(dict(
            storm=storm, iso=iso, move=mv,
            phi500=sh500["wn1_phi"], phi200=sh200["wn1_phi"],
            ellipt500=sh500["ellipt"], ellipt200=sh200["ellipt"],
            amp500=sh500["wn1_amp"], amp200=sh200["wn1_amp"],
            err500=ang_diff(sh500["wn1_phi"], mv),
            err200=ang_diff(sh200["wn1_phi"], mv),
            wind=r["wind"],
        ))

    def stat(rows, label, key='err'):
        if not rows:
            print(f'  {label}: n=0')
            return
        e = [x[f'{key}'] for x in rows]
        print(f'  {label}: n={len(e):>3}  mean={np.mean(e):5.1f}  median={np.median(e):5.1f}  <45={np.mean([x<45 for x in e])*100:.0f}%')

    print('=== 500 hPa vs 200 hPa WN1 對比（前瞻 t→t+12h）===\n')
    stat(results, '500 hPa 全部', 'err500')
    stat(results, '200 hPa 全部', 'err200')
    print()

    good500 = [x for x in results if x['ellipt500'] <= 0.4]
    good200 = [x for x in results if x['ellipt200'] <= 0.4]
    stat(good500, '500 ellipt<=0.4', 'err500')
    stat(good200, '200 ellipt<=0.4', 'err200')
    print()

    from scipy import stats as st
    dphi = np.array([ang_diff(x['phi500'], x['phi200']) for x in results])
    print(f'兩層相位差: mean={np.mean(dphi):.1f}°  median={np.median(dphi):.1f}°')
    print(f'  <30°: {np.mean(dphi<30)*100:.0f}%  <60°: {np.mean(dphi<60)*100:.0f}%  >90°: {np.mean(dphi>90)*100:.0f}%')
    print()

    better500 = sum(1 for x in results if x['err500'] < x['err200'])
    better200 = sum(1 for x in results if x['err200'] < x['err500'])
    print(f'500 較好: {better500}/{len(results)} ({better500/len(results)*100:.0f}%)')
    print(f'200 較好: {better200}/{len(results)} ({better200/len(results)*100:.0f}%)')
    print()

    oracle = [min(x['err500'], x['err200']) for x in results]
    print(f'Oracle: mean={np.mean(oracle):.1f}  <45={np.mean([x<45 for x in oracle])*100:.0f}%')
    mixed = []
    for x in results:
        if x['ellipt500'] <= x['ellipt200']:
            mixed.append(x['err500'])
        else:
            mixed.append(x['err200'])
    print(f'混合（揀 ellipt 較細嗰層）: mean={np.mean(mixed):.1f}  <45={np.mean([x<45 for x in mixed])*100:.0f}%')
    print()

    print('--- 每颱風 500 vs 200 ---')
    for s in sorted({x['storm'] for x in results}):
        sub = [x for x in results if x['storm'] == s]
        e500 = [x['err500'] for x in sub]
        e200 = [x['err200'] for x in sub]
        print(f'  {s}: n={len(sub):>2}  500={np.mean(e500):5.1f}°  200={np.mean(e200):5.1f}°')

    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f'\nsaved: {OUT}')


if __name__ == "__main__":
    main()
