#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌀 WN1 影子 UQ 門檻驗證（Shape-Transition Uncertainty Quantification）
=======================================================================
承接 shape_transition_analysis_report.md：
WN1 相位係「純偶極形狀」嘅算子 — 形狀偏離偶極（ellipt 高）或形狀轉換中（Δe/Δφ 大）時相位不可靠。

驗證 UQ 門檻效果：
  ellipt ≤ 0.6（形狀近純偶極）
  Δφ ≤ 60°/12h（相位穩定）
  Δe ≤ 0.5/12h（形狀轉換唔劇烈）

預期：剔除壞形狀樣本後，精度由 31.8° 升到 ~24°，<45° 比例 75%→87%
"""
import xarray as xr
import numpy as np
import json
from pathlib import Path

TRACKS = Path("/app/working/workspaces/tygtDc/projects/typhoon-dh-curl/era5_downloader/tracks_6typhoon.json")
DATA_ROOT = Path("/app/working/workspaces/tygtDc/projects/typhoon-dh-curl/era5_6typhoon")
ROWS_JSON = Path("/app/working/workspaces/tygtDc/projects/typhoon-dh-curl/wn1_forward_validation.json")
OUT_JSON = Path("/app/working/workspaces/tygtDc/projects/typhoon-dh-curl/wn1_uq_shape.json")

R_WN1 = (1.0, 5.0)

# ── UQ 門檻 ──
ELLIPT_MAX = 0.6   # 形狀偏離純偶極上限
DPHI_MAX = 60.0    # WN1 相位跳躍上限（°/12h）
DELLIPT_MAX = 0.5  # 形狀轉換率上限（/12h）


def shape_at(f, clat, clon):
    """200 hPa 徑向風環帶傅立葉分解：WN1 相位 + ellipt（WN2/WN1）"""
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
    li200 = np.argmin(np.abs(levs - 200))
    u = ds["u"].values[0, li200]; v = ds["v"].values[0, li200]
    m = (r_deg >= R_WN1[0]) & (r_deg <= R_WN1[1])
    th = np.radians(theta[m])
    uu = u[m]; vv = v[m]
    te = np.radians(90.0 - theta[m])
    u_r = uu * np.cos(te) + vv * np.sin(te)
    A1 = 2.0*np.mean(u_r*np.cos(th)); B1 = 2.0*np.mean(u_r*np.sin(th))
    wn1_r = np.hypot(A1, B1); wn1_phi = (np.degrees(np.arctan2(B1, A1))+360) % 360
    A2 = 2.0*np.mean(u_r*np.cos(2*th)); B2 = 2.0*np.mean(u_r*np.sin(2*th))
    wn2_r = np.hypot(A2, B2)
    ellipt = wn2_r / max(wn1_r, 1e-9)
    ds.close()
    return dict(wn1_phi=wn1_phi, ellipt=ellipt)


def find_point_at(track_points, target_iso):
    from datetime import datetime
    target = datetime.strptime(target_iso, "%Y-%m-%d %H:%M")
    best, best_delta = None, None
    for p in track_points:
        iso = p["iso"]
        if not iso:
            continue
        try:
            t = datetime.strptime(iso[:16], "%Y-%m-%d %H:%M")
        except ValueError:
            continue
        delta = abs((t - target).total_seconds())
        if best_delta is None or delta < best_delta:
            best_delta, best = delta, p
    return best


def ang_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def main():
    with open(TRACKS) as f:
        tracks = json.load(f)
    with open(ROWS_JSON) as f:
        rows = json.load(f)

    results = []
    for r in rows:
        storm, iso = r["storm"], r["iso"]
        fname = f"{storm}_{iso[:10].replace('-','')}{iso[11:13]}.nc"
        fpath = DATA_ROOT / storm / fname
        if not fpath.exists():
            continue
        c = find_point_at(tracks[storm]["points"], iso)
        if c is None:
            continue
        sh = shape_at(fpath, c["lat"], c["lon"])

        # 前 12h 形狀（計轉換率）
        from datetime import datetime, timedelta
        t0 = datetime.strptime(iso, "%Y-%m-%d %H:%M")
        t_m = t0 - timedelta(hours=12)
        f_m = DATA_ROOT / storm / f"{storm}_{t_m.strftime('%Y%m%d%H')}.nc"
        dphi, dell = None, None
        c_m = find_point_at(tracks[storm]["points"], t_m.strftime("%Y-%m-%d %H:%M"))
        if f_m.exists() and c_m is not None:
            sh_m = shape_at(f_m, c_m["lat"], c_m["lon"])
            dphi = ang_diff(sh["wn1_phi"], sh_m["wn1_phi"])
            dell = abs(sh["ellipt"] - sh_m["ellipt"])

        err = ang_diff(r["wn1"], r["move_fwd"])
        # UQ 判定
        flags = []
        if sh["ellipt"] > ELLIPT_MAX:
            flags.append("ellipt")
        if dphi is not None and dphi > DPHI_MAX:
            flags.append("dphi")
        if dell is not None and dell > DELLIPT_MAX:
            flags.append("dell")
        results.append(dict(storm=storm, iso=iso, err=err, ellipt=round(sh["ellipt"], 3),
                            dphi=None if dphi is None else round(dphi, 1),
                            dell=None if dell is None else round(dell, 3),
                            flags=flags, trustworthy=len(flags) == 0))

    # 統計
    all_e = [x["err"] for x in results]
    good = [x for x in results if x["trustworthy"]]
    bad = [x for x in results if not x["trustworthy"]]
    good_e = [x["err"] for x in good]
    bad_e = [x["err"] for x in bad]

    def stat(e, label):
        print(f"  {label}: n={len(e):>3}  mean={np.mean(e):5.1f}°  median={np.median(e):5.1f}°  <45°={np.mean([x<45 for x in e])*100:.0f}%")

    print("═══ WN1 UQ 門檻驗證 ═══")
    print(f"門檻: ellipt≤{ELLIPT_MAX}, Δφ≤{DPHI_MAX}°, Δe≤{DELLIPT_MAX}")
    stat(all_e, "全部")
    stat(good_e, "✅ 可信（pass 全部）")
    stat(bad_e, "❌ 剔除（fail ≥1）")
    print(f"\n剔除比例: {len(bad)}/{len(results)} = {len(bad)/len(results)*100:.0f}%")

    # 每颱風
    print("\n每颱風（可信樣本）:")
    for s in sorted({x["storm"] for x in results}):
        sub = [x for x in good if x["storm"] == s]
        if sub:
            e = [x["err"] for x in sub]
            print(f"  {s}: {len(sub):>2}  mean={np.mean(e):5.1f}°  <45°={np.mean([x<45 for x in e])*100:.0f}%")

    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n💾 saved: {OUT_JSON}")


if __name__ == "__main__":
    main()
