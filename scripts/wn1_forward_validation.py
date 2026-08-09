#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌀 WN1 影子前瞻驗證 Pipeline（反擊 GEMINI 第 1、2 點）
=====================================================
需要：MKP 電腦已下載嘅完整 ERA5 序列（每 12h 一檔，含颱風形成前/消散後背景場）
      era5_6typhoon/{STORM}/{STORM}_{YYYYMMDDHH}.nc

設計（用「t 時刻相位 → t+12h 移動」做真正前瞻）：
1. 對每個時間步 t：200 hPa r1-5° 徑向風 WN1 相位 φ(t)
2. 移動方向：track 上 t → t+12h 嘅方位角 move(t→t+12)
3. 診斷對照：φ(t) vs move(t−12→t)（同期）
4. 前瞻測試：φ(t) vs move(t→t+12)（提前 12h）
5. 逆風局分層：移動方向變化率 |Δmove| 大嘅時間步（轉向/滯留）單獨統計
6. VWS 分層：高 VWS vs 低 VWS 時間步分開統計
"""
import xarray as xr
import numpy as np
import json
from pathlib import Path

TRACKS = Path("/app/working/workspaces/tygtDc/projects/typhoon-dh-curl/era5_downloader/tracks_6typhoon.json")
# 由 MKP 電腦 copy 過嚟嘅數據根目錄
DATA_ROOT = Path("/app/working/workspaces/tygtDc/projects/typhoon-dh-curl/era5_6typhoon")

STORMS = {
    "HATO":     dict(peak="2017-08-23 00:00", hours_min=-216, hours_max=143),
    "MANGKHUT": dict(peak="2018-09-11 12:00", hours_min=-150, hours_max=185),
    "SAOLA":    dict(peak="2023-08-30 00:00", hours_min=-198, hours_max=161),
    "MERANTI":  dict(peak="2016-09-13 12:00", hours_min=-162, hours_max=173),
    "HAGIBIS":  dict(peak="2019-10-07 12:00", hours_min=-816, hours_max=647),
    "GONI":     dict(peak="2020-10-31 18:00", hours_min=-666, hours_max=797),
}

R_WN1 = (1.0, 5.0)
R_ENV = (8.0, 12.0)

# ── 樣本質量門檻 ──
# 背景場（颱風未成形/已消散）WN1 相位係 noise，必須用強度 filter
WIND_MIN = 40.0   # kt（颱風級別 TS 下限 ~34kt，取 40kt 確保結構成熟）
AMP_MIN = 1.0     # m/s（WN1 振幅門檻，太低 = 相位無意義）


def load_tracks():
    with open(TRACKS) as f:
        return json.load(f)


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


def bearing_deg(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    return (np.degrees(np.arctan2(x, y)) + 360) % 360


def ang_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def signed_offset(phi, ref):
    d = (phi - ref) % 360
    return d if d <= 180 else d - 360


def wn1_and_vws(ds, clat, clon):
    """單一檔案：200 hPa 徑向風 WN1 相位 + 環境 VWS 方向"""
    lat = ds["latitude"].values
    lon = ds["longitude"].values
    LON, LAT = np.meshgrid(lon, lat)
    dlat = (LAT - clat) * np.pi / 180.0
    dlon = (LON - clon) * np.pi / 180.0
    r_deg = np.degrees(np.arccos(np.clip(
        np.sin(clat * np.pi / 180) * np.sin(LAT * np.pi / 180) +
        np.cos(clat * np.pi / 180) * np.cos(LAT * np.pi / 180) * np.cos(dlon), -1, 1)))
    theta = np.degrees(np.arctan2(
        np.sin(dlon) * np.cos(LAT * np.pi / 180),
        np.cos(clat * np.pi / 180) * np.sin(LAT * np.pi / 180) -
        np.sin(clat * np.pi / 180) * np.cos(LAT * np.pi / 180) * np.cos(dlon)))
    theta = (theta + 360) % 360

    levs = ds["pressure_level"].values
    li200 = np.argmin(np.abs(levs - 200))
    li850 = np.argmin(np.abs(levs - 850))
    u = ds["u"].values[0]; v = ds["v"].values[0]

    # WN1（200 hPa 徑向風）
    m = (r_deg >= R_WN1[0]) & (r_deg <= R_WN1[1])
    th = np.radians(theta[m])
    uu = u[li200][m]; vv = v[li200][m]
    te = np.radians(90.0 - theta[m])
    u_r = uu * np.cos(te) + vv * np.sin(te)
    A1 = 2.0 * np.mean(u_r * np.cos(th))
    B1 = 2.0 * np.mean(u_r * np.sin(th))
    wn1 = (np.degrees(np.arctan2(B1, A1)) + 360) % 360
    wn1_amp = float(np.hypot(A1, B1))

    # VWS（環境環帶 850→200）
    me = (r_deg >= R_ENV[0]) & (r_deg <= R_ENV[1])
    du = u[li200][me].mean() - u[li850][me].mean()
    dv = v[li200][me].mean() - v[li850][me].mean()
    vws = (np.degrees(np.arctan2(du, dv)) + 360) % 360
    vws_mag = float(np.hypot(du, dv))
    return wn1, wn1_amp, vws, vws_mag


def main():
    tracks = load_tracks()
    all_rows = []

    for storm, cfg in STORMS.items():
        storm_dir = DATA_ROOT / storm
        if not storm_dir.exists():
            print(f"⏭️ {storm}: 未見 {storm_dir}")
            continue
        files = sorted(storm_dir.glob(f"{storm}_*.nc"))
        if not files:
            print(f"⏭️ {storm}: 冇檔案")
            continue
        print(f"🌀 {storm}: {len(files)} 檔")

        track_points = tracks[storm]["points"]
        for f in files:
            # 由檔名解析時間：{STORM}_{YYYYMMDDHH}.nc
            iso = f.stem.replace(f"{storm}_", "")
            tstr = f"{iso[:4]}-{iso[4:6]}-{iso[6:8]} {iso[8:10]}:00"
            center = find_point_at(track_points, tstr)
            if center is None:
                continue
            clat, clon = center["lat"], center["lon"]

            # 強度門檻：背景場（未成形/已消散）WN1 相位係 noise
            wind_raw = center.get("wind", "")
            try:
                wind = float(wind_raw)
            except (ValueError, TypeError):
                wind = None
            if wind is None or wind < WIND_MIN:
                continue

            # 移動方向：t → t+12h（前瞻目標）
            from datetime import datetime, timedelta
            t0 = datetime.strptime(tstr, "%Y-%m-%d %H:%M")
            p_fwd = find_point_at(track_points, (t0 + timedelta(hours=12)).strftime("%Y-%m-%d %H:%M"))
            p_bwd = find_point_at(track_points, (t0 - timedelta(hours=12)).strftime("%Y-%m-%d %H:%M"))
            if p_fwd is None or p_bwd is None:
                continue
            move_fwd = bearing_deg(clat, clon, p_fwd["lat"], p_fwd["lon"])  # 前瞻
            move_bwd = bearing_deg(p_bwd["lat"], p_bwd["lon"], clat, clon)  # 同期（過去→現在）

            try:
                ds = xr.open_dataset(f)
                wn1, wn1_amp, vws, vws_mag = wn1_and_vws(ds, clat, clon)
                ds.close()
            except Exception as e:
                print(f"  ⚠️ {f.name}: {e}")
                continue
            if wn1_amp < AMP_MIN:
                continue

            # 移動方向變化率（轉向度數/12h）
            turn = ang_diff(move_bwd, move_fwd)
            all_rows.append(dict(storm=storm, iso=tstr, wn1=wn1, wn1_amp=wn1_amp,
                                 vws=vws, vws_mag=vws_mag, wind=wind,
                                 move_fwd=move_fwd, move_bwd=move_bwd, turn=turn))

    if not all_rows:
        print("❌ 冇數據。請將 MKP 電腦嘅 era5_6typhoon/ 目錄 copy 到", DATA_ROOT)
        return

    n = len(all_rows)
    print(f"\n總樣本: {n}")

    # ── 統計 ──
    def circ_stats(rows, label):
        if not rows:
            print(f"  {label}: 冇樣本")
            return
        arr = np.array(rows)
        m = arr.mean(); md = np.median(arr); lo = sum(arr < 45) / len(arr)
        print(f"  {label}: n={len(arr):>4}  mean={m:5.1f}°  median={md:5.1f}°  <45°={lo*100:.0f}%")

    print("\n═══ 1. 同期 vs 前瞻 ═══")
    diag = [ang_diff(r["wn1"], r["move_bwd"]) for r in all_rows]   # φ(t) vs 移動(t-12→t)
    fore = [ang_diff(r["wn1"], r["move_fwd"]) for r in all_rows]   # φ(t) vs 移動(t→t+12)
    circ_stats(diag, "診斷 φ(t) vs 過去→現在移動")
    circ_stats(fore, "前瞻 φ(t) vs 未來 12h 移動  ← 真正預測價值")

    print("\n═══ 2. 逆風局分層（轉向度/12h）═══")
    straight = [ang_diff(r["wn1"], r["move_fwd"]) for r in all_rows if r["turn"] <= 20]
    turning = [ang_diff(r["wn1"], r["move_fwd"]) for r in all_rows if r["turn"] > 20]
    sharp = [ang_diff(r["wn1"], r["move_fwd"]) for r in all_rows if r["turn"] > 40]
    circ_stats(straight, "直行段 (Δmove≤20°)")
    circ_stats(turning, "轉向段 (Δmove>20°)")
    circ_stats(sharp, "急轉段 (Δmove>40°)  ← 真正考驗")

    print("\n═══ 3. VWS 分層 ═══")
    low_vws = [ang_diff(r["wn1"], r["move_fwd"]) for r in all_rows if r["vws_mag"] < 10]
    high_vws = [ang_diff(r["wn1"], r["move_fwd"]) for r in all_rows if r["vws_mag"] >= 10]
    circ_stats(low_vws, "低 VWS (<10 m/s)")
    circ_stats(high_vws, "高 VWS (≥10 m/s)  ← 風切污染風險最高")

    print("\n═══ 4. 隨機基準 ═══")
    rng = np.random.default_rng(42)
    nulls = []
    for _ in range(2000):
        ph = rng.uniform(0, 360, len(fore)); mv = rng.uniform(0, 360, len(fore))
        nulls.append(np.mean([min(abs(p - m) % 360, 360 - abs(p - m) % 360) for p, m in zip(ph, mv)]))
    nulls = np.array(nulls)
    print(f"  隨機期望 mean={nulls.mean():.1f}°  p1={np.percentile(nulls, 1):.1f}°")

    out = Path("/app/working/workspaces/tygtDc/projects/typhoon-dh-curl/wn1_forward_validation.json")
    # float32 → float 才可 JSON 序列化
    def _conv(o):
        if isinstance(o, dict):
            return {k: _conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_conv(v) for v in o]
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        return o
    with open(out, "w", encoding="utf-8") as f:
        json.dump(_conv(all_rows), f, indent=2, ensure_ascii=False)
    print(f"\n💾 saved: {out}")


if __name__ == "__main__":
    main()
