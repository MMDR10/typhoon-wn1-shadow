#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌀 Typhoon Monitor — WN1 路徑 + dH_curl 強度 + 鞍點環 一體化追蹤
=================================================================
每 6h 由 GitHub Actions 觸發，對每個活躍颱風同時收集三類數據：

  1. **路徑（WN1 Shadow）** — 500 hPa 徑向風 wavenumber-1 相位 → 移動方向
  2. **強度（dH_curl）** — 850 hPa 渦度 H_shell − H_core（core=5° / shell=8°）
     🔴 負數 → 組織化結構（已發展颱風）；🟢 正數 → 發散結構
  3. **鞍點環（Saddle Ring）** — 850 hPa ζ 場 Morse 分類，眼牆環帶 (1.5°, 3.0°)
     內鞍點點集 box-count 維數 → D_fold（負 = 比隨機更凝聚 = 環狀結構證據）

輸出：
  - typhoon_history.json  — 每颱風 × 每 cycle 一條（冪等），三類數據同列
  - typhoon_latest.md     — 最新 summary

用法：
  python typhoon_monitor.py                    # 自動：最新 GFS + 活躍颱風
  python typhoon_monitor.py --grib /tmp/g.grib2
  python typhoon_monitor.py --storms 'X:1,2'   # 覆寫颱風位置（測試用）
  python typhoon_monitor.py --cycle 20260811,00
"""
import argparse
import json
import re
import sys
import time
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

HISTORY_PATH = "typhoon_history.json"
SUMMARY_PATH = "typhoon_latest.md"
CYCLOCANE_HOME = "https://www.cyclocane.com/"

# 鞍點環參數（沿用 phase15 設定）
EW_IN = 1.5          # 眼牆環帶內半徑（°）
EW_OUT = 3.0         # 眼牆環帶外半徑（°）
WINDOW_DEG = 6.0     # 鞍點分析窗口（中心 ±6°，夠 cover 環帶 + 對照）
N_NULL = 5           # null 次數（live 版減到 5，慳時間）

# dH_curl 參數（沿用 dolphin v8）
CORE_DEG = 5.0
SHELL_DEG = 8.0


# ─── 0. IGRF-13 地磁場（簡化 n=1,2,3，2020 係數） ──────────────
IGRF_COEFS = {
    'g10': -29404.8, 'g11': -1450.9, 'h11': 4652.5,
    'g20': -2499.6, 'g21': 2982.0, 'h21': -2991.6,
    'g22': 1677.0, 'h22': -734.6,
    'g30': 1363.2, 'g31': -2381.2, 'h31': -82.1,
    'g32': 1236.2, 'h32': 241.9, 'g33': 525.7, 'h33': -543.4,
}

def igrf_field(lat, lon, alt_km=0):
    """IGRF-13 簡化版（n=1,2,3）：返回 total field F (nT) + declination D"""
    theta = np.radians(90 - lat)
    phi = np.radians(lon)
    a = 6371.2
    r = a + alt_km
    Br = Btheta = Bphi = 0.0

    # n=1
    Br += 2 * (a/r)**3 * (IGRF_COEFS['g10'] * np.cos(theta) +
                           IGRF_COEFS['g11'] * np.sin(theta) * np.cos(phi) +
                           IGRF_COEFS['h11'] * np.sin(theta) * np.sin(phi))
    Btheta += -(a/r)**3 * (-IGRF_COEFS['g10'] * np.sin(theta) +
                            IGRF_COEFS['g11'] * np.cos(theta) * np.cos(phi) +
                            IGRF_COEFS['h11'] * np.cos(theta) * np.sin(phi))
    Bphi += -(a/r)**3 * (IGRF_COEFS['g11'] * np.sin(phi) -
                          IGRF_COEFS['h11'] * np.cos(phi)) / np.sin(theta)
    # n=2
    P20 = 0.5 * (3 * np.cos(theta)**2 - 1)
    P21 = 3 * np.sin(theta) * np.cos(theta)
    P22 = 3 * np.sin(theta)**2
    Br += 3 * (a/r)**4 * (IGRF_COEFS['g20'] * P20 +
                           (IGRF_COEFS['g21'] * np.cos(phi) + IGRF_COEFS['h21'] * np.sin(phi)) * P21 +
                           (IGRF_COEFS['g22'] * np.cos(2*phi) + IGRF_COEFS['h22'] * np.sin(2*phi)) * P22)
    Btheta += -(a/r)**4 * (IGRF_COEFS['g20'] * (1.5*np.sin(2*theta)) +
                            IGRF_COEFS['g21'] * np.cos(phi) * (3*np.cos(2*theta) - 1) * 0.5 +
                            IGRF_COEFS['h21'] * np.sin(phi) * (3*np.cos(2*theta) - 1) * 0.5 +
                            IGRF_COEFS['g22'] * np.cos(2*phi) * (3*np.sin(2*theta)) +
                            IGRF_COEFS['h22'] * np.sin(2*phi) * (3*np.sin(2*theta)))
    Bphi += -(a/r)**4 * (IGRF_COEFS['g21'] * (-np.sin(phi)) * (3*np.sin(theta)*np.cos(theta)) +
                          IGRF_COEFS['h21'] * np.cos(phi) * (3*np.sin(theta)*np.cos(theta)) +
                          IGRF_COEFS['g22'] * (-2*np.sin(2*phi)) * (3*np.sin(theta)**2) +
                          IGRF_COEFS['h22'] * (2*np.cos(2*phi)) * (3*np.sin(theta)**2))
    # n=3（只 g30，足夠總強度近似）
    P30 = 0.5 * (5 * np.cos(theta)**3 - 3 * np.cos(theta))
    Br += 4 * (a/r)**5 * IGRF_COEFS['g30'] * P30

    Bx = -Btheta; By = -Bphi; Bz = -Br
    F = float(np.sqrt(Bx**2 + By**2 + Bz**2))
    D = float(np.degrees(np.arctan2(By, Bx)))
    return F, D


def classify_uq(ellipt, amp, mag_F, lat):
    """
    WN1 UQ 機制 v3（2026-08-15 confound test 修正版）：
      🟢🟢 高信心：ellipt ≤ 0.4 + amp ≥ 門檻（經 confound test 驗證）
      🟡 中信心：ellipt ≤ 0.4（任何 amp）
      🔴 低信心：ellipt > 0.4（尤其 amp < 7 = Q4 區域，30% 準確度）
    緯度依賴門檻：低緯 <20° 用 Amp≥10、中緯 20-30° 用 Amp≥7、高緯 ≥30° 用 Amp≥5
    ⚠️ v3 變更：移除「Mag≥35000 升級」——confound test 證明地磁場 r(amp,F|lat)
       ≈ 0（200hPa -0.114 p=0.21；500hPa 反號 -0.315），磁場係緯度偽相關，
       「地磁場強度調制」撤回。緯度門檻保留（partial r≈0.45 獨立顯著）。
    """
    amp_thresh = 10.0
    if lat >= 20 and lat < 30:
        amp_thresh = 7.0
    elif lat >= 30:
        amp_thresh = 5.0

    if ellipt <= 0.4:
        if amp >= amp_thresh:
            return "VERY_HIGH", amp_thresh
        return "MEDIUM", amp_thresh
    return "LOW", amp_thresh


# ─── 1. 活躍颱風列表（cyclocane 主頁） ──────────────────────────
def fetch_active_storms():
    """攞 NW Pacific 活躍颱風：name, slug, 風速, 是否 final advisory"""
    req = urllib.request.Request(CYCLOCANE_HOME, headers={"User-Agent": "Mozilla/5.0"})
    html = urllib.request.urlopen(req, timeout=30).read().decode("utf-8", errors="replace")
    m = re.search(r"Northwest Pacific Storms</h3>(.*?)<h2>Tropical Disturbances", html, re.S)
    if not m:
        return []
    seg = m.group(1)
    storms = []
    for li in re.finditer(r"<li>(.*?)</li>", seg, re.S):
        block = li.group(1)
        link = re.search(r"href='/([a-z\-]+)-storm-tracker/'", block)
        if not link:
            continue
        slug = link.group(1)
        name_m = re.search(r">\s*([A-Z][A-Z\- ]+?)\s*</span>", block)
        name = name_m.group(1).strip() if name_m else slug.upper()
        final = "final advisory" in block.lower()
        wind_m = re.search(r"Current Wind:\s*([\d]+)\s*knots", block)
        wind = int(wind_m.group(1)) if wind_m else None
        storms.append({
            "name": name.replace("TYPHOON ", "").replace("TROPICAL STORM ", "")
                       .replace("TROPICAL DEPRESSION ", ""),
            "slug": slug,
            "wind_kts": wind,
            "final": final,
        })
    return [s for s in storms if not s["final"]]


# ─── 2. 颱風最新位置（cyclocane tracker 頁內嵌 JTWC advisory） ──
def fetch_storm_position(slug):
    """攞 position: (lat, lon), storm_id, advisory_time, wind_kts"""
    url = f"https://www.cyclocane.com/{slug}-storm-tracker"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    html = urllib.request.urlopen(req, timeout=30).read().decode("utf-8", errors="replace")
    adv = re.search(r"WTPN\d+\s+PGTW\s+(\d{6})", html)
    adv_time = adv.group(1) if adv else None
    sid = re.search(r"(?:TROPICAL STORM|TYPHOON|SUPER TYPHOON|TROPICAL DEPRESSION)\s+(\d+W)", html)
    storm_id = sid.group(1) if sid else None
    pos = None
    wp = re.search(r"WARNING POSITION:.*?NEAR\s+([\d.]+)([NS])\s+([\d.]+)([EW])", html, re.S)
    if wp:
        lat = float(wp.group(1)) * (1 if wp.group(2) == "N" else -1)
        lon = float(wp.group(3)) * (1 if wp.group(4) == "E" else -1)
        pos = (lat, lon)
    if pos is None:
        near = re.search(r"NEAR\s+([\d.]+)([NS])\s+([\d.]+)([EW])", html)
        if near:
            lat = float(near.group(1)) * (1 if near.group(2) == "N" else -1)
            lon = float(near.group(3)) * (1 if near.group(4) == "E" else -1)
            pos = (lat, lon)
    wind = None
    wm = re.search(r"MAX SUSTAINED WINDS\s*-\s*(\d+)\s*KT", html)
    if wm:
        wind = int(wm.group(1))
    return {"position": pos, "storm_id": storm_id, "advisory_time": adv_time, "wind_kts": wind}


# ─── 3. GFS 下載（500 + 850 兩層一次過） ────────────────────────
def download_gfs_latest(out_path="/tmp/gfs_latest.grib2", cycle_hint=None):
    """下載最新 GFS 分析場（500 + 850 hPa U/V 一次過）。

    cycle_hint: 指定 cycle（今日 + 昨日）。
    冇 hint：由而家 UTC 最近嘅 6h cycle 開始向後試（最多 4 個），
    保證攞到「最新已出」嘅分析場。
    """
    utc = datetime.now(timezone.utc)
    candidates = []
    if cycle_hint:
        for day_offset in [0, -1]:
            day = (utc + timedelta(days=day_offset)).strftime("%Y%m%d")
            candidates.append((day, cycle_hint))
    else:
        nearest = (utc.hour // 6) * 6
        for back in range(4):
            t = (utc.replace(hour=0, minute=0, second=0, microsecond=0)
                 + timedelta(hours=nearest - back * 6))
            candidates.append((t.strftime("%Y%m%d"), f"{t.hour:02d}"))
    for day, cycle in candidates:
        url = (f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl"
               f"?file=gfs.t{cycle}z.pgrb2.0p25.f000"
               f"&lev_850_mb=on&lev_500_mb=on&var_UGRD=on&var_VGRD=on"
               f"&dir=%2Fgfs.{day}%2F{cycle}%2Fatmos")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                data = r.read()
            Path(out_path).write_bytes(data)
            print(f"✅ GFS {day} {cycle}Z → {out_path} ({len(data)/1e6:.1f} MB)")
            return out_path, day, cycle
        except Exception as e:
            print(f"  ⏳ {day} {cycle}Z: {str(e)[:70]}")
    raise RuntimeError("GFS 下載失敗")


def load_gfs(path, level=500):
    """讀特定 level（500 或 850）u/v + lat/lon。"""
    import xarray as xr
    ds = xr.open_dataset(path, engine="cfgrib",
                         backend_kwargs={"filter_by_keys": {"typeOfLevel": "isobaricInhPa",
                                                            "level": level}})
    lat = ds["latitude"].values
    lon = ds["longitude"].values
    u = ds["u"].values
    v = ds["v"].values
    valid = str(ds["valid_time"].values)
    ds.close()
    return lat, lon, u, v, valid


# ─── 4. WN1 路徑（500 hPa 徑向風相位） ───────────────────────────
def shape_at(lat, lon, u, v, clat, clon, rmin=1.0, rmax=5.0):
    """徑向風傅立葉分解：WN1 相位 + ellipt（用度數網格）"""
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
    m = (r_deg >= rmin) & (r_deg <= rmax)
    if m.sum() < 10:
        raise ValueError("太少網格點")
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
    return dict(wn1_phi=wn1_phi, ellipt=ellipt, wn1_amp=wn1_r,
                asym=wn1_r/max(abs(wn0_r), 1e-9))


# ─── 5. dH_curl 強度（850 hPa 渦度） ─────────────────────────────
def compute_vorticity(u, v, lon, lat):
    """ζ = dv/dx − du/dy（s⁻¹），用 m 單位網格間距。"""
    R = 6371000.0
    mid = np.mean(lat) * np.pi / 180.0
    dlon_m = np.mean(np.diff(lon)) * np.pi / 180.0 * R * np.cos(mid)
    dlat_m = np.mean(np.diff(lat)) * np.pi / 180.0 * R
    return np.gradient(v, dlon_m, axis=1) - np.gradient(u, dlat_m, axis=0)


def compute_dh_curl(zeta, lat, lon, clat, clon, core_deg=CORE_DEG, shell_deg=SHELL_DEG):
    """dH_curl = H_shell(ζ) − H_core(ζ)。core/shell 用度數半徑。"""
    LON, LAT = np.meshgrid(lon, lat)
    dlat = np.abs(LAT - clat)
    dlon = np.abs(((LON - clon + 180) % 360) - 180)
    dist = np.sqrt(dlat**2 + dlon**2)
    core = dist <= core_deg
    shell = (dist > core_deg) & (dist <= shell_deg)
    Hc = float(zeta[core].mean())
    Hs = float(zeta[shell].mean())
    return Hs - Hc, Hc, Hs


def classify_mode(dh):
    if dh < -0.3e-4:
        return "collapse"
    if dh < -0.1e-4:
        return "organized"
    if dh > 0.3e-4:
        return "developing"
    return "neutral"


# ─── 6. 鞍點環（850 hPa ζ 場 Morse 分類 + box-count） ───────────
def hessian_2d(f, spacing=(1.0, 1.0)):
    """2D Hessian via double gradient. f: (ny, nx). Returns (fxx, fyy, fxy)."""
    gy, gx = np.gradient(f, *spacing)
    gyy, gxy = np.gradient(gy, *spacing)
    gyx, gxx = np.gradient(gx, *spacing)
    fxy = (gxy + gyx) / 2.0
    return gxx, gyy, fxy


def morse_classify(field, lat, lon, lat_range, lon_range, parabolic_frac=0.10):
    """Morse classification via Hessian eigenvalues。
    peak: λ1,λ2<0 / trough: λ1,λ2>0 / saddle: λ1·λ2<0 / parabolic。
    """
    lat_mask = (lat >= lat_range[0]) & (lat <= lat_range[1])
    lon_mask = (lon >= lon_range[0]) & (lon <= lon_range[1])
    crop = field[np.ix_(lat_mask, lon_mask)]
    lat_crop = lat[lat_mask]
    lon_crop = lon[lon_mask]
    dlat = np.abs(np.diff(lat_crop).mean())
    dlon = np.abs(np.diff(lon_crop).mean())
    fxx, fyy, fxy = hessian_2d(crop, spacing=(dlat, dlon))
    trace = fxx + fyy
    det = fxx * fyy - fxy ** 2
    disc = np.maximum(trace ** 2 - 4 * det, 0)
    l1 = (trace + np.sqrt(disc)) / 2.0
    l2 = (trace - np.sqrt(disc)) / 2.0
    adet = np.abs(det)
    thr = np.percentile(adet, parabolic_frac * 100)
    parabolic = adet < thr
    peak = (~parabolic) & (l1 < 0) & (l2 < 0)
    trough = (~parabolic) & (l1 > 0) & (l2 > 0)
    saddle = (~parabolic) & (l1 * l2 < 0)
    return {"peak": peak, "trough": trough, "saddle": saddle, "parabolic": parabolic,
            "lat_mask": lat_mask, "lon_mask": lon_mask}


def box_count_dim(points, n_boxes_list=None):
    """2D box-counting。points: (N, 2)。返回 (D, r2, intercept, counts, sizes)。"""
    if n_boxes_list is None:
        n_boxes_list = [2, 4, 8, 16, 32, 64, 128]
    pts = np.asarray(points, dtype=float)
    if pts.shape[0] < 10:
        return None
    pmin = pts.min(axis=0)
    pmax = pts.max(axis=0)
    span = (pmax - pmin)
    span[span == 0] = 1.0
    pts_n = (pts - pmin) / span
    sizes = []
    counts = []
    for nb in n_boxes_list:
        idx = np.floor(pts_n * nb).astype(int)
        idx = np.clip(idx, 0, nb - 1)
        unique = len(np.unique(idx, axis=0))
        counts.append(unique)
        sizes.append(1.0 / nb)
    sizes = np.array(sizes)
    counts = np.array(counts, dtype=float)
    A = np.vstack([np.log(1.0 / sizes), np.ones_like(sizes)]).T
    coef, res, *_ = np.linalg.lstsq(A, np.log(counts), rcond=None)
    D = coef[0]
    pred = A @ coef
    ss_res = np.sum((np.log(counts) - pred) ** 2)
    ss_tot = np.sum((np.log(counts) - np.log(counts).mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    return D, r2, coef[1], counts, sizes


def sample_uniform_annulus(n, r_in=EW_IN, r_out=EW_OUT, seed=0):
    """均勻喺 2D 環帶 (r_in, r_out) 內抽 n 點。"""
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0, 2 * np.pi, n)
    r = np.sqrt(rng.uniform(r_in ** 2, r_out ** 2, n))
    return np.stack([r * np.cos(theta), r * np.sin(theta)], axis=1)


def sample_circle(n, r=(EW_IN + EW_OUT) / 2.0, seed=0):
    """均勻喺 1D 圓周抽 n 點。"""
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0, 2 * np.pi, n)
    return np.stack([r * np.cos(theta), r * np.sin(theta)], axis=1)


def saddle_ring_analysis(zeta, lat, lon, clat, clon):
    """眼牆環帶 (EW_IN, EW_OUT) 內鞍點點集 → box-count D_fold。

    回傳 dict（n<30 時只報 n，唔計維數——沿用 phase15 門檻）。
    """
    # 中心窗口（減計算量）
    lat_mask = (lat >= clat - WINDOW_DEG) & (lat <= clat + WINDOW_DEG)
    lon_mask = (lon >= clon - WINDOW_DEG) & (lon <= clon + WINDOW_DEG)
    if lat_mask.sum() < 10 or lon_mask.sum() < 10:
        return {"error": "窗口太小"}
    lat_range = (float(lat[lat_mask].min()), float(lat[lat_mask].max()))
    lon_range = (float(lon[lon_mask].min()), float(lon[lon_mask].max()))
    M = morse_classify(zeta, lat, lon, lat_range, lon_range)
    lm = M["lat_mask"]; lom = M["lon_mask"]
    lat_crop = lat[lm]; lon_crop = lon[lom]

    # 眼牆環帶 mask（中心歸零，度數座標）
    LON, LAT = np.meshgrid(lon_crop, lat_crop)
    dlat = np.abs(LAT - clat)
    dlon = np.abs(((LON - clon + 180) % 360) - 180)
    dist = np.sqrt(dlat**2 + dlon**2)
    ew = (dist >= EW_IN) & (dist < EW_OUT)
    saddles_ew = ew & M["saddle"]
    ys, xs = np.where(saddles_ew)
    n_sad = len(ys)
    out = {"n_saddle": int(n_sad), "window_deg": WINDOW_DEG}
    if n_sad < 30:
        out["error"] = "n<30"
        return out

    dx = ((lon_crop[xs] - clon + 180) % 360) - 180
    dy = lat_crop[ys] - clat
    pts = np.stack([dx, dy], axis=1)

    obs_dim = box_count_dim(pts)
    if obs_dim is None:
        out["error"] = "box-count 失敗"
        return out
    obs_d = obs_dim[0]

    # null A: 同 n 隨機 2D 環帶填充
    null_dims_fill = []
    for s in range(N_NULL):
        null_pts = sample_uniform_annulus(n_sad, r_in=EW_IN, r_out=EW_OUT, seed=s)
        d = box_count_dim(null_pts)
        if d:
            null_dims_fill.append(d[0])
    fill_mean = float(np.mean(null_dims_fill)) if null_dims_fill else float('nan')

    # null B: 同 n 1D 圓周（完美環）
    null_dims_1d = []
    for s in range(N_NULL):
        null_pts = sample_circle(n_sad, seed=s)
        d = box_count_dim(null_pts)
        if d:
            null_dims_1d.append(d[0])
    circle_mean = float(np.mean(null_dims_1d)) if null_dims_1d else float('nan')
    circle_std = float(np.std(null_dims_1d)) if len(null_dims_1d) > 1 else float('nan')

    D_fold = obs_d - fill_mean
    z_vs_circle = (obs_d - circle_mean) / circle_std if circle_std and not np.isnan(circle_std) else float('nan')

    # 全部環帶點嘅維數（對照：2D 填充應該 ≈2）
    idx_all = np.where(ew)
    ew_all_pts = np.stack([
        ((lon_crop[idx_all[1]] - clon + 180) % 360) - 180,
        lat_crop[idx_all[0]] - clat,
    ], axis=1)
    all_dim = box_count_dim(ew_all_pts)

    out.update({
        "D_obs": round(obs_d, 3),
        "D_fill_null": round(fill_mean, 3) if not np.isnan(fill_mean) else None,
        "D_circle_null": round(circle_mean, 3) if not np.isnan(circle_mean) else None,
        "D_all_points": round(all_dim[0], 3) if all_dim else None,
        "D_fold": round(D_fold, 3) if not np.isnan(D_fold) else None,
        "z_vs_circle": round(z_vs_circle, 2) if not np.isnan(z_vs_circle) else None,
    })
    return out


# ─── 7. History 累積（冪等） ───────────────────────────────────
def load_history():
    if Path(HISTORY_PATH).exists():
        return json.loads(Path(HISTORY_PATH).read_text())
    return {"records": []}


def save_history(hist):
    Path(HISTORY_PATH).write_text(json.dumps(hist, indent=2, ensure_ascii=False))


def bearing_name(deg):
    names = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
             "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    return names[int((deg + 11.25) // 22.5) % 16]


# ─── 8. Main ────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grib", default="", help="指定 grib 檔")
    ap.add_argument("--storms", default="", help="覆寫: NAME:lat,lon;...")
    ap.add_argument("--cycle", default="", help="指定 GFS cycle: YYYYMMDD,HH")
    ap.add_argument("--dry-run", action="store_true", help="唔寫 history")
    args = ap.parse_args()

    # Storms
    if args.storms:
        storms = {}
        for part in args.storms.split(";"):
            name, pos = part.split(":")
            la, lo = map(float, pos.split(","))
            storms[name] = {"position": (la, lo)}
    else:
        active = fetch_active_storms()
        if not active:
            print("❌ 冇活躍颱風（NW Pacific）— 完成")
            return
        storms = {}
        for s in active:
            try:
                info = fetch_storm_position(s["slug"])
            except Exception as e:
                print(f"  ⚠️ {s['name']}: 位置攞失敗 {e}")
                continue
            if not info["position"]:
                print(f"  ⚠️ {s['name']}: 冇位置，跳過")
                continue
            storms[s["name"]] = {"position": info["position"],
                                 "storm_id": info["storm_id"],
                                 "advisory_time": info["advisory_time"],
                                 "wind_kts": info["wind_kts"] or s["wind_kts"]}
        if not storms:
            print("❌ 全部颱風位置攞失敗 — 完成")
            return

    # GFS
    if args.grib:
        path = args.grib
        day, cycle = "", ""
        print(f"📂 用 {path}")
    else:
        hint = None
        if args.cycle:
            d, c = args.cycle.split(",")
            hint = c
        path, day, cycle = download_gfs_latest(cycle_hint=hint)

    # 讀兩層
    lat5, lon5, u5, v5, valid5 = load_gfs(path, level=500)
    lat8, lon8, u8, v8, valid8 = load_gfs(path, level=850)
    valid_str = valid5[:16]
    print(f"🕐 GFS valid_time: {valid_str}\n")

    # 850 hPa 渦度（強度 + 鞍點共用）
    zeta = compute_vorticity(u8, v8, lon8, lat8)

    hist = load_history()
    existing = {(r["storm"], r["cycle"]) for r in hist["records"]}
    new_records = []

    print("=" * 78)
    print(f"🌀 Typhoon Monitor — GFS {valid_str}")
    print("=" * 78)

    for name, info in storms.items():
        clat, clon = info["position"]
        print(f"\n▸ {name}  ({clat:.1f}N, {clon:.1f}E)"
              + (f"  [{info.get('storm_id','')}]" if info.get("storm_id") else ""))

        # 1) WN1 路徑
        try:
            s = shape_at(lat5, lon5, u5, v5, clat, clon)
            phi = s["wn1_phi"]; el = s["ellipt"]; amp = s["wn1_amp"]
            # 地磁場（IGRF-13 簡化，500 hPa ≈ 5.5 km）+ UQ v3
            mag_F, mag_D = igrf_field(clat, clon, alt_km=5.5)
            uq_level, amp_thresh = classify_uq(el, amp, mag_F, clat)
            uq_emoji = {"VERY_HIGH": "🟢🟢", "HIGH": "🟢", "MEDIUM": "🟡", "LOW": "🔴"}[uq_level]
            print(f"  路徑: WN1 相位 {phi:5.1f}° ({bearing_name(phi)})  ellipt={el:.2f}"
                  + f"  amp={amp:.1f} (門檻 {amp_thresh:.0f})  UQ={uq_emoji} {uq_level}")
        except Exception as e:
            print(f"  路徑: ❌ {e}")
            s = None
            mag_F = None
            uq_level = None
            amp_thresh = None

        # 2) dH_curl 強度
        try:
            dh, Hc, Hs = compute_dh_curl(zeta, lat8, lon8, clat, clon)
            mode = classify_mode(dh)
            print(f"  強度: dH_curl = {dh:.3e}  [{mode}]  (Hc={Hc:.2e}, Hs={Hs:.2e})")
        except Exception as e:
            print(f"  強度: ❌ {e}")
            dh, Hc, Hs, mode = None, None, None, None

        # 3) 鞍點環
        try:
            sr = saddle_ring_analysis(zeta, lat8, lon8, clat, clon)
            if "error" in sr:
                print(f"  鞍點: n_saddle={sr['n_saddle']} ({sr['error']})")
            else:
                print(f"  鞍點: n={sr['n_saddle']}  D_obs={sr['D_obs']}  "
                      f"D_fold={sr['D_fold']}  z_1D={sr['z_vs_circle']}")
        except Exception as e:
            print(f"  鞍點: ❌ {e}")
            sr = None

        rec = {
            "storm": name,
            "storm_id": info.get("storm_id"),
            "cycle": valid_str,
            "center_lat": round(clat, 2),
            "center_lon": round(clon, 2),
            # 路徑
            "wn1_phi": round(phi, 1) if s else None,
            "ellipt": round(el, 3) if s else None,
            "wn1_amp": round(s["wn1_amp"], 2) if s else None,
            "asym": round(s["asym"], 2) if s else None,
            # UQ v3（2026-08-15 confound test 修正：移除磁場條件）
            "mag_F_nT": round(mag_F, 0) if mag_F is not None else None,
            "uq_level": uq_level,
            "amp_thresh": amp_thresh,
            # 強度
            "dh_curl": round(dh, 8) if dh is not None else None,
            "H_core": round(Hc, 8) if Hc is not None else None,
            "H_shell": round(Hs, 8) if Hs is not None else None,
            "dh_mode": mode,
            # 鞍點環
            "saddle": sr,
            # 元數據
            "advisory_time": info.get("advisory_time"),
            "wind_kts": info.get("wind_kts"),
            "recorded_at": datetime.now(timezone.utc).isoformat(),
        }
        key = (name, valid_str)
        if key in existing:
            print(f"  ↪ {name} {valid_str} 已存在，跳過")
            continue
        new_records.append(rec)

    if new_records:
        hist["records"].extend(new_records)
        hist["records"].sort(key=lambda r: (r["storm"], r["cycle"]))
        if not args.dry_run:
            save_history(hist)
        print(f"\n✅ 新增 {len(new_records)} 條記錄 → {HISTORY_PATH}（總 {len(hist['records'])} 條）")
    else:
        print("\nℹ️ 冇新記錄（全部已存在）")

    # Summary
    if new_records and not args.dry_run:
        lines = ["# 🌀 Typhoon Monitor 最新追蹤", "",
                 f"**GFS {valid_str}** — 自動更新 "
                 f"({datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC)", "",
                 "| 颱風 | 位置 | WN1 相位 | ellipt | UQ | dH_curl | 強度模式 | 鞍點 n | D_fold |",
                 "|------|------|---------|--------|-----|---------|---------|--------|--------|"]
        uq_emoji = {"VERY_HIGH": "🟢🟢", "HIGH": "🟢", "MEDIUM": "🟡", "LOW": "🔴"}
        for r in new_records:
            s = r["saddle"] or {}
            dfold = s.get("D_fold", "—")
            uq = r.get("uq_level") or "—"
            lines.append(f"| {r['storm']} | {r['center_lat']}N {r['center_lon']}E | "
                         f"{r['wn1_phi']}° | {r['ellipt']} | "
                         f"{uq_emoji.get(uq, '—')} {uq} | "
                         f"{r['dh_curl']:.2e} | {r['dh_mode']} | "
                         f"{s.get('n_saddle','—')} | {dfold} |")
        lines += ["", "**UQ 機制 v3**（2026-08-15 confound test 修正）：🟢🟢 Very High = ellipt≤0.4 + Amp≥門檻"
                      "（<20°→10 / 20-30°→7 / ≥30°→5）；🟡 Medium = ellipt≤0.4；🔴 Low = ellipt>0.4。"
                      "（v3 移除 Mag≥35k 條件：地磁場係緯度偽相關，partial r≈0）", ""]
        Path(SUMMARY_PATH).write_text("\n".join(lines) + "\n")
        print(f"✅ Summary → {SUMMARY_PATH}")


if __name__ == '__main__':
    main()
