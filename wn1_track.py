#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🌀 WN1 Shadow Auto-Track — 活躍颱風 500 hPa WN1 相位前瞻自動追蹤
=================================================================
每 6h 由 GitHub Actions 觸發：
  1. 從 cyclocane 攞活躍颱風列表（NW Pacific，排除 final advisory）
  2. 對每個颱風攞最新位置（JTWC advisory 內嵌）
  3. 下載最新 GFS 500 hPa 分析場（NOMADS 0.25°）
  4. 計算 WN1 相位 + ellipt（UQ）
  5. 寫入 wn1_history.json（每個 cycle 每颱風一條，冪等）
  6. 輸出 summary

用法：
  python wn1_track.py                    # 自動：最新 GFS + 活躍颱風
  python wn1_track.py --grib /tmp/g.grib2
  python wn1_track.py --storms 'X:1,2'   # 覆寫颱風位置（測試用）
  python wn1_track.py --cycle 20260809,12
"""
import argparse
import json
import re
import sys
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

HISTORY_PATH = "wn1_history.json"
SUMMARY_PATH = "wn1_latest.md"
CYCLOCANE_HOME = "https://www.cyclocane.com/"

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

    # 搵 "Northwest Pacific Storms" 段
    m = re.search(r"Northwest Pacific Storms</h3>(.*?)<h2>Tropical Disturbances", html, re.S)
    if not m:
        return []
    seg = m.group(1)

    storms = []
    # 每個颱風 block: <li> <a href='/slug-storm-tracker/'> <span>NAME</span> </a> ... final advisory? ... wind
    for li in re.finditer(r"<li>(.*?)</li>", seg, re.S):
        block = li.group(1)
        link = re.search(r"href='/([a-z\-]+)-storm-tracker/'", block)
        if not link:
            continue
        slug = link.group(1)
        name_m = re.search(r">\s*([A-Z][A-Z\- ]+?)\s*</span>", block)
        name = name_m.group(1).strip() if name_m else slug.upper()
        final = "final advisory" in block.lower()
        wind_m = re.search(r"Wind Speed:\s*<[^>]*>([\d]+)\s*knots", block)
        wind = int(wind_m.group(1)) if wind_m else None
        storms.append({
            "name": name.replace("TYPHOON ", "").replace("TROPICAL STORM ", ""),
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

    # WTPN advisory header: WTPN33 PGTW 091500
    adv = re.search(r"WTPN\d+\s+PGTW\s+(\d{6})", html)
    adv_time = adv.group(1) if adv else None

    # storm id: TROPICAL STORM 14W (CHAN-HOM)
    sid = re.search(r"(?:TROPICAL STORM|TYPHOON|SUPER TYPHOON|TROPICAL DEPRESSION)\s+(\d+W)", html)
    storm_id = sid.group(1) if sid else None

    # 最新位置 — 用 "WARNING POSITION" 或第一個 NEAR
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

    # 風速
    wind = None
    wm = re.search(r"MAX SUSTAINED WINDS\s*-\s*(\d+)\s*KT", html)
    if wm:
        wind = int(wm.group(1))
    return {"position": pos, "storm_id": storm_id, "advisory_time": adv_time, "wind_kts": wind}


# ─── 3. GFS 下載 ───────────────────────────────────────────────
def download_gfs_latest(out_path="/tmp/gfs_latest_500.grib2", cycle_hint=None):
    from datetime import datetime as dt
    utc = dt.utcnow()
    cycles = [cycle_hint] if cycle_hint else ["12", "06", "00", "18"]
    for cycle in cycles:
        for day_offset in [0, -1]:
            day = (utc + timedelta(days=day_offset)).strftime("%Y%m%d")
            url = (f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl"
                   f"?file=gfs.t{cycle}z.pgrb2.0p25.f000"
                   f"&lev_500_mb=on&var_UGRD=on&var_VGRD=on"
                   f"&dir=%2Fgfs.{day}%2F{cycle}%2Fatmos")
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=60) as r:
                    data = r.read()
                Path(out_path).write_bytes(data)
                print(f"✅ GFS {day} {cycle}Z → {out_path} ({len(data)/1e6:.1f} MB)")
                return out_path, day, cycle
            except Exception as e:
                print(f"  ⏳ {day} {cycle}Z: {str(e)[:70]}")
    raise RuntimeError("GFS 下載失敗")


def load_gfs(path):
    import cfgrib
    ds = cfgrib.open_dataset(path)
    lat = ds["latitude"].values
    lon = ds["longitude"].values
    u = ds["u"].values
    v = ds["v"].values
    valid = str(ds["valid_time"].values)
    ds.close()
    return lat, lon, u, v, valid


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


def bearing_name(deg):
    names = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
             "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    return names[int((deg + 11.25) // 22.5) % 16]


# ─── 4. History 累積（冪等） ───────────────────────────────────
def load_history():
    if Path(HISTORY_PATH).exists():
        return json.loads(Path(HISTORY_PATH).read_text())
    return {"records": []}


def save_history(hist):
    Path(HISTORY_PATH).write_text(json.dumps(hist, indent=2, ensure_ascii=False))


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

    lat, lon, u, v, valid = load_gfs(path)
    valid_str = valid[:16]
    print(f"🕐 GFS valid_time: {valid_str}\n")

    hist = load_history()
    existing = {(r["storm"], r["cycle"]) for r in hist["records"]}
    new_records = []

    print(f"{'='*72}")
    print(f"🌀 WN1 Shadow Auto-Track — GFS {valid_str}")
    print(f"{'='*72}\n")

    for name, info in storms.items():
        clat, clon = info["position"]
        try:
            s = shape_at(lat, lon, u, v, clat, clon)
        except Exception as e:
            print(f"  {name}: ❌ {e}")
            continue
        phi = s["wn1_phi"]; el = s["ellipt"]; amp = s["wn1_amp"]
        # 地磁場（IGRF-13 簡化，500 hPa ≈ 5.5 km）
        mag_F, mag_D = igrf_field(clat, clon, alt_km=5.5)
        uq_level, amp_thresh = classify_uq(el, amp, mag_F, clat)
        uq_emoji = {"VERY_HIGH": "🟢🟢", "HIGH": "🟢", "MEDIUM": "🟡", "LOW": "🔴"}[uq_level]
        uq_txt = {"VERY_HIGH": "ellipt≤0.4+Amp≥門檻",
                  "HIGH": "ellipt≤0.4+Amp≥門檻",
                  "MEDIUM": "ellipt≤0.4",
                  "LOW": "ellipt>0.4"}[uq_level]
        print(f"  {name}  ({clat:.1f}N, {clon:.1f}E)"
              + (f"  [{info.get('storm_id','')}]" if info.get("storm_id") else ""))
        print(f"    500 hPa WN1 相位 = {phi:5.1f}° ({bearing_name(phi)})")
        print(f"    ellipt = {el:.2f}  |  WN1 amp = {amp:.1f} m/s (門檻 {amp_thresh:.0f})  |  Mag F = {mag_F:.0f} nT")
        print(f"    UQ = {uq_emoji} {uq_level} ({uq_txt})")
        print()

        rec = {
            "storm": name,
            "storm_id": info.get("storm_id"),
            "cycle": valid_str,
            "center_lat": round(clat, 2),
            "center_lon": round(clon, 2),
            "wn1_phi": round(phi, 1),
            "ellipt": round(el, 3),
            "wn1_amp": round(amp, 2),
            "asym": round(s["asym"], 2),
            "mag_F_nT": round(mag_F, 0),
            "uq_level": uq_level,
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
        print(f"✅ 新增 {len(new_records)} 條記錄 → {HISTORY_PATH}（總 {len(hist['records'])} 條）")
    else:
        print("ℹ️ 冇新記錄（全部已存在）")

    # Summary
    if new_records and not args.dry_run:
        lines = ["# 🌀 WN1 Shadow 最新追蹤", "",
                 f"**GFS {valid_str}** — 自動更新 "
                 f"({datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC)", "",
                 "| 颱風 | 位置 | WN1 相位 | 方位 | ellipt | amp | MagF | UQ |",
                 "|------|------|---------|------|--------|-----|------|-----|"]
        uq_emoji = {"VERY_HIGH": "🟢🟢", "HIGH": "🟢", "MEDIUM": "🟡", "LOW": "🔴"}
        for r in new_records:
            lines.append(f"| {r['storm']} | {r['center_lat']}N {r['center_lon']}E | "
                         f"{r['wn1_phi']}° | {bearing_name(r['wn1_phi'])} | "
                         f"{r['ellipt']} | {r['wn1_amp']} | {r['mag_F_nT']} | "
                         f"{uq_emoji[r['uq_level']]} {r['uq_level']} |")
        lines += ["", "**UQ 機制 v3**（2026-08-15 confound test 修正）：🟢🟢 Very High = ellipt≤0.4 + Amp≥門檻；"
                  "🟡 Medium = ellipt≤0.4；🔴 Low = ellipt>0.4。"
                  "門檻按緯度：<20° Amp≥10、20-30° Amp≥7、≥30° Amp≥5。"
                  "（v3 移除 Mag≥35k 條件：地磁場係緯度偽相關，partial r≈0）", ""]
        Path(SUMMARY_PATH).write_text("\n".join(lines) + "\n")
        print(f"✅ Summary → {SUMMARY_PATH}")


if __name__ == "__main__":
    main()
