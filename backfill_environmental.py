#!/usr/bin/env python3
"""backfill_environmental.py — 歷史記錄回填環境側（steering/遠環/分離）

策略：typhoon_history.json 入面沒有 steer_dir 嘅記錄，用 NOMADS filter
只拉 500hPa u/v（3.1MB/cycle），用記錄中心位置計四訊號，更新入 history。

⚠️ 只支援 NOMADS filter 保留範圍（約 8/27 起）；更早要 AWS 全檔（500MB）跳過。
"""
import sys, json, urllib.request, urllib.error
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from typhoon_monitor import load_gfs, steering_vec, ang_diff, bearing_name

HISTORY_PATH = Path(__file__).parent / "typhoon_history.json"
TMP = Path("/tmp/gfs_hist.grib2")


def download_500(day, cycle):
    url = (f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl"
           f"?file=gfs.t{cycle}z.pgrb2.0p25.f000"
           f"&lev_500_mb=on&var_UGRD=on&var_VGRD=on"
           f"&dir=%2Fgfs.{day}%2F{cycle}%2Fatmos")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read()
    TMP.write_bytes(data)
    return TMP


def main(limit=None, dry=False):
    hist = json.loads(HISTORY_PATH.read_text())
    recs = hist["records"]
    targets = [r for r in recs if r.get("steer_dir") is None]
    print(f"待回填: {len(targets)} 條 (總 {len(recs)})")
    done = skipped = fail = 0
    for r in targets:
        cyc = r["cycle"]  # e.g. 2026-08-27T18:00
        day = cyc[:10].replace("-", "")
        hour = cyc[11:13]
        try:
            path = download_500(day, hour)
            lat5, lon5, u5, v5, valid = load_gfs(path, level=500)
            clat, clon = r["center_lat"], r["center_lon"]
            steer_spd, steer_dir = steering_vec(lat5, lon5, u5, v5, clat, clon, 8.0, 12.0)
            far_spd, far_dir = steering_vec(lat5, lon5, u5, v5, clat, clon, 12.0, 18.0)
            near_spd, near_dir = steering_vec(lat5, lon5, u5, v5, clat, clon, 3.0, 6.0)
            if steer_dir is None:
                raise ValueError("steering 網格不足")
            phi = r.get("wn1_phi")
            sep = ang_diff(phi, steer_dir) if phi is not None else None
            r["steer_spd"] = round(steer_spd, 2)
            r["steer_dir"] = round(steer_dir, 1)
            r["near_spd"] = round(near_spd, 2) if near_spd is not None else None
            r["near_dir"] = round(near_dir, 1) if near_dir is not None else None
            r["far_spd"] = round(far_spd, 2) if far_spd is not None else None
            r["far_dir"] = round(far_dir, 1) if far_dir is not None else None
            r["wn1_steer_sep"] = round(sep, 1) if sep is not None else None
            print(f"✅ {r['storm']} {r['cycle']}: steer {steer_dir:.1f}° ({bearing_name(steer_dir)}) "
                  f"{steer_spd:.1f} m/s | far {far_dir:.1f}° | sep {sep:+.1f}" if sep is not None
                  else f"✅ {r['storm']} {r['cycle']}: steer {steer_dir:.1f}° | far {far_dir:.1f}° | (無 WN1)")
            done += 1
        except urllib.error.HTTPError as e:
            print(f"❌ {r['storm']} {r['cycle']}: HTTP {e.code} — 跳過")
            skipped += 1
        except Exception as e:
            print(f"❌ {r['storm']} {r['cycle']}: {e}")
            fail += 1
        if limit and done >= limit:
            break
    if done and not dry:
        hist["records"] = recs
        HISTORY_PATH.write_text(json.dumps(hist, ensure_ascii=False, indent=2) + "\n")
        print(f"\n✅ 已更新 {done} 條 → {HISTORY_PATH}")
    else:
        print(f"\nℹ️ 完成: {done} 成功 / {skipped} 跳過(HTTP) / {fail} 失敗 (dry={dry})")
    return done, skipped, fail


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--dry", action="store_true")
    args = ap.parse_args()
    main(limit=args.limit, dry=args.dry)