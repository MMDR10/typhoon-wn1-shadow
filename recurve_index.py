#!/usr/bin/env python3
"""recurve_index.py — 「回馬槍指數」監測工具（operational）

由 2026-09-01 之 PILANDOK 診斷 + 88 條歷史回填 + 非線性驗證總結出嚟嘅規則：

核心結論（9d/§10）：
  - 主訊號 = steering 方向已轉非西向（環境唔再支持西行）——線性已夠
  - 輔訊號 = 遠環 NE 繞行、WN1−Steer 分離（唔係獨立訊號！MOKE/NARRA 反例）、結構退化
  - 滯後效應（新）：環境先轉、移動後追 ≥6h（BANG-LANG +84° 個案、lag-1 負相關）
     → 「steering 已轉 + 移動未跟 = 回馬槍計分」嘅核心

usage:
  python recurve_index.py --storm PILANDOK        # 單颱風最新一次計分 + 解析
  python recurve_index.py --all                   # 全部已測颱風計分表
  python recurve_index.py --storm SAUDEL --cycles 3   # 睇最近 3 段

規則（weighted sum, max 10）：
  S1 steering 已轉（3 pts）   — 向西支持度 w = cos(steer−270°) < 0.5
  S2 遠環 NE 繞行（2 pts）   — far 指向 NE 象限 cos(far−45°) 高
  S3 分離度（2 pts）         — |sep|/60 clump
  S4 結構退化（1 pts）       — ellipt + amp 退化
  S5 滯後確認（2 pts）       — 有前段：steer 已轉但移動未跟；冇前段：0
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HISTORY_PATH = Path(__file__).parent / "typhoon_history.json"

W = {"S1": 3.0, "S2": 2.0, "S3": 2.0, "S4": 1.0, "S5": 2.0}   # 權重
LEVELS = [(7.0, "HIGH"), (4.0, "MEDIUM"), (0.0, "LOW")]


def clamp(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))


def ang_diff(a, b):
    return ((a - b + 180.0) % 360.0) - 180.0


def bearing(lat1, lon1, lat2, lon2):
    dlon = math.radians(lon2 - lon1)
    x = math.sin(dlon) * math.cos(math.radians(lat2))
    y = (math.cos(math.radians(lat1)) * math.sin(math.radians(lat2))
         - math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(dlon))
    return (math.degrees(math.atan2(x, y)) + 360) % 360


def load_records():
    data = json.loads(HISTORY_PATH.read_text())
    return data["records"]


def by_storm_sorted(records):
    out = {}
    for r in records:
        out.setdefault(r["storm"], []).append(r)
    for s in out:
        out[s].sort(key=lambda r: r["cycle"])
    return out


def compute_recurve_index(rec, prev=None):
    """計單一 record 嘅回馬槍指數。

    prev: 同一颱風前一段 record（有前段先計到 S5 滯後）。
    回傳 dict(score, level, signals{...}, weights, raw)。
    """
    sig = {}

    # S1: steering 已指向「回馬槍方向」（主訊號）
    # 定義：回馬槍典型方向 = NE (45°)。steer 對 NE 對齊度：cos(steer−45°)
    #   steer=45° (NE) → 1.0  ✅ 已轉
    #   steer=0° (N) 或 90° (E) → 0.71（偏北/偏東都算偏向回馬槍）
    #   steer=152° (SSE) → 0（NARRA 環境穩定東南，唔算已轉）✅ 排除假陽性
    #   steer=256~272° (W) → 0（MOKE 環境一直西，冇轉向）✅ 排除
    steer = rec.get("steer_dir")
    if steer is None:
        return None
    sig["S1_steering_aligned_NE"] = clamp(math.cos(math.radians(steer - 45.0)))

    # S2: 遠環 NE 繞行
    far = rec.get("far_dir")
    if far is not None:
        sig["S2_far_NE"] = clamp(math.cos(math.radians(far - 45.0)))
    else:
        sig["S2_far_NE"] = 0.0

    # S3: 分離度（輔助，唔獨立）
    sep = rec.get("wn1_steer_sep")
    if sep is not None:
        sig["S3_separation"] = clamp(abs(sep) / 60.0)
    else:
        sig["S3_separation"] = 0.0

    # S4: 結構退化（ellipt 高 + amp 低）
    ell = rec.get("ellipt")
    amp = rec.get("wn1_amp")
    decay = 0.0
    if ell is not None:
        decay += 0.5 if ell > 0.6 else 0.25 if ell > 0.4 else 0.0
    if amp is not None:
        decay += 0.5 if amp < 7.0 else 0.25 if amp < 9.0 else 0.0
    sig["S4_decay"] = clamp(decay)

    # S5: 滯後確認（環境已轉 + 移動未跟）
    lag = 0.0
    if prev is not None:
        # 環境轉向幅度（前段→現段）
        steer_turn = ang_diff(steer, prev.get("steer_dir"))
        # 實際移動方向改變（由兩段位置變化計）
        p1, p2 = prev, rec
        if all(p1.get(k) is not None for k in ("center_lat", "center_lon")) and \
           all(p2.get(k) is not None for k in ("center_lat", "center_lon")):
            dist = math.hypot(p2["center_lat"] - p1["center_lat"],
                              p2["center_lon"] - p1["center_lon"])
            if dist >= 0.3:   # 有實質移動先計
                move_dir = bearing(p1["center_lat"], p1["center_lon"],
                                   p2["center_lat"], p2["center_lon"])
                # 移動有冇開始向 steering 靠攏？
                align = abs(ang_diff(move_dir, steer))
                # 環境明顯轉咗 >15°，但移動同 steering 仍然差 >45°
                if abs(steer_turn) > 15.0 and align > 45.0:
                    lag = 1.0
                else:
                    lag = 0.3   # 有移動但未確認滯後 → 淺層
    sig["S5_lag"] = lag

    # 加權總分（signals key = 訊號名，直接逐項）
    score_w = (W["S1"] * sig["S1_steering_aligned_NE"] + W["S2"] * sig["S2_far_NE"]
               + W["S3"] * sig["S3_separation"] + W["S4"] * sig["S4_decay"]
               + W["S5"] * sig["S5_lag"])
    level = next(lv for thr, lv in LEVELS if score_w >= thr)
    return {"score": round(score_w, 2), "level": level, "signals": sig,
            "weights": W, "raw": {"steer": steer, "far": far, "sep": sep,
                                  "ellipt": ell, "amp": amp, "steer_turn": None}}


def storm_recurves(records, cycles=None):
    """全部颱風：last N cycles 嘅指數序列（計分用完整序列抓 prev，顯示先截尾）"""
    storms = by_storm_sorted(records)
    out = {}
    for s, rs in storms.items():
        seq = []
        for i, r in enumerate(rs):
            prev = rs[i - 1] if i > 0 else None
            idx = compute_recurve_index(r, prev)
            if idx is None:
                continue
            idx["cycle"] = r["cycle"]
            seq.append(idx)
        if cycles is not None and len(seq) > cycles:
            seq = seq[-cycles:]
        if seq:
            out[s] = seq
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--storm", default="")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--cycles", type=int, default=1)
    args = ap.parse_args()

    recs = load_records()
    recs = [r for r in recs if r.get("steer_dir") is not None]
    if not recs:
        print("❌ 未有環境側記錄（先跑 monitor 或 backfill）"); return

    curves = storm_recurves(recs, cycles=args.cycles)

    if args.storm:
        if args.storm not in curves:
            print(f"❌ 冇 {args.storm} 嘅環境側記錄")
            sys.exit(1)
        print(f"=== {args.storm} — 最近 {args.cycles} 段「回馬槍指數」 ===\n")
        for it in curves[args.storm]:
            s = it["signals"]
            print(f"  {it['cycle'][:16]}  指數={it['score']:.2f}/10  [🟢LOW 🟡MEDIUM 🔴HIGH]→ {it['level']}")
            print(f"      S1 steering_aligned_NE={s['S1_steering_aligned_NE']:.2f} "
                  f"(steer={it['raw']['steer']}°)  | S2 far_NE={s['S2_far_NE']:.2f} "
                  f"(far={it['raw']['far']}°)")
            print(f"      S3 sep={s['S3_separation']:.2f} (|sep|={abs(it['raw']['sep']) if it['raw']['sep'] is not None else 0:.0f}°)  "
                  f"| S4 decay={s['S4_decay']:.2f} (ell={it['raw']['ellipt']}, amp={it['raw']['amp']})  "
                  f"| S5 lag={s['S5_lag']:.2f}\n")
    elif args.all:
        print("=== 全部颱風「回馬槍指數」（最近一段） ===\n")
        print(f"{'颱風':<16} {'cycle':<18} {'指數':<7} {'level':<8} S1  S2  S3  S4  S5")
        print("-" * 70)
        for s in sorted(curves):
            last = curves[s][-1]
            sig = last["signals"]
            print(f"{s:<16} {last['cycle'][:16]:<18} {last['score']:<7.2f} {last['level']:<8} "
                  f"{sig['S1_steering_aligned_NE']:.2f} {sig['S2_far_NE']:.2f} {sig['S3_separation']:.2f} "
                  f"{sig['S4_decay']:.2f} {sig['S5_lag']:.2f}")
    else:
        ap.print_help()


if __name__ == "__main__":
    main()