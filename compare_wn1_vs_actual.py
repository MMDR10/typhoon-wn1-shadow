#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WN1 相位 vs 實際移動方向 live 驗證（GitHub 自動化累積數據）
- 讀 wn1_history.json（舊格式，8/9-8/11）+ typhoon_history.json（新格式，8/11-8/12）
- 對每個颱風：按 cycle 排序，用連續有位置變化嘅點計實際移動方位角（大圓起點方位）
- 同 WN1 相位 φ 比對（φ 係 GFS 分析場測嘅移動方向預測）
- ellipt≤0.4 = UQ 可信門檻
"""
import json, math
from datetime import datetime

def bearing(lat1, lon1, lat2, lon2):
    """大圓初始方位角（度，0=北，順時針）"""
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    return (math.degrees(math.atan2(x, y)) + 360) % 360

def ang_diff(a, b):
    d = abs(a - b) % 360
    return d if d <= 180 else 360 - d

def dist_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2-lat1, lon2-lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return 2*R*math.asin(math.sqrt(a))

# ---- 讀兩份 history ----
wn1 = json.load(open('/app/working/workspaces/tygtDc/wn1-shadow/wn1_history.json'))['records']
ty = json.load(open('/tmp/typhoon_history_latest.json'))['records']

# 合併（新格式優先，用 cycle 去重）
records = {}
for r in ty:
    key = (r['storm'].upper(), r['cycle'])
    records[key] = dict(r)
for r in wn1:
    key = (r['storm'].upper(), r['cycle'])
    if key not in records:
        records[key] = dict(r)

# ---- 按颱風分組、按時間排序 ----
storms = {}
for (storm, cycle), r in records.items():
    storms.setdefault(storm, []).append(r)
for s in storms:
    storms[s].sort(key=lambda r: r['cycle'])

# ---- 對比：t 時刻 phi 預測 vs t→t+dt 實際移動 ----
print(f"{'颱風':<14s}{'cycle':<18s}{'φ預測':>7s}{'ellipt':>7s}{'UQ':>4s}{'實際方位':>8s}{'誤差':>6s}{'位移km':>7s}")
print('-'*78)
rows = []
for storm in sorted(storms):
    recs = storms[storm]
    for i in range(len(recs)-1):
        t0, t1 = recs[i], recs[i+1]
        # 位置一樣就跳過（advisory 未更新，冇移動資訊）
        if abs(t0['center_lat']-t1['center_lat']) < 1e-9 and abs(t0['center_lon']-t1['center_lon']) < 1e-9:
            continue
        actual = bearing(t0['center_lat'], t0['center_lon'], t1['center_lat'], t1['center_lon'])
        err = ang_diff(t0['wn1_phi'], actual)
        km = dist_km(t0['center_lat'], t0['center_lon'], t1['center_lat'], t1['center_lon'])
        uq = '✅' if t0['ellipt'] <= 0.4 else '⚠️'
        name = storm[:13]
        print(f"{name:<14s}{t0['cycle'][5:16]:<18s}{t0['wn1_phi']:7.1f}{t0['ellipt']:7.2f}{uq:>4s}{actual:8.1f}{err:6.1f}{km:7.0f}")
        rows.append(dict(storm=storm, cycle=t0['cycle'], phi=t0['wn1_phi'], ellipt=t0['ellipt'],
                         actual=actual, err=err, km=km, amp=t0.get('wn1_amp'),
                         lat0=t0['center_lat'], lon0=t0['center_lon'],
                         lat1=t1['center_lat'], lon1=t1['center_lon']))

# ---- 統計 ----
print()
if rows:
    all_err = [r['err'] for r in rows]
    uq_rows = [r for r in rows if r['ellipt'] <= 0.4]
    print(f"總比對樣本: {len(rows)}  (每樣本 = 用 t 時刻 φ 預測 t→t+6h 實際移動)")
    print(f"全部: mean 誤差 {sum(all_err)/len(all_err):.1f}°, median {sorted(all_err)[len(all_err)//2]:.1f}°, <45° 佔 {sum(1 for e in all_err if e<45)/len(all_err)*100:.0f}%, <90° 佔 {sum(1 for e in all_err if e<90)/len(all_err)*100:.0f}%")
    if uq_rows:
        uq_err = [r['err'] for r in uq_rows]
        print(f"UQ 可信 (ellipt≤0.4, n={len(uq_rows)}): mean {sum(uq_err)/len(uq_err):.1f}°, median {sorted(uq_err)[len(uq_err)//2]:.1f}°, <45° 佔 {sum(1 for e in uq_err if e<45)/len(uq_err)*100:.0f}%, <90° 佔 {sum(1 for e in uq_err if e<90)/len(uq_err)*100:.0f}%")
    # 每颱風 mean
    by_storm = {}
    for r in rows:
        by_storm.setdefault(r['storm'], []).append(r['err'])
    print()
    for s, es in sorted(by_storm.items()):
        print(f"  {s}: n={len(es)}, mean 誤差 {sum(es)/len(es):.1f}°")
