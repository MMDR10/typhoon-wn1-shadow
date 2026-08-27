#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WN1 影子 backtest — 2026-08-27 全量版
「預算」= t 時刻 GFS 500hPa WN1 相位 φ（預測移動方向）
「測測」= t→t+dt 實際中心位移方位角（大圓起始方位）
改进 vs 舊版 compare_wn1_vs_actual.py：
  1. 名稱正規化（HURRICANE LALA == LALA；去 TYPHOON/TS/HURRICANE 前綴）
  2. 按 (storm, cycle) 去重（舊 wn1_history + 新 typhoon_history 重疊）
  3. 只配對 dt ∈ [3h, 12h] 嘅連續記錄（跳過 gap / 重複 cycle）
  4. 分層統計：Amp、ellipt、速度、UQ 複合門檻（8/21 建議版）
"""
import json, math, re
from datetime import datetime

def bearing(lat1, lon1, lat2, lon2):
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

def norm_name(s):
    s = s.upper().strip()
    s = re.sub(r'^(TYPHOON|HURRICANE|TROPICAL STORM|TS|TD|PTC|SUPER TYPHOON)\s+', '', s)
    s = s.replace('POTENTIAL TROPICAL CYCLE ONE', 'ONE-C')
    s = re.sub(r'\s+', ' ', s)
    return s

def parse_cycle(c):
    return datetime.fromisoformat(c)

records = {}
for fname in ['wn1_history.json', 'typhoon_history.json']:
    data = json.load(open(fname))
    recs = data['records'] if isinstance(data, dict) else data
    for r in recs:
        if r.get('wn1_phi') is None or r.get('center_lat') is None:
            continue
        key = (norm_name(r['storm']), r['cycle'])
        # 新格式（有 wn1_amp）優先覆蓋舊格式
        if key not in records or r.get('wn1_amp') is not None:
            records[key] = r

storms = {}
for (storm, cycle), r in records.items():
    storms.setdefault(storm, []).append(r)
for s in storms:
    storms[s].sort(key=lambda r: parse_cycle(r['cycle']))

rows = []
for storm in sorted(storms):
    recs = storms[storm]
    for i in range(len(recs)-1):
        t0, t1 = recs[i], recs[i+1]
        dt = (parse_cycle(t1['cycle']) - parse_cycle(t0['cycle'])).total_seconds()/3600
        if not (3 <= dt <= 12):
            continue
        if abs(t0['center_lat']-t1['center_lat']) < 1e-9 and abs(t0['center_lon']-t1['center_lon']) < 1e-9:
            continue  # 位置未更新
        km = dist_km(t0['center_lat'], t0['center_lon'], t1['center_lat'], t1['center_lon'])
        if km < 15:
            continue  # 幾乎靜止，方向無意義
        actual = bearing(t0['center_lat'], t0['center_lon'], t1['center_lat'], t1['center_lon'])
        err = ang_diff(t0['wn1_phi'], actual)
        rows.append(dict(storm=storm, cycle=t0['cycle'], dt=dt, phi=t0['wn1_phi'],
                         ellipt=t0.get('ellipt'), amp=t0.get('wn1_amp'),
                         wind=t0.get('wind_kts'), actual=actual, err=err,
                         km=km, speed=km/dt))

def stats(rs, label):
    if not rs:
        print(f"{label:<38s} n=0")
        return
    es = sorted(r['err'] for r in rs)
    n = len(es)
    mean = sum(es)/n
    med = es[n//2]
    h45 = sum(1 for e in es if e < 45)/n*100
    h90 = sum(1 for e in es if e < 90)/n*100
    print(f"{label:<38s} n={n:<4d} mean={mean:5.1f}° med={med:5.1f}° <45°={h45:3.0f}% <90°={h90:3.0f}%")

print(f"{'颱風':<12s}{'cycle':<13s}{'φ預':>6s}{'實測':>7s}{'誤差':>6s}{'Amp':>6s}{'ell':>5s}{'km/6h':>7s}")
print('-'*64)
for r in sorted(rows, key=lambda x: (x['storm'], x['cycle'])):
    amp = f"{r['amp']:.1f}" if r['amp'] is not None else "  —"
    ell = f"{r['ellipt']:.2f}" if r['ellipt'] is not None else "  —"
    print(f"{r['storm'][:11]:<12s}{r['cycle'][5:16]:<13s}{r['phi']:6.1f}{r['actual']:7.1f}{r['err']:6.1f}{amp:>6s}{ell:>5s}{r['km']/r['dt']*6:7.0f}")

print()
print('=== 分層統計 ===')
stats(rows, '全部')
stats([r for r in rows if r['amp'] is not None and r['amp'] >= 7], 'Amp>=7 (8/21 最強 UQ)')
stats([r for r in rows if r['amp'] is not None and r['amp'] >= 10], 'Amp>=10')
stats([r for r in rows if r['amp'] is not None and r['amp'] < 7], 'Amp<7')
stats([r for r in rows if r['ellipt'] is not None and r['ellipt'] <= 0.4], 'ellipt<=0.4 (舊門檻)')
stats([r for r in rows if r['ellipt'] is not None and 0.4 < r['ellipt'] <= 0.6], 'ellipt 0.4-0.6')
stats([r for r in rows if r['ellipt'] is not None and r['ellipt'] > 0.6], 'ellipt>0.6')
stats([r for r in rows if r['amp'] is not None and r['ellipt'] is not None
       and r['amp'] >= 7 and r['ellipt'] <= 0.6], '複合: Amp>=7 + ellipt<=0.6')
stats([r for r in rows if r['speed'] >= 25], '快移 >=150km/6h')
stats([r for r in rows if 13 <= r['speed'] < 25], '中速 80-150km/6h')
stats([r for r in rows if r['speed'] < 13], '慢移 <80km/6h')

print()
print('=== 每颱風 ===')
by_storm = {}
for r in rows:
    by_storm.setdefault(r['storm'], []).append(r)
for s, rs in sorted(by_storm.items(), key=lambda x: sum(r['err'] for r in x[1])/len(x[1])):
    es = [r['err'] for r in rs]
    amps = [r['amp'] for r in rs if r['amp'] is not None]
    print(f"  {s:<12s} n={len(rs):<3d} mean={sum(es)/len(es):5.1f}°  worst={max(es):5.1f}°  amp={min(amps) if amps else 0:.0f}-{max(amps) if amps else 0:.0f}")

json.dump(rows, open('output/backtest_20260827_rows.json', 'w'), ensure_ascii=False, indent=1)
print(f"\nsaved {len(rows)} rows -> output/backtest_20260827_rows.json")
