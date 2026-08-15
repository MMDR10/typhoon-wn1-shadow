#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""WN1 live 驗證：穩健性分析 + 視覺化
1. 排除慢移動樣本（位移<150km）睇誤差會唔會改善
2. 每颱風相位軌跡 vs 實際移動方向圖
"""
import json, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

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

wn1 = json.load(open('/app/working/workspaces/tygtDc/wn1-shadow/wn1_history.json'))['records']
ty = json.load(open('/tmp/typhoon_history_latest.json'))['records']
records = {}
for r in ty: records[(r['storm'].upper(), r['cycle'])] = dict(r)
for r in wn1: records.setdefault((r['storm'].upper(), r['cycle']), dict(r))
storms = {}
for (storm, cycle), r in records.items():
    storms.setdefault(storm, []).append(r)
for s in storms: storms[s].sort(key=lambda r: r['cycle'])

rows = []
for storm in sorted(storms):
    recs = storms[storm]
    for i in range(len(recs)-1):
        t0, t1 = recs[i], recs[i+1]
        if abs(t0['center_lat']-t1['center_lat']) < 1e-9 and abs(t0['center_lon']-t1['center_lon']) < 1e-9:
            continue
        actual = bearing(t0['center_lat'], t0['center_lon'], t1['center_lat'], t1['center_lon'])
        err = ang_diff(t0['wn1_phi'], actual)
        km = dist_km(t0['center_lat'], t0['center_lon'], t1['center_lat'], t1['center_lon'])
        rows.append(dict(storm=storm, cycle=t0['cycle'], phi=t0['wn1_phi'], ellipt=t0['ellipt'],
                         actual=actual, err=err, km=km, amp=t0.get('wn1_amp'),
                         lat0=t0['center_lat'], lon0=t0['center_lon'],
                         lat1=t1['center_lat'], lon1=t1['center_lon']))

print("=== 位移速度分層（6h 位移 → 速度） ===")
for thr in [100, 150, 200, 250]:
    sub = [r for r in rows if r['km'] >= thr]
    if sub:
        errs = [r['err'] for r in sub]
        print(f"位移≥{thr}km (n={len(sub)}): mean {sum(errs)/len(errs):.1f}°, median {sorted(errs)[len(errs)//2]:.1f}°, <45° {sum(1 for e in errs if e<45)/len(errs)*100:.0f}%")
slow = [r for r in rows if r['km'] < 150]
fast = [r for r in rows if r['km'] >= 150]
for label, sub in [("慢移動(<150km)", slow), ("快移動(≥150km)", fast)]:
    if sub:
        errs = [r['err'] for r in sub]
        print(f"{label} (n={len(sub)}): mean {sum(errs)/len(errs):.1f}°, median {sorted(errs)[len(errs)//2]:.1f}°")

print()
print("=== 排除 CHAN-HOM + TD FIFTEEN（慢/弱颱風） ===")
sub = [r for r in rows if r['storm'] in ('PEILOU', 'NANGKA')]
errs = [r['err'] for r in sub]
print(f"PEILOU+NANGKA (n={len(sub)}): mean {sum(errs)/len(errs):.1f}°, median {sorted(errs)[len(errs)//2]:.1f}°, <45° {sum(1 for e in errs if e<45)/len(errs)*100:.0f}%, <90° {sum(1 for e in errs if e<90)/len(errs)*100:.0f}%")

print()
print("=== amp 分層 ===")
for thr in [4, 6, 8, 10]:
    sub = [r for r in rows if r['amp'] and r['amp'] >= thr]
    if sub:
        errs = [r['err'] for r in sub]
        print(f"amp≥{thr} (n={len(sub)}): mean {sum(errs)/len(errs):.1f}°")

# ============ 圖：相位軌跡 vs 實際移動 ============
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
for ax, storm in zip(axes, ['CHAN-HOM', 'PEILOU', 'NANGKA']):
    if storm not in storms: continue
    recs = storms[storm]
    times = [r['cycle'][5:16] for r in recs]
    phis = [r['wn1_phi'] for r in recs]
    # 實際移動方向：每段 t0->t1
    acts = []
    atimes = []
    for i in range(len(recs)-1):
        t0, t1 = recs[i], recs[i+1]
        if abs(t0['center_lat']-t1['center_lat']) < 1e-9 and abs(t0['center_lon']-t1['center_lon']) < 1e-9:
            continue
        acts.append(bearing(t0['center_lat'], t0['center_lon'], t1['center_lat'], t1['center_lon']))
        atimes.append(t0['cycle'][5:16])
    x = np.arange(len(times))
    ax.plot(x, phis, 'o-', color='tab:blue', label='WN1 相位 φ', lw=2)
    # 實際移動點畫喺 t0 位置
    ax2x = [times.index(at) for at in atimes]
    ax.plot(ax2x, acts, 's--', color='tab:red', label='實際移動方向', lw=1.5)
    for i, (xi, a) in enumerate(zip(ax2x, acts)):
        e = ang_diff(phis[xi], a)
        ax.annotate(f"{e:.0f}°", (xi, a), textcoords="offset points", xytext=(0, -14), fontsize=8, color='darkred')
    # ellipt 標註
    for i, r in enumerate(recs):
        ax.scatter([i], [phis[i]], c='lime' if r['ellipt']<=0.4 else 'orange', s=60, zorder=5, edgecolors='k', linewidths=0.5)
    ax.set_title(f"{storm}  (green=ellipt<=0.4 UQ ok, orange=warn)\nactual motion=red dashed (err deg), phi=blue solid")
    ax.set_ylabel('Direction (deg, 0=N)')
    ax.set_ylim(0, 360)
    ax.set_xticks(x)
    ax.set_xticklabels(times, rotation=45, fontsize=7)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig('/app/working/workspaces/tygtDc/wn1-shadow/wn1_live_vs_actual.png', dpi=110)
print()
print("圖已存: wn1_live_vs_actual.png")
