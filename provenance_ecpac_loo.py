#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐樣本溯源：印出 E/C Pac Amp>=7 全部樣本、raw δ、LOO bias、校正後誤差
→ 證明 62.2° / 16.9° 每個數字由邊一嚿數據計出。"""
import json, math
import numpy as np

rows = json.load(open('output/backtest_20260827_rows.json'))
def cd(a, b): return (a - b + 180) % 360 - 180

EC = {'LALA', 'TWO-C', 'POTENTIAL TROPICAL CYCLONE ONE'}
sel = [r for r in rows if r['storm'] in EC and r['amp'] is not None and r['amp'] >= 7]
sel.sort(key=lambda r: r['cycle'])
print(f"入選樣本 n={len(sel)}（條件：storm ∈ E/C Pacific、amp≥7）")
deltas = [cd(r['phi'], r['actual']) for r in sel]
print(f"\n{'storm':<8s}{'cycle':<14s}{'φ':>7s}{'actual':>8s}{'δ=φ−act':>9s}{'amp':>6s}")
for r, d in zip(sel, deltas):
    print(f"{r['storm']:<8s}{r['cycle']:<14s}{r['phi']:7.1f}{r['actual']:8.1f}{d:+9.1f}{r['amp']:6.1f}")
raw = [abs(d) for d in deltas]
print(f"\n[1] 未校正：med|δ| = {np.median(raw):.1f}°  mean = {np.mean(raw):.1f}°   ← 「62.2°」來源")
# LOO: each sample corrected by circ-mean bias of the OTHER 4
def cmean(ds):
    t = np.radians(ds); mu = math.degrees(math.atan2(np.sin(t).mean(), np.cos(t).mean()))
    return (mu + 180) % 360 - 180
print("\n[2] LOO 校正逐樣本（bias 用其餘 4 個樣本嘅 circular mean）：")
print(f"{'sample':<22s}{'LOO bias':>9s}{'校後|δ|':>9s}")
corr = []
for i in range(len(sel)):
    others = [deltas[j] for j in range(len(sel)) if j != i]
    b = cmean(others)
    e = abs(cd(sel[i]['phi'] - b, sel[i]['actual']))
    corr.append(e)
    print(f"{sel[i]['storm']+' '+sel[i]['cycle'][5:16]:<22s}{b:+9.1f}{e:9.1f}")
print(f"\n[3] 校正後：med = {np.median(corr):.1f}°  mean = {np.mean(corr):.1f}°   ← 「16.9°」來源")
print(f"[4] 樣本時間跨度：{sel[0]['cycle']} → {sel[-1]['cycle']}，系統數 = {len(set(s['storm'] for s in sel))}")
json.dump({'samples': sel, 'raw_abs': raw, 'corrected_abs': corr},
          open('output/provenance_ecpac_amp7_loo.json', 'w'), indent=1)
print("[5] saved -> output/provenance_ecpac_amp7_loo.json")
