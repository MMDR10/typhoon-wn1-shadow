#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WN1 backtest n=128 非線性驗算（§8.5 鐵律合規版）
問題：琴晚版用「角度差的算術 mean + <45° 二值計數」驗算非線性特徵 φ —— 屬線性化操作
本版改用：
  1. 有符號循環誤差 δ = circ_diff(φ, actual) → circular mean + Rayleigh R（向量強度）
  2. Jammalamadaka circular correlation（φ vs actual 兩角度變量嘅同構相關）
  3. dCor on circular distances（循環距離矩陣）+ shuffle null p-value
"""
import json, math
import numpy as np

rows = json.load(open('output/backtest_20260827_rows.json'))
phi = np.array([r['phi'] for r in rows])
act = np.array([r['actual'] for r in rows])
amp = np.array([r['amp'] if r['amp'] is not None else np.nan for r in rows])
ell = np.array([r['ellipt'] if r['ellipt'] is not None else np.nan for r in rows])
spd = np.array([r['speed'] for r in rows])
n = len(rows)

def circ_diff(a, b):
    """signed a-b in (-180,180]"""
    return (a - b + 180) % 360 - 180

def rayleigh(deltas_deg):
    """circular mean + R (0..1 vector strength) + p (large-sample)"""
    t = np.radians(deltas_deg)
    C, S = t.cos().mean() if hasattr(t,'cos') else (np.cos(t).mean(), np.sin(t).mean())
    C, S = np.cos(t).mean(), np.sin(t).mean()
    R = math.hypot(C, S)
    cmean = math.degrees(math.atan2(S, C)) % 360
    # Rayleigh p (Mardia & Jupp): 1 - exp(-n R^2) (approx, R small)
    p = math.exp(-n * R * R)  # upper-tail approx
    return cmean, R, p

def circ_corr(a_deg, b_deg):
    """Jammalamadaka circular correlation between two angles"""
    ta, tb = np.radians(a_deg), np.radians(b_deg)
    ma, mb = np.arctan2(np.sin(ta).mean(), np.cos(ta).mean()), np.arctan2(np.sin(tb).mean(), np.cos(tb).mean())
    num = np.sum(np.sin(ta - ma) * np.sin(tb - mb))
    den = math.sqrt(np.sum(np.sin(ta - ma)**2) * np.sum(np.sin(tb - mb)**2))
    return num / den if den > 0 else float('nan')

def dcor(x, y):
    """distance correlation on 1-D arrays"""
    ax = np.abs(x[:, None] - x[None, :]); ay = np.abs(y[:, None] - y[None, :])
    A = ax - ax.mean(0)[None, :] - ax.mean(1)[:, None] + ax.mean()
    B = ay - ay.mean(0)[None, :] - ay.mean(1)[:, None] + ay.mean()
    dcov2 = (A * B).mean(); dvx = (A * A).mean(); dvy = (B * B).mean()
    if dvx <= 0 or dvy <= 0: return 0.0
    return math.sqrt(max(dcov2, 0) / math.sqrt(dvx * dvy))

def circ_dist(a, b):
    return np.abs(circ_diff(a[:, None] if np.ndim(a) else a[None,:] , b[None,:] if np.ndim(b)!=1 else b[:,None]))

def dcor_circ(a_deg, b_deg):
    """dCor using circular distance matrices"""
    D1 = np.abs(circ_diff(a_deg[:, None], a_deg[None, :]))
    D2 = np.abs(circ_diff(b_deg[:, None], b_deg[None, :]))
    A = D1 - D1.mean(0)[None,:] - D1.mean(1)[:,None] + D1.mean()
    B = D2 - D2.mean(0)[None,:] - D2.mean(1)[:,None] + D2.mean()
    dcov2 = (A*B).mean(); dvx=(A*A).mean(); dvy=(B*B).mean()
    return math.sqrt(max(dcov2,0)/math.sqrt(dvx*dvy)) if dvx>0 and dvy>0 else 0.0

rng = np.random.default_rng(42)
def shuffle_p(fn, a, b, nperm=1000):
    obs = fn(a, b)
    cnt = sum(1 for _ in range(nperm) if fn(a, rng.permutation(b)) >= obs)
    return obs, (cnt+1)/(nperm+1)

def report(mask, label):
    k = mask.sum()
    if k < 5: print(f"{label:<30s} n={k} (too few)"); return
    d = circ_diff(phi[mask], act[mask])
    cmean, R, p = rayleigh(d)
    cd = np.abs(d)
    cc = circ_corr(phi[mask], act[mask])
    dc, pd_ = shuffle_p(dcor_circ, phi[mask], act[mask])
    print(f"{label:<30s} n={k:<4d} circMeanErr={abs(cmean) if cmean<=180 else 360-cmean:5.1f}° R={R:.2f}(p={p:.1e}) "
          f"circCorr={cc:+.2f} dCor_circ={dc:.2f}(p={pd_:.3f}) med|δ|={np.median(cd):4.1f}°")

print(f"=== WN1 預算 vs 實測：非線性（循環）驗算 n={n} ===")
report(np.ones(n, bool), '全部')
report(~np.isnan(amp) & (amp>=7), 'Amp≥7')
report(~np.isnan(amp) & (amp>=10), 'Amp≥10')
report(~np.isnan(amp) & (amp<7), 'Amp<7')
report((~np.isnan(ell)) & (ell<=0.4) & (amp>=7), 'ellipt≤0.4 + Amp≥7')
report(spd>=25, '快移≥150km/6h')
report(spd<13, '慢移<80km/6h')
# E/C pacific
ecp = np.array([r['storm'] for r in rows])
mask_ec = np.isin(ecp, ['LALA','TWO-C','POTENTIAL TROPICAL CYCLONE ONE'])
report(mask_ec, 'E/C Pacific')
report(~mask_ec, 'NW Pacific')

# 對照：線性化統計（琴晚版）
errs = np.abs(circ_diff(phi, act))
hi = ~np.isnan(amp) & (amp>=7)
print()
print("=== 對照：線性化統計（8/27 琴晚版） ===")
print(f"全部:   算術mean={errs.mean():.1f}° median={np.median(errs):.1f}°")
print(f"Amp≥7:  算術mean={errs[hi].mean():.1f}° median={np.median(errs[hi]):.1f}°")
print(f"Amp<7:  算術mean={errs[~hi].mean():.1f}° median={np.median(errs[~hi]):.1f}°")
# 雙峰檢查：Amp<7 組角度差分佈
import collections
h,edges = np.histogram(errs[~hi], bins=[0,30,60,90,120,150,180])
print("Amp<7 誤差分佈直方圖 [0-30,30-60,60-90,90-120,120-150,150-180]:", h.tolist())
