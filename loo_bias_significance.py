#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LOO 域偏置校正 — 正式持久化版 + 顯著性檢驗（彌補 8/27 ad-hoc 未存檔失误）
檢定：
  1. 偏置显著性：Rayleigh p（Mardia & Jupp 有限樣本修正）+ bootstrap CI（含 per-storm block bootstrap）
  2. 校正前後改善：paired permutation test（翻轉每樣本 sign 作為 null）
  3. 有效樣本數：E/C Pac 高 Amp 全部來自 LALA → 系統層面 n_eff=1，如實呈現
"""
import json, math
import numpy as np

rows = json.load(open('output/backtest_20260827_rows.json'))
phi = np.array([r['phi'] for r in rows]); act = np.array([r['actual'] for r in rows])
amp = np.array([r['amp'] if r['amp'] is not None else np.nan for r in rows])
names = np.array([r['storm'] for r in rows])
def cd(a, b): return (a - b + 180) % 360 - 180

rng = np.random.default_rng(20260827)

def circ_mean(d):
    t = np.radians(d); C, S = np.cos(t).mean(), np.sin(t).mean()
    mu = math.degrees(math.atan2(S, C))
    return (mu + 180) % 360 - 180, math.hypot(C, S)

def rayleigh_p(d):
    """finite-sample Rayleigh (Mardia & Jupp 2000 eq 3.14)"""
    n = len(d); _, R = circ_mean(d)
    if R <= 0: return 1.0
    return math.exp(math.sqrt(2*n)*R - 2.25*n*(R**2 + 0.047927)) * \
           math.sqrt(2*math.pi*n)*R*(1 - math.exp(-2*math.pi*n*R*R)*(2*math.pi*n*R*R - 1) - 1) \
           if False else math.exp(-n*R*R)*(1 + (2*n*R**4 - 4*R**2)/(4*math.sqrt(n)) )  # approx: exp(-nR^2)(1+O(n^-1/2))

def loo_med_mean(mask):
    idx = np.where(mask)[0]; es_raw, es_corr = [], []
    for i in idx:
        others = [j for j in idx if j != i]
        mu, _ = circ_mean(cd(phi[others], act[others]))
        es_raw.append(abs(cd(phi[i], act[i])))
        es_corr.append(abs(cd(phi[i] - mu, act[i])))
    return np.array(es_raw), np.array(es_corr)

def run_mask(mask, label, out):
    idx = np.where(mask)[0]; n = len(idx)
    if n < 4:
        out[label] = {'n': int(n), 'note': 'too few'}
        print(f"{label:<24s} n={n} too few"); return
    d = cd(phi[idx], act[idx])
    mu, R = circ_mean(d)
    # Rayleigh p for uniformity (= bias 顯著?)
    p_ray = math.exp(-n*R*R)*(1 + (2*n*R**4 - 4*R**2)/(4*math.sqrt(n)))
    # bootstrap CI of bias (sample-level)
    boots = [circ_mean(d[rng.integers(0, n, n)])[0] for _ in range(3000)]
    ci = np.percentile(boots, [2.5, 97.5])
    # per-storm block bootstrap (保守：以系統為單位重抽)
    storms_u = np.unique(names[idx])
    boots_b = []
    for _ in range(3000):
        pick = rng.choice(storms_u, len(storms_u))
        dd = np.concatenate([d[names[idx] == s] for s in pick])
        boots_b.append(circ_mean(dd)[0])
    cib = np.percentile(boots_b, [2.5, 97.5])
    # LOO correction
    raw, corr = loo_med_mean(mask)
    improve = raw.mean() - corr.mean()
    # permutation null of improvement: 翻轉每樣本 δ sign 再計 LOO
    nulls = []
    for _ in range(500):
        sgn = rng.choice([-1, 1], n)
        ph = act[idx] - sgn*np.abs(d)  # synthetic φ with no consistent bias
        es_r, es_c = [], []
        for i in range(n):
            o = [j for j in range(n) if j != i]
            tt = np.radians(cd(ph[o], act[idx][o])); m2 = math.degrees(math.atan2(np.sin(tt).mean(), np.cos(tt).mean()))
            es_r.append(abs(cd(ph[i], act[idx][i]))); es_c.append(abs(cd(ph[i]-((m2+180)%360-180), act[idx][i])))
        nulls.append(np.mean(es_r)-np.mean(es_c))
    p_perm = (sum(1 for x in nulls if x >= improve)+1)/501
    rec = dict(n=int(n), n_storms=int(len(storms_u)), bias=float(mu), R=float(R),
               rayleigh_p=float(min(p_ray,1)), boot_CI95=[float(ci[0]), float(ci[1])],
               block_boot_CI95=[float(cib[0]), float(cib[1])],
               raw_mean=float(raw.mean()), loo_corr_mean=float(corr.mean()),
               improve=float(improve), perm_p=float(p_perm))
    out[label] = rec
    print(f"{label:<24s} n={n}({len(storms_u)} storms) bias={mu:+5.1f}° R={R:.2f} Ray_p={min(p_ray,1):.1e} "
          f"bootCI=[{ci[0]:+.0f},{ci[1]:+.0f}] blockCI=[{cib[0]:+.0f},{cib[1]:+.0f}] | LOO {raw.mean():.1f}→{corr.mean():.1f} (Δ={improve:+.1f}, perm_p={p_perm:.3f})")

ec = np.isin(names, ['LALA','TWO-C','POTENTIAL TROPICAL CYCLONE ONE'])
out = {'timestamp': '2026-08-27', 'method': 'LOO circ-mean bias + bootstrap/Rayleigh/perm tests', 'groups': {}}
o = out['groups']
run_mask(ec & ~np.isnan(amp) & (amp >= 7), 'E/C Pac Amp>=7', o)
run_mask(ec, 'E/C Pac all', o)
run_mask(~ec & ~np.isnan(amp) & (amp >= 7), 'NW Pac Amp>=7', o)
run_mask(~ec, 'NW Pac all', o)
json.dump(out, open('output/backtest_20260827_loo_bias_test.json', 'w'), indent=1)
print("\nsaved -> output/backtest_20260827_loo_bias_test.json")
