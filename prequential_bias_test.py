#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prequential（online，向後睇）域偏置校正測試 — 模擬真實部署
之前 LOO 有 cross-sample 資訊（用未來樣本嘅 bias 校過去）→ 唔算真正上線提升。
呢度嚴格按 cycle 時間排序：評估第 i 個樣本時，bias 只用同一域、時間上早過佢嘅樣本
（warmup >= 6 先開始評估）。呢先係「準度有冇提升」嘅誠實答案。
"""
import json, math
import numpy as np

rows = json.load(open('output/backtest_20260827_rows.json'))
rows.sort(key=lambda r: r['cycle'])
def cd(a, b): return (a - b + 180) % 360 - 180
EC = {'LALA', 'TWO-C', 'POTENTIAL TROPICAL CYCLONE ONE'}

def circ_mean(d):
    t = np.radians(d)
    mu = math.degrees(math.atan2(np.sin(t).mean(), np.cos(t).mean()))
    return (mu + 180) % 360 - 180

WARMUP = 6
for domain in ['NW Pacific', 'E/C Pacific', 'ALL']:
    hist_d, raw_e, cor_e, n_eval = [], [], [], 0
    eval_start_idx = None
    for k, r in enumerate(rows):
        is_ec = r['storm'] in EC
        dom = 'E/C Pacific' if is_ec else 'NW Pacific'
        if domain != 'ALL' and dom != domain: continue
        # 只有 Amp>=7 樣本先進入「評估流」；warmup 樣本仍要餵入 bias 流? 
        # 嚴謹起見：bias 流同評估流都用 Amp>=7（部署時 UQ 門先決）
        if r['amp'] is None or r['amp'] < 7: continue
        d = cd(r['phi'], r['actual'])
        if len(hist_d) >= WARMUP:
            bias = circ_mean(np.array(hist_d))
            raw_e.append(abs(d)); cor_e.append(abs(cd(r['phi'] - bias, r['actual'])))
        hist_d.append(d)
    if len(raw_e) < 5:
        print(f"{domain:<12s} Amp>=7 prequential: 評估樣本不足 ({len(raw_e)})"); continue
    raw_e, cor_e = np.array(raw_e), np.array(cor_e)
    # paired sign-flip test on improvement diffs
    diffs = raw_e - cor_e
    obs = diffs.mean()
    rng = np.random.default_rng(7)
    nulls = [np.mean(diffs * rng.choice([-1, 1], len(diffs))) for _ in range(2000)]
    p = (sum(1 for x in nulls if x >= obs) + 1) / 2001
    print(f"{domain:<12s} Amp>=7 prequential (warmup={WARMUP}): eval n={len(raw_e)}  "
          f"未校正 mean={raw_e.mean():5.1f}° med={np.median(raw_e):5.1f}° | 校正後 mean={cor_e.mean():5.1f}° med={np.median(cor_e):5.1f}° | Δ={obs:+.1f}° (paired sign-flip p={p:.3f})")
    # hit-rate improvement
    h_r = (raw_e < 45).mean() * 100; h_c = (cor_e < 45).mean() * 100
    print(f"{'':<12s} <45° hit: {h_r:.0f}% → {h_c:.0f}%")
    json.dump(dict(domain=domain, n_eval=len(raw_e), warmup=WARMUP,
                   raw_mean=float(raw_e.mean()), corr_mean=float(cor_e.mean()),
                   raw_med=float(np.median(raw_e)), corr_med=float(np.median(cor_e)),
                   raw_hit45=float(h_r), corr_hit45=float(h_c)),
              open(f'output/prequential_{domain.replace("/","_").replace(" ","")}.json', 'w'), indent=1)
