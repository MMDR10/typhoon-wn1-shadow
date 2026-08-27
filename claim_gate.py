#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
claim_gate.py — 提升宣稱嘅 deterministic 守門員（8/27 MKP 事件驅動）
規則（全部機械可判，唔靠模型自覺）：
  G1 evidence_file_exists : claim 引用嘅 evidence JSON 必須存在
  G2 file_precedes_claim  : evidence mtime <= claim 時間戳（唔准事後補檔）
  G3 n_eff                : 有效獨立樣本 = 系統數（唔係觀測數）；effect claim 要 n_eff>=3
  G4 metric_consistency   : headline 用嘅 statistic（mean/median）要喺 evidence 記錄到
  G5 sign_or_bias_claim   : 「偏置存在」要 Rayleigh p<0.05；「校後提升」要 prequential eval n>0
輸出：PASS / FAIL + 逐條 reasons。任何一條 FAIL → 成個宣稱唔准入報告。
"""
import json, os, sys
from datetime import datetime

def gate(claim):
    """claim = dict(
        label, effect_str,            # 例如 "45° improvement (median 62.2->16.9)"
        evidence_file, claim_ts,      # ISO
        n_obs, n_systems,
        statistic,                     # "median" | "mean"
        bias_p=None, prequential_eval_n=None)"""
    r = {}
    ev = claim.get('evidence_file')
    r['G1_evidence_file_exists'] = bool(ev) and os.path.exists(ev)
    if r['G1_evidence_file_exists']:
        m = datetime.fromtimestamp(os.path.getmtime(ev))
        try: cts = datetime.fromisoformat(claim['claim_ts'])
        except Exception: cts = datetime.max
        r['G2_file_precedes_claim'] = m <= cts
    else:
        r['G2_file_precedes_claim'] = False
    r['G3_n_eff'] = claim.get('n_systems', 0) >= 3   # 宣稱層面：獨立系統先至算數
    if r['G1_evidence_file_exists']:
        try:
            evd = json.load(open(ev))
            blob = json.dumps(evd)
            r['G4_metric_consistency'] = claim['statistic'] in blob
        except Exception:
            r['G4_metric_consistency'] = False
    else:
        r['G4_metric_consistency'] = False
    checks = []
    if claim.get('bias_p') is not None: checks.append(claim['bias_p'] < 0.05)
    if claim.get('prequential_eval_n') is not None: checks.append(claim['prequential_eval_n'] > 0)
    r['G5_statistical_support'] = all(checks) if checks else False
    ok = all(r.values())
    return ok, r

if __name__ == '__main__':
    # 負樣本：8/27 我犯過嘅 +45° 宣稱（當時狀態還原）
    negative = dict(label="E/C Pac +45° improvement (median 62.2→16.9)",
        evidence_file='output/provenance_ecpac_amp7_loo.json',  # 呢檔係事後先寫！
        claim_ts='2026-08-27T14:45:00',          # 宣稱時間（heredoc 嗰輪）
        n_obs=5, n_systems=1, statistic='median',
        bias_p=0.015, prequential_eval_n=0)      # prequential 未跑過 → 當時無此數 = None→fail
    # 正樣本：NW Pac prequential（宣稱 median → 要指向真有 median 嘅 evidence 檔）
    positive = dict(label="NW Pac Amp≥7 prequential med 25.2→22.7 (Δ=−1.9°, n=30 obs, 6 systems)",
        evidence_file='output/prequential_NWPacific.json',
        claim_ts='2099-01-01T00:00:00',          # 事後所有引用都合法
        n_obs=36, n_systems=6, statistic='median',
        bias_p=1.5e-12, prequential_eval_n=30)
    for c in (negative, positive):
        ok, reasons = gate(c)
        print(f"{'✅ PASS' if ok else '🚫 FAIL'}  {c['label']}")
        for k, v in reasons.items(): print(f"      {k}: {v}")
