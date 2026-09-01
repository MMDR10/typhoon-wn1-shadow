#!/usr/bin/env python3
"""recurve_nonlinear_validation.py — 回馬槍診斷非線性驗證

MKP 問題：「回馬槍診斷」用嘅係線性方法（sep 角度差、steering 平均），
框架核心關係（WN1↔移動）已知係非線性（dCor 0.541 vs Pearson -0.258，8/22）。
用 88 條歷史回填記錄，驗證「steering 轉向 → 颱風移動轉向」嘅關係
喺線性 vs 非線性（dCor / MI）下表現。

對照 8/27 餘震線教訓：唔預設立場話「一定有非線性信號」，
要實際測——線性主導就係線性主導，如實報告。
"""
import sys, json, math
from pathlib import Path
import numpy as np
from scipy.spatial import distance

sys.path.insert(0, str(Path(__file__).parent))


def mutual_info(x, y, n_bins=15):
    """Mutual information (bits) — 非線性一般依賴"""
    x, y = np.array(x, dtype=float), np.array(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 10 or np.std(x) < 1e-10 or np.std(y) < 1e-10:
        return np.nan
    bins_x = np.linspace(np.nanmin(x), np.nanmax(x), n_bins + 1)
    bins_y = np.linspace(np.nanmin(y), np.nanmax(y), n_bins + 1)
    idx_x = np.clip(np.digitize(x, bins_x) - 1, 0, n_bins - 1)
    idx_y = np.clip(np.digitize(y, bins_y) - 1, 0, n_bins - 1)
    joint = np.zeros((n_bins, n_bins))
    for ix, iy in zip(idx_x, idx_y):
        joint[iy, ix] += 1
    joint /= joint.sum()
    px = joint.sum(axis=0)
    py = joint.sum(axis=1)
    mi = 0.0
    for i in range(n_bins):
        for j in range(n_bins):
            if joint[j, i] > 0 and px[i] > 0 and py[j] > 0:
                mi += joint[j, i] * np.log2(joint[j, i] / (px[i] * py[j]))
    return mi


def dcor(x, y):
    """Distance correlation (Székely et al. 2009) — 捕捉任何依賴（線性+非線性）"""
    x, y = np.array(x, dtype=float), np.array(y, dtype=float)
    if len(x) < 10 or len(y) < 10:
        return np.nan
    n = min(len(x), len(y))
    x, y = x[:n], y[:n]
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 10 or np.std(x) < 1e-10 or np.std(y) < 1e-10:
        return np.nan
    n = len(x)
    x_c = x - x.mean()
    y_c = y - y.mean()
    a = distance.pdist(np.column_stack([x_c, np.zeros(n)]))
    b = distance.pdist(np.column_stack([np.zeros(n), y_c]))
    A = distance.squareform(a)
    B = distance.squareform(b)
    A_bar = A.mean()
    B_bar = B.mean()
    dcov2 = (A * B).mean() - A_bar * B_bar
    if dcov2 < 0:
        dcov2 = 0
    dvar_x = max(0, (A**2).mean() - A_bar**2)
    dvar_y = max(0, (B**2).mean() - B_bar**2)
    denom = np.sqrt(dvar_x) * np.sqrt(dvar_y)
    if denom < 1e-12:
        return np.nan
    return float(np.clip(np.sqrt(dcov2 / denom), 0, 1))

HISTORY_PATH = Path(__file__).parent / "typhoon_history.json"


def ang_diff(a, b):
    """有號角差 (a-b) in (-180, 180]"""
    return ((a - b + 180.0) % 360.0) - 180.0


def bearing(lat1, lon1, lat2, lon2):
    dlon = math.radians(lon2 - lon1)
    x = math.sin(dlon) * math.cos(math.radians(lat2))
    y = (math.cos(math.radians(lat1)) * math.sin(math.radians(lat2))
         - math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(dlon))
    return (math.degrees(math.atan2(x, y)) + 360) % 360


def build_samples():
    """做「轉向事件」樣本：每段 6h 之間，steer 轉向 vs 實際移動轉向。

    回傳 list of dict：each = {steer_turn, move_turn, |sep at start|, ellipt, steer_spd}
    steer_turn = 環境引導方向 6h 內變化（有號，正 = 順時針轉向 NE 方向）
    move_turn  = 實際移動方向 6h 內變化
    """
    data = json.loads(HISTORY_PATH.read_text())
    recs = [r for r in data["records"] if r.get("steer_dir") is not None]
    # 按 storm + cycle 排序
    by_storm = {}
    for r in recs:
        by_storm.setdefault(r["storm"], []).append(r)
    samples = []
    for storm, rs in by_storm.items():
        rs.sort(key=lambda r: r["cycle"])
        for i in range(1, len(rs)):
            a, b = rs[i - 1], rs[i]
            if None in (a.get("steer_dir"), b.get("steer_dir"), a.get("center_lat"),
                        a.get("center_lon"), b.get("center_lat"), b.get("center_lon")):
                continue
            # 移動距離太細（meander 停滯）→ 轉向角無意義
            dist = math.hypot(b["center_lat"] - a["center_lat"],
                              b["center_lon"] - a["center_lon"])
            if dist < 0.5:  # <0.5° 移動（~55km/6h 以下）跳過
                continue
            move_dir = bearing(a["center_lat"], a["center_lon"],
                               b["center_lat"], b["center_lon"])
            steer_turn = ang_diff(b["steer_dir"], a["steer_dir"])   # 環境轉向
            move_turn = ang_diff(move_dir, 0)  # 冇用，下面計相對
            # 移動方向相對於「環境」嘅變化：上一段移動 vs 今段移動
            # 直接用 move_dir 連續體：同 storm 內上一段移動方向
            if i >= 2:
                c = rs[i - 2]
                prev_dist = math.hypot(a["center_lat"] - c["center_lat"],
                                       a["center_lon"] - c["center_lon"])
                if prev_dist >= 0.5:
                    prev_move = bearing(c["center_lat"], c["center_lon"],
                                        a["center_lat"], a["center_lon"])
                    samples.append({
                        "storm": storm,
                        "cycle_b": b["cycle"],
                        "steer_turn": steer_turn,
                        "move_turn": ang_diff(move_dir, prev_move),
                        "sep_at_a": a.get("wn1_steer_sep"),
                        "sep_at_b": b.get("wn1_steer_sep"),
                        "ellipt": b.get("ellipt"),
                        "steer_spd_b": b.get("steer_spd"),
                        "steer_dir_b": b.get("steer_dir"),
                        "move_dir_b": move_dir,
                    })
    return samples


def perm_pvalue(x, y, stat_fn, n_perm=500, seed=42):
    """非參數 permutation p-value"""
    rng = np.random.default_rng(seed)
    obs = stat_fn(x, y)
    count = 0
    for _ in range(n_perm):
        yp = rng.permutation(y)
        if stat_fn(x, yp) >= obs:
            count += 1
    return obs, (count + 1) / (n_perm + 1)


def main():
    samples = build_samples()
    print(f"轉向事件樣本: {len(samples)} 段（6h 間隔，移動>0.5°）\n")

    # === 分析 1: steer_turn ↔ move_turn ===
    st = np.array([s["steer_turn"] for s in samples], dtype=float)
    mt = np.array([s["move_turn"] for s in samples], dtype=float)
    # 移除 NaN
    mask = np.isfinite(st) & np.isfinite(mt)
    st, mt = st[mask], mt[mask]
    if len(st) < 10:
        print("❌ 樣本太少"); return

    print("── 分析 1: 環境轉向 (steer_turn) ↔ 颱風移動轉向 (move_turn) ──")
    pearson = np.corrcoef(st, mt)[0, 1]
    dcor_val = dcor(st, mt)
    mi_val = mutual_info(st, mt)
    print(f"  線性  Pearson r   = {pearson:+.3f}")
    print(f"  非線性 dCor       = {dcor_val:.3f}")
    print(f"  非線性 MI         = {mi_val:.3f} bits")
    dcor_p = perm_pvalue(st, mt, dcor, n_perm=300)
    mi_p = perm_pvalue(st, mt, mutual_info, n_perm=300)
    pear_p = perm_pvalue(st, mt, lambda x, y: np.corrcoef(x, y)[0, 1], n_perm=300)
    print(f"  Pearson p={pear_p[1]:.3f} | dCor p={dcor_p[1]:.3f} | MI p={mi_p[1]:.3f}")

    # === 分析 2: |sep| ↔ move_turn（分離度大 ⇒ 移動偏離環境？） ===
    print("\n── 分析 2: |sep| (開始時) ↔ 移動轉向幅度 |move_turn| ──")
    sep = np.array([abs(s["sep_at_a"]) if s["sep_at_a"] is not None else np.nan
                    for s in samples], dtype=float)
    mt_abs = np.abs(mt)
    # 用 samples 原始序
    st_full = np.array([s["steer_turn"] for s in samples], dtype=float)
    mask2 = np.isfinite(sep) & np.isfinite(mt_abs)
    if mask2.sum() >= 10:
        pear2 = np.corrcoef(sep[mask2], mt_abs[mask2])[0, 1]
        dcor2 = dcor(sep[mask2], mt_abs[mask2])
        mi2 = mutual_info(sep[mask2], mt_abs[mask2])
        print(f"  線性  Pearson r   = {pear2:+.3f}")
        print(f"  非線性 dCor       = {dcor2:.3f}")
        print(f"  非線性 MI         = {mi2:.3f} bits")
        dcor2_p = perm_pvalue(sep[mask2], mt_abs[mask2], dcor, n_perm=300)
        print(f"  dCor p={dcor2_p[1]:.3f}")

    # === 分析 3: steer 方向類別 → 移動轉向（環境已轉 NE 時移動轉向大？） ===
    print("\n── 分析 3: 「steer 已轉 NE」子集內移動轉向 vs 冇轉子集 ──")
    turned = [(s["move_turn"], s.get("sep_at_b"), s.get("steer_dir_b"))
              for s in samples if s.get("steer_dir_b") is not None]
    ne_group = [t for t in turned if 0 <= t[2] <= 90]       # steer 指向 N~E
    w_group = [t for t in turned if 180 < t[2] <= 270]      # steer 指向 S~W（或 270-360）
    w2_group = [t for t in turned if 270 < t[2] <= 360]     # W~N
    other = [t for t in turned if t not in ne_group and t not in w_group and t not in w2_group]
    ne_turns = np.array([t[0] for t in ne_group])
    w_turns = np.array([t[0] for t in w_group]) if w_group else np.array([])
    w2_turns = np.array([t[0] for t in w2_group])
    print(f"  steer N~E (0-90°):  n={len(ne_turns)}  移動轉向 mean={np.mean(ne_turns):+.1f}°  std={np.std(ne_turns):.1f}°")
    if len(w2_turns):
        print(f"  steer W~N (270-360°): n={len(w2_turns)}  移動轉向 mean={np.mean(w2_turns):+.1f}°  std={np.std(w2_turns):.1f}°")
    if len(ne_turns) >= 3 and len(w2_turns) >= 3:
        # Mann-Whitney
        try:
            from scipy.stats import mannwhitneyu
            u, p = mannwhitneyu(ne_turns, w2_turns, alternative="two-sided")
            print(f"  Mann-Whitney: u={u:.1f} p={p:.3f}")
        except Exception as e:
            print(f"  Mann-Whitney: ❌ {e}")
    else:
        print("  Mann-Whitney: 樣本太少跳過")

    # === 分析 4: 三個「關鍵案例」逐個睇（PILANDOK / MOKE 反例 / BANG-LANG 正例） ===
    print("\n── 分析 4: 關鍵案例序列 ──")
    for s in samples:
        if s["storm"] in ("PILANDOK", "MOKE", "BANG-LANG", "ETAU", "SAUDEL"):
            print(f"  {s['storm']} {s['cycle_b'][:16]}  steer_turn={s['steer_turn']:+.1f}°  "
                  f"move_turn={s['move_turn']:+.1f}°  |sep_a|={abs(s['sep_at_a']) if s['sep_at_a'] is not None else 0:4.1f}  "
                  f"steer_dir_b={s['steer_dir_b']:.0f}°")

    print("\n── 分析 5: 時間延遲效應（環境轉向 6h 後移動先轉？）──")
    # 用「同時刻 steer_turn」對「滯後 1 段嘅 move_turn」
    # samples 已按 storm+cycle 排序，順序即時間序
    lags = {}
    for i in range(len(samples) - 1):
        cur = samples[i]
        nxt = samples[i + 1]
        if cur["storm"] != nxt["storm"]:
            continue
        # 確保 nxt 係 cur 嘅下一段（cycle 相鄰）
        if nxt["cycle_b"] <= cur["cycle_b"]:
            continue
        key = "lag1"
        lags.setdefault(key, []).append((cur["steer_turn"], nxt["move_turn"]))
    print(f"  lag-1 樣本: {len(lags.get('lag1', []))} 對")
    for key, pairs in lags.items():
        xs = np.array([p[0] for p in pairs], dtype=float)
        ys = np.array([p[1] for p in pairs], dtype=float)
        mask = np.isfinite(xs) & np.isfinite(ys)
        xs, ys = xs[mask], ys[mask]
        if len(xs) < 10:
            print(f"  {key}: 樣本太少 ({len(xs)})")
            continue
        pear = np.corrcoef(xs, ys)[0, 1]
        dc = dcor(xs, ys)
        mi = mutual_info(xs, ys)
        dcp = perm_pvalue(xs, ys, dcor, n_perm=300)
        print(f"  {key}: Pearson r={pear:+.3f} | dCor={dc:.3f} (p={dcp[1]:.3f}) | MI={mi:.3f} bits")

    # 摘要
    print("\n── 判讀 ──")
    dcor_ratio = dcor_val / max(abs(pearson), 1e-9) if abs(pearson) > 0.05 else float('inf')
    print(f"  dCor/Pearson 比 = {dcor_ratio:.1f}x  "
          f"({'非線性主導' if dcor_ratio > 1.5 else '線性主導' if dcor_ratio < 0.8 else '混合'})")


if __name__ == "__main__":
    main()