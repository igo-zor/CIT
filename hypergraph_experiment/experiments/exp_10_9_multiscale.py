"""
§10.9 多尺度時間視窗：固定微觀歷史（可多條彙總），對不同聚合窗口 w 計算宏觀量與窗口間 JS。

論文《約束世界論 30》§10.9.2（五）：宏觀有效連通度 C_eff 可用多種操作化；本模組固定採
「各窗口位置取子序列末態之超邊數 |E|」作為 C_eff^{(w)}(ℓ) 之代理，再對其序列計算
方差與平台摘要（與正文「正文只需固定一種」一致）。
"""

from __future__ import annotations

import math
import random
import statistics
from typing import Any, Dict, List, Sequence  # Any 用於宏觀標籤鍵

from hypergraph_experiment.ch10_paper_presets import CH10_9_BASELINE
from hypergraph_experiment.core import (
    SIGNATURES,
    HypergraphConfig,
    is_connected_2section,
    run_trajectory,
    sample_candidates_and_filter,
)
from hypergraph_experiment.output_rounding import round_floats_for_output
from hypergraph_experiment.progress_types import ProgressCallback
from hypergraph_experiment.refinement import js_divergence_bits
from hypergraph_experiment.time_series_metrics import (
    mean_edge_turnover_rate,
    plateau_length_max_abs_diff,
    run_length_stats_for_labels,
)

# 多條微觀歷史時軌道種子步進（避免與採樣種子簇聚）
_HIST_TRAJ_SEED_STRIDE: int = 7919
# 論文 N_seed 次重跑時之基底種子步進
_SEED_BATCH_STRIDE: int = 500_017


def _signature_fn(sig_obs: str):
    """依鍵取得解析簽名函式；無效時回退 medium。"""
    if sig_obs in SIGNATURES:
        return SIGNATURES[sig_obs]
    return SIGNATURES["medium"]


def _macro_label(c: HypergraphConfig, r: int, sig_obs: str) -> tuple[Any, ...]:
    """宏觀標籤：σ_obs 簽名 + 2-section 連通旗標 + 超邊數摘要。"""
    sig = _signature_fn(sig_obs)(c)
    conn = 1 if is_connected_2section(c) else 0
    return (sig, conn, len(c.hyperedges))


def _window_macro_series(
    traj: Sequence[HypergraphConfig],
    w: int,
    r: int,
    delta_t: int,
    sig_obs: str,
) -> List[tuple[Any, ...]]:
    """
    對每個起點 ℓ（步進 Δt）取窗口 [ℓ,ℓ+w) 內配置之宏觀標籤，末態加上窗口內平均邊數四捨五入。

    Args:
        delta_t: 論文 §10.9.4 之聚合步長 Δt（宏觀觀測每次沿時間前進之步數）；小於 1 時視為 1。
    """
    T = len(traj)
    out: List[tuple[Any, ...]] = []
    if w <= 0 or T < w:
        return out
    step = max(1, int(delta_t))
    for ell in range(0, T - w + 1, step):
        chunk = traj[ell : ell + w]
        avg_m = round(statistics.mean(len(c.hyperedges) for c in chunk))
        out.append(_macro_label(chunk[-1], r, sig_obs) + (avg_m,))
    return out


def _ceff_series(
    traj: Sequence[HypergraphConfig],
    w: int,
    delta_t: int,
) -> List[int]:
    """與宏觀標籤序列相同之 ℓ 網格下，C_eff 代理（末態超邊數）。"""
    T = len(traj)
    if w <= 0 or T < w:
        return []
    step = max(1, int(delta_t))
    return [
        len(traj[min(ell + w - 1, len(traj) - 1)].hyperedges)
        for ell in range(0, T - w + 1, step)
    ]


def _per_window_for_trajectory(
    traj: Sequence[HypergraphConfig],
    *,
    window_sizes: Sequence[int],
    r: int,
    sig_obs: str,
    delta_t: int,
    epsilon_plat: float,
) -> tuple[float, List[Dict[str, Any]]]:
    """單條軌跡：回傳 R_edge_bar 與各 w 之指標列（與舊版單歷史相容）。"""
    r_edge_bar = float(mean_edge_turnover_rate(traj))

    def _hist_for_w(wv: int) -> tuple[Dict[Any, int], int, List[tuple[Any, ...]]]:
        labs = _window_macro_series(traj, int(wv), r, delta_t, sig_obs)
        ct: Dict[Any, int] = {}
        for lb in labs:
            ct[lb] = ct.get(lb, 0) + 1
        return ct, sum(ct.values()), labs

    ws_sorted = sorted(set(int(x) for x in window_sizes))
    if 1 not in ws_sorted:
        ws_sorted = [1] + [x for x in ws_sorted if x != 1]
        ws_sorted.sort()

    ct1, tot1, _labs1 = _hist_for_w(1)
    hist_w1: Dict[Any, float] = {k: v / tot1 for k, v in ct1.items()} if tot1 else {}

    rows: List[Dict[str, Any]] = []
    for w in sorted(set(int(x) for x in window_sizes)):
        ct, tot, labs_w = _hist_for_w(w)
        pvec = [c / tot for c in ct.values()] if tot else []
        h_macro = -sum(p * math.log(p, 2) for p in pvec if p > 0)
        ceff = _ceff_series(traj, w, delta_t)
        var_c = statistics.pvariance(ceff) if len(ceff) > 1 else 0.0
        tau_stats = run_length_stats_for_labels(labs_w)
        eps = float(epsilon_plat)
        plat_c = plateau_length_max_abs_diff([float(x) for x in ceff], epsilon=eps)

        js_vs_w1 = None
        if w != 1 and tot and hist_w1:
            hist_w = {k: v / tot for k, v in ct.items()}
            keys = sorted(set(hist_w1.keys()) | set(hist_w.keys()), key=lambda x: str(x))
            p1 = [hist_w1.get(k, 0.0) for k in keys]
            p2 = [hist_w.get(k, 0.0) for k in keys]
            js_vs_w1 = js_divergence_bits(p1, p2)

        rows.append(
            {
                "w": w,
                "H_macro_bits": round(h_macro, 6),
                "Var_C_eff_edges": round(var_c, 6),
                "n_window_positions": len(ceff),
                "JS_vs_w1_bits": None if js_vs_w1 is None else round(js_vs_w1, 6),
                "R_edge_bar": round(r_edge_bar, 6),
                "tau_unit_max": tau_stats.get("tau_unit_max"),
                "tau_unit_mean": tau_stats.get("tau_unit_mean"),
                "L_plat_max_len_Ceff": plat_c.get("max_plateau_length"),
                "L_plat_num_seg_Ceff": plat_c.get("num_plateau_segments"),
            }
        )

    return round(r_edge_bar, 6), rows


def _mean_pstdev(values: List[float | None]) -> tuple[float | None, float | None]:
    """略過 None 後計算平均與母體標準差。"""
    nums = [float(v) for v in values if v is not None and not isinstance(v, bool)]
    if not nums:
        return None, None
    mu = float(statistics.mean(nums))
    sig = float(statistics.pstdev(nums)) if len(nums) > 1 else 0.0
    return mu, sig


def _aggregate_per_window_rows(
    list_of_row_lists: List[List[Dict[str, Any]]],
) -> List[Dict[str, Any]]:
    """將多份 per_window（同序 w）對應欄位做算術平均；可選 *_std 於 N>1 時填入。"""
    if not list_of_row_lists:
        return []
    n = len(list_of_row_lists)
    # 依 w 分組
    by_w: Dict[int, List[Dict[str, Any]]] = {}
    for rows in list_of_row_lists:
        for row in rows:
            wv = int(row["w"])
            by_w.setdefault(wv, []).append(row)

    out: List[Dict[str, Any]] = []
    for wv in sorted(by_w.keys()):
        group = by_w[wv]
        agg: Dict[str, Any] = {"w": wv}

        def _col_mean_std(key: str, round_to: int | None = 6) -> None:
            raw = [r.get(key) for r in group]
            mu, sig = _mean_pstdev(raw)
            if mu is not None and round_to is not None:
                agg[key] = round(mu, round_to)
            else:
                agg[key] = mu
            if n > 1 and sig is not None:
                agg[f"{key}_std"] = round(sig, 6) if round_to is not None else sig

        _col_mean_std("H_macro_bits")
        _col_mean_std("Var_C_eff_edges")
        pos_vals = [r.get("n_window_positions") for r in group]
        pmu, _ = _mean_pstdev([float(x) if x is not None else None for x in pos_vals])
        agg["n_window_positions"] = int(round(pmu)) if pmu is not None else None

        js_vals: List[float | None] = []
        for r in group:
            j = r.get("JS_vs_w1_bits")
            js_vals.append(float(j) if j is not None else None)
        jmu, jsig = _mean_pstdev(js_vals)
        if wv == 1:
            agg["JS_vs_w1_bits"] = None
        else:
            agg["JS_vs_w1_bits"] = None if jmu is None else round(jmu, 6)
            if n > 1 and jsig is not None and wv != 1:
                agg["JS_vs_w1_bits_std"] = round(jsig, 6)

        _col_mean_std("R_edge_bar")
        for tk in ("tau_unit_max", "tau_unit_mean"):
            raw = [r.get(tk) for r in group]
            mu, sig = _mean_pstdev([float(x) if x is not None else None for x in raw])
            agg[tk] = None if mu is None else round(mu, 6)
            if n > 1 and sig is not None:
                agg[f"{tk}_std"] = round(sig, 6)
        for tk in ("L_plat_max_len_Ceff", "L_plat_num_seg_Ceff"):
            raw = [r.get(tk) for r in group]
            mu, sig = _mean_pstdev([float(x) if x is not None else None for x in raw])
            agg[tk] = None if mu is None else int(round(mu))
            if n > 1 and sig is not None:
                agg[f"{tk}_std"] = round(sig, 6)

        out.append(agg)
    return out


def run_experiment_10_9(
    *,
    n: int | None = None,
    max_edge_size: int | None = None,
    max_edges: int | None = None,
    sample_limit: int | None = None,
    seed: int = 3,
    connected: bool | None = None,
    max_degree: int | None = None,
    forbid_pair_triangles: bool = False,
    steps: int | None = None,
    window_sizes: Sequence[int] = (1, 3, 6),
    r: int | None = None,
    sig_obs: str | None = None,
    m_trial: int | None = None,
    delta_t: int | None = None,
    epsilon_plat: float | None = None,
    n_hist: int | None = None,
    n_seed: int | None = None,
    progress: ProgressCallback = None,
) -> dict[str, Any]:
    """
    固定一條或多條合法演化軌跡（可重跑 seed 批次），對各時間窗口 w 建立宏觀型別直方圖，並與 w=1 比較 JS。

    Args:
        n: 節點數；None 時取 CH10_9_BASELINE。
        max_edge_size: 最大超邊階數 k_max。
        max_edges: 超邊數上限 m_max。
        sample_limit: 候選採樣上限。
        seed: 偽隨機基底種子。
        connected: 是否要求 2-section 連通（域型）。
        max_degree: 頂點最大度數 d_max；None 時由基準帶入。
        forbid_pair_triangles: 是否禁止二元三角形 motif。
        steps: 軌跡步數 T（軌跡長度為 T+1 含起點）。
        window_sizes: 時間聚合寬度集合；應包含 1 作為最細參照（未包含時會自動補 1 供 JS）。
        r: 空間鄰域深度（傳入宏觀標籤）。
        sig_obs: 解析簽名 weak／medium／strong（論文 σ_obs）；None 時由基準帶入。
        m_trial: 每步候選更新數 M_trial（§10.9.6）；None 時由基準帶入。
        delta_t: 宏觀觀測時間步進 Δt（§10.9.4）；None 時由基準帶入。
        epsilon_plat: 平台判定閾 ε_plat（§10.9.4）；None 時由基準帶入。
        n_hist: 微觀歷史條數 N_hist（§10.9.6）；None 時由基準帶入。
        n_seed: 每組參數重跑次數 N_seed（§10.9.6）；1 表示單一批次；None 時由基準帶入。

    Returns:
        ``parameters``、``per_window``（各 w 已對 N_hist／N_seed 平均過）、``R_edge_bar``、
        ``trajectory_length``（最後一條軌跡長度）等。
    """
    _b = CH10_9_BASELINE
    n_i = int(_b["n"]) if n is None else int(n)
    k_max = int(_b["k_max"]) if max_edge_size is None else int(max_edge_size)
    m_max = int(_b["m_max"]) if max_edges is None else int(max_edges)
    slimit = int(_b["sample_limit"]) if sample_limit is None else int(sample_limit)
    conn = bool(_b["connected"]) if connected is None else bool(connected)
    d_max = int(_b["d_max"]) if max_degree is None else int(max_degree)
    steps_i = int(_b["steps"]) if steps is None else int(steps)
    r_i = int(_b["r"]) if r is None else int(r)
    sig = str(_b["sig_obs"]) if sig_obs is None else str(sig_obs)
    mt = int(_b["m_trial"]) if m_trial is None else int(m_trial)
    dt = int(_b["delta_t"]) if delta_t is None else int(delta_t)
    eps = float(_b["epsilon_plat"]) if epsilon_plat is None else float(epsilon_plat)
    n_h = int(_b["n_hist"]) if n_hist is None else int(n_hist)
    n_s = int(_b["n_seed_runs"]) if n_seed is None else int(n_seed)
    n_h = max(1, n_h)
    n_s = max(1, n_s)

    parameters: dict[str, Any] = {
        "n": n_i,
        "max_edge_size": k_max,
        "max_edges": m_max,
        "sample_limit": slimit,
        "seed": seed,
        "connected": conn,
        "max_degree": d_max,
        "forbid_pair_triangles": forbid_pair_triangles,
        "steps": steps_i,
        "window_sizes": list(window_sizes),
        "r": r_i,
        "sig_obs": sig,
        "m_trial": mt,
        "delta_t": dt,
        "epsilon_plat": eps,
        "n_hist": n_h,
        "n_seed_runs": n_s,
    }

    total_ticks = max(1, n_s * (1 + n_h))
    tick = 0

    def _report(msg: str) -> None:
        nonlocal tick
        if progress:
            progress(tick, total_ticks, msg)
        tick += 1

    per_seed_aggregates: List[List[Dict[str, Any]]] = []
    r_bars_all: List[float] = []
    last_traj_len = 0

    for s_run in range(n_s):
        seed_batch = int(seed) + s_run * _SEED_BATCH_STRIDE
        _report(f"§10.9 採樣合法域（批次 {s_run + 1}/{n_s}）…")
        _, cfgs = sample_candidates_and_filter(
            n=n_i,
            max_edge_size=k_max,
            max_edges=m_max,
            sample_limit=slimit,
            seed=seed_batch,
            connected=conn,
            max_degree=d_max,
            forbid_pair_triangles=forbid_pair_triangles,
        )
        if not cfgs:
            return round_floats_for_output(
                {
                    "experiment": "10.9",
                    "error": "無合法配置",
                    "parameters": parameters,
                }
            )

        cfg_list = list(cfgs)
        cfg_set = set(cfgs)
        hist_row_lists: List[List[Dict[str, Any]]] = []

        for h in range(n_h):
            _report(f"§10.9 軌跡 {h + 1}/{n_h}（批次 {s_run + 1}/{n_s}）…")
            rng_s = random.Random(seed_batch + 17 + h * 1337)
            start = rng_s.choice(cfg_list)
            traj = run_trajectory(
                start,
                steps=steps_i,
                seed=seed_batch + 1 + h * _HIST_TRAJ_SEED_STRIDE,
                max_edge_size=k_max,
                max_edges=m_max,
                connected_required=conn,
                max_degree=d_max,
                forbid_pair_triangles=forbid_pair_triangles,
                m_trial=max(1, mt),
                allowed_configs=cfg_set,
            )
            last_traj_len = len(traj)
            r_b, rows_one = _per_window_for_trajectory(
                traj,
                window_sizes=window_sizes,
                r=r_i,
                sig_obs=sig,
                delta_t=dt,
                epsilon_plat=eps,
            )
            r_bars_all.append(r_b)
            hist_row_lists.append(rows_one)

        per_seed_aggregates.append(_aggregate_per_window_rows(hist_row_lists))

    if progress:
        progress(total_ticks, total_ticks, "§10.9 彙總完成")

    if len(per_seed_aggregates) == 1:
        final_rows = per_seed_aggregates[0]
    else:
        # 對各 seed 之「已縱向平均」per_window 再橫向平均
        final_rows = _aggregate_per_window_rows(per_seed_aggregates)

    r_edge_mean = float(statistics.mean(r_bars_all)) if r_bars_all else 0.0

    return round_floats_for_output(
        {
            "experiment": "10.9",
            "parameters": parameters,
            "trajectory_length": last_traj_len,
            "R_edge_bar": round(r_edge_mean, 6),
            "per_window": final_rows,
            "aggregated_N_hist": n_h,
            "aggregated_n_seed_runs": n_s,
        }
    )


def section_10_9_output_parameters_df(result: dict[str, Any]) -> "pd.DataFrame":
    r"""
    組裝《約束世界論 30》§10.9.5「輸出參數」對照表（與 §10.8.5 表頭同形）。

    Args:
        result: ``run_experiment_10_9`` 之回傳字典。

    Returns:
        五欄 ``DataFrame``：論文小節、輸出參數、論文記號、數值、論文語義摘要。
    """
    import pandas as pd

    cols = ("論文小節", "輸出參數", "論文記號", "數值", "論文語義摘要")

    if result.get("error"):
        return pd.DataFrame(
            [
                {
                    "論文小節": "—",
                    "輸出參數": "執行狀態",
                    "論文記號": "—",
                    "數值": "—",
                    "論文語義摘要": str(result.get("error")),
                }
            ],
            columns=list(cols),
        )

    def _fmt(v: object) -> str:
        if v is None:
            return "—"
        if isinstance(v, float):
            return f"{v:.6g}"
        return str(v)

    pw = result.get("per_window")
    if not isinstance(pw, list) or not pw:
        return pd.DataFrame(
            [
                {
                    "論文小節": "—",
                    "輸出參數": "無法組表",
                    "論文記號": "—",
                    "數值": "—",
                    "論文語義摘要": "結果缺少 per_window 列表。",
                }
            ],
            columns=list(cols),
        )

    r_bar = result.get("R_edge_bar")
    n_hist = int(result.get("aggregated_N_hist") or 1)
    n_seed_b = int(result.get("aggregated_n_seed_runs") or 1)
    agg_note = ""
    if n_hist > 1 or n_seed_b > 1:
        agg_note = f"（已對 N_hist={n_hist}、N_seed_runs={n_seed_b} 平均）"

    def _by_w_summary(key: str, *, int_vals: bool = False) -> str:
        parts: List[str] = []
        for row in sorted(pw, key=lambda d: int(d["w"])):
            r_w = row.get("w")
            val = row.get(key)
            if int_vals and val is not None:
                val = int(round(float(val)))
            std_k = f"{key}_std"
            if row.get(std_k) is not None:
                parts.append(f"w={r_w}:{_fmt(val)}±{_fmt(row.get(std_k))}")
            else:
                parts.append(f"w={r_w}:{_fmt(val)}")
        return "；".join(parts)

    rows_out: List[dict[str, str]] = [
        {
            "論文小節": "10.9.5（一）",
            "輸出參數": f"微觀邊周轉率平均{agg_note}",
            "論文記號": r"$\bar R_{\mathrm{edge}}$",
            "數值": _fmt(r_bar),
            "論文語義摘要": "整段軌跡邊周轉率之平均；多歷史／多批次時取算術平均。",
        },
        {
            "論文小節": "10.9.5（二）",
            "輸出參數": "宏觀有效連通度（時間序列）",
            "論文記號": r"$C_{\mathrm{eff}}^{(w)}(\ell)$",
            "數值": "見各 w 之 per_window／完整軌跡 JSON",
            "論文語義摘要": "程式固定以「窗口末態超邊數 |E|」作為 C_eff 代理；序列用於 Var 與平台摘要。",
        },
        {
            "論文小節": "10.9.5（三）",
            "輸出參數": "宏觀連通度方差",
            "論文記號": r"$\mathrm{Var}\!\bigl(C_{\mathrm{eff}}^{(w)}\bigr)$",
            "數值": _by_w_summary("Var_C_eff_edges"),
            "論文語義摘要": "同上 C_eff 代理序列之母體方差（單條軌跡內）；多歷史時為跨軌跡平均後之欄位。",
        },
        {
            "論文小節": "10.9.5（四）",
            "輸出參數": "宏觀型別熵",
            "論文記號": r"$H_{\mathrm{macro}}^{(w)}$",
            "數值": _by_w_summary("H_macro_bits"),
            "論文語義摘要": "窗口宏觀型別標籤分布之熵（bits）。",
        },
        {
            "論文小節": "10.9.5（五）",
            "輸出參數": "單元持續時間（最大／平均）",
            "論文記號": r"$\tau_{\mathrm{unit}}^{(w)}$",
            "數值": "max: "
            + _by_w_summary("tau_unit_max", int_vals=False)
            + "｜mean: "
            + _by_w_summary("tau_unit_mean", int_vals=False),
            "論文語義摘要": "宏觀標籤運行長度之簡化統計（tau_unit_max／tau_unit_mean）。",
        },
        {
            "論文小節": "10.9.5（六）",
            "輸出參數": "窗口間分布差異（對 w=1）",
            "論文記號": r"$\mathrm{JS}_w$",
            "數值": _by_w_summary("JS_vs_w1_bits"),
            "論文語義摘要": r"宏觀型別分布與最細窗口 w=1 之 Jensen–Shannon（bits）；w=1 列為 —。",
        },
        {
            "論文小節": "10.9.5（七）",
            "輸出參數": "宏觀平台（C_eff 序列）",
            "論文記號": r"$L_{\mathrm{plat}}^{(w)}$",
            "數值": "最長平台長度: "
            + _by_w_summary("L_plat_max_len_Ceff", int_vals=True)
            + "｜平台段數: "
            + _by_w_summary("L_plat_num_seg_Ceff", int_vals=True),
            "論文語義摘要": "對 C_eff 代理序列以 ε_plat 做平台偵測之最大段長與段數。",
        },
    ]
    return pd.DataFrame(rows_out, columns=list(cols))
