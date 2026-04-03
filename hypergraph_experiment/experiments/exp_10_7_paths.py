"""
§10.7 補強：同一合法配置域與更新核下，以兩組解析觀測 (簽名,δ) 追蹤終態類別，比較終端分布。
"""

from __future__ import annotations

import random
from collections import Counter
from typing import Any, Dict, List, Sequence, Tuple

from hypergraph_experiment.core import (
    HypergraphConfig,
    entropy,
    run_trajectory,
    sample_candidates_and_filter,
)
from hypergraph_experiment.output_rounding import round_floats_for_output
from hypergraph_experiment.progress_types import ProgressCallback
from hypergraph_experiment.refinement import js_divergence_bits, partition_analytic_units


def _terminal_histogram(
    configs: Sequence[HypergraphConfig],
    traj_ends: Sequence[HypergraphConfig],
    sig: str,
    delta: int,
) -> Tuple[List[float], Dict[HypergraphConfig, int]]:
    """終點配置在給定解析層上之類別索引直方圖（與配置域上該層分割對齊）。"""
    classes, idx = partition_analytic_units(configs, sig, delta)
    k = len(classes)
    hist = [0.0] * k
    for c in traj_ends:
        if c in idx:
            hist[idx[c]] += 1.0
    tot = sum(hist)
    if tot > 0:
        hist = [h / tot for h in hist]
    return hist, idx


def run_experiment_10_7(
    *,
    n: int = 5,
    max_edge_size: int = 3,
    max_edges: int = 4,
    sample_limit: int = 800,
    seed: int = 7,
    connected: bool = False,
    max_degree: int = 4,
    forbid_pair_triangles: bool = False,
    runs: int = 40,
    steps: int = 60,
    sig_path_a: str = "weak",
    delta_a: int = 2,
    sig_path_b: str = "strong",
    delta_b: int = 0,
    progress: ProgressCallback = None,
) -> dict[str, Any]:
    """
    雙解析路徑終端比較：軌道終點在兩種觀測層下之經驗類別分布，對齊維度後計算 JS 與熵差。

    Args:
        sig_path_a, delta_a: 第一條解析觀測層。
        sig_path_b, delta_b: 第二條解析觀測層。

    Returns:
        含 ``metrics.js_terminal_bits``、``H_terminal_a_bits``、``H_terminal_b_bits`` 等。
    """
    if progress:
        progress(0, max(1, runs), "§10.7 抽樣配置…")
    _, configs = sample_candidates_and_filter(
        n=n,
        max_edge_size=max_edge_size,
        max_edges=max_edges,
        sample_limit=sample_limit,
        seed=seed,
        connected=connected,
        max_degree=max_degree,
        forbid_pair_triangles=forbid_pair_triangles,
    )
    if not configs:
        return round_floats_for_output(
            {
                "experiment": "10.7",
                "error": "可採用配置為空，請放寬域型條件。",
                "parameters": {
                    "n": n,
                    "sample_limit": sample_limit,
                    "seed": seed,
                },
            }
        )

    cfg_set = set(configs)
    rng = random.Random(seed)
    ends: List[HypergraphConfig] = []
    for r in range(runs):
        if progress:
            progress(r + 1, runs, f"§10.7 軌道 {r + 1}/{runs}")
        start = rng.choice(configs)
        traj = run_trajectory(
            start,
            steps=steps,
            seed=seed + r + 11,
            max_edge_size=max_edge_size,
            max_edges=max_edges,
            connected_required=connected,
            max_degree=max_degree,
            forbid_pair_triangles=forbid_pair_triangles,
            allowed_configs=cfg_set,
        )
        ends.append(traj[-1])

    ha, _ = _terminal_histogram(configs, ends, sig_path_a, delta_a)
    hb, _ = _terminal_histogram(configs, ends, sig_path_b, delta_b)
    m = max(len(ha), len(hb))
    pa = ha + [0.0] * (m - len(ha))
    pb = hb + [0.0] * (m - len(hb))
    js = js_divergence_bits(pa, pb)
    e_a = entropy(ha, base=2.0) if ha else 0.0
    e_b = entropy(hb, base=2.0) if hb else 0.0

    agree = 0
    _, idx_a = partition_analytic_units(configs, sig_path_a, delta_a)
    _, idx_b = partition_analytic_units(configs, sig_path_b, delta_b)
    for c in ends:
        agree += 1 if idx_a.get(c) == idx_b.get(c) else 0

    return round_floats_for_output(
        {
            "experiment": "10.7",
            "parameters": {
                "n": n,
                "max_edge_size": max_edge_size,
                "max_edges": max_edges,
                "sample_limit": sample_limit,
                "seed": seed,
                "runs": runs,
                "steps": steps,
                "sig_path_a": sig_path_a,
                "delta_a": delta_a,
                "sig_path_b": sig_path_b,
                "delta_b": delta_b,
            },
            "metrics": {
                "js_terminal_bits": None if js is None else round(js, 8),
                "H_terminal_a_bits": round(e_a, 8),
                "H_terminal_b_bits": round(e_b, 8),
                "entropy_abs_diff_terminal": round(abs(e_a - e_b), 8),
                "terminal_class_agree_rate": round(agree / max(1, len(ends)), 8),
                "|S|_path_a": len(ha),
                "|S|_path_b": len(hb),
            },
        }
    )
