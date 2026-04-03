"""
論文第十章表 10-2、10-3：域型約束梯子與解析簽名層次比較（批次指標）。

不依賴 Streamlit，供 CLI／網頁共用。
"""

from __future__ import annotations

import math
import random
from collections import Counter
from typing import Any, Dict, List, Sequence, Set, Tuple

from hypergraph_experiment.core import (
    SIGNATURES,
    VIOLATION_DEGREE_EXCESS,
    VIOLATION_DISCONNECTED,
    VIOLATION_EDGE_SIZE_BAD,
    VIOLATION_PAIR_TRIANGLE_FORBIDDEN,
    VIOLATION_TOO_MANY_EDGES,
    HypergraphConfig,
    all_possible_hyperedges,
    analyze_static,
    domain_constraint_violation_primary,
    filter_configs,
    powerset_limited,
    subsample_obs_configs,
)
from hypergraph_experiment.output_rounding import round_floats_for_output
from hypergraph_experiment.progress_types import ProgressCallback

# 主因統計順序（與 ``domain_constraint_violation_primary`` 檢查順序一致，§10.2.7 第二點）
PRIMARY_VIOLATION_CODES: Tuple[str, ...] = (
    VIOLATION_TOO_MANY_EDGES,
    VIOLATION_EDGE_SIZE_BAD,
    VIOLATION_DEGREE_EXCESS,
    VIOLATION_DISCONNECTED,
    VIOLATION_PAIR_TRIANGLE_FORBIDDEN,
)

# 批次扁平化時 ``level_i_<tail>`` 之尾碼（對應 ``build_ch10_column_name_map``）
TABLE_10_2_VIOL_BATCH_TAIL_ALL: Dict[str, str] = {
    VIOLATION_TOO_MANY_EDGES: "viol_all_too_many_edges",
    VIOLATION_EDGE_SIZE_BAD: "viol_all_edge_size_bad",
    VIOLATION_DEGREE_EXCESS: "viol_all_degree_excess",
    VIOLATION_DISCONNECTED: "viol_all_disconnected",
    VIOLATION_PAIR_TRIANGLE_FORBIDDEN: "viol_all_pair_triangle",
}
TABLE_10_2_VIOL_BATCH_TAIL_NEW: Dict[str, str] = {
    VIOLATION_TOO_MANY_EDGES: "viol_new_too_many_edges",
    VIOLATION_EDGE_SIZE_BAD: "viol_new_edge_size_bad",
    VIOLATION_DEGREE_EXCESS: "viol_new_degree_excess",
    VIOLATION_DISCONNECTED: "viol_new_disconnected",
    VIOLATION_PAIR_TRIANGLE_FORBIDDEN: "viol_new_pair_triangle",
}


def _violation_zh(code: str) -> str:
    """違規主因碼對應之短中文標籤（用於表 10-2 欄名）。"""
    return {
        VIOLATION_TOO_MANY_EDGES: "超邊數超限",
        VIOLATION_EDGE_SIZE_BAD: "超邊大小不符",
        VIOLATION_DEGREE_EXCESS: "度數超限",
        VIOLATION_DISCONNECTED: "連通未滿足",
        VIOLATION_PAIR_TRIANGLE_FORBIDDEN: "禁二元三角",
    }[code]


def _count_primary_violations(
    configs: Sequence[HypergraphConfig],
    *,
    max_edge_size: int,
    max_edges: int,
    connected_required: bool,
    max_degree: int | None,
    forbid_pair_triangles: bool,
) -> Counter[str]:
    """對候選序列逐筆標記首個未通過之域型條件，並計數。"""
    ctr: Counter[str] = Counter()
    for c in configs:
        r = domain_constraint_violation_primary(
            c,
            max_edge_size=max_edge_size,
            max_edges=max_edges,
            connected_required=connected_required,
            max_degree=max_degree,
            forbid_pair_triangles=forbid_pair_triangles,
        )
        if r is not None:
            ctr[r] += 1
    return ctr


def _violation_concentration(ctr: Counter[str]) -> float | None:
    """主因集中度：占比最高之主因類別所占比例（介於 0 與 1；無排除時為 None）。"""
    total = int(sum(ctr.values()))
    if total <= 0:
        return None
    return round(max(ctr.values()) / total, 6)


def exhaustive_candidate_count(
    n: int,
    max_edge_size: int,
    max_edges: int,
    min_edge_size: int = 2,
) -> Tuple[int, int]:
    """
    在 ``sample_limit=0`` 時，枚舉候選超圖之**精確**筆數與單層候選超邊數 M。

    Returns:
        (總候選數, M)，其中總候選數 = Σ_{r=0}^{min(m_max,M)} C(M,r)。
    """
    vertices = tuple(range(1, n + 1))
    m_edges = len(all_possible_hyperedges(vertices, min_edge_size, max_edge_size))
    max_r = min(max_edges, m_edges)
    total = sum(math.comb(m_edges, r) for r in range(0, max_r + 1))
    return total, m_edges


def _generate_candidates_for_domain_ladder(
    *,
    n: int,
    min_edge_size: int,
    max_edge_size: int,
    max_edges: int,
    sample_limit: int,
    seed: int,
    num_seeds: int,
) -> List[HypergraphConfig]:
    """依 k_min/k_max 與多種子設定產生候選配置（去重後）。"""
    vertices = tuple(range(1, n + 1))
    candidate_edges = all_possible_hyperedges(vertices, min_edge_size, max_edge_size)
    if sample_limit in (None, 0):
        configs: List[HypergraphConfig] = []
        for edge_subset in powerset_limited(candidate_edges, max_edges):
            configs.append(HypergraphConfig(vertices=vertices, hyperedges=frozenset(edge_subset)))
        return configs

    target = max(1, int(sample_limit))
    configs_set: Set[HypergraphConfig] = set()
    for i in range(max(1, num_seeds)):
        if len(configs_set) >= target:
            break
        rng = random.Random(seed + i * 1009)
        attempts = 0
        max_attempts = max(target * 20, 1000)
        while len(configs_set) < target and attempts < max_attempts:
            attempts += 1
            m = rng.randint(0, max_edges)
            chosen = rng.sample(candidate_edges, k=min(m, len(candidate_edges)))
            configs_set.add(HypergraphConfig(vertices=vertices, hyperedges=frozenset(chosen)))
    return list(configs_set)


def table_10_2_domain_ladder(
    *,
    n: int,
    min_edge_size: int,
    max_edge_size: int,
    max_edges: int,
    sample_limit: int,
    seed: int,
    num_seeds: int,
    levels: Sequence[Dict[str, Any]],
    progress: ProgressCallback = None,
) -> Tuple[List[Dict[str, Any]], int]:
    """
    表 10-2：固定候選產生方式，對嵌套域型約束逐層篩選，輸出 |Cfg| 與比例。

    Args:
        levels: 由粗到細之層級，每層鍵包含 ``label``、``connected``、``max_degree``、
            ``forbid_pair_triangles``（與 ``filter_configs`` 一致）。

    Returns:
        (列資料, |候選|)。
    """
    random.seed(seed)
    if progress:
        progress(0, max(3, len(levels) + 1), "表 10-2：產生候選…")
    candidates = _generate_candidates_for_domain_ladder(
        n=n,
        min_edge_size=min_edge_size,
        max_edge_size=max_edge_size,
        max_edges=max_edges,
        sample_limit=sample_limit,
        seed=seed,
        num_seeds=num_seeds,
    )
    n_cand = len(candidates)
    rows: List[Dict[str, Any]] = []
    prev_n: int | None = None
    prev_cfg_set: Set[HypergraphConfig] | None = None
    chain_subset_ok = True
    tot_lv = max(1, len(levels))
    for li, spec in enumerate(levels):
        if progress:
            progress(li + 1, tot_lv + 1, f"表 10-2：篩選層級 {li + 1}/{tot_lv}")
        cfg = filter_configs(
            candidates,
            max_edge_size=max_edge_size,
            max_edges=max_edges,
            connected_required=bool(spec["connected"]),
            max_degree=int(spec["max_degree"]) if spec.get("max_degree") is not None else None,
            forbid_pair_triangles=bool(spec["forbid_pair_triangles"]),
            progress=None,
        )
        cfg_set = set(cfg)
        sz = len(cfg_set)
        n_forbid_i = max(0, n_cand - sz)
        n_forbid_delta = max(0, (len(prev_cfg_set) - sz) if prev_cfg_set is not None else n_forbid_i)
        subset_prev = cfg_set.issubset(prev_cfg_set) if prev_cfg_set is not None else True
        chain_subset_ok = chain_subset_ok and subset_prev
        rho_i_to_prev = (sz / len(prev_cfg_set)) if prev_cfg_set else None
        candidate_ratio = (sz / n_cand) if n_cand else 0.0
        md = int(spec["max_degree"]) if spec.get("max_degree") is not None else None
        kwargs_dom = dict(
            max_edge_size=max_edge_size,
            max_edges=max_edges,
            connected_required=bool(spec["connected"]),
            max_degree=md,
            forbid_pair_triangles=bool(spec["forbid_pair_triangles"]),
        )
        excluded_all = [c for c in candidates if c not in cfg_set]
        ctr_all = _count_primary_violations(excluded_all, **kwargs_dom)
        if prev_cfg_set is not None:
            excluded_new = [c for c in prev_cfg_set if c not in cfg_set]
        else:
            excluded_new = excluded_all
        ctr_new = _count_primary_violations(excluded_new, **kwargs_dom)
        row: Dict[str, Any] = {
            "域型約束層級": spec.get("label", ""),
            "合法配置數": sz,
            "佔候選比例": candidate_ratio,
            "相對候選收縮率": 1.0 - candidate_ratio,
            "累計排除數": n_forbid_i,
            "本層新增排除數": n_forbid_delta,
            "是否為上一層子集": subset_prev,
            "鏈式子集成立": chain_subset_ok,
            "相對上一層保留率": rho_i_to_prev,
            "違規主因集中度": _violation_concentration(ctr_all),
            "本層新增違規主因集中度": _violation_concentration(ctr_new),
        }
        for code in PRIMARY_VIOLATION_CODES:
            zh = _violation_zh(code)
            row[f"該層排除主因{zh}筆數"] = int(ctr_all.get(code, 0))
            row[f"本層新增排除主因{zh}筆數"] = int(ctr_new.get(code, 0))
        if prev_n is not None and prev_n > 0:
            prev_ratio = sz / prev_n
            row["相對上一層保留率"] = prev_ratio
            row["相對上一層收縮率"] = 1.0 - prev_ratio
        prev_n = sz
        prev_cfg_set = cfg_set
        rows.append(row)
    if progress:
        progress(tot_lv + 1, tot_lv + 1, "表 10-2 完成")
    return round_floats_for_output(rows), n_cand


def table_10_3_signature_comparison(
    configs: Sequence[HypergraphConfig],
    delta: int,
    *,
    s_min: int = 0,
    n_cfg: int | None = None,
    seed: int = 7,
    progress: ProgressCallback = None,
) -> List[Dict[str, Any]]:
    """
    表 10-3：在與 §10.3 單次實驗相同之觀測集上，對 weak／medium／strong 比較解析指標。

    Args:
        configs: 可採用配置全集（域型篩選後）。
        delta: 簽名距離閾值 δ。
        s_min: 重疊率鄰域最小支持。
        n_cfg: 論文 :math:`N_{cfg}`；``None`` 時以全部 ``configs`` 為觀測集。
        seed: 與 ``run_full_experiment`` / ``subsample_obs_configs`` 一致之基底種子。
        progress: 可選進度回呼。
    """
    obs, _req, _n_act, _notice = subsample_obs_configs(list(configs), n_cfg, seed=seed)
    rows: List[Dict[str, Any]] = []
    names = sorted(SIGNATURES.keys())
    tn = max(1, len(names))
    for ni, name in enumerate(names):
        if progress:
            progress(ni, tn, f"表 10-3：簽名 {name}")
        a = analyze_static(obs, name, delta, s_min=s_min)
        u = a.get("compression_ratio_U")
        rows.append(
            {
                "解析簽名": name,
                "解析單元數": a["num_equivalence_classes"],
                "平均單元大小": round(float(a.get("avg_class_size", 0.0)), 6),
                "解析壓縮比": (float(u) if u is not None and u != math.inf else None),
                "重疊率": round(float(a.get("overlap_rate", 0.0)), 6),
                "相容孤立率": round(float(a.get("isol_rate_compat_graph", 0.0)), 6),
                "傳遞違反率": round(float(a.get("transitivity_violation_rate", 0.0)), 6),
                "解析熵位元": round(float(a.get("entropy_bits", 0.0)), 6),
            }
        )
    if progress:
        progress(tn, tn, "表 10-3 完成")
    return round_floats_for_output(rows)


TABLE_10_2_VIOLATION_COUNT_COLUMNS: List[Tuple[str, str, str]] = [
    (c, f"該層排除主因{_violation_zh(c)}筆數", f"本層新增排除主因{_violation_zh(c)}筆數")
    for c in PRIMARY_VIOLATION_CODES
]
