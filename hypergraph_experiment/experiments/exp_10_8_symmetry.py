"""
§10.8 對稱初態、微擾與局部型別擴張：N_type、H_type、Var(d)、ρ_type。
"""

from __future__ import annotations

import math
import random
import statistics
from collections import Counter, deque
from typing import Any, Dict, List, Sequence, Set, Tuple

from hypergraph_experiment.core import (
    SIGNATURES,
    HypergraphConfig,
    all_legal_successors,
    degree_sequence,
    satisfies_domain_constraints,
    two_section_adjacency,
)
from hypergraph_experiment.output_rounding import round_floats_for_output
from hypergraph_experiment.progress_types import ProgressCallback


def _make_symmetric_ring(n: int, m_target: int, *, max_edge_size: int, rng: random.Random) -> HypergraphConfig:
    """建近似對稱之二元環，再隨機補三元超邊至 m_target。"""
    verts = tuple(range(1, n + 1))
    edges: Set[frozenset] = set()
    for i in range(n):
        a = verts[i]
        b = verts[(i + 1) % n]
        edges.add(frozenset({a, b}))
    pool3 = [
        frozenset({verts[i], verts[(i + 1) % n], verts[(i + 2) % n]}) for i in range(n)
    ]
    while len(edges) < m_target:
        if pool3:
            edges.add(rng.choice(pool3))
        else:
            break
    return HypergraphConfig(vertices=verts, hyperedges=frozenset(edges))


def _perturb_edges(
    c: HypergraphConfig,
    eta: float,
    *,
    max_edge_size: int,
    max_edges: int,
    rng: random.Random,
) -> HypergraphConfig:
    """以比例 eta 嘗試重接／替換少許超邊（保持 |E| 近似）。"""
    el = list(c.hyperedges)
    if not el:
        return c
    k = max(1, int(round(len(el) * eta)))
    new_e = set(c.hyperedges)
    verts = list(c.vertices)
    for _ in range(k):
        if not new_e:
            break
        rem = rng.choice(list(new_e))
        new_e.remove(rem)
        if len(verts) >= 2:
            a, b = rng.sample(verts, 2)
            if max_edge_size >= 3 and rng.random() < 0.35:
                c3 = rng.choice(verts)
                new_e.add(frozenset({a, b, c3}))
            else:
                new_e.add(frozenset({a, b}))
    cand = HypergraphConfig(c.vertices, frozenset(new_e))
    if not satisfies_domain_constraints(
        cand,
        max_edge_size=max_edge_size,
        max_edges=max_edges,
        connected_required=False,
        max_degree=max(2, max_edges),
        forbid_pair_triangles=False,
    ):
        return c
    return cand


def _local_type_signature(
    c: HypergraphConfig,
    v: int,
    r: int,
) -> Tuple[Any, ...]:
    """以 2-section 圖上 r 步鄰域之度數多重集為局部型別（簡化版 σ_obs）。"""
    adj = two_section_adjacency(c)
    seen = {v}
    frontier = {v}
    for _ in range(r):
        nxt: Set[int] = set()
        for u in frontier:
            nxt |= adj[u]
        frontier = nxt - seen
        seen |= frontier
    degs = sorted(len(adj[u]) for u in seen)
    return tuple(degs)


def _iso_ball_canonical_signature(c: HypergraphConfig, v: int, r: int) -> Tuple[Any, ...]:
    """
    N_iso 之有限代理：以 r 步球上之誘導超圖規範簽名（度序列 + 邊型樣式）。

    與論文「球面同構類」不同：此處不判定完整同構，僅供比較對稱破缺之實務指標。
    """
    adj = two_section_adjacency(c)
    seen: Set[int] = {v}
    frontier: Set[int] = {v}
    for _ in range(r):
        nxt: Set[int] = set()
        for u in frontier:
            nxt |= adj[u]
        frontier = nxt - seen
        seen |= frontier
    ball = seen
    sub_edges: List[Tuple[int, ...]] = []
    for e in c.hyperedges:
        if e <= ball:
            sub_edges.append(tuple(sorted(e)))
    sub_pat = tuple(sorted(sub_edges, key=lambda x: (len(x), x)))
    deg_sub: Counter[int] = Counter()
    for e in c.hyperedges:
        if e <= ball:
            for u in e:
                deg_sub[u] += 1
    degs = tuple(sorted(deg_sub.get(u, 0) for u in sorted(ball)))
    return (degs, sub_pat)


def _mean_shortest_distance_2section(c: HypergraphConfig) -> Tuple[float, float]:
    """
    A_reach：2-section 上所有可達頂點對之無權最短路平均長度，以及可達對占全體無序對比例。

    Returns:
        (平均距離, 可達對比例)。若全無可達對則平均為 float('inf')。
    """
    adj = two_section_adjacency(c)
    verts = list(c.vertices)
    n = len(verts)
    if n < 2:
        return 0.0, 1.0
    total_dist = 0
    n_pairs = 0
    for i, s in enumerate(verts):
        dist: Dict[int, int] = {s: 0}
        q: deque[int] = deque([s])
        while q:
            u = q.popleft()
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    q.append(w)
        for t in verts[i + 1 :]:
            if t in dist:
                total_dist += dist[t]
                n_pairs += 1
    denom = n * (n - 1) // 2
    frac = n_pairs / denom if denom else 0.0
    mean_d = total_dist / n_pairs if n_pairs else float("inf")
    return mean_d, frac


def _sample_random_legal(
    n: int,
    m: int,
    *,
    max_edge_size: int,
    seed: int,
    rng: random.Random,
) -> HypergraphConfig:
    """隨機合法構型（簡便採用對稱環再強擾動）。"""
    c0 = _make_symmetric_ring(n, m, max_edge_size=max_edge_size, rng=rng)
    return _perturb_edges(c0, 0.6, max_edge_size=max_edge_size, max_edges=m, rng=rng)


def run_experiment_10_8(
    *,
    n: int = 12,
    m: int = 18,
    max_edge_size: int = 3,
    init_family: str = "sym",
    eta: float = 0.1,
    T_sb: int = 20,
    r: int = 2,
    sig_obs: str = "medium",
    n_samples: int = 15,
    seed: int = 0,
    dynamics_steps: int = 0,
    progress: ProgressCallback | None = None,
) -> dict[str, Any]:
    """
    比較對稱初態 / 微擾 / 隨機構型之局部型別數與型別熵。

    Args:
        init_family: ``sym`` | ``pert`` | ``rand``。
        eta, T_sb: 微擾比例與（可選）後續動力學步數占位；目前 dynamcis_steps 由參數 dynamics_steps 控制額外演化。
        r: 局部鄰域深度。
        sig_obs: 用於整體簽名對照之 weak/medium/strong。

    Returns:
        聚合指標與每樣本之列資料 ``per_series``。
    """
    rng = random.Random(seed)
    series: List[dict[str, Any]] = []

    for si in range(n_samples):
        if progress:
            progress(si + 1, n_samples, f"§10.8 樣本 {si + 1}/{n_samples}")
        if init_family == "sym":
            c = _make_symmetric_ring(n, m, max_edge_size=max_edge_size, rng=rng)
        elif init_family == "pert":
            c0 = _make_symmetric_ring(n, m, max_edge_size=max_edge_size, rng=rng)
            c = _perturb_edges(c0, eta, max_edge_size=max_edge_size, max_edges=m, rng=rng)
        else:
            c = _sample_random_legal(n, m, max_edge_size=max_edge_size, seed=seed + si, rng=rng)

        for _ in range(dynamics_steps):
            succ = all_legal_successors(
                c,
                max_edge_size=max_edge_size,
                max_edges=len(c.hyperedges) + 2,
                connected_required=False,
                max_degree=None,
                forbid_pair_triangles=False,
            )
            if succ:
                c = rng.choice(succ)

        lt: Counter[Tuple[Any, ...]] = Counter()
        iso_sigs: Set[Tuple[Any, ...]] = set()
        for v in c.vertices:
            lt[_local_type_signature(c, v, r)] += 1
            iso_sigs.add(_iso_ball_canonical_signature(c, v, r))
        n_type = len(lt)
        n_iso = len(iso_sigs)
        totv = sum(lt.values())
        probs = [cnt / totv for cnt in lt.values()]
        h_type = -sum(p * math.log(p, 2) for p in probs if p > 0)
        deg = degree_sequence(c)
        var_d = statistics.pvariance(deg) if len(deg) > 1 else 0.0
        sig_global = SIGNATURES[sig_obs](c) if sig_obs in SIGNATURES else ()
        a_d, a_frac = _mean_shortest_distance_2section(c)
        a_reach = float(a_d) if math.isfinite(a_d) else None

        series.append(
            {
                "sample": si,
                "N_type": n_type,
                "N_iso": n_iso,
                "H_type_bits": round(h_type, 6),
                "Var_degree": round(var_d, 6),
                "A_reach_mean_dist": None if a_reach is None else round(a_reach, 6),
                "A_reach_pair_fraction": round(a_frac, 6),
                "sig_obs_repr": str(sig_global)[:200],
            }
        )

    n_sym = statistics.mean(s["N_type"] for s in series) if series else 0.0
    n_iso_m = statistics.mean(s["N_iso"] for s in series) if series else 0.0
    a_ok = [s["A_reach_mean_dist"] for s in series if s.get("A_reach_mean_dist") is not None]
    return round_floats_for_output(
        {
            "experiment": "10.8",
            "parameters": {
                "n": n,
                "m": m,
                "init_family": init_family,
                "eta": eta,
                "T_sb": T_sb,
                "r": r,
                "sig_obs": sig_obs,
                "n_samples": n_samples,
                "dynamics_steps": dynamics_steps,
                "seed": seed,
            },
            "metrics": {
                "N_type_mean": round(n_sym, 6),
                "N_iso_mean": round(n_iso_m, 6),
                "H_type_mean": round(statistics.mean(s["H_type_bits"] for s in series), 6)
                if series
                else 0.0,
                "Var_degree_mean": round(statistics.mean(s["Var_degree"] for s in series), 6)
                if series
                else 0.0,
                "A_reach_mean": round(statistics.mean(a_ok), 6) if a_ok else None,
                "A_reach_pair_fraction_mean": round(
                    statistics.mean(s["A_reach_pair_fraction"] for s in series), 6
                )
                if series
                else None,
            },
            "per_sample": series,
        }
    )


def rho_type_expansion(sym_mean: float, pert_mean: float, eps: float = 1e-6) -> float:
    """ρ_type = (N_pert - N_sym) / (N_sym + ε)。"""
    return (pert_mean - sym_mean) / (sym_mean + eps)


# 三臂並跑時各臂種子間隔，避免 RNG 流重疊（與樣本內迴圈無關之顯式分離）。
_THREE_ARM_SEED_STRIDE: int = 10_000_000


def run_experiment_10_8_three_arm(
    *,
    n: int = 12,
    m: int = 18,
    max_edge_size: int = 3,
    eta: float = 0.1,
    T_sb: int = 20,
    r: int = 2,
    sig_obs: str = "medium",
    n_samples: int = 15,
    seed: int = 0,
    dynamics_steps: int = 0,
    progress: ProgressCallback | None = None,
) -> dict[str, Any]:
    """
    在相同超參數下連續執行 sym、pert、rand 三類初態族，便於並列對照。

    論文記號 **ρ_type** 仍僅定義為 pert 相對 sym 之型別擴張；rand 另以
    ``rho_type_rand_over_sym``（同公式型態之輔助量）呈現，不取代 ρ_type 語意。

    Args:
        n, m, max_edge_size: 超圖規模與邊大小上限。
        eta: 微擾強度（pert 臂；rand 臂內部仍依既有採樣邏輯）。
        T_sb: 與單臂實驗一致之占位參數（寫入 parameters）。
        r: 局部鄰域深度。
        sig_obs: 整體觀測簽名 weak／medium／strong。
        n_samples: 每臂重複樣本數。
        seed: 基底種子；三臂分別使用 ``seed``、``seed + stride``、``seed + 2*stride``。
        dynamics_steps: 每臂演化步數。
        progress: 可選進度回呼，總進度為 ``3 * n_samples`` 步。

    Returns:
        含 ``by_family``（三份單臂結果）、``comparison``（並列均值與 ρ）、
        以及共用 ``parameters`` 之字典。
    """
    families: Tuple[str, ...] = ("sym", "pert", "rand")
    by_family: dict[str, dict[str, Any]] = {}

    for arm_idx, fam in enumerate(families):
        arm_seed = seed + arm_idx * _THREE_ARM_SEED_STRIDE

        def _wrap_progress(
            cur: int,
            total: int,
            msg: str,
            *,
            _ai: int = arm_idx,
            _fam: str = fam,
        ) -> None:
            if progress is None:
                return
            base = _ai * n_samples
            progress(base + cur, max(1, 3 * n_samples), f"§10.8 [{_fam}] {msg}")

        by_family[fam] = run_experiment_10_8(
            n=n,
            m=m,
            max_edge_size=max_edge_size,
            init_family=fam,
            eta=eta,
            T_sb=T_sb,
            r=r,
            sig_obs=sig_obs,
            n_samples=n_samples,
            seed=arm_seed,
            dynamics_steps=dynamics_steps,
            progress=_wrap_progress if progress else None,
        )

    m_sym = by_family["sym"]["metrics"]
    m_pert = by_family["pert"]["metrics"]
    m_rand = by_family["rand"]["metrics"]
    n_sym = float(m_sym.get("N_type_mean", 0.0))
    n_pert = float(m_pert.get("N_type_mean", 0.0))
    n_rand = float(m_rand.get("N_type_mean", 0.0))
    rho_type = rho_type_expansion(n_sym, n_pert)
    rho_rand_over_sym = rho_type_expansion(n_sym, n_rand)

    comparison: dict[str, Any] = {
        "N_type_mean_sym": m_sym.get("N_type_mean"),
        "N_type_mean_pert": m_pert.get("N_type_mean"),
        "N_type_mean_rand": m_rand.get("N_type_mean"),
        "H_type_mean_sym": m_sym.get("H_type_mean"),
        "H_type_mean_pert": m_pert.get("H_type_mean"),
        "H_type_mean_rand": m_rand.get("H_type_mean"),
        "N_iso_mean_sym": m_sym.get("N_iso_mean"),
        "N_iso_mean_pert": m_pert.get("N_iso_mean"),
        "N_iso_mean_rand": m_rand.get("N_iso_mean"),
        "Var_degree_mean_sym": m_sym.get("Var_degree_mean"),
        "Var_degree_mean_pert": m_pert.get("Var_degree_mean"),
        "Var_degree_mean_rand": m_rand.get("Var_degree_mean"),
        "A_reach_mean_sym": m_sym.get("A_reach_mean"),
        "A_reach_mean_pert": m_pert.get("A_reach_mean"),
        "A_reach_mean_rand": m_rand.get("A_reach_mean"),
        "rho_type_pert_over_sym": round(rho_type, 6),
        "rho_type_rand_over_sym": round(rho_rand_over_sym, 6),
    }

    return round_floats_for_output(
        {
            "experiment": "10.8_three_arm",
            "parameters": {
                "n": n,
                "m": m,
                "max_edge_size": max_edge_size,
                "eta": eta,
                "T_sb": T_sb,
                "r": r,
                "sig_obs": sig_obs,
                "n_samples": n_samples,
                "dynamics_steps": dynamics_steps,
                "seed_base": seed,
                "seed_stride": _THREE_ARM_SEED_STRIDE,
                "seeds_by_family": {
                    "sym": seed,
                    "pert": seed + _THREE_ARM_SEED_STRIDE,
                    "rand": seed + 2 * _THREE_ARM_SEED_STRIDE,
                },
            },
            "comparison": comparison,
            "by_family": by_family,
        }
    )


def section_10_8_output_parameters_df(result: dict[str, Any]) -> "pd.DataFrame":
    r"""
    組裝《約束世界論 30》§10.8.5「輸出參數」之主表（與扁平化寬表並列，供匯出／複製）。

    列舉（一）$N_{\mathrm{type}}$、（二）$H_{\mathrm{type}}$、（三）$N_{\mathrm{iso}}^{(r)}$、
    （四）$\mathrm{Var}(d)$、（五）$A_{\mathrm{reach}}$，並附（六）$\rho_{\mathrm{type}}$ 的計算提示。

    Note:
        單臂結果中（六）$\rho_{\mathrm{type}}$ 僅提示公式；三臂結果（``experiment=="10.8_three_arm"``）
        則直接填入 pert∥sym 之數值，並另列 rand∥sym 輔助列。

    Args:
        result: ``run_experiment_10_8`` 或 ``run_experiment_10_8_three_arm`` 之回傳字典。

    Returns:
        五欄 ``DataFrame``：論文小節、輸出參數、論文記號、數值、論文語義摘要。
    """
    import pandas as pd

    cols = ("論文小節", "輸出參數", "論文記號", "數值", "論文語義摘要")

    def _fmt(v: object) -> str:
        if v is None:
            return "—"
        if isinstance(v, float):
            return f"{v:.6g}"
        return str(v)

    # 三臂並跑：以 comparison 組表，ρ_type 列為論文定義之 pert／sym。
    if result.get("experiment") == "10.8_three_arm":
        comp = result.get("comparison")
        if not isinstance(comp, dict):
            return pd.DataFrame(
                [
                    {
                        "論文小節": "—",
                        "輸出參數": "無法組表",
                        "論文記號": "—",
                        "數值": "—",
                        "論文語義摘要": "三臂結果缺少 comparison 字典。",
                    }
                ],
                columns=list(cols),
            )

        def _trip(k_sym: str, k_pert: str, k_rand: str) -> str:
            a, b, c = comp.get(k_sym), comp.get(k_pert), comp.get(k_rand)
            return f"sym={_fmt(a)}；pert={_fmt(b)}；rand={_fmt(c)}"

        rho_p = comp.get("rho_type_pert_over_sym")
        rho_r = comp.get("rho_type_rand_over_sym")

        rows_3: List[dict[str, str]] = [
            {
                "論文小節": "10.8.5（一）",
                "輸出參數": "局部型別數（均值；三臂）",
                "論文記號": r"$N_{\mathrm{type}}$",
                "數值": _trip("N_type_mean_sym", "N_type_mean_pert", "N_type_mean_rand"),
                "論文語義摘要": "同參數下 sym／pert／rand 三類初態族各跑 n_samples 後之 N_type 均值並列。",
            },
            {
                "論文小節": "10.8.5（二）",
                "輸出參數": "型別熵（均值；三臂；bits）",
                "論文記號": r"$H_{\mathrm{type}}$",
                "數值": _trip("H_type_mean_sym", "H_type_mean_pert", "H_type_mean_rand"),
                "論文語義摘要": "三臂並列之型別熵均值。",
            },
            {
                "論文小節": "10.8.5（三）",
                "輸出參數": "鄰域同構類數（均值；三臂；代理）",
                "論文記號": r"$N_{\mathrm{iso}}^{(r)}$",
                "數值": _trip("N_iso_mean_sym", "N_iso_mean_pert", "N_iso_mean_rand"),
                "論文語義摘要": "三臂並列之 N_iso 代理指標。",
            },
            {
                "論文小節": "10.8.5（四）",
                "輸出參數": "度數方差（均值；三臂）",
                "論文記號": r"$\mathrm{Var}(d)$",
                "數值": _trip(
                    "Var_degree_mean_sym",
                    "Var_degree_mean_pert",
                    "Var_degree_mean_rand",
                ),
                "論文語義摘要": "三臂並列之度數方差均值。",
            },
            {
                "論文小節": "10.8.5（五）",
                "輸出參數": "平均可達性（三臂）",
                "論文記號": r"$A_{\mathrm{reach}}$",
                "數值": _trip(
                    "A_reach_mean_sym",
                    "A_reach_mean_pert",
                    "A_reach_mean_rand",
                ),
                "論文語義摘要": "三臂並列之 A_reach；若某臂無可達對則可能為 —。",
            },
            {
                "論文小節": "10.8.5（六）",
                "輸出參數": "型別擴張率（論文 ρ_type）",
                "論文記號": r"$\rho_{\mathrm{type}}$",
                "數值": _fmt(rho_p),
                "論文語義摘要": r"$\rho_{\mathrm{type}}=(N_{\mathrm{type}}^{(\mathrm{pert})}-N_{\mathrm{type}}^{(\mathrm{sym})})/(N_{\mathrm{type}}^{(\mathrm{sym})}+\varepsilon)$；由上列三臂之 sym／pert 均值代入。",
            },
            {
                "論文小節": "10.8.5（輔助）",
                "輸出參數": "rand 相對 sym 之型別擴張（非論文 ρ_type）",
                "論文記號": r"$\rho_{\mathrm{rand}/\mathrm{sym}}$（輔助）",
                "數值": _fmt(rho_r),
                "論文語義摘要": "與 ρ_type 同形之對照量：(N_type_rand−N_type_sym)/(N_type_sym+ε)；論文主線仍以 pert∥sym 為 ρ_type。",
            },
        ]
        return pd.DataFrame(rows_3, columns=list(cols))

    metrics = result.get("metrics")
    if not isinstance(metrics, dict):
        return pd.DataFrame(
            [
                {
                    "論文小節": "—",
                    "輸出參數": "無法組表",
                    "論文記號": "—",
                    "數值": "—",
                    "論文語義摘要": "結果缺少 metrics 字典。",
                }
            ],
            columns=list(cols),
        )
    if metrics.get("error"):
        return pd.DataFrame(
            [
                {
                    "論文小節": "—",
                    "輸出參數": "執行狀態",
                    "論文記號": "—",
                    "數值": "—",
                    "論文語義摘要": str(metrics["error"]),
                }
            ],
            columns=list(cols),
        )

    n_type = metrics.get("N_type_mean")
    h_type = metrics.get("H_type_mean")
    n_iso = metrics.get("N_iso_mean")
    var_d = metrics.get("Var_degree_mean")
    a_reach = metrics.get("A_reach_mean")

    rows: List[dict[str, str]] = [
        {
            "論文小節": "10.8.5（一）",
            "輸出參數": "局部型別數（均值）",
            "論文記號": r"$N_{\mathrm{type}}$",
            "數值": _fmt(n_type),
            "論文語義摘要": "在固定 r 與 σ_obs 下，節點局部鄰域型別之總數；此處為樣本平均。",
        },
        {
            "論文小節": "10.8.5（二）",
            "輸出參數": "型別熵（均值；bits）",
            "論文記號": r"$H_{\mathrm{type}}$",
            "數值": _fmt(h_type),
            "論文語義摘要": "局部型別分布之熵；越大表示可分辨性越高。",
        },
        {
            "論文小節": "10.8.5（三）",
            "輸出參數": "鄰域同構類數（均值；代理）",
            "論文記號": r"$N_{\mathrm{iso}}^{(r)}$",
            "數值": _fmt(n_iso),
            "論文語義摘要": "r 步球之同構類數；本實作採規範簽名代理，非完整同構判定。",
        },
        {
            "論文小節": "10.8.5（四）",
            "輸出參數": "度數方差（均值）",
            "論文記號": r"$\mathrm{Var}(d)$",
            "數值": _fmt(var_d),
            "論文語義摘要": "節點度數分布之方差；對稱初態通常較低。",
        },
        {
            "論文小節": "10.8.5（五）",
            "輸出參數": "平均可達性（2-section 平均最短距離）",
            "論文記號": r"$A_{\mathrm{reach}}$",
            "數值": _fmt(a_reach),
            "論文語義摘要": "以 2-section 上可達頂點對之平均距離操作化；越小表示越緊密可達。",
        },
        {
            "論文小節": "10.8.5（六）",
            "輸出參數": "型別擴張率（提示）",
            "論文記號": r"$\rho_{\mathrm{type}}$",
            "數值": "需 sym 與 pert 兩次實驗",
            "論文語義摘要": r"$\rho_{\mathrm{type}}=(N_{\mathrm{type}}^{(\mathrm{pert})}-N_{\mathrm{type}}^{(\mathrm{sym})})/(N_{\mathrm{type}}^{(\mathrm{sym})}+\varepsilon)$；請分別執行 init_family=sym 與 pert 後再計算。",
        },
    ]
    return pd.DataFrame(rows, columns=list(cols))
