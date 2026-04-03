#!/usr/bin/env python3
"""
第十章有限超圖 toy model：域型約束、解析鄰域、等價穩定化、靜態指標與簡易動力學。

模組刻意以標準函式庫為主，供 CLI 與網頁介面共用；圖表繪製請於上層以 matplotlib 處理。
"""

from __future__ import annotations

import json
import math
import random
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from typing import Any, Callable, Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from hypergraph_experiment.output_rounding import round_floats_for_output
from hypergraph_experiment.progress_types import ProgressCallback
from hypergraph_experiment.time_series_metrics import (
    plateau_length_max_abs_diff,
    state_change_fraction_traj,
)

Vertex = int
Hyperedge = FrozenSet[Vertex]


@dataclass(frozen=True)
class HypergraphConfig:
    """有限超圖配置 c = (V, E)，其中每條超邊為頂點之非空子集。"""

    vertices: Tuple[Vertex, ...]
    hyperedges: FrozenSet[Hyperedge]

    def to_dict(self) -> dict[str, Any]:
        """轉成可 JSON 序列化之字典（頂點與超邊皆以排序後列表表示）。"""
        return {
            "vertices": list(self.vertices),
            "hyperedges": [
                sorted(e) for e in sorted(self.hyperedges, key=lambda x: (len(x), tuple(sorted(x))))
            ],
        }


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def all_possible_hyperedges(vertices: Sequence[Vertex], min_size: int, max_size: int) -> List[Hyperedge]:
    edges: List[Hyperedge] = []
    for k in range(min_size, max_size + 1):
        edges.extend(frozenset(c) for c in combinations(vertices, k))
    return edges


def powerset_limited(items: Sequence[Hyperedge], max_subset_size: int) -> Iterable[Tuple[Hyperedge, ...]]:
    for r in range(0, max_subset_size + 1):
        yield from combinations(items, r)


def two_section_adjacency(c: HypergraphConfig) -> Dict[Vertex, Set[Vertex]]:
    adj: Dict[Vertex, Set[Vertex]] = {v: set() for v in c.vertices}
    for e in c.hyperedges:
        for a in e:
            for b in e:
                if a != b:
                    adj[a].add(b)
    return adj


def is_connected_2section(c: HypergraphConfig) -> bool:
    if not c.vertices:
        return True
    adj = two_section_adjacency(c)
    start = c.vertices[0]
    stack = [start]
    seen: Set[Vertex] = set()
    while stack:
        x = stack.pop()
        if x in seen:
            continue
        seen.add(x)
        stack.extend(adj[x] - seen)
    return len(seen) == len(c.vertices)


def degree_sequence(c: HypergraphConfig) -> Tuple[int, ...]:
    deg = Counter()
    for e in c.hyperedges:
        for v in e:
            deg[v] += 1
    return tuple(sorted(deg.get(v, 0) for v in c.vertices))


def edge_size_multiset(c: HypergraphConfig) -> Tuple[int, ...]:
    return tuple(sorted(len(e) for e in c.hyperedges))


def pair_adjacency_signature(c: HypergraphConfig) -> Tuple[int, ...]:
    # 2-section 圖之成對鄰接編碼（依頂點順序）
    adj = two_section_adjacency(c)
    vals: List[int] = []
    verts = list(c.vertices)
    for i, a in enumerate(verts):
        for b in verts[i + 1 :]:
            vals.append(1 if b in adj[a] else 0)
    return tuple(vals)


def motif_counts(c: HypergraphConfig) -> Tuple[int, int]:
    # (二元超邊數, 三元超邊數)
    n2 = sum(1 for e in c.hyperedges if len(e) == 2)
    n3 = sum(1 for e in c.hyperedges if len(e) == 3)
    return (n2, n3)


def canonical_edge_pattern(c: HypergraphConfig) -> Tuple[Tuple[int, ...], ...]:
    return tuple(sorted((tuple(sorted(e)) for e in c.hyperedges), key=lambda x: (len(x), x)))


# -----------------------------------------------------------------------------
# Domain constraints
# -----------------------------------------------------------------------------


def forbidden_triangle_of_pairs(c: HypergraphConfig) -> bool:
    """若存在僅由 2-超邊構成的三角形則回傳 True。"""
    pair_edges = {tuple(sorted(e)) for e in c.hyperedges if len(e) == 2}
    verts = list(c.vertices)
    for a, b, d in combinations(verts, 3):
        if (a, b) in pair_edges and (a, d) in pair_edges and (b, d) in pair_edges:
            return True
    return False


VIOLATION_TOO_MANY_EDGES = "too_many_edges"
VIOLATION_EDGE_SIZE_BAD = "edge_size_bad"
VIOLATION_DEGREE_EXCESS = "degree_excess"
VIOLATION_DISCONNECTED = "disconnected"
VIOLATION_PAIR_TRIANGLE_FORBIDDEN = "pair_triangle_forbidden"


def domain_constraint_violation_primary(
    c: HypergraphConfig,
    *,
    max_edge_size: int,
    max_edges: int,
    connected_required: bool,
    max_degree: int | None,
    forbid_pair_triangles: bool,
) -> str | None:
    """
    若配置未通過 ``satisfies_domain_constraints``，回傳**第一個**未通過條件對應之原因碼；
    通過則回傳 ``None``。順序與 ``satisfies_domain_constraints`` 一致，供 §10.2.7「違規主因」統計。

    Returns:
        以下之一：超邊數、超邊大小、度數、二階連通、禁二元三角；通過時為 ``None``。
    """
    if len(c.hyperedges) > max_edges:
        return VIOLATION_TOO_MANY_EDGES

    deg: Counter[Vertex] = Counter()
    for e in c.hyperedges:
        if len(e) == 0 or len(e) > max_edge_size:
            return VIOLATION_EDGE_SIZE_BAD
        for v in e:
            deg[v] += 1
            if max_degree is not None and deg[v] > max_degree:
                return VIOLATION_DEGREE_EXCESS

    if connected_required and not is_connected_2section(c):
        return VIOLATION_DISCONNECTED

    if forbid_pair_triangles and forbidden_triangle_of_pairs(c):
        return VIOLATION_PAIR_TRIANGLE_FORBIDDEN

    return None


def satisfies_domain_constraints(
    c: HypergraphConfig,
    *,
    max_edge_size: int,
    max_edges: int,
    connected_required: bool,
    max_degree: int | None,
    forbid_pair_triangles: bool,
) -> bool:
    """檢查配置是否滿足給定之域型約束 Λ_dom。"""
    return (
        domain_constraint_violation_primary(
            c,
            max_edge_size=max_edge_size,
            max_edges=max_edges,
            connected_required=connected_required,
            max_degree=max_degree,
            forbid_pair_triangles=forbid_pair_triangles,
        )
        is None
    )


# -----------------------------------------------------------------------------
# Signatures / tolerance / equivalence
# -----------------------------------------------------------------------------


def signature_weak(c: HypergraphConfig):
    return edge_size_multiset(c)


def signature_medium(c: HypergraphConfig):
    return (edge_size_multiset(c), degree_sequence(c))


def signature_strong(c: HypergraphConfig):
    return (edge_size_multiset(c), degree_sequence(c), pair_adjacency_signature(c), motif_counts(c))


SIGNATURES: Dict[str, Callable[[HypergraphConfig], object]] = {
    "weak": signature_weak,
    "medium": signature_medium,
    "strong": signature_strong,
}


def flatten_numeric_structure(x: object) -> List[int]:
    out: List[int] = []
    if isinstance(x, int):
        out.append(x)
    elif isinstance(x, (list, tuple)):
        for item in x:
            out.extend(flatten_numeric_structure(item))
    else:
        out.append(abs(hash(repr(x))) % 10_000)
    return out


def signature_distance(sig1: object, sig2: object) -> int:
    a = flatten_numeric_structure(sig1)
    b = flatten_numeric_structure(sig2)
    m = max(len(a), len(b))
    a = a + [0] * (m - len(a))
    b = b + [0] * (m - len(b))
    return sum(abs(x - y) for x, y in zip(a, b))


def approx(
    c1: HypergraphConfig,
    c2: HypergraphConfig,
    signature_fn: Callable[[HypergraphConfig], object],
    delta: int,
) -> bool:
    return signature_distance(signature_fn(c1), signature_fn(c2)) <= delta


class UnionFind:
    """並查集：將相容鄰域之聯通閉包凝聚為等價類。"""

    def __init__(self, items: Iterable[HypergraphConfig]) -> None:
        self.parent = {x: x for x in items}

    def find(self, x: HypergraphConfig) -> HypergraphConfig:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: HypergraphConfig, b: HypergraphConfig) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def _padded_flat_signature_vectors(
    cfg_list: Sequence[HypergraphConfig],
    signature_fn: Callable[[HypergraphConfig], object],
) -> tuple[Dict[HypergraphConfig, List[int]], int]:
    """
    將各配置之簽名展開為等長之整數向量（右側零填充），供 L1 距離與 δ=0 分桶共用。

    Returns:
        (padded, M)：每配置對應 padding 後向量、全域長度 M。
    """
    sigs = {c: signature_fn(c) for c in cfg_list}
    flats = {c: flatten_numeric_structure(sigs[c]) for c in cfg_list}
    M = max((len(flats[c]) for c in cfg_list), default=0)
    padded: Dict[HypergraphConfig, List[int]] = {
        c: flats[c] + [0] * (M - len(flats[c])) for c in cfg_list
    }
    return padded, M


def _l1_padded_leq_delta(
    pc: List[int], pd: List[int], M: int, delta: int
) -> bool:
    """兩個已等長 padding 向量之 L1 距離是否 ≤ δ（超過時提前結束）。"""
    s = 0
    for i in range(M):
        s += abs(pc[i] - pd[i])
        if s > delta:
            return False
    return True


def tolerance_classes_and_neighborhoods(
    configs: Sequence[HypergraphConfig],
    signature_fn: Callable[[HypergraphConfig], object],
    delta: int,
) -> Tuple[List[Set[HypergraphConfig]], Dict[HypergraphConfig, Set[HypergraphConfig]]]:
    """
    建立 δ-鄰域表與對應之傳遞閉包等價類（與舊版 ``build_tolerance_neighborhoods`` +
    ``stabilize_to_equivalence`` 語義一致）。

    複雜度：
        - **δ=0**：O(n·L)（分桶；L 為展平簽名長度），避免 O(N²) 全對全比較。
        - **δ>0**：仍為 O(N²·L)（需列舉可相容對），但不再對每一對重複展開簽名。

    Args:
        configs: 觀測或候選配置序列。
        signature_fn: 簽名映射。
        delta: 整數閾值（與 ``signature_distance`` 之尺度一致）。

    Returns:
        (classes, neighborhoods)；``neighborhoods[c]=T(c)``。
    """
    cfg_list = list(configs)
    if not cfg_list:
        return [], {}
    padded, M = _padded_flat_signature_vectors(cfg_list, signature_fn)
    neighborhoods: Dict[HypergraphConfig, Set[HypergraphConfig]]

    if delta == 0:
        buckets: Dict[Tuple[int, ...], Set[HypergraphConfig]] = defaultdict(set)
        for c in cfg_list:
            buckets[tuple(padded[c])].add(c)
        neighborhoods = {c: buckets[tuple(padded[c])] for c in cfg_list}
        classes = list(buckets.values())
        return classes, neighborhoods

    neighborhoods = {}
    for c in cfg_list:
        pc = padded[c]
        neigh: Set[HypergraphConfig] = set()
        for d in cfg_list:
            if _l1_padded_leq_delta(pc, padded[d], M, delta):
                neigh.add(d)
        neighborhoods[c] = neigh

    uf = UnionFind(cfg_list)
    for a, neigh in neighborhoods.items():
        for b in neigh:
            uf.union(a, b)
    class_map: Dict[HypergraphConfig, Set[HypergraphConfig]] = defaultdict(set)
    for c in cfg_list:
        class_map[uf.find(c)].add(c)
    classes = list(class_map.values())
    return classes, neighborhoods


def tolerance_equivalence_classes(
    configs: Sequence[HypergraphConfig],
    signature_fn: Callable[[HypergraphConfig], object],
    delta: int,
) -> List[Set[HypergraphConfig]]:
    """
    僅回傳等價類列表；內部共用 ``tolerance_classes_and_neighborhoods`` 之最佳化（δ=0 為 O(n) 級）。

    Args:
        configs: 配置序列。
        signature_fn: 簽名映射。
        delta: 距離閾值。

    Returns:
        非空之解析單元（集合）列表。
    """
    cls, _ = tolerance_classes_and_neighborhoods(configs, signature_fn, delta)
    return cls


def build_tolerance_neighborhoods(
    configs: Sequence[HypergraphConfig],
    signature_fn: Callable[[HypergraphConfig], object],
    delta: int,
) -> Dict[HypergraphConfig, Set[HypergraphConfig]]:
    _, nh = tolerance_classes_and_neighborhoods(configs, signature_fn, delta)
    return nh


def stabilize_to_equivalence(
    neighborhoods: Dict[HypergraphConfig, Set[HypergraphConfig]],
) -> List[Set[HypergraphConfig]]:
    """
    由 δ-鄰域無向圖求聯通分量（並查集）。

    若鄰域表係以 **δ=0** 之分桶建構（同一輻條之多個節點共用同一 ``set`` 物件），
    則直接還原唯一集合物件為等價類，避免 O(N²) 次 ``union``。
    """
    if not neighborhoods:
        return []
    n = len(neighborhoods)
    unique_ids = len({id(s) for s in neighborhoods.values()})
    if unique_ids < n:
        return list({id(s): s for s in neighborhoods.values()}.values())
    uf = UnionFind(neighborhoods.keys())
    for a, neigh in neighborhoods.items():
        for b in neigh:
            uf.union(a, b)
    classes: Dict[HypergraphConfig, Set[HypergraphConfig]] = defaultdict(set)
    for c in neighborhoods.keys():
        classes[uf.find(c)].add(c)
    return list(classes.values())


# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------


def overlap_rate(
    neighborhoods: Dict[HypergraphConfig, Set[HypergraphConfig]],
    *,
    s_min: int = 0,
) -> float:
    """
    鄰域對之 Jaccard 風格重疊率平均。

    Args:
        neighborhoods: 各配置之 δ-鄰域 T(c)。
        s_min: 僅當 |T(c)| >= s_min 之配置參與配對平均；0 表示不過濾（向後相容）。
    """
    keys = [k for k in neighborhoods.keys() if len(neighborhoods[k]) >= s_min]
    if len(keys) < 2:
        return 0.0
    vals: List[float] = []
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            inter = len(neighborhoods[a] & neighborhoods[b])
            union = len(neighborhoods[a] | neighborhoods[b])
            vals.append((inter / union) if union else 0.0)
    return statistics.mean(vals) if vals else 0.0


def compatibility_isolated_rate(
    configs: Sequence[HypergraphConfig],
    signature_fn: Callable[[HypergraphConfig], object],
    delta: int,
) -> float:
    """
    R_iso 操作化：在「不同實例、簽名距離 ≤ δ」之無向關係圖上，沒有此類鄰邊之節點占比。

    亦即僅與自身距離為 0 之點（在有限樣本中不與其他配置 δ-相容）視為近似孤立。
    """
    if not configs:
        return 0.0
    sigs = {c: signature_fn(c) for c in configs}
    n_iso = 0
    for c in configs:
        n_others = sum(
            1 for d in configs if d is not c and signature_distance(sigs[c], sigs[d]) <= delta
        )
        if n_others == 0:
            n_iso += 1
    return n_iso / len(configs)


def transitivity_violation_rate(
    configs: Sequence[HypergraphConfig],
    signature_fn: Callable[[HypergraphConfig], object],
    delta: int,
    max_triplets: int = 20_000,
) -> float:
    sigs = {c: signature_fn(c) for c in configs}
    triplets_checked = 0
    violations = 0
    n = len(configs)
    if n < 3:
        return 0.0

    indices = list(range(n))
    random.shuffle(indices)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                a, b, d = configs[indices[i]], configs[indices[j]], configs[indices[k]]
                ab = signature_distance(sigs[a], sigs[b]) <= delta
                bd = signature_distance(sigs[b], sigs[d]) <= delta
                ad = signature_distance(sigs[a], sigs[d]) <= delta
                if ab and bd and not ad:
                    violations += 1
                triplets_checked += 1
                if triplets_checked >= max_triplets:
                    return violations / triplets_checked
    return violations / triplets_checked if triplets_checked else 0.0


def class_distribution(
    classes: Sequence[Set[HypergraphConfig]],
    weight_fn: Callable[[HypergraphConfig], float] | None = None,
) -> List[float]:
    if weight_fn is None:
        weight_fn = lambda _c: 1.0
    weights = [sum(weight_fn(c) for c in cls) for cls in classes]
    total = sum(weights)
    if total <= 0:
        return []
    return [w / total for w in weights]


def entropy(probs: Sequence[float], base: float = 2.0) -> float:
    if base <= 0 or base == 1:
        raise ValueError("對數底須為正且不得為 1。")
    return -sum(p * (math.log(p) / math.log(base)) for p in probs if p > 0)


# -----------------------------------------------------------------------------
# Enumeration / sampling
# -----------------------------------------------------------------------------


def generate_candidate_configs(
    n: int,
    *,
    max_edge_size: int,
    max_edges: int,
    sample_limit: int | None,
    seed: int,
    progress: ProgressCallback = None,
) -> List[HypergraphConfig]:
    """產生候選超圖；``progress(done,total,msg)`` 於抽樣時回報，枚舉模式僅報起訖。"""
    vertices = tuple(range(1, n + 1))
    candidate_edges = all_possible_hyperedges(vertices, 2, max_edge_size)

    if sample_limit in (None, 0):
        if progress:
            progress(0, 1, "完整枚舉候選超圖（可能耗時／耗記憶體）…")
        configs: List[HypergraphConfig] = []
        for edge_subset in powerset_limited(candidate_edges, max_edges):
            configs.append(HypergraphConfig(vertices=vertices, hyperedges=frozenset(edge_subset)))
            if progress and len(configs) % 50_000 == 0:
                progress(0, 1, f"枚舉中… 已得 {len(configs)} 筆候選")
        if progress:
            progress(1, 1, f"枚舉完成，|候選|={len(configs)}")
        return configs

    rng = random.Random(seed)
    configs_set: Set[HypergraphConfig] = set()
    attempts = 0
    max_attempts = max(sample_limit * 20, 1000)
    tgt = max(1, int(sample_limit))
    if progress:
        progress(0, tgt, "隨機抽樣候選超圖…")
    while len(configs_set) < sample_limit and attempts < max_attempts:
        attempts += 1
        m = rng.randint(0, max_edges)
        chosen = rng.sample(candidate_edges, k=min(m, len(candidate_edges)))
        configs_set.add(HypergraphConfig(vertices=vertices, hyperedges=frozenset(chosen)))
        if progress and len(configs_set) % max(1, tgt // 25) == 0:
            progress(min(len(configs_set), tgt), tgt, f"候選 {len(configs_set)}/{tgt}")
    if progress:
        progress(tgt, tgt, f"抽樣完成，|候選|={len(configs_set)}")
    return list(configs_set)


def filter_configs(
    configs: Sequence[HypergraphConfig],
    *,
    max_edge_size: int,
    max_edges: int,
    connected_required: bool,
    max_degree: int | None,
    forbid_pair_triangles: bool,
    progress: ProgressCallback = None,
) -> List[HypergraphConfig]:
    """套用 Λ_dom；可選 ``progress`` 於長串列時回報篩選進度。"""
    out: List[HypergraphConfig] = []
    seq = list(configs)
    total = max(1, len(seq))
    step = max(1, total // 80)
    for i, c in enumerate(seq):
        if satisfies_domain_constraints(
            c,
            max_edge_size=max_edge_size,
            max_edges=max_edges,
            connected_required=connected_required,
            max_degree=max_degree,
            forbid_pair_triangles=forbid_pair_triangles,
        ):
            out.append(c)
        if progress and (i + 1) % step == 0:
            progress(i + 1, total, "域型約束篩選…")
    if progress:
        progress(total, total, "域型篩選完成")
    return out


# -----------------------------------------------------------------------------
# Dynamics
# -----------------------------------------------------------------------------


def all_legal_successors(
    c: HypergraphConfig,
    *,
    max_edge_size: int,
    max_edges: int,
    connected_required: bool,
    max_degree: int | None,
    forbid_pair_triangles: bool,
) -> List[HypergraphConfig]:
    vertices = c.vertices
    candidate_edges = all_possible_hyperedges(vertices, 2, max_edge_size)
    current_edges = set(c.hyperedges)
    successors: Set[HypergraphConfig] = set()

    for e in candidate_edges:
        if e not in current_edges:
            new_edges = frozenset(current_edges | {e})
            cand = HypergraphConfig(vertices, new_edges)
            if satisfies_domain_constraints(
                cand,
                max_edge_size=max_edge_size,
                max_edges=max_edges,
                connected_required=connected_required,
                max_degree=max_degree,
                forbid_pair_triangles=forbid_pair_triangles,
            ):
                successors.add(cand)

    for e in list(current_edges):
        new_edges = set(current_edges)
        new_edges.remove(e)
        cand = HypergraphConfig(vertices, frozenset(new_edges))
        if satisfies_domain_constraints(
            cand,
            max_edge_size=max_edge_size,
            max_edges=max_edges,
            connected_required=connected_required,
            max_degree=max_degree,
            forbid_pair_triangles=forbid_pair_triangles,
        ):
            successors.add(cand)

    for old_e in list(current_edges):
        reduced = set(current_edges)
        reduced.remove(old_e)
        for new_e in candidate_edges:
            if new_e in reduced:
                continue
            cand = HypergraphConfig(vertices, frozenset(reduced | {new_e}))
            if satisfies_domain_constraints(
                cand,
                max_edge_size=max_edge_size,
                max_edges=max_edges,
                connected_required=connected_required,
                max_degree=max_degree,
                forbid_pair_triangles=forbid_pair_triangles,
            ):
                successors.add(cand)

    successors.discard(c)
    return list(successors)


def choose_successor(
    c: HypergraphConfig,
    rng: random.Random,
    *,
    max_edge_size: int,
    max_edges: int,
    connected_required: bool,
    max_degree: int | None,
    forbid_pair_triangles: bool,
    allowed_configs: Set[HypergraphConfig] | None = None,
) -> HypergraphConfig:
    successors = all_legal_successors(
        c,
        max_edge_size=max_edge_size,
        max_edges=max_edges,
        connected_required=connected_required,
        max_degree=max_degree,
        forbid_pair_triangles=forbid_pair_triangles,
    )
    if allowed_configs is not None:
        successors = [s for s in successors if s in allowed_configs]
    if not successors:
        return c

    def score(cfg: HypergraphConfig) -> float:
        deg = degree_sequence(cfg)
        variance = statistics.pvariance(deg) if len(deg) > 1 else 0.0
        n3 = sum(1 for e in cfg.hyperedges if len(e) == 3)
        return -variance + 0.25 * n3

    weights = [math.exp(score(s)) for s in successors]
    total = sum(weights)
    r = rng.random() * total
    acc = 0.0
    for s, w in zip(successors, weights):
        acc += w
        if acc >= r:
            return s
    return successors[-1]


def _score_successor(cfg: HypergraphConfig) -> float:
    """後繼候選之權重分數（供軌道選擇用；非論文參數）。"""
    deg = degree_sequence(cfg)
    variance = statistics.pvariance(deg) if len(deg) > 1 else 0.0
    n3 = sum(1 for e in cfg.hyperedges if len(e) == 3)
    return -variance + 0.25 * n3


def _softmax_choice(
    rng: random.Random, candidates: Sequence[HypergraphConfig]
) -> HypergraphConfig:
    """對候選套用 softmax 權重後抽樣（與 choose_successor 一致的偏好）。"""
    if not candidates:
        raise ValueError("candidates 不可為空。")
    weights = [math.exp(_score_successor(s)) for s in candidates]
    total = sum(weights)
    r = rng.random() * total
    acc = 0.0
    for s, w in zip(candidates, weights):
        acc += w
        if acc >= r:
            return s
    return candidates[-1]


def _propose_local_update(
    current: HypergraphConfig,
    rng: random.Random,
    *,
    candidate_edges: Sequence[FrozenSet[int]],
    max_edges: int,
) -> HypergraphConfig:
    """
    產生一個「候選局部變動」之新配置（不先過濾是否合法）。

    依論文 §10.7.2（四）之語意，這裡的候選更新先提出，再由合法性檢查決定是否可接受。
    """
    cur = set(current.hyperedges)
    if not candidate_edges:
        return current

    # 依目前邊數狀態調整操作集合，避免必然失敗
    can_add = len(cur) < int(max_edges) and len(cur) < len(candidate_edges)
    can_remove = len(cur) > 0
    can_replace = can_remove and len(candidate_edges) > 0
    ops: list[str] = []
    if can_add:
        ops.append("add")
    if can_remove:
        ops.append("remove")
    if can_replace:
        ops.append("replace")
    if not ops:
        return current

    op = rng.choice(ops)
    if op == "add":
        # 隨機加入一條不存在之候選邊
        tries = 0
        while tries < 8:
            tries += 1
            e = rng.choice(candidate_edges)
            if e not in cur:
                cur.add(e)
                break
        return HypergraphConfig(current.vertices, frozenset(cur))
    if op == "remove":
        e = rng.choice(list(cur))
        cur.remove(e)
        return HypergraphConfig(current.vertices, frozenset(cur))

    # replace：先移除一條既有邊，再加入一條不存在之候選邊
    old_e = rng.choice(list(cur))
    cur.remove(old_e)
    tries = 0
    while tries < 12:
        tries += 1
        new_e = rng.choice(candidate_edges)
        if new_e not in cur:
            cur.add(new_e)
            break
    return HypergraphConfig(current.vertices, frozenset(cur))


def choose_successor_m_trial(
    c: HypergraphConfig,
    rng: random.Random,
    *,
    m_trial: int,
    max_edge_size: int,
    max_edges: int,
    connected_required: bool,
    max_degree: int | None,
    forbid_pair_triangles: bool,
    allowed_configs: Set[HypergraphConfig] | None = None,
) -> tuple[HypergraphConfig, float]:
    """
    單步動力學：依論文 §10.7.2（四）提出 M_trial 個候選更新，回傳下一狀態與 r_adm^{(ℓ)}。

    Note:
        - 本函式與既有 choose_successor 的差異在於「先提案、再檢查、最後在可接受者中抽樣」。
        - r_adm^{(ℓ)} 以「本步 admissible 候選數 / M_trial」操作化。
    """
    m = int(max(1, m_trial))
    vertices = c.vertices
    candidate_edges = all_possible_hyperedges(vertices, 2, max_edge_size)

    admissible: list[HypergraphConfig] = []
    adm_count = 0
    for _ in range(m):
        cand = _propose_local_update(c, rng, candidate_edges=candidate_edges, max_edges=max_edges)
        ok = satisfies_domain_constraints(
            cand,
            max_edge_size=max_edge_size,
            max_edges=max_edges,
            connected_required=connected_required,
            max_degree=max_degree,
            forbid_pair_triangles=forbid_pair_triangles,
        )
        if ok and allowed_configs is not None and cand not in allowed_configs:
            ok = False
        if ok:
            adm_count += 1
            admissible.append(cand)

    r_adm = adm_count / m if m > 0 else 0.0
    if not admissible:
        return c, r_adm
    # 在 admissible 候選中抽下一狀態（沿用既有偏好）
    return _softmax_choice(rng, admissible), r_adm


def run_trajectory(
    start: HypergraphConfig,
    *,
    steps: int,
    seed: int,
    max_edge_size: int,
    max_edges: int,
    connected_required: bool,
    max_degree: int | None,
    forbid_pair_triangles: bool,
    m_trial: int = 1,
    collect_r_adm: bool = False,
    allowed_configs: Set[HypergraphConfig] | None = None,
) -> Any:
    """
    產生一條長度為 T 的軌跡（含起點，共 T+1 筆）。

    Args:
        m_trial: 論文 §10.7.2（四）之每步候選更新數 M_trial；1 表示退化為舊版「每步取一個後繼」。
        collect_r_adm: 若 True，回傳 (traj, r_adm_series)；否則僅回傳 traj。
    """
    rng = random.Random(seed)
    traj = [start]
    r_adm_series: List[float] = []
    current = start
    for _ in range(steps):
        if int(m_trial) <= 1 and not collect_r_adm:
            current = choose_successor(
                current,
                rng,
                max_edge_size=max_edge_size,
                max_edges=max_edges,
                connected_required=connected_required,
                max_degree=max_degree,
                forbid_pair_triangles=forbid_pair_triangles,
                allowed_configs=allowed_configs,
            )
        else:
            current, r_adm = choose_successor_m_trial(
                current,
                rng,
                m_trial=int(max(1, m_trial)),
                max_edge_size=max_edge_size,
                max_edges=max_edges,
                connected_required=connected_required,
                max_degree=max_degree,
                forbid_pair_triangles=forbid_pair_triangles,
                allowed_configs=allowed_configs,
            )
            r_adm_series.append(float(r_adm))
        traj.append(current)
    if collect_r_adm:
        # r_adm_series 長度為 T（每步一值）
        return traj, r_adm_series
    return traj


# -----------------------------------------------------------------------------
# Analysis routines
# -----------------------------------------------------------------------------


def analyze_static(
    configs: Sequence[HypergraphConfig],
    signature_name: str,
    delta: int,
    *,
    s_min: int = 0,
) -> dict[str, Any]:
    """
    靜態解析指標（§10.3）：鄰域、重疊、傳違反、壓縮比、熵與 R_iso。

    Args:
        configs: 觀測集合 :math:`\\mathcal C_{obs}`（通常為自 :math:`\\mathrm{Cfg}_\\Lambda` 抽取之
            :math:`N_{cfg}` 筆；未抽樣時即全體可採用配置）。
        signature_name: weak / medium / strong。
        delta: 整數簽名距離閾值 δ（與 ``signature_distance`` 之離散尺度一致）。
        s_min: 重疊率計算時僅納入 |T(c)| >= s_min 之節點；0 為預設（不過濾）。

    Note:
        ``compression_ratio_U`` 取論文 :math:`U_\\Lambda = |S_{\\Lambda,\\delta}| / N_{cfg}`（此處
        :math:`N_{cfg}=``len(configs)``）；非舊版「平均壓縮倍數」:math:`N/|S|`。
    """
    signature_fn = SIGNATURES[signature_name]
    classes, neighborhoods = tolerance_classes_and_neighborhoods(configs, signature_fn, delta)
    probs = class_distribution(classes)

    avg_neighborhood = statistics.mean(len(v) for v in neighborhoods.values()) if neighborhoods else 0.0
    avg_class_size = statistics.mean(len(cls) for cls in classes) if classes else 0.0
    n_obs = len(configs)
    k_classes = len(classes)
    # 論文 §10.3.5：U_Λ = |S| / N_cfg（此處 N_cfg 即觀測集大小 n_obs）
    U = (k_classes / n_obs) if n_obs > 0 else 0.0

    return round_floats_for_output(
        {
            "num_configs": len(configs),
            "signature": signature_name,
            "delta": delta,
            "s_min": int(s_min),
            "avg_neighborhood_size": avg_neighborhood,
            "overlap_rate": overlap_rate(neighborhoods, s_min=s_min),
            "isol_rate_compat_graph": compatibility_isolated_rate(configs, signature_fn, delta),
            "transitivity_violation_rate": transitivity_violation_rate(configs, signature_fn, delta),
            "num_equivalence_classes": len(classes),
            "avg_class_size": avg_class_size,
            "compression_ratio_U": U,
            "entropy_bits": entropy(probs, base=2.0),
            "largest_classes": sorted((len(cls) for cls in classes), reverse=True)[:10],
            "sample_neighborhood_sizes": sorted((len(v) for v in neighborhoods.values()), reverse=True)[
                :10
            ],
        }
    )


def make_class_index(classes: Sequence[Set[HypergraphConfig]]) -> Dict[HypergraphConfig, int]:
    idx: Dict[HypergraphConfig, int] = {}
    for i, cls in enumerate(classes):
        for c in cls:
            idx[c] = i
    return idx


def time_distributions_from_trajectories(
    trajectories: Sequence[Sequence[HypergraphConfig]],
    classes: Sequence[Set[HypergraphConfig]],
) -> List[List[float]]:
    if not trajectories:
        return []
    idx = make_class_index(classes)
    T = len(trajectories[0])
    K = len(classes)
    out: List[List[float]] = []
    for t in range(T):
        counts = [0] * K
        for traj in trajectories:
            counts[idx[traj[t]]] += 1
        total = sum(counts)
        out.append([c / total for c in counts])
    return out


def analyze_dynamics(
    configs: Sequence[HypergraphConfig],
    *,
    signature_name: str,
    delta: int,
    runs: int,
    steps: int,
    seed: int,
    max_edge_size: int,
    max_edges: int,
    connected_required: bool,
    max_degree: int | None,
    forbid_pair_triangles: bool,
    m_trial: int = 1,
    w_h: int = 1,
    w_a: int = 20,
    p_max: int = 20,
    epsilon_plat: float = 0.02,
    progress: ProgressCallback = None,
) -> dict[str, Any]:
    if not configs:
        return {"error": "尚無可採用之配置，無法執行動力學。"}

    signature_fn = SIGNATURES[signature_name]
    classes, neighborhoods = tolerance_classes_and_neighborhoods(configs, signature_fn, delta)

    rng = random.Random(seed)
    config_set = set(configs)
    starts = [rng.choice(list(configs)) for _ in range(runs)]
    trajectories: List[List[HypergraphConfig]] = []
    r_adm_series_list: List[List[float]] = []
    total_runs = max(1, runs)
    for i, start in enumerate(starts):
        if progress:
            progress(i, total_runs, f"動力學軌道 {i + 1}/{total_runs}")
        traj, r_adm = run_trajectory(
            start,
            steps=steps,
            seed=seed + i + 1,
            max_edge_size=max_edge_size,
            max_edges=max_edges,
            connected_required=connected_required,
            max_degree=max_degree,
            forbid_pair_triangles=forbid_pair_triangles,
            m_trial=int(max(1, m_trial)),
            collect_r_adm=True,
            allowed_configs=config_set,
        )
        trajectories.append(traj)
        r_adm_series_list.append(r_adm)
    if progress:
        progress(total_runs, total_runs, "軌道模擬完成，彙總熵序列…")

    # 各時間步之類別計數（runs × 時間步），同時供 H_Λ 與窗口熵使用
    class_idx = make_class_index(classes)
    T = len(trajectories[0]) if trajectories else 0
    K = len(classes)
    time_counts: List[List[int]] = [[0] * K for _ in range(T)]
    labels_per_run: List[List[int]] = []
    for traj in trajectories:
        labels = [class_idx[c] for c in traj]
        labels_per_run.append(labels)
        for t, lab in enumerate(labels):
            time_counts[t][lab] += 1

    time_probs = []
    for t in range(T):
        tot = sum(time_counts[t])
        time_probs.append([c / tot for c in time_counts[t]] if tot > 0 else [0.0] * K)
    time_entropy = [entropy(p, base=2.0) for p in time_probs]

    # 窗口熵 H_{Λ,w_H}^{(ℓ)}：對 pooled counts 做滑動窗口加總
    w_h_eff = int(max(1, w_h))
    time_entropy_wH: List[float] = []
    if T > 0 and w_h_eff == 1:
        time_entropy_wH = list(time_entropy)
    elif T > 0:
        win_counts = [0] * K
        for t in range(T):
            # 加入當前 t
            ct = time_counts[t]
            for j in range(K):
                win_counts[j] += ct[j]
            # 移除離窗 t-w_h
            out_t = t - w_h_eff
            if out_t >= 0:
                ct_out = time_counts[out_t]
                for j in range(K):
                    win_counts[j] -= ct_out[j]
            win_len = min(w_h_eff, t + 1)
            tot = runs * win_len
            probs = [c / tot for c in win_counts] if tot > 0 else [0.0] * K
            time_entropy_wH.append(entropy(probs, base=2.0))

    step_fracs = [state_change_fraction_traj(tr) for tr in trajectories]
    plat = plateau_length_max_abs_diff(time_entropy, epsilon=float(epsilon_plat))
    plat_wH = (
        plateau_length_max_abs_diff(time_entropy_wH, epsilon=float(epsilon_plat))
        if time_entropy_wH
        else {"max_plateau_length": None, "num_plateau_segments": 0, "mean_plateau_length": None}
    )

    # r_adm：每步 admissible/M_trial；彙總為 per-run mean、以及各步 mean
    r_adm_per_run_mean: List[float] = [
        (statistics.mean(s) if s else 0.0) for s in r_adm_series_list
    ]
    r_adm_mean = statistics.mean(r_adm_per_run_mean) if r_adm_per_run_mean else 0.0
    r_adm_time_series_mean: List[float] = []
    if r_adm_series_list:
        for t in range(steps):
            vals = [rs[t] for rs in r_adm_series_list if t < len(rs)]
            r_adm_time_series_mean.append(statistics.mean(vals) if vals else 0.0)

    branch_sizes = [
        len(
            [
                x
                for x in all_legal_successors(
                    s,
                    max_edge_size=max_edge_size,
                    max_edges=max_edges,
                    connected_required=connected_required,
                    max_degree=max_degree,
                    forbid_pair_triangles=forbid_pair_triangles,
                )
                if x in config_set
            ]
        )
        for s in starts[: min(50, len(starts))]
    ]

    return_times: List[int] = []
    terminal_classes = Counter()
    for traj in trajectories:
        start_class = class_idx[traj[0]]
        terminal_classes[class_idx[traj[-1]]] += 1
        rt = None
        for t in range(1, len(traj)):
            if class_idx[traj[t]] == start_class:
                rt = t
                break
        if rt is not None:
            return_times.append(rt)

    # §10.7.5（五）可達狀態數：各軌跡上相異解析單元（類別索引）個數，再對軌跡彙總
    n_reach_per_run: List[int] = [
        len({class_idx[c] for c in traj}) for traj in trajectories
    ]

    # 週期與吸引子進入時刻：以解析類別標籤序列操作化（§10.7.2（八））
    w_a_eff = int(max(1, w_a))
    p_max_eff = int(max(1, p_max))
    p_cycle_per_run: List[int | None] = []
    ell_a_cycle_per_run: List[int | None] = []
    for labels in labels_per_run:
        best_p: int | None = None
        best_t: int | None = None
        # 窗口長度至少 2 才有意義；否則視為無法判定週期
        if len(labels) >= 2 and w_a_eff >= 2:
            max_t = len(labels) - 1
            for p in range(1, min(p_max_eff, max_t) + 1):
                found = False
                # t + p + w_a - 1 <= max_t
                for t0 in range(0, max_t - p - (w_a_eff - 1) + 1):
                    ok = True
                    for i in range(w_a_eff):
                        if labels[t0 + i] != labels[t0 + p + i]:
                            ok = False
                            break
                    if ok:
                        best_p = p
                        best_t = t0
                        found = True
                        break
                if found:
                    break
        p_cycle_per_run.append(best_p)
        ell_a_cycle_per_run.append(best_t)

    # 以熵序列判定平台吸引子進入時刻（連續 w_A-1 步皆 |ΔH|<=ε）
    def _ell_a_plat(series: Sequence[float]) -> int | None:
        if len(series) < 2 or w_a_eff < 2:
            return None
        for t0 in range(0, len(series) - w_a_eff + 1):
            ok = True
            for i in range(w_a_eff - 1):
                if abs(float(series[t0 + i + 1]) - float(series[t0 + i])) > float(epsilon_plat):
                    ok = False
                    break
            if ok:
                return t0
        return None

    ell_a_plat = _ell_a_plat(time_entropy_wH if time_entropy_wH else time_entropy)

    # 綜合 ell_A：取最早進入週期或平台者（以單一值供主表；同時輸出來源）
    ell_a_per_run: List[int | None] = []
    for t_cycle in ell_a_cycle_per_run:
        if ell_a_plat is None and t_cycle is None:
            ell_a_per_run.append(None)
        elif ell_a_plat is None:
            ell_a_per_run.append(t_cycle)
        elif t_cycle is None:
            ell_a_per_run.append(ell_a_plat)
        else:
            ell_a_per_run.append(min(int(ell_a_plat), int(t_cycle)))

    def _summ_int(xs: Sequence[int | None]) -> dict[str, float | int | None]:
        vals = [int(x) for x in xs if x is not None]
        if not vals:
            return {"mean": None, "min": None, "max": None}
        return {
            "mean": round(statistics.mean(vals), 6),
            "min": int(min(vals)),
            "max": int(max(vals)),
        }

    p_cycle_summary = _summ_int(p_cycle_per_run)
    ell_a_summary = _summ_int(ell_a_per_run)

    # largest_terminal_basins：將 Counter 鍵轉成可 JSON 之 int（類別索引）
    basin_list = [(int(k), int(v)) for k, v in terminal_classes.most_common(10)]

    return round_floats_for_output(
        {
            "num_configs": len(configs),
            "num_equivalence_classes": len(classes),
            "runs": runs,
            "steps": steps,
            "avg_branching_factor": statistics.mean(branch_sizes) if branch_sizes else 0.0,
            "avg_return_time": statistics.mean(return_times) if return_times else None,
            "num_distinct_terminal_classes": len(terminal_classes),
            "largest_terminal_basins": basin_list,
            "entropy_time_series": time_entropy,
            "entropy_time_series_wH": time_entropy_wH,
            "n_reach_per_run": n_reach_per_run,
            "n_reach_mean": (
                statistics.mean(n_reach_per_run) if n_reach_per_run else None
            ),
            "n_reach_max": max(n_reach_per_run) if n_reach_per_run else None,
            "n_reach_min": min(n_reach_per_run) if n_reach_per_run else None,
            "legal_update_step_fraction_mean": (
                statistics.mean(step_fracs) if step_fracs else 0.0
            ),
            "r_adm_time_series_mean": r_adm_time_series_mean,
            "r_adm_per_run_mean": r_adm_per_run_mean,
            "r_adm_mean": r_adm_mean,
            "m_trial": int(max(1, m_trial)),
            "w_h": int(max(1, w_h)),
            "w_a": int(max(1, w_a)),
            "p_max": int(max(1, p_max)),
            "p_cycle_per_run": p_cycle_per_run,
            "ell_a_cycle_per_run": ell_a_cycle_per_run,
            "ell_a_plat": ell_a_plat,
            "p_cycle_summary": p_cycle_summary,
            "ell_a_summary": ell_a_summary,
            "epsilon_plat": float(epsilon_plat),
            "entropy_summary": {
                "start": time_entropy[0] if time_entropy else None,
                "end": time_entropy[-1] if time_entropy else None,
                "max": max(time_entropy) if time_entropy else None,
                "min": min(time_entropy) if time_entropy else None,
                "mean": statistics.mean(time_entropy) if time_entropy else None,
                "plateau_max_length": plat.get("max_plateau_length"),
                "plateau_num_segments": plat.get("num_plateau_segments"),
                "plateau_mean_length": plat.get("mean_plateau_length"),
                "epsilon_plat": float(epsilon_plat),
            },
            "entropy_summary_wH": {
                "start": time_entropy_wH[0] if time_entropy_wH else None,
                "end": time_entropy_wH[-1] if time_entropy_wH else None,
                "max": max(time_entropy_wH) if time_entropy_wH else None,
                "min": min(time_entropy_wH) if time_entropy_wH else None,
                "mean": statistics.mean(time_entropy_wH) if time_entropy_wH else None,
                "plateau_max_length": plat_wH.get("max_plateau_length"),
                "plateau_num_segments": plat_wH.get("num_plateau_segments"),
                "plateau_mean_length": plat_wH.get("mean_plateau_length"),
                "epsilon_plat": float(epsilon_plat),
            },
        }
    )


# -----------------------------------------------------------------------------
# 完整實驗執行（CLI / GUI 共用）
# -----------------------------------------------------------------------------


def sample_candidates_and_filter(
    *,
    n: int,
    max_edge_size: int,
    max_edges: int,
    sample_limit: int,
    seed: int,
    connected: bool,
    max_degree: int | None,
    forbid_pair_triangles: bool,
    progress: ProgressCallback = None,
) -> Tuple[List[HypergraphConfig], List[HypergraphConfig]]:
    """
    產生候選超圖集合並套用域型約束，回傳 (候選列表, 可採用配置列表)。

    與 ``run_full_experiment`` 前半段邏輯相同，供 Streamlit 表 10-2／10-3 或外部工具重複使用。
    """
    random.seed(seed)
    candidates = generate_candidate_configs(
        n,
        max_edge_size=max_edge_size,
        max_edges=max_edges,
        sample_limit=sample_limit,
        seed=seed,
        progress=progress,
    )
    configs = filter_configs(
        candidates,
        max_edge_size=max_edge_size,
        max_edges=max_edges,
        connected_required=connected,
        max_degree=max_degree,
        forbid_pair_triangles=forbid_pair_triangles,
        progress=progress,
    )
    return candidates, configs


# 與候選母集抽樣 seed 區隔，供 §10.3 觀測集 \\mathcal C_obs 可重現子樣本
_OBS_SUBSAMPLE_SEED_OFFSET = 10_303
# §10.5 合法配置子樣本（論文 N_cfg）專用偏移，與 §10.3 解耦
OBS_SUBSAMPLE_RNG_OFFSET_10_5 = 10_505


def subsample_obs_configs(
    configs: Sequence[HypergraphConfig],
    n_cfg: Optional[int],
    *,
    seed: int,
    rng_chain_offset: int = _OBS_SUBSAMPLE_SEED_OFFSET,
) -> Tuple[List[HypergraphConfig], Optional[int], int, Optional[str]]:
    """
    論文 §10.3.2：自可採用域 :math:`\\mathrm{Cfg}_\\Lambda` 得到觀測集 :math:`\\mathcal C_{obs}`。

    Args:
        configs: 域型篩選後之可採用配置序列。
        n_cfg: 欲抽取之配置數 :math:`N_{cfg}`；``None`` 表示以**全部**可採用配置為觀測集（向後相容）。
        seed: 偽隨機基底，與候選抽樣共用同一專案慣例但加偏移以避免耦合。
        rng_chain_offset: 加在 ``seed`` 上之鏈節偏移；§10.5 請傳 ``OBS_SUBSAMPLE_RNG_OFFSET_10_5``，避免與 §10.3 共用序列。

    Returns:
        (觀測集列表, 請求之 N_cfg 或 None, 實際觀測筆數, 若可採用數不足請求時之繁中提示否则 None)。
    """
    cfg_list = list(configs)
    if not cfg_list:
        return [], (int(n_cfg) if n_cfg is not None else None), 0, None
    if n_cfg is None:
        return cfg_list, None, len(cfg_list), None
    n_req = max(1, int(n_cfg))
    n_avail = len(cfg_list)
    if n_avail >= n_req:
        rng = random.Random(int(seed) + int(rng_chain_offset))
        obs = rng.sample(cfg_list, n_req)
        return obs, n_req, n_req, None
    notice = (
        f"可採用配置僅 {n_avail} 筆，低於請求之 N_cfg={n_req}；"
        "已以全部可採用集合作為觀測集。請加大候選採樣上限或放寬域型條件。"
    )
    return cfg_list, n_req, n_avail, notice


def _aggregate_static_repetitions(rep_rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """
    匯總 §10.3 靜態重抽結果（N_rep）。

    目前僅對主要標量指標輸出平均與標準差；其餘結構欄位（如 largest_classes）
    由各 repetition 明細保留，避免語義混淆。
    """
    scalar_keys: Tuple[str, ...] = (
        "overlap_rate",
        "transitivity_violation_rate",
        "isol_rate_compat_graph",
        "num_equivalence_classes",
        "avg_class_size",
        "compression_ratio_U",
        "entropy_bits",
        "avg_neighborhood_size",
        "num_configs",
    )
    out: dict[str, Any] = {"n_rep_effective": len(rep_rows)}
    for key in scalar_keys:
        vals: List[float] = []
        for row in rep_rows:
            v = row.get(key)
            if isinstance(v, (int, float)):
                vals.append(float(v))
        if not vals:
            continue
        out[f"{key}_mean"] = statistics.mean(vals)
        out[f"{key}_std"] = statistics.pstdev(vals) if len(vals) > 1 else 0.0
        out[f"{key}_min"] = min(vals)
        out[f"{key}_max"] = max(vals)
    return round_floats_for_output(out)


def _get_nested_value(d: dict[str, Any], path: Sequence[str]) -> Any:
    """安全取得巢狀字典值；任一路徑不存在則回傳 None。"""
    cur: Any = d
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur


def _aggregate_dynamics_seed_runs(seed_rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """
    匯總 §10.7 動力學多 seed 重跑結果（N_seed）。

    註：此處僅匯總「標量」指標；序列（如 entropy_time_series）仍由各 seed 明細保留，
    以免混合後失去語義。
    """
    out: dict[str, Any] = {"n_seed_effective": len(seed_rows)}

    # (輸出鍵, 取值路徑)；路徑可指向巢狀欄位
    scalar_paths: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
        ("r_adm_mean", ("r_adm_mean",)),
        ("legal_update_step_fraction_mean", ("legal_update_step_fraction_mean",)),
        ("n_reach_mean", ("n_reach_mean",)),
        ("entropy_mean", ("entropy_summary", "mean")),
        ("entropy_wH_mean", ("entropy_summary_wH", "mean")),
        ("plateau_max_length", ("entropy_summary", "plateau_max_length")),
        ("plateau_wH_max_length", ("entropy_summary_wH", "plateau_max_length")),
        ("p_cycle_mean", ("p_cycle_summary", "mean")),
        ("ell_A_mean", ("ell_a_summary", "mean")),
    )

    for out_key, path in scalar_paths:
        vals: List[float] = []
        for row in seed_rows:
            v = _get_nested_value(row, path)
            if isinstance(v, (int, float)):
                vals.append(float(v))
        if not vals:
            continue
        out[f"{out_key}_mean"] = statistics.mean(vals)
        out[f"{out_key}_std"] = statistics.pstdev(vals) if len(vals) > 1 else 0.0
        out[f"{out_key}_min"] = min(vals)
        out[f"{out_key}_max"] = max(vals)
    return round_floats_for_output(out)


def run_full_experiment(
    *,
    mode: str = "static",
    n: int = 5,
    max_edge_size: int = 3,
    max_edges: int = 4,
    max_degree: int = 4,
    connected: bool = False,
    forbid_pair_triangles: bool = False,
    sample_limit: int = 2000,
    n_cfg: Optional[int] = None,
    n_rep: int = 1,
    signature: str = "medium",
    delta: int = 0,
    s_min: int = 0,
    epsilon_plat: float = 0.02,
    runs: int = 50,
    steps: int = 30,
    m_trial: int = 1,
    w_h: int = 1,
    w_a: int = 20,
    p_max: int = 20,
    n_seed_107: int = 1,
    seed: int = 7,
    show_sample_configs: int = 3,
    refinement_enabled: bool = False,
    refine_coarse_signature: str = "weak",
    refine_coarse_delta: int = 3,
    refine_fine_signature: str = "medium",
    refine_fine_delta: int = 0,
    refine_kernel: str = "uniform",
    refine_compare_chains: bool = True,
    refine_coarse_sample_size: int | None = None,
    refine_fine_sample_size: int | None = None,
    refine_fiber_sample_size: int | None = None,
    progress: ProgressCallback = None,
) -> dict[str, Any]:
    """
    執行與 `experiment.py` CLI 相同流程之完整實驗，回傳可 JSON 序列化之結果字典。

    Args:
        mode: ``\"static\"`` 或 ``\"dynamics\"``。
        n: 節點數 |V|。
        max_edge_size: 超邊最大容量 k_max。
        max_edges: 超邊數上限 m_max。
        max_degree: 參與超邊數上限 d_max（各頂點）。
        connected: 是否要求 2-section 連通。
        forbid_pair_triangles: 是否禁止僅由二元邊構成之三角形。
        sample_limit: 候選配置抽樣數；0 表示在可行時完整枚舉幂集子集。
        n_cfg: 論文 :math:`N_{cfg}`，自可採用域抽樣之**觀測配置數**；``None`` 表示以全部可採用配置分析。
        n_rep: 論文 :math:`N_{rep}`，觀測集重抽次數（僅 ``static`` 模式；>=1）。
        signature: 解析簽名層次 ``weak`` / ``medium`` / ``strong``。
        delta: 簽名距離閾值 δ（整數；與 ``signature_distance`` 一致）。
        s_min: §10.3 重疊率可選過濾，僅 |T(c)| >= s_min 參與；0 表示不過濾。
        epsilon_plat: §10.7 熵時間序列平台判定 |ΔH| 閾值。
        runs: 動力學軌道條數。
        steps: 每條軌道時間步數 T。
        m_trial: §10.7 每步候選更新數 M_trial（提案—過濾—抽樣）。
        w_h: §10.7 熵滑動窗口寬度 w_H；1 表示逐步熵。
        w_a: §10.7 吸引子／週期檢測窗口 w_A（以解析標籤或熵差分操作化）。
        p_max: §10.7 週期搜尋上限 P_max。
        n_seed_107: §10.7 建議之同參數重跑次數 N_seed（穩健性）；1 表示不重跑。
        seed: 隨機種子。
        show_sample_configs: 輸出中保留之範例超圖筆數。
        refinement_enabled: 是否附加 §10.4 細化／投影／條件核分析（僅 ``static`` 模式）。
        refine_coarse_signature, refine_coarse_delta: 粗解析 Λ。
        refine_fine_signature, refine_fine_delta: 細解析 Λ′。
        refine_kernel: ``uniform`` 或 ``proportional``（纖維條件核）。
        refine_compare_chains: 是否計算 §10.4.1 雙路徑終端 JS／熵差。
        refine_coarse_sample_size: §10.4 粗層樣本數（None 表示不額外截取）。
        refine_fine_sample_size: §10.4 細層樣本數（None 表示不額外截取）。
        refine_fiber_sample_size: §10.4 每粗單元纖維樣本上限（None 表示不截斷）。
        progress: 可選進度回呼 ``(done,total,message)``（外層約分四相：候選、分析、細化、完成）。

    Returns:
        含 ``parameters``、``num_candidates``、``num_admissible_configs``、
        ``num_obs_configs``、``sample_configs``、``analysis`` 等；細化時另含 ``refinement_10_4``。
    """
    candidates, configs = sample_candidates_and_filter(
        n=n,
        max_edge_size=max_edge_size,
        max_edges=max_edges,
        sample_limit=sample_limit,
        seed=seed,
        connected=connected,
        max_degree=max_degree,
        forbid_pair_triangles=forbid_pair_triangles,
        progress=progress,
    )

    obs, n_req, n_obs_actual, obs_notice = subsample_obs_configs(configs, n_cfg, seed=seed)

    if progress:
        progress(1, 4, "解析／動力學分析…")

    result: dict[str, Any] = {
        "parameters": {
            "mode": mode,
            "n": n,
            "max_edge_size": max_edge_size,
            "max_edges": max_edges,
            "max_degree": max_degree,
            "connected": connected,
            "forbid_pair_triangles": forbid_pair_triangles,
            "sample_limit": sample_limit,
            "n_cfg": n_cfg,
            "n_rep": int(max(1, n_rep)),
            "signature": signature,
            "delta": delta,
            "s_min": s_min,
            "epsilon_plat": epsilon_plat,
            "runs": runs,
            "steps": steps,
            "m_trial": int(m_trial),
            "w_h": int(w_h),
            "w_a": int(w_a),
            "p_max": int(p_max),
            "n_seed_107": int(n_seed_107),
            "seed": seed,
            "refinement_enabled": refinement_enabled,
            "refine_coarse_signature": refine_coarse_signature,
            "refine_coarse_delta": refine_coarse_delta,
            "refine_fine_signature": refine_fine_signature,
            "refine_fine_delta": refine_fine_delta,
            "refine_kernel": refine_kernel,
            "refine_compare_chains": refine_compare_chains,
            "refine_coarse_sample_size": refine_coarse_sample_size,
            "refine_fine_sample_size": refine_fine_sample_size,
            "refine_fiber_sample_size": refine_fiber_sample_size,
        },
        "num_candidates": len(candidates),
        "num_admissible_configs": len(configs),
        "num_obs_configs": n_obs_actual,
        "n_cfg_requested": n_req,
        "sample_configs": [c.to_dict() for c in obs[:show_sample_configs]],
    }
    if obs_notice:
        result["n_cfg_notice"] = obs_notice

    n_rep_eff = int(max(1, n_rep))
    if mode == "static":
        if n_rep_eff == 1:
            result["analysis"] = analyze_static(obs, signature, delta, s_min=s_min)
        else:
            rep_rows: List[dict[str, Any]] = []
            rep_meta: List[dict[str, Any]] = []
            # 使用固定偏移規則派生各次重抽種子，確保可重現且與候選生成解耦。
            for rep_idx in range(n_rep_eff):
                rep_seed = int(seed) + rep_idx * 1009
                rep_obs, _rep_req, rep_obs_n, _rep_notice = subsample_obs_configs(
                    configs, n_cfg, seed=rep_seed
                )
                rep_analysis = analyze_static(rep_obs, signature, delta, s_min=s_min)
                rep_rows.append(rep_analysis)
                rep_meta.append(
                    {
                        "rep_index": rep_idx,
                        "rep_seed": rep_seed,
                        "num_obs_configs": rep_obs_n,
                    }
                )
            # 向後相容：analysis 仍提供單次結構（取第 0 次），圖表不需改動。
            result["analysis"] = rep_rows[0] if rep_rows else {"error": "N_rep 重抽無結果。"}
            result["analysis_repetitions"] = [
                {"meta": meta, "analysis": row} for meta, row in zip(rep_meta, rep_rows)
            ]
            result["analysis_rep_summary"] = _aggregate_static_repetitions(rep_rows)
            result["analysis_rep_summary"]["seed_rule"] = "seed + rep_index * 1009"
            obs_counts = [m["num_obs_configs"] for m in rep_meta]
            if obs_counts:
                result["num_obs_configs_min"] = min(obs_counts)
                result["num_obs_configs_max"] = max(obs_counts)
    else:
        n_seed_eff = int(max(1, n_seed_107))
        if n_seed_eff == 1:
            result["analysis"] = analyze_dynamics(
                obs,
                signature_name=signature,
                delta=delta,
                runs=runs,
                steps=steps,
                seed=seed,
                max_edge_size=max_edge_size,
                max_edges=max_edges,
                connected_required=connected,
                max_degree=max_degree,
                forbid_pair_triangles=forbid_pair_triangles,
                m_trial=int(m_trial),
                w_h=int(w_h),
                w_a=int(w_a),
                p_max=int(p_max),
                epsilon_plat=epsilon_plat,
                # 與 run_full_experiment 外層進度刻度分離，避免子回呼覆寫相位；單獨呼叫 analyze_dynamics 仍可傳 progress。
                progress=None,
            )
        else:
            seed_rows: List[dict[str, Any]] = []
            seed_meta: List[dict[str, Any]] = []
            for seed_idx in range(n_seed_eff):
                run_seed = int(seed) + seed_idx * 1009
                an = analyze_dynamics(
                    obs,
                    signature_name=signature,
                    delta=delta,
                    runs=runs,
                    steps=steps,
                    seed=run_seed,
                    max_edge_size=max_edge_size,
                    max_edges=max_edges,
                    connected_required=connected,
                    max_degree=max_degree,
                    forbid_pair_triangles=forbid_pair_triangles,
                    m_trial=int(m_trial),
                    w_h=int(w_h),
                    w_a=int(w_a),
                    p_max=int(p_max),
                    epsilon_plat=epsilon_plat,
                    progress=None,
                )
                seed_rows.append(an)
                seed_meta.append(
                    {
                        "seed_index": seed_idx,
                        "seed": run_seed,
                    }
                )
            # 向後相容：analysis 仍提供第 0 次結構（圖表與表格不需改動）
            result["analysis"] = seed_rows[0] if seed_rows else {"error": "N_seed 重跑無結果。"}
            result["analysis_seeds"] = [
                {"meta": meta, "analysis": row} for meta, row in zip(seed_meta, seed_rows)
            ]
            result["analysis_seed_summary"] = _aggregate_dynamics_seed_runs(seed_rows)
            result["analysis_seed_summary"]["seed_rule"] = "seed + seed_index * 1009"

    if progress:
        progress(2, 4, "主要分析完成")

    if refinement_enabled and mode == "static" and configs:
        from hypergraph_experiment.refinement import analyze_section_10_4_bundle

        km = refine_kernel if refine_kernel in ("uniform", "proportional") else "uniform"
        if progress:
            progress(3, 4, "§10.4 細化／纖維分析…")
        result["refinement_10_4"] = analyze_section_10_4_bundle(
            configs,
            coarse_sig=refine_coarse_signature,
            coarse_delta=refine_coarse_delta,
            fine_sig=refine_fine_signature,
            fine_delta=refine_fine_delta,
            kernel_mode=km,
            compare_chains=refine_compare_chains,
            coarse_sample_size=refine_coarse_sample_size,
            fine_sample_size=refine_fine_sample_size,
            sample_seed=int(seed),
            max_fiber_size=refine_fiber_sample_size,
        )

    if progress:
        progress(4, 4, "實驗完成")

    return round_floats_for_output(result)


def result_to_json(result: dict[str, Any], *, indent: int = 2) -> str:
    """將實驗結果序列化為 UTF-8 JSON 字串（供 CLI 列印或寫檔）。"""
    return json.dumps(round_floats_for_output(result), ensure_ascii=False, indent=indent)
