"""
第十章 §10.4：解析細化 Λ⊆Λ′、投影 π_{Λ′→Λ}、纖維條件核 K_s、推前一致性與雙路徑細化鏈比較。

對應論文 Definition 10.4-A～D、Proposition 10.4-C 與 §10.4.1 之有限離散實作。
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Literal, Sequence, Set, Tuple, TypedDict

from hypergraph_experiment.output_rounding import round_floats_for_output
from hypergraph_experiment.core import (
    HypergraphConfig,
    SIGNATURES,
    tolerance_equivalence_classes,
    entropy,
    make_class_index,
    subsample_obs_configs,
)

KernelMode = Literal["uniform", "proportional"]

# ---------------------------------------------------------------------------
# 實驗 B（§10.4.7）：A/B 子步驟以「映射類型 + 參數」表徵，對齊論文 9.6-C 之 R_A、R_B 語義；
# 下方操作化對應 core.py 之簽名特徵族（非直接輸入 weak/medium/strong 當作子步驟名稱）。
# ---------------------------------------------------------------------------

SubstepKind = Literal["edge_scale_motif", "degree_split", "adjacency_motif_fine"]


class RefinementSubstepSpec(TypedDict):
    """
    單一子步驟之有限操作化規格。

    Attributes:
        kind: 映射種類——超邊尺度、度序列細分、或鄰接+動機之較精細特徵。
        delta: 相容閾值（與 §10.3 整數 δ 及 signature_distance 一致）。
    """

    kind: SubstepKind
    delta: int


def substep_spec_to_signature_delta(spec: RefinementSubstepSpec) -> Tuple[str, int]:
    """
    將子步驟規格轉成內部分割所用之 (signature_name, delta)。

    對應關係（toy 固定）：
    - edge_scale_motif → weak（僅邊階多重集，對應動機/尺度面向之較粗辨識）
    - degree_split → medium（邊階 + 度數序列）
    - adjacency_motif_fine → strong（鄰接 + motif_counts 等）
    """
    m: Dict[SubstepKind, str] = {
        "edge_scale_motif": "weak",
        "degree_split": "medium",
        "adjacency_motif_fine": "strong",
    }
    kind = spec["kind"]
    if kind not in m:
        raise ValueError(f"未知的子步驟種類: {kind!r}")
    return m[kind], int(spec["delta"])


def validate_substep_spec(spec: Any, *, label: str) -> RefinementSubstepSpec:
    """
    驗證實驗 B 子步驟規格是否完整且語義正確。

    Args:
        spec: 使用者輸入的子步驟規格（預期為包含 kind 與 delta 的字典）。
        label: 錯誤訊息中的參數名稱（step_a 或 step_b）。

    Returns:
        正規化後的 RefinementSubstepSpec（delta 轉為 int）。
    """
    if not isinstance(spec, dict):
        raise ValueError(f"{label} 必須是 dict，格式為 {{'kind': <SubstepKind>, 'delta': <int>}}。")
    if "kind" not in spec or "delta" not in spec:
        raise ValueError(f"{label} 缺少必要欄位，需同時包含 'kind' 與 'delta'。")
    kind = spec["kind"]
    delta = spec["delta"]
    if kind not in {"edge_scale_motif", "degree_split", "adjacency_motif_fine"}:
        raise ValueError(
            f"{label}.kind 無效：{kind!r}。可用值為 edge_scale_motif、degree_split、adjacency_motif_fine。"
        )
    try:
        delta_int = int(delta)
    except (TypeError, ValueError) as e:
        raise ValueError(f"{label}.delta 必須可轉為整數，收到：{delta!r}") from e
    return {"kind": kind, "delta": delta_int}


def preset_layer_to_substep_spec(layer_sig: str, delta: int) -> RefinementSubstepSpec:
    """
    將主表 preset 之粗/細層標籤（weak／medium／strong）轉為實驗 B 預設子步驟規格。

    僅供 UI 預設對齊論文 10.4.6 固定鏈；實驗 B 仍可手動覆寫 kind與 δ。
    """
    m = {"weak": "edge_scale_motif", "medium": "degree_split", "strong": "adjacency_motif_fine"}
    if layer_sig not in m:
        raise ValueError(f"無法由 preset 層轉子步驟：{layer_sig!r}")
    return {"kind": m[layer_sig], "delta": int(delta)}


def partition_analytic_units(
    configs: Sequence[HypergraphConfig],
    signature_name: str,
    delta: int,
) -> Tuple[List[Set[HypergraphConfig]], Dict[HypergraphConfig, int]]:
    """
    在給定簽名與 δ 下，建立解析等價類（即 S_Λ 之元素為配置子集）與配置到類別索引之映射 q_Λ。

    Args:
        configs: 可採用配置集合 Cfg_Λ。
        signature_name: weak / medium / strong。
        delta: 簽名距離閾值。

    Returns:
        (classes, config_to_index)，其中 classes[k] 為第 k 個解析單元所含之配置集合。
    """
    if not configs:
        return [], {}
    fn = SIGNATURES[signature_name]
    classes = tolerance_equivalence_classes(configs, fn, delta)
    idx = make_class_index(classes)
    return classes, idx


def refinement_holds(
    classes_fine: Sequence[Set[HypergraphConfig]],
    config_to_coarse: Dict[HypergraphConfig, int],
) -> bool:
    """
    檢查細層分割是否為粗層分割之細化：每個細解析單元必須完全落在某個粗解析單元內。

    對應論文中 Λ⊆Λ′ 時，每個 s′ 所屬之粗類型唯一。
    """
    for fcls in classes_fine:
        if not fcls:
            continue
        rep = next(iter(fcls))
        i0 = config_to_coarse[rep]
        for c in fcls:
            if config_to_coarse[c] != i0:
                return False
    return True


def projection_pi_fine_to_coarse(
    classes_fine: Sequence[Set[HypergraphConfig]],
    config_to_coarse: Dict[HypergraphConfig, int],
) -> List[int]:
    """
    建構 π_{Λ′→Λ}：對每個細層類別索引 j，回傳其投影之粗層索引 π(j)。

    實作上取該細類中任一配置 c，令 π(j)=q_Λ(c)；在細化成立時與代表元選取無關。
    """
    pi_map: List[int] = []
    for fcls in classes_fine:
        if not fcls:
            pi_map.append(-1)
            continue
        pi_map.append(config_to_coarse[next(iter(fcls))])
    return pi_map


def fibers_from_pi(pi: Sequence[int], n_coarse: int) -> Dict[int, List[int]]:
    """
    對每個粗層單元 s∈S_Λ，列出纖維 π^{-1}(s) 內之細層索引集合。

    Args:
        pi: π(j) 列表，j 為細層類別索引。
        n_coarse: |S_Λ|（粗層類別數）。

    Returns:
        字典 s -> [j1, j2, ...]。
    """
    fib: Dict[int, List[int]] = {s: [] for s in range(n_coarse)}
    for j, s in enumerate(pi):
        if 0 <= s < n_coarse:
            fib[s].append(j)
    return fib


def kernel_stability_summary(
    kernels: Dict[int, Dict[int, float]],
    *,
    mass_eps: float = 1e-12,
) -> Dict[str, Any]:
    """
    纖維核穩定度摘要：各 K_s 之有效非零質量比例與熵（§10.4 觀察補助）。

    Args:
        kernels: 粗索引 s → 細索引 j → 非負質量（通常為機率）。
        mass_eps: 視為「有效非零」之門檻。

    Returns:
        平均非零比例、平均熵（bit）、纖維數等可 JSON 欄位。
    """
    if not kernels:
        return {"kernel_n_fibers": 0}
    nnz_ratios: List[float] = []
    entropies: List[float] = []
    for kdict in kernels.values():
        vals = [float(v) for v in kdict.values()]
        if not vals:
            continue
        nnz = sum(1 for v in vals if v > mass_eps)
        nnz_ratios.append(nnz / len(vals))
        tot = sum(vals)
        if tot <= 0:
            entropies.append(0.0)
            continue
        pn = [v / tot for v in vals]
        entropies.append(entropy(pn, base=2.0))
    return {
        "kernel_n_fibers": len(kernels),
        "kernel_mean_nnz_ratio": sum(nnz_ratios) / len(nnz_ratios) if nnz_ratios else 0.0,
        "kernel_mean_entropy_bits": sum(entropies) / len(entropies) if entropies else 0.0,
    }


def build_kernels(
    fibers: Dict[int, List[int]],
    classes_fine: Sequence[Set[HypergraphConfig]],
    mode: KernelMode,
    max_fiber_size: int | None = None,
) -> Dict[int, Dict[int, float]]:
    """
    在每個粗單元 s 之纖維上建條件核 K_s（細層索引 j→機率）。

    - uniform：纖維內均匀。
    - proportional：與細層類別大小 |C_j| 成正比。
    """
    kernels: Dict[int, Dict[int, float]] = {}
    for s, js in fibers.items():
        if not js:
            continue
        js_use = list(js)
        if max_fiber_size is not None and max_fiber_size > 0 and len(js_use) > max_fiber_size:
            # 纖維樣本數上限：固定取排序後前 K 個，確保可重現。
            js_use = sorted(js_use)[: int(max_fiber_size)]
        if mode == "uniform":
            w = 1.0 / len(js_use)
            kernels[s] = {j: w for j in js_use}
        else:
            sizes = [max(1, len(classes_fine[j])) for j in js_use]
            tot = float(sum(sizes))
            kernels[s] = {j: sizes[k] / tot for k, j in enumerate(js_use)}
    return kernels


def induce_fine_distribution(
    p_coarse: Sequence[float],
    pi: Sequence[int],
    kernels: Dict[int, Dict[int, float]],
) -> List[float]:
    """由 p_Λ 與 {K_s} 得到 p_{Λ′}。"""
    n_f = len(pi)
    p_fine = [0.0] * n_f
    for j in range(n_f):
        s = pi[j]
        if s < 0 or s >= len(p_coarse):
            continue
        kdict = kernels.get(s)
        if not kdict:
            continue
        w = kdict.get(j, 0.0)
        p_fine[j] = float(p_coarse[s]) * w
    return p_fine


def pushforward_from_fine(
    p_fine: Sequence[float],
    fibers: Dict[int, List[int]],
) -> List[float]:
    """由 p_{Λ′} 在纖維上求和得到推前 (π_* p_{Λ′})_s。"""
    n_c = max(fibers.keys(), default=-1) + 1
    recovered = [0.0] * n_c
    for s, js in fibers.items():
        if s < 0 or s >= n_c:
            continue
        for j in js:
            if j < len(p_fine):
                recovered[s] += float(p_fine[j])
    return recovered


def max_pushforward_error(p_coarse: Sequence[float], recovered: Sequence[float]) -> float:
    """dTV 型：逐座標最大絕對差作為推前誤差上界之簡化度量。"""
    if not p_coarse:
        return 0.0
    m = min(len(p_coarse), len(recovered))
    return max(abs(float(p_coarse[i]) - float(recovered[i])) for i in range(m)) if m else 0.0


def uniform_distribution_on_configs(classes_coarse: Sequence[Set[HypergraphConfig]]) -> List[float]:
    """在配置層均勻時，粗層誘導邊際 p_Λ(s) = |C_s| / |Cfg|。"""
    sizes = [len(c) for c in classes_coarse]
    tot = sum(sizes)
    if tot <= 0:
        return []
    return [s / tot for s in sizes]


def _kl_divergence_bits(p: Sequence[float], q: Sequence[float]) -> float:
    """以 2 為底之 KL(p||q)。"""
    s = 0.0
    for pi, qi in zip(p, q):
        if pi <= 0.0:
            continue
        if qi <= 0.0:
            return float("inf")
        s += pi * math.log(pi / qi) / math.log(2.0)
    return s


def js_divergence_bits(p: Sequence[float], q: Sequence[float]) -> float | None:
    """
    Jensen–Shannon 散度（比特），用於 §10.4.1 終端分布差異。

    若支撐不交疊導致無法定義，回傳 None。
    """
    sp, sq = sum(p), sum(q)
    if sp <= 0 or sq <= 0:
        return None
    pn = [float(x) / sp for x in p]
    qn = [float(x) / sq for x in q]
    if len(pn) != len(qn):
        return None
    m = [0.5 * (a + b) for a, b in zip(pn, qn)]
    if any(x <= 0.0 for x in m):
        return None
    return 0.5 * _kl_divergence_bits(pn, m) + 0.5 * _kl_divergence_bits(qn, m)


def apply_refinement_step(
    p_coarse: Sequence[float],
    classes_coarse: Sequence[Set[HypergraphConfig]],
    classes_fine: Sequence[Set[HypergraphConfig]],
    kernel_mode: KernelMode,
    max_fiber_size: int | None = None,
) -> Tuple[List[float], Dict[str, Any]]:
    """
    單步 Λ→Λ′：由 p_Λ 與 {K_s} 得到 p_{Λ′}，並回傳可序列化之中間結果。

    Returns:
        (p_fine, detail) 其中 detail 含 pi、fibers、kernels、推前誤差、細化是否成立。
    """
    config_to_coarse = make_class_index(classes_coarse)
    valid = refinement_holds(classes_fine, config_to_coarse)
    pi = projection_pi_fine_to_coarse(classes_fine, config_to_coarse)
    n_c = len(classes_coarse)
    fibers = fibers_from_pi(pi, n_c)
    kernels = build_kernels(fibers, classes_fine, kernel_mode, max_fiber_size=max_fiber_size)
    p_fine = induce_fine_distribution(p_coarse, pi, kernels)
    recovered = pushforward_from_fine(p_fine, fibers)
    err = max_pushforward_error(p_coarse, recovered)

    detail: Dict[str, Any] = {
        "refinement_valid": valid,
        "pi_fine_to_coarse": pi,
        "fibers": {str(s): js for s, js in fibers.items()},
        "kernels_K_s": {
            str(s): {str(j): round(kernels[s][j], 12) for j in sorted(kernels[s].keys())}
            for s in sorted(kernels.keys())
            if s in kernels
        },
        "kernel_stability": kernel_stability_summary(kernels),
        "fiber_sample_size": (int(max_fiber_size) if max_fiber_size is not None else None),
        "pushforward_max_error": err,
        "entropy_coarse_bits": entropy(p_coarse, base=2.0) if p_coarse else 0.0,
        "entropy_fine_bits": entropy(p_fine, base=2.0) if p_fine else 0.0,
    }
    return p_fine, detail


def run_refinement_chain(
    configs: Sequence[HypergraphConfig],
    levels: Sequence[Tuple[str, int]],
    kernel_mode: KernelMode,
) -> Tuple[List[float], List[Dict[str, Any]], List[List[Set[HypergraphConfig]]]]:
    """
    沿細化鏈 Λ=Λ₀→Λ₁→…→Λ_m 逐步條件化：levels[0] 由粗到細依序。

    初值取配置層均勻在 S_{Λ₀} 上之誘導邊際 p_{Λ₀}。
    """
    if not levels:
        return [], [], []
    partitions = [partition_analytic_units(configs, sig, d)[0] for sig, d in levels]
    p = uniform_distribution_on_configs(partitions[0])
    details: List[Dict[str, Any]] = []
    for i in range(len(levels) - 1):
        p, d = apply_refinement_step(p, partitions[i], partitions[i + 1], kernel_mode)
        details.append(d)
    return p, details, partitions


def analyze_refinement_pair(
    configs: Sequence[HypergraphConfig],
    coarse_sig: str,
    coarse_delta: int,
    fine_sig: str,
    fine_delta: int,
    kernel_mode: KernelMode,
    *,
    coarse_sample_size: int | None = None,
    fine_sample_size: int | None = None,
    sample_seed: int = 7,
    max_fiber_size: int | None = None,
) -> Dict[str, Any]:
    """單步細化：回傳 p_fine、明細與粗／細熵。"""
    cfg_all = list(configs)
    coarse_cfgs, _rq_c, _n_c, _nt_c = subsample_obs_configs(
        cfg_all, coarse_sample_size, seed=int(sample_seed) + 401
    )
    # 細層樣本需為粗層樣本子集，避免投影映射缺鍵。
    fine_cfgs, _rq_f, _n_f, _nt_f = subsample_obs_configs(
        coarse_cfgs, fine_sample_size, seed=int(sample_seed) + 809
    )
    classes_c, _ = partition_analytic_units(coarse_cfgs, coarse_sig, coarse_delta)
    classes_f, _ = partition_analytic_units(fine_cfgs, fine_sig, fine_delta)
    p0 = uniform_distribution_on_configs(classes_c)
    p1, detail = apply_refinement_step(
        p0, classes_c, classes_f, kernel_mode, max_fiber_size=max_fiber_size
    )
    detail["coarse_sample_size_used"] = len(coarse_cfgs)
    detail["fine_sample_size_used"] = len(fine_cfgs)
    return {"p_coarse": p0, "p_fine": p1, "detail": detail}


def _joint_terminal_from_two_steps(
    configs: Sequence[HypergraphConfig],
    sig_a: str,
    delta_a: int,
    sig_b: str,
    delta_b: int,
) -> tuple[
    list[set[HypergraphConfig]],
    dict[HypergraphConfig, int],
    dict[HypergraphConfig, int],
    list[int],
    list[int],
]:
    """
    由兩個子步驟 A/B 生成共同終端單元（交叉細分），確保 A→B 與 B→A 比較在同一終端空間。
    """
    classes_a, idx_a = partition_analytic_units(configs, sig_a, delta_a)
    classes_b, idx_b = partition_analytic_units(configs, sig_b, delta_b)
    pair_to_term: dict[tuple[int, int], int] = {}
    term_classes: list[set[HypergraphConfig]] = []
    pi_term_to_a: list[int] = []
    pi_term_to_b: list[int] = []
    for cfg in configs:
        ia = idx_a[cfg]
        ib = idx_b[cfg]
        k = (ia, ib)
        t = pair_to_term.get(k)
        if t is None:
            t = len(term_classes)
            pair_to_term[k] = t
            term_classes.append(set())
            pi_term_to_a.append(ia)
            pi_term_to_b.append(ib)
        term_classes[t].add(cfg)
    return term_classes, idx_a, idx_b, pi_term_to_a, pi_term_to_b


def analyze_section_10_4_bundle(
    configs: Sequence[HypergraphConfig],
    *,
    coarse_sig: str,
    coarse_delta: int,
    fine_sig: str,
    fine_delta: int,
    kernel_mode: KernelMode,
    compare_chains: bool,
    coarse_sample_size: int | None = None,
    fine_sample_size: int | None = None,
    sample_seed: int = 7,
    max_fiber_size: int | None = None,
) -> Dict[str, Any]:
    """
    §10.4 套件：單步 Λ→Λ′ 推前檢查，以及可選之雙路徑終端比較（細化順序 A→B 與 B→A）。

    雙路徑實作：在三層簽名 (粗,中,細) 固定為 weak→medium→strong 時，
    比較 (weak→medium)→strong 與 weak→(medium→strong) 之終端細層鍵分布（對齊後 JS）。
    """
    out: Dict[str, Any] = {}
    try:
        bundle = analyze_refinement_pair(
            configs,
            coarse_sig,
            coarse_delta,
            fine_sig,
            fine_delta,
            kernel_mode,
            coarse_sample_size=coarse_sample_size,
            fine_sample_size=fine_sample_size,
            sample_seed=sample_seed,
            max_fiber_size=max_fiber_size,
        )
        d = bundle["detail"]
        out["single_step_Lambda_to_Lambda_prime"] = {
            "pushforward_max_error": d.get("pushforward_max_error"),
            "H_p_Lambda_bits": d.get("entropy_coarse_bits"),
            "H_p_Lambda_prime_bits": d.get("entropy_fine_bits"),
            "refinement_valid": d.get("refinement_valid"),
            "detail": d,
        }
    except Exception as e:
        out["single_step_Lambda_to_Lambda_prime"] = {"error": str(e)}

    if compare_chains and len(configs) > 0:
        try:
            cmp_out = compare_ordered_refinement_paths(
                configs,
                step_a=preset_layer_to_substep_spec(coarse_sig, coarse_delta),
                step_b=preset_layer_to_substep_spec(fine_sig, fine_delta),
                kernel_mode=kernel_mode,
                coarse_sample_size=coarse_sample_size,
                fine_sample_size=fine_sample_size,
                sample_seed=sample_seed,
                max_fiber_size=max_fiber_size,
            )
            out["two_path_ab_ba_10_4_7"] = {
                "js_divergence_bits_terminal": cmp_out.get("js_divergence_bits_terminal_ab_ba"),
                "entropy_abs_diff_terminal": cmp_out.get("entropy_abs_diff_terminal_ab_ba"),
                "A_to_B": cmp_out.get("A_to_B"),
                "B_to_A": cmp_out.get("B_to_A"),
            }
        except Exception as e:
            out["two_path_ab_ba_10_4_7"] = {"error": str(e)}
    return round_floats_for_output(out)


def compare_ordered_refinement_paths(
    configs: Sequence[HypergraphConfig],
    *,
    step_a: RefinementSubstepSpec | None = None,
    step_b: RefinementSubstepSpec | None = None,
    kernel_mode: KernelMode,
    coarse_sample_size: int | None = None,
    fine_sample_size: int | None = None,
    sample_seed: int = 7,
    max_fiber_size: int | None = None,
    **legacy_kwargs: Any,
) -> Dict[str, Any]:
    """
    進階比較：同一組輸入上計算 A→B 與 B→A 兩條細化鏈之摘要。

    Args:
        step_a: 第一子步驟 R_A 之操作化規格（映射種類 + δ）。
        step_b: 第二子步驟 R_B 之操作化規格。

    Returns:
        每條鏈的推前誤差、熵、纖維核穩定度，以及終端 JS／熵差摘要。
    """
    if legacy_kwargs:
        legacy_keys = sorted(str(k) for k in legacy_kwargs.keys())
        raise ValueError(
            "compare_ordered_refinement_paths 不再接受舊欄位參數。"
            f"收到：{legacy_keys}。請改用 step_a 與 step_b（格式：{{'kind': ..., 'delta': ...}}）。"
        )
    if step_a is None or step_b is None:
        raise ValueError(
            "compare_ordered_refinement_paths 需要 step_a 與 step_b。"
            "請傳入格式：step_a={'kind': ..., 'delta': ...}, step_b={'kind': ..., 'delta': ...}。"
        )
    step_a = validate_substep_spec(step_a, label="step_a")
    step_b = validate_substep_spec(step_b, label="step_b")
    out: Dict[str, Any] = {}
    sig_a, delta_a = substep_spec_to_signature_delta(step_a)
    sig_b, delta_b = substep_spec_to_signature_delta(step_b)
    cfg_all = list(configs)
    cfg_a, _rq_a, _n_a, _nt_a = subsample_obs_configs(
        cfg_all, coarse_sample_size, seed=int(sample_seed) + 1301
    )
    cfg_obs, _rq_b, _n_b, _nt_b = subsample_obs_configs(
        cfg_a, fine_sample_size, seed=int(sample_seed) + 1709
    )
    if not cfg_obs:
        return round_floats_for_output(
            {
                "A_to_B": {"path_key": "A→B", "step_a": dict(step_a), "step_b": dict(step_b), "error": "empty_configs"},
                "B_to_A": {"path_key": "B→A", "step_a": dict(step_a), "step_b": dict(step_b), "error": "empty_configs"},
                "js_divergence_bits_terminal_ab_ba": None,
                "entropy_abs_diff_terminal_ab_ba": 0.0,
            }
        )

    term_classes, idx_a, idx_b, pi_t2a, pi_t2b = _joint_terminal_from_two_steps(
        cfg_obs, sig_a, delta_a, sig_b, delta_b
    )
    classes_a, _ = partition_analytic_units(cfg_obs, sig_a, delta_a)
    classes_b, _ = partition_analytic_units(cfg_obs, sig_b, delta_b)
    p_a = uniform_distribution_on_configs(classes_a)
    p_b = uniform_distribution_on_configs(classes_b)

    fibers_a = fibers_from_pi(pi_t2a, len(classes_a))
    fibers_b = fibers_from_pi(pi_t2b, len(classes_b))
    kernels_a = build_kernels(fibers_a, term_classes, kernel_mode, max_fiber_size=max_fiber_size)
    kernels_b = build_kernels(fibers_b, term_classes, kernel_mode, max_fiber_size=max_fiber_size)

    p_ab = induce_fine_distribution(p_a, pi_t2a, kernels_a)
    p_ba = induce_fine_distribution(p_b, pi_t2b, kernels_b)
    rec_a = pushforward_from_fine(p_ab, fibers_a)
    rec_b = pushforward_from_fine(p_ba, fibers_b)
    err_ab = max_pushforward_error(p_a, rec_a)
    err_ba = max_pushforward_error(p_b, rec_b)
    js_cross = js_divergence_bits(p_ab, p_ba)
    h_ab = entropy(p_ab, base=2.0) if p_ab else 0.0
    h_ba = entropy(p_ba, base=2.0) if p_ba else 0.0
    out["A_to_B"] = {
        "path_key": "A→B",
        "step_a": dict(step_a),
        "step_b": dict(step_b),
        "pushforward_max_error": err_ab,
        "entropy_fine_bits": h_ab,
        "kernel_stability": kernel_stability_summary(kernels_a),
        "terminal_dim": len(p_ab),
    }
    out["B_to_A"] = {
        "path_key": "B→A",
        "step_a": dict(step_a),
        "step_b": dict(step_b),
        "pushforward_max_error": err_ba,
        "entropy_fine_bits": h_ba,
        "kernel_stability": kernel_stability_summary(kernels_b),
        "terminal_dim": len(p_ba),
    }
    out["js_divergence_bits_terminal_ab_ba"] = None if js_cross is None else round(js_cross, 8)
    out["entropy_abs_diff_terminal_ab_ba"] = round(abs(float(h_ab) - float(h_ba)), 8)
    return round_floats_for_output(out)
