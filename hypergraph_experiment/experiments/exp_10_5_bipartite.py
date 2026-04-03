"""
§10.5 解析不可分解性與非可分機率：二分載體、跨區塊邊、ρ_irred、D_sep、I(A;B)。

「可分」之操作化：同一解析單元內，若兩配置之 (Sig_A, Sig_B) 相同則 Sig_AB 必須相同；
違反者視為解析不可分解單元。詳見《約束世界論 30》§10.5.2（四）。
"""

from __future__ import annotations

import math
import random
from collections import Counter
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, cast

import pandas as pd

from hypergraph_experiment.core import (
    OBS_SUBSAMPLE_RNG_OFFSET_10_5,
    SIGNATURES,
    HypergraphConfig,
    all_possible_hyperedges,
    satisfies_domain_constraints,
    subsample_obs_configs,
    tolerance_equivalence_classes,
)
from hypergraph_experiment.output_rounding import round_floats_for_output
from hypergraph_experiment.progress_types import ProgressCallback
from hypergraph_experiment.refinement import js_divergence_bits


def _restriction(c: HypergraphConfig, verts: Set[int]) -> HypergraphConfig:
    """限制在頂點子集上之導出子超圖（僅保留完全含於 verts 之超邊）。"""
    sub_v = tuple(sorted(v for v in c.vertices if v in verts))
    es = frozenset(e for e in c.hyperedges if e.issubset(verts))
    return HypergraphConfig(vertices=sub_v, hyperedges=es)


def _sig_ab_summary(c: HypergraphConfig, n_a: int) -> Tuple[Any, ...]:
    """整體簽名 Sig_AB：強簽名 + 跨區塊邊摘要。"""
    va = frozenset(range(1, n_a + 1))
    vb = frozenset(v for v in c.vertices if v > n_a)
    cross = [e for e in c.hyperedges if not (e.issubset(va) or e.issubset(vb))]
    cross_part = (len(cross), tuple(sorted(len(e) for e in cross)))
    s_full = SIGNATURES["strong"](c)
    return (s_full, cross_part)


def _classify_configs_entropic(
    configs: Sequence[HypergraphConfig], delta: int
) -> List[Set[HypergraphConfig]]:
    """以整體強簽名 + dist 建立與 core 相容之解析單元（近似 q_Λ^ent）。"""

    def sig_ent(c: HypergraphConfig) -> object:
        return SIGNATURES["strong"](c)

    # δ=0 時 core 內改為分桶 O(n)；δ>0 仍為 O(n²·L) 級
    return tolerance_equivalence_classes(list(configs), sig_ent, delta)


def _joint_ab_keys(configs: Sequence[HypergraphConfig], n_a: int) -> Dict[Tuple[Any, ...], int]:
    """(Sig_A, Sig_B) 聯合計數。"""
    joint: Counter[Tuple[Any, ...]] = Counter()
    va = frozenset(range(1, n_a + 1))
    for c in configs:
        vb = frozenset(v for v in c.vertices if v > n_a)
        ca = _restriction(c, va)
        cb = _restriction(c, vb)
        ka = SIGNATURES["medium"](ca)
        kb = SIGNATURES["medium"](cb)
        joint[(ka, kb)] += 1
    return dict(joint)


def _tv_distance(p: List[float], q: List[float]) -> float:
    return 0.5 * sum(abs(a - b) for a, b in zip(p, q))


def _entropy_bits(probs: Sequence[float]) -> float:
    return -sum(p * math.log(p, 2) for p in probs if p > 0)


def _mutual_intrinsic_bits(joint_counts: Dict[Tuple[Any, ...], int], total: int) -> float:
    """由 (Sig_A,Sig_B) 聯合計數求 I(A;B)（bits）。"""
    if total <= 0:
        return 0.0
    p_ab = [cnt / total for cnt in joint_counts.values()]
    h_ab = _entropy_bits(p_ab)
    pa: Counter[Any] = Counter()
    pb: Counter[Any] = Counter()
    for (ka, kb), cnt in joint_counts.items():
        pa[ka] += cnt
        pb[kb] += cnt
    h_a = _entropy_bits([c / total for c in pa.values()])
    h_b = _entropy_bits([c / total for c in pb.values()])
    return max(0.0, h_a + h_b - h_ab)


def section_10_5_output_parameters_df(result: dict[str, Any]) -> pd.DataFrame:
    """
    組裝《約束世界論 30》§10.5.5「輸出參數」之主表（與扁平化寬表並列，供匯出／複製）。

    列舉（一）ρ_irred 與其分子分母、（二）ρ_cross、（三）D_sep（TV 為正文可固定之主距離；JS 表列為補充）、（四）I(A;B)。
    若無合法可採用配置或 ``metrics`` 含錯誤鍵，則回傳單列提醒表。

    Args:
        result: ``run_experiment_10_5`` 之回傳字典。

    Returns:
        五欄 ``DataFrame``：論文小節、輸出參數、論文記號、數值、論文語義摘要。
    """
    cols = ("論文小節", "輸出參數", "論文記號", "數值", "論文語義摘要")
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

    n_cls = result.get("num_classes")
    n_irred = metrics.get("N_irred")
    rho_ir = metrics.get("rho_irred")
    rho_c = metrics.get("rho_cross_mean")
    d_tv = metrics.get("D_sep_total_variation")
    d_js = metrics.get("D_sep_JS_bits")
    i_ab = metrics.get("I_A_B_bits")

    def _fmt(v: object) -> str:
        if v is None:
            return "—"
        if isinstance(v, float):
            return f"{v:.6g}"
        return str(v)

    rows: List[dict[str, str]] = [
        {
            "論文小節": "10.5.5（一）",
            "輸出參數": "解析不可分解比例",
            "論文記號": r"$\rho_{\mathrm{irred}}=N_{\mathrm{irred}}/|S_\Lambda|$",
            "數值": _fmt(rho_ir),
            "論文語義摘要": "整體解析型別中無法被局部直積吸收的單元占比。",
        },
        {
            "論文小節": "10.5.5（一）",
            "輸出參數": "不可分解單元數（分子）",
            "論文記號": r"$N_{\mathrm{irred}}$",
            "數值": _fmt(n_irred),
            "論文語義摘要": "解析不可分解單元之個數。",
        },
        {
            "論文小節": "10.5.5（一）",
            "輸出參數": "解析單元數（分母）/觀測集",
            "論文記號": r"$|S_\Lambda|$／$|\mathcal C_{obs}|$",
            "數值": _fmt(n_cls),
            "論文語義摘要": "觀測集上之解析單元個數（程式以 $N_{cfg}$ 抽樣後之合法樣本做統計）。",
        },
        {
            "論文小節": "10.5.5（二）",
            "輸出參數": "跨區塊耦合密度",
            "論文記號": r"$\rho_{\mathrm{cross}}=m_{\mathrm{cross}}/m$",
            "數值": _fmt(rho_c),
            "論文語義摘要": "結構驅動量；實作取各解析單元代表配置之 m_cross/m 再算術平均。",
        },
        {
            "論文小節": "10.5.5（三）",
            "輸出參數": "非可分偏離量（總變差，正文主選）",
            "論文記號": r"$D_{\mathrm{sep}}$（TV）",
            "數值": _fmt(d_tv),
            "論文語義摘要": r"§10.5.5（三）：$D(p,\;p_A\otimes p_B)$ 可固定為 TV；本列為主要報告欄。",
        },
        {
            "論文小節": "10.5.5（三）",
            "輸出參數": "非可分偏離量（JS，程式補充）",
            "論文記號": r"$D_{\mathrm{sep}}$（JS）",
            "數值": _fmt(d_js),
            "論文語義摘要": "論文正文建議只固定一種距離；JS 僅供對照（bits），非必報主表欄。",
        },
        {
            "論文小節": "10.5.5（四）",
            "輸出參數": "區塊互資訊",
            "論文記號": r"$I(A;B)$",
            "數值": _fmt(i_ab),
            "論文語義摘要": "由 (Sig_A,Sig_B) 聯合計數求得之互資訊（bits）。",
        },
    ]
    return pd.DataFrame(rows, columns=list(cols))


def generate_bipartite_candidates(
    n_a: int,
    n_b: int,
    m_edges: int,
    *,
    max_edge_size: int,
    k_min: int,
    seed: int,
    alpha_cross: float,
    sample_limit: int,
    progress: ProgressCallback = None,
) -> List[HypergraphConfig]:
    """
    在 V_A ⊔ V_B 上隨機產生 m_edges 條超邊；以 alpha_cross 控制跨區塊邊比例傾向。
    """
    n = n_a + n_b
    vertices = tuple(range(1, n + 1))
    va = frozenset(range(1, n_a + 1))
    vb = frozenset(range(n_a + 1, n + 1))
    pool = all_possible_hyperedges(vertices, k_min, max_edge_size)
    pool_local = [e for e in pool if e.issubset(va) or e.issubset(vb)]
    pool_cross = [e for e in pool if not (e.issubset(va) or e.issubset(vb))]
    rng = random.Random(seed)
    out: List[HypergraphConfig] = []
    tgt = max(1, sample_limit)
    for trial in range(sample_limit):
        if progress and trial % max(1, tgt // 20) == 0:
            progress(min(trial + 1, tgt), tgt, f"§10.5 候選 {trial + 1}/{tgt}")
        edges_set: Set[frozenset] = set()
        attempts = 0
        while len(edges_set) < m_edges and attempts < m_edges * 30:
            attempts += 1
            if rng.random() < alpha_cross and pool_cross:
                e = rng.choice(pool_cross)
            elif pool_local:
                e = rng.choice(pool_local)
            else:
                e = rng.choice(pool)
            edges_set.add(e)
        out.append(HypergraphConfig(vertices=vertices, hyperedges=frozenset(edges_set)))
    if progress:
        progress(tgt, tgt, "§10.5 候選產生完畢")
    return out


def run_experiment_10_5(
    *,
    n_a: int = 6,
    n_b: int = 6,
    m_edges: int = 16,
    max_edge_size: int = 3,
    k_min: int = 2,
    alpha_cross: float = 0.3,
    sample_limit: int = 500,
    n_cfg: Optional[int] = None,
    seed: int = 42,
    delta_ent: int = 0,
    connected: bool = False,
    max_degree: int | None = 6,
    forbid_pair_triangles: bool = False,
    progress: ProgressCallback = None,
) -> dict[str, Any]:
    """
    執行 §10.5 批次統計：跨區塊密度、不可分解比例、D_sep（TV）、I(A;B)。

    Args:
        n_a, n_b: 二分兩側節點數。
        m_edges: 每個候選超邊總數。
        alpha_cross: 生成時採跨區邊之機率傾向。
        sample_limit: 候選樣本數 N_cand（過濾前，論文 §10.5.10）。
        n_cfg: 過濾後用於解析統計之觀測筆數（論文 $N_{\mathrm{cfg}}$）；``None`` 則以全部合法配置為觀測集（向後相容）。
        delta_ent: 解析映射 q_Λ^ent 所使用之簽名距離閾值（整數刻度，與 core 一致）。
            **δ=0**（預設）時等價類分桶為 O(N_adm) 級；**δ>0** 時鄰域建構為 O(N_adm²·L)，
            L 為展平簽名長度——大樣本請優先維持 δ=0 或降低 N_adm。

    Returns:
        可 JSON 序列化之結果字典，含 ``parameters``、``metrics``、``num_admissible_filtered``、
        ``num_obs_configs``、``num_admissible``（同觀測集大小，向後相容）等。
    """
    n = n_a + n_b
    candidates = generate_bipartite_candidates(
        n_a,
        n_b,
        m_edges,
        max_edge_size=max_edge_size,
        k_min=k_min,
        seed=seed,
        alpha_cross=alpha_cross,
        sample_limit=sample_limit,
        progress=progress,
    )
    if progress:
        progress(0, 2, "§10.5 域型過濾…")
    adm_configs = [
        c
        for c in candidates
        if satisfies_domain_constraints(
            c,
            max_edge_size=max_edge_size,
            max_edges=m_edges,
            connected_required=connected,
            max_degree=max_degree,
            forbid_pair_triangles=forbid_pair_triangles,
        )
    ]
    if progress:
        progress(2, 2, "§10.5 解析分類…")

    if not adm_configs:
        return {
            "experiment": "10.5",
            "parameters": {
                "n_a": n_a,
                "n_b": n_b,
                "m_edges": m_edges,
                "k_min": k_min,
                "k_max": max_edge_size,
                "alpha_cross": alpha_cross,
                "sample_limit": sample_limit,
                "n_cfg": n_cfg,
                "seed": seed,
                "delta_ent": delta_ent,
            },
            "num_admissible_filtered": 0,
            "num_obs_configs": 0,
            "num_admissible": 0,
            "n_cfg_requested": int(n_cfg) if n_cfg is not None else None,
            "n_cfg_notice": None,
            "metrics": {"error": "無合法配置，請放寬域型條件或增加樣本數。"},
        }

    obs_cfgs, n_req, n_obs, n_notice = subsample_obs_configs(
        adm_configs,
        n_cfg,
        seed=seed,
        rng_chain_offset=OBS_SUBSAMPLE_RNG_OFFSET_10_5,
    )
    configs = obs_cfgs

    classes = _classify_configs_entropic(configs, delta_ent)
    n_classes = len(classes)
    va = frozenset(range(1, n_a + 1))
    vb = frozenset(range(n_a + 1, n + 1))

    n_irred = 0
    rho_cross_list: List[float] = []
    for cls in classes:
        sig_ab_by_pair: Dict[Tuple[Any, Any], Set[Tuple[Any, ...]]] = {}
        for c in cls:
            ca = _restriction(c, va)
            cb = _restriction(c, vb)
            ka = SIGNATURES["medium"](ca)
            kb = SIGNATURES["medium"](cb)
            sab = _sig_ab_summary(c, n_a)
            sig_ab_by_pair.setdefault((ka, kb), set()).add(sab)
        is_irred = any(len(v) > 1 for v in sig_ab_by_pair.values())
        if is_irred:
            n_irred += 1
        rep = next(iter(cls))
        cross_e = [e for e in rep.hyperedges if not (e.issubset(va) or e.issubset(vb))]
        rho_cross_list.append(len(cross_e) / max(1, len(rep.hyperedges)))

    joint = _joint_ab_keys(configs, n_a)
    keys_sorted = sorted(joint.keys(), key=lambda t: (str(t[0]), str(t[1])))
    n_joint = len(keys_sorted)
    total = float(len(configs))
    p_joint = [joint[k] / total for k in keys_sorted]
    pa_full = Counter()
    pb_full = Counter()
    for k, cnt in joint.items():
        pa_full[k[0]] += cnt
        pb_full[k[1]] += cnt
    pa_aligned = [pa_full[k[0]] / total for k in keys_sorted]
    pb_aligned = [pb_full[k[1]] / total for k in keys_sorted]
    p_prod = [pa_aligned[i] * pb_aligned[i] for i in range(n_joint)]
    d_sep_tv = _tv_distance(p_joint, p_prod)
    js_sep = js_divergence_bits(p_joint, p_prod)
    ia_b = _mutual_intrinsic_bits(dict(joint), int(total))

    rho_irred = n_irred / max(1, n_classes)
    rho_cross_mean = sum(rho_cross_list) / max(1, len(rho_cross_list))

    return round_floats_for_output(
        {
            "experiment": "10.5",
            "parameters": {
                "n_a": n_a,
                "n_b": n_b,
                "m_edges": m_edges,
                "k_min": k_min,
                "k_max": max_edge_size,
                "alpha_cross": alpha_cross,
                "sample_limit": sample_limit,
                "n_cfg": n_cfg,
                "seed": seed,
                "delta_ent": delta_ent,
                "connected": connected,
                "max_degree": max_degree,
                "forbid_pair_triangles": forbid_pair_triangles,
            },
            "num_admissible_filtered": len(adm_configs),
            "num_obs_configs": n_obs,
            "num_admissible": len(configs),
            "n_cfg_requested": n_req,
            "n_cfg_notice": n_notice,
            "num_classes": n_classes,
            "metrics": {
                "rho_irred": round(rho_irred, 6),
                "N_irred": n_irred,
                "rho_cross_mean": round(rho_cross_mean, 6),
                "D_sep_total_variation": round(d_sep_tv, 6),
                "D_sep_JS_bits": None if js_sep is None else round(js_sep, 6),
                "I_A_B_bits": round(ia_b, 6),
            },
        }
    )


def section_10_5_batch_sweep_dataframe(batch_runs: Sequence[dict[str, Any]]) -> pd.DataFrame:
    """
    將 ``run_batch_per_run_rows`` 產出之 §10.5 批次列表彙整為單一 DataFrame，供掃描軸拆線圖使用。

    僅收錄無例外、且 ``metrics`` 不含 ``error`` 之成功列；欄位含 ``run_index``、
    參數表之原始鍵（來自 ``param_row``），以及論文 §10.5.5 對應之指標與診斷欄。

    Args:
        batch_runs: 各元素須具 ``run_index``、``error``、``result``、``param_row`` 等鍵，
            與 ``hypergraph_experiment.streamlit_common.run_batch_per_run_rows`` 回傳格式一致。

    Returns:
        可能為空之 DataFrame；指標欄名為繪圖友善之英文鍵（如 ``rho_irred``、``D_sep_TV``）。
    """
    rows: List[dict[str, Any]] = []
    for item in batch_runs:
        if not isinstance(item, dict):
            continue
        if item.get("error"):
            continue
        res = item.get("result")
        if not isinstance(res, dict):
            continue
        metrics = res.get("metrics")
        if not isinstance(metrics, dict) or metrics.get("error"):
            continue
        pr_raw = item.get("param_row")
        pr = cast(Dict[str, Any], pr_raw) if isinstance(pr_raw, dict) else {}
        row: dict[str, Any] = {"run_index": item.get("run_index")}
        for k, v in pr.items():
            row[str(k)] = v
        row["num_admissible_filtered"] = res.get("num_admissible_filtered")
        row["num_obs_configs"] = res.get("num_obs_configs")
        row["num_classes"] = res.get("num_classes")
        row["rho_irred"] = metrics.get("rho_irred")
        row["rho_cross_mean"] = metrics.get("rho_cross_mean")
        row["D_sep_TV"] = metrics.get("D_sep_total_variation")
        row["D_sep_JS_bits"] = metrics.get("D_sep_JS_bits")
        row["I_A_B_bits"] = metrics.get("I_A_B_bits")
        row["N_irred"] = metrics.get("N_irred")
        rows.append(row)
    return pd.DataFrame(rows)


def section_10_5_sweep_axis_candidates(df: pd.DataFrame) -> List[str]:
    """
    自彙整表挑出「至少兩筆有效值、且跨列有變化」之數值欄，作為掃描拆線圖橫軸候選。

    排除明顯屬輸出指標或序號之欄位，避免誤選。

    Args:
        df: :func:`section_10_5_batch_sweep_dataframe` 之回傳值。

    Returns:
        欄名列表（保留參數表原始鍵名，如中文「跨區傾向」）。
    """
    exclude = {
        "run_index",
        "num_admissible_filtered",
        "num_obs_configs",
        "num_classes",
        "rho_irred",
        "rho_cross_mean",
        "D_sep_TV",
        "D_sep_JS_bits",
        "I_A_B_bits",
        "N_irred",
    }
    out: List[str] = []
    for col in df.columns:
        if col in exclude:
            continue
        ser = pd.to_numeric(df[col], errors="coerce")
        if int(ser.notna().sum()) < 2:
            continue
        if int(ser.nunique(dropna=True)) < 2:
            continue
        out.append(str(col))
    return out
