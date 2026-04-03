"""
§10.6 局部解析族與拼合障礙：循環 parity 主例、χ_glue、ρ_glue、r_min、ρ_val。

在 GF(2) 上將「變數為節點 0..n-1 之位元、約束為 x_i⊕x_j=b」之系統作為全域延拓判準。
"""

from __future__ import annotations

import itertools
import random
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd

from hypergraph_experiment.output_rounding import round_floats_for_output
from hypergraph_experiment.progress_types import ProgressCallback

# 單一約束：(i, j, b) 代表 x_i ⊕ x_j = b（i,j 為 0-based 節點索引）
XorConstraint = Tuple[int, int, int]


def _gf2_add(a: int, b: int) -> int:
    return a ^ b


def system_globally_satisfiable(n: int, equations: Sequence[XorConstraint]) -> bool:
    """以並查集維護 xor 差分；若出現矛盾則不可滿足。"""
    parent: List[int] = list(range(n))
    diff_to_root: List[int] = [0] * n

    def find(x: int) -> Tuple[int, int]:
        if parent[x] != x:
            r, d = find(parent[x])
            parent[x] = r
            diff_to_root[x] = _gf2_add(diff_to_root[x], d)
        return parent[x], diff_to_root[x]

    def union(i: int, j: int, b: int) -> bool:
        ri, di = find(i)
        rj, dj = find(j)
        if ri == rj:
            need = _gf2_add(_gf2_add(di, dj), b)
            return need == 0
        parent[rj] = ri
        diff_to_root[rj] = _gf2_add(di, _gf2_add(dj, b))
        return True

    for i, j, b in equations:
        if not union(i, j, b):
            return False
    return True


def minimal_conflict_subset_size(n: int, equations: Sequence[XorConstraint]) -> Optional[int]:
    """若不可滿足，回傳導致矛盾之最小子約束數；可滿足則 None。"""
    eqs = list(equations)
    m = len(eqs)
    if system_globally_satisfiable(n, eqs):
        return None
    for k in range(1, m + 1):
        for subset in itertools.combinations(range(m), k):
            sub = [eqs[t] for t in subset]
            if not system_globally_satisfiable(n, sub):
                return k
    return m


def canonical_cycle_obstruction() -> Tuple[int, List[XorConstraint]]:
    """
    論文 §10.6.2 四窗口主例：x0⊕x1=0, x1⊕x2=0, x2⊕x3=0, x3⊕x0=1。
    局部兩兩可相容，整體不可滿足。
    """
    n = 4
    eqs = [(0, 1, 0), (1, 2, 0), (2, 3, 0), (3, 0, 1)]
    return n, eqs


def validate_cyclic_window_params(n_nodes: int, M: int, w: int, eta: int) -> Optional[str]:
    """
    檢查 §10.6.4 循環 M 視窗幾何是否可在 n 個節點上嵌入。

    Returns:
        若非法則為繁中錯誤字串，否則 ``None``。
    """
    if n_nodes < 2:
        return "節點數 n 至少為 2。"
    if M < 2:
        return "上下文數 M 至少為 2（循環視窗族）。"
    if w < 2:
        return "窗口大小 w 至少為 2。"
    if eta < 1 or eta >= w:
        return "交疊 η 需滿足 1 ≤ η < w（§10.6.4（三）與（二））。"
    step = w - eta
    if step < 1:
        return "需滿足 w > η。"
    ring_len = M * step
    if ring_len > n_nodes:
        return (
            f"循環嵌入長度 M(w−η)={ring_len} 超過節點數 n={n_nodes}；"
            "請增大 n 或降低 M／w 或提高 η。"
        )
    if w == 2 and eta == 1 and M == 2:
        return (
            "當 w=2、η=1 時，M=2 會使兩視窗對應同一直線段之反覆；論文 §10.6.4 建議 M≥4。"
        )
    return None


def minimal_conflict_context_subset_size(
    n: int,
    context_groups: Sequence[Sequence[XorConstraint]],
) -> Optional[int]:
    """
    以「上下文」為單位（每組對應一個視窗之局部方程組）求最小矛盾子族大小。

    對應論文 §10.6.5（三）r_min 之「上下文子族」操作化（非單條 XOR）。
    """
    groups: List[List[XorConstraint]] = [list(g) for g in context_groups]
    flat_all: List[XorConstraint] = [e for g in groups for e in g]
    if not flat_all:
        return None
    if system_globally_satisfiable(n, flat_all):
        return None
    mctx = len(groups)
    for k in range(1, mctx + 1):
        for subset in itertools.combinations(range(mctx), k):
            merged: List[XorConstraint] = [e for t in subset for e in groups[t]]
            if not system_globally_satisfiable(n, merged):
                return k
    return mctx


def sample_cyclic_window_xor_groups(
    rng: random.Random,
    n_nodes: int,
    M: int,
    w: int,
    eta: int,
    *,
    mode: str,
    t_loc: int,
) -> List[List[XorConstraint]]:
    """
    產生 §10.6.6（一）（二）之循環交疊視窗：每窗 w 個節點、相鄰窗交疊 η，
    每窗內為樹狀 parity 鏈（w−1 條 XOR）。

    Args:
        rng: 隨機源。
        n_nodes: 全域節點數 n。
        M, w, eta: 論文 §10.6.4 之上下文數、窗口大小、交疊大小。
        mode: ``obstruction`` 或 ``satisfiable``。
        t_loc: 局部條件型別數 T_loc（目前僅 parity；用於位元抽樣之輕微摺疊，避免假掃描完全等價）。

    Returns:
        長度 M 之串列，每元素為該上下文內之 XOR 方程列表。
    """
    err = validate_cyclic_window_params(n_nodes, M, w, eta)
    if err:
        raise ValueError(err)
    step = w - eta
    ring_len = M * step
    cycle_vertices = rng.sample(range(n_nodes), k=ring_len)
    bit_mix = (int(t_loc) - 2) & 1

    groups: List[List[XorConstraint]] = []
    for k in range(M):
        start = k * step
        verts = [cycle_vertices[(start + j) % ring_len] for j in range(w)]
        local: List[XorConstraint] = []
        for j in range(w - 1):
            b_bit = rng.randint(0, 1) ^ bit_mix
            local.append((verts[j], verts[j + 1], b_bit))
        groups.append(local)

    # 正文 parity 主線 w=2：每窗一條邊，可閉合環上 XOR 使整體可滿足（對照 mode）
    if mode == "satisfiable" and w == 2 and M >= 2:
        xor_sum = 0
        for k in range(M - 1):
            xor_sum ^= groups[k][0][2]
        u, v, _ = groups[M - 1][0]
        groups[M - 1][0] = (u, v, xor_sum)

    flat = [e for g in groups for e in g]
    if mode == "obstruction":
        # 在仍可滿足時逐條翻轉 b，直到整體不可滿足
        _obs_done = False
        for gi in range(len(groups)):
            for ej in range(len(groups[gi])):
                if not system_globally_satisfiable(n_nodes, flat):
                    _obs_done = True
                    break
                u, v, b0 = groups[gi][ej]
                groups[gi][ej] = (u, v, b0 ^ 1)
                flat = [e for g in groups for e in g]
            if _obs_done:
                break
        flat = [e for g in groups for e in g]
        # 路徑複雜時仍可能可滿足：後備以 x0⊕x0=1 強制矛盾（理論上極少觸發）
        if system_globally_satisfiable(n_nodes, flat):
            groups = [list(g) for g in groups]
            groups.append([(0, 0, 1)])

    return groups


def section_10_6_output_parameters_df(result: dict[str, Any]) -> pd.DataFrame:
    r"""
    組裝《約束世界論 30》§10.6.5「輸出參數」主表（與扁平化寬表並列，供匯出／複製）。

    列舉（一）$\chi_{\mathrm{glue}}$ 之批量計數語意、（二）$\rho_{\mathrm{glue}}$、
    （三）$r_{\min}$ 及其分布摘要、（四）$\rho_{\mathrm{val}}$。
    若 ``experiment`` 為 ``10.6_canonical``，則改為單例主例之三列簡表。

    Args:
        result: ``run_experiment_10_6`` 或 ``run_canonical_demo_10_6`` 之回傳字典。

    Returns:
        五欄 ``DataFrame``：論文小節、輸出參數、論文記號、數值、論文語義摘要。
    """
    cols = ("論文小節", "輸出參數", "論文記號", "數值", "論文語義摘要")
    exp = result.get("experiment")

    if exp == "10.6_canonical":
        chi = result.get("chi_glue")
        rm = result.get("r_min")
        sat = result.get("globally_satisfiable")
        rows_c: List[dict[str, str]] = [
            {
                "論文小節": "10.6.5（一）",
                "輸出參數": "不可延拓指標（論文主例單族）",
                "論文記號": r"$\chi_{\mathrm{glue}}$",
                "數值": str(int(chi)) if chi is not None else "—",
                "論文語義摘要": "循環 parity 主例：1 表示不可全域延拓（拼合障礙）。",
            },
            {
                "論文小節": "10.6.5（三）",
                "輸出參數": "最小衝突核大小",
                "論文記號": r"$r_{\min}$",
                "數值": str(rm) if rm is not None else "—",
                "論文語義摘要": "為導致全域矛盾之最小（上下文）子族大小；主例應為 4。",
            },
            {
                "論文小節": "10.6.5（四）",
                "輸出參數": "共同賦值可行（單例對應）",
                "論文記號": r"$\rho_{\mathrm{val}}$ 之 0/1 對照",
                "數值": ("1" if sat is True else "0" if sat is False else "—"),
                "論文語義摘要": r"單一例上：1 表全域共同賦值存在；與 $\chi_{\mathrm{glue}}$ 互斥。",
            },
        ]
        return pd.DataFrame(rows_c, columns=list(cols))

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

    def _fmt(v: object) -> str:
        if v is None:
            return "—"
        if isinstance(v, float):
            return f"{v:.6g}"
        return str(v)

    params = result.get("parameters") or {}
    n_req_raw = params.get("n_ctx")
    n_eff = result.get("n_ctx_effective")
    sum_chi = result.get("sum_chi_glue")
    if n_eff is None and n_req_raw is not None:
        n_eff = int(n_req_raw)
    if sum_chi is None and n_eff is not None and metrics.get("rho_glue") is not None:
        sum_chi = int(round(float(metrics["rho_glue"]) * int(n_eff)))
    n_eff_i = int(n_eff) if n_eff is not None else 0
    sum_chi_i = int(sum_chi) if sum_chi is not None else 0

    rho_g = metrics.get("rho_glue")
    rho_v = metrics.get("rho_val")
    n_val_m = metrics.get("n_val")
    if n_val_m is None and rho_v is not None and n_eff_i > 0:
        n_val_m = int(round(float(rho_v) * n_eff_i))
    r_mean = metrics.get("r_min_mean")
    r_hist = metrics.get("r_min_hist")
    hist_str = "—"
    if isinstance(r_hist, dict) and r_hist:
        hist_str = "；".join(f"{k}→{v}" for k, v in sorted(r_hist.items(), key=lambda t: str(t[0])))

    chi_count_str: str
    if n_eff_i > 0:
        chi_count_str = f"{sum_chi_i}／{n_eff_i}（χ=1 的樣本數／有效樣本數）"
    else:
        chi_count_str = "—（有效樣本數為 0，請提高 n_nodes 或檢查嵌入條件）"

    rows: List[dict[str, str]] = [
        {
            "論文小節": "10.6.5（一）",
            "輸出參數": "不可延拓指標（逐樣本計數）",
            "論文記號": r"$\chi_{\mathrm{glue}}^{(i)}\in\{0,1\}$",
            "數值": chi_count_str,
            "論文語義摘要": "§10.6.5（一）：每個上下文族樣本之不可延拓指示；右欄給出 χ=1 之筆數與分母。",
        },
        {
            "論文小節": "10.6.5（二）",
            "輸出參數": "拼合障礙比例",
            "論文記號": r"$\rho_{\mathrm{glue}}=\frac{1}{N_{\mathrm{ctx}}}\sum_i \chi_{\mathrm{glue}}^{(i)}$",
            "數值": _fmt(rho_g),
            "論文語義摘要": "批量平均，與（一）之計數比一致（若分母為有效樣本數則以程式實際分母為準）。",
        },
        {
            "論文小節": "10.6.5（三）",
            "輸出參數": "最小衝突核大小（樣本平均）",
            "論文記號": r"$r_{\min}$",
            "數值": _fmt(r_mean),
            "論文語義摘要": "§10.6.5（三）：以「上下文」為單位之最小矛盾子族大小（程式操作化）；僅 χ=1 樣本有值。",
        },
        {
            "論文小節": "10.6.5（三）",
            "輸出參數": "最小衝突核大小（分布，程式補充）",
            "論文記號": r"$r_{\min}$ 直方",
            "數值": hist_str,
            "論文語義摘要": "子族大小 k 之次數；便於觀察矛盾是否多由極小核引起。",
        },
        {
            "論文小節": "10.6.5（四）",
            "輸出參數": "共同賦值可行樣本數",
            "論文記號": r"$N_{\mathrm{val}}$",
            "數值": _fmt(n_val_m),
            "論文語義摘要": "分子：存在全域共同賦值之樣本數；分母請對照有效樣本數或請求之 N_ctx。",
        },
        {
            "論文小節": "10.6.5（四）",
            "輸出參數": "共同賦值可行率",
            "論文記號": r"$\rho_{\mathrm{val}}=N_{\mathrm{val}}/N_{\mathrm{ctx}}$",
            "數值": _fmt(rho_v),
            "論文語義摘要": "§10.6.5（四）：ρ_val 之分母於程式為有效樣本數（通常等於請求之 N_ctx）。",
        },
    ]
    return pd.DataFrame(rows, columns=list(cols))


def run_experiment_10_6(
    *,
    n_ctx: int = 200,
    n_nodes: int = 8,
    mode: str = "obstruction",
    seed: int = 1,
    M: int = 4,
    w: int = 2,
    eta: int = 1,
    t_loc: int = 2,
    n_search: int = 5000,
    progress: ProgressCallback = None,
) -> dict[str, Any]:
    r"""
    §10.6 批次：依循環 M 視窗 parity 系統估計 $\rho_{\mathrm{glue}}$、$\rho_{\mathrm{val}}$ 等。

    Args:
        n_ctx: 論文 $N_{\mathrm{ctx}}$，上下文族樣本數。
        n_nodes: 論文 $n$，位元變數個數（節點 0..n-1）。
        mode: ``obstruction`` 在嵌入後翻轉 XOR 使整體不可滿足（對照 §10.6.2）；
            ``satisfiable`` 不強制注入矛盾。
        seed: 偽隨機種子。
        M: §10.6.4（一）上下文數。
        w: §10.6.4（二）窗口大小。
        eta: §10.6.4（三）交疊大小 $\eta$。
        t_loc: §10.6.4（四）$T_{\mathrm{loc}}$；目前僅 parity，此參數小幅改變位元抽樣摺疊。
        n_search: §10.6.6 建議表 $N_{\mathrm{search}}$；GF(2) 合併為多項式時間，
            本欄保留與論文對齊並寫入輸出（未作硬性截斷步數）。

    Returns:
        含完整 ``parameters``、``n_ctx_effective``、``sum_chi_glue``、``metrics.n_val`` 等。
    """
    geo_err = validate_cyclic_window_params(n_nodes, int(M), int(w), int(eta))
    t_loc_i = int(max(2, min(4, t_loc)))
    base_params: Dict[str, Any] = {
        "n_ctx": int(n_ctx),
        "n_nodes": int(n_nodes),
        "mode": str(mode),
        "seed": int(seed),
        "M": int(M),
        "w_ctx": int(w),
        "eta_ctx": int(eta),
        "T_loc": t_loc_i,
        "n_search": int(n_search),
    }
    if geo_err:
        return round_floats_for_output(
            {
                "experiment": "10.6",
                "parameters": base_params,
                "n_ctx_requested": int(n_ctx),
                "n_ctx_effective": 0,
                "sum_chi_glue": 0,
                "metrics": {"error": geo_err},
                "canonical_demo": {"chi_glue": 1, "globally_sat": False, "r_min": 4},
            }
        )

    rng = random.Random(int(seed))
    chi_list: List[int] = []
    r_mins: List[int] = []
    val_ok: List[int] = []

    for t in range(n_ctx):
        if progress and t % max(1, n_ctx // 25) == 0:
            progress(t + 1, n_ctx, f"§10.6 樣本 {t + 1}/{n_ctx}")
        try:
            groups = sample_cyclic_window_xor_groups(
                rng,
                int(n_nodes),
                int(M),
                int(w),
                int(eta),
                mode=str(mode),
                t_loc=t_loc_i,
            )
        except ValueError as e:
            return round_floats_for_output(
                {
                    "experiment": "10.6",
                    "parameters": base_params,
                    "n_ctx_requested": int(n_ctx),
                    "n_ctx_effective": len(chi_list),
                    "sum_chi_glue": int(sum(chi_list)),
                    "metrics": {"error": str(e)},
                    "canonical_demo": {"chi_glue": 1, "globally_sat": False, "r_min": 4},
                }
            )
        flat = [e for g in groups for e in g]
        sat = system_globally_satisfiable(int(n_nodes), flat)
        chi = 0 if sat else 1
        chi_list.append(chi)
        if chi == 1:
            rm_ctx = minimal_conflict_context_subset_size(int(n_nodes), groups)
            if rm_ctx is not None:
                r_mins.append(rm_ctx)
        val_ok.append(1 if sat else 0)

    n_ctx_effective = len(chi_list)
    sum_chi_glue = sum(chi_list)
    sum_n_val = sum(val_ok)
    rho_glue = sum_chi_glue / max(1, n_ctx_effective)
    rho_val = sum_n_val / max(1, len(val_ok))

    return round_floats_for_output(
        {
            "experiment": "10.6",
            "parameters": base_params,
            "n_ctx_requested": int(n_ctx),
            "n_ctx_effective": int(n_ctx_effective),
            "sum_chi_glue": int(sum_chi_glue),
            "metrics": {
                "rho_glue": round(rho_glue, 6),
                "rho_val": round(rho_val, 6),
                "n_val": int(sum_n_val),
                "r_min_mean": round(sum(r_mins) / len(r_mins), 4) if r_mins else None,
                "r_min_hist": dict((str(k), r_mins.count(k)) for k in sorted(set(r_mins))) if r_mins else {},
            },
            "canonical_demo": {"chi_glue": 1, "globally_sat": False, "r_min": int(M)},
        }
    )


def run_canonical_demo_10_6() -> dict[str, Any]:
    """回傳論文最小例之 χ_glue、可滿足性與 r_min。"""
    n, eqs = canonical_cycle_obstruction()
    sat = system_globally_satisfiable(n, eqs)
    rm = minimal_conflict_subset_size(n, eqs)
    return round_floats_for_output(
        {
            "experiment": "10.6_canonical",
            "n": n,
            "equations": list(eqs),
            "chi_glue": 0 if sat else 1,
            "globally_satisfiable": sat,
            "r_min": rm,
        }
    )
