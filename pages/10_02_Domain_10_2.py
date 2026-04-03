"""
§10.2 配置域與域型約束收縮（表 10-2）。
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

from hypergraph_experiment.ch10_symbol_glossary import (
    L10_CAND_ENUM,
    L10_CONNECTED,
    L10_D_MAX,
    L10_FORBID_TRI,
    L10_K_MIN,
    L10_K_MAX,
    L10_M_MAX,
    L10_N,
    L10_N_SEED_BATCHES,
    L10_SAMPLE_LIMIT,
    L10_SEED,
    render_section_glossary,
)
from hypergraph_experiment.paper_tables import exhaustive_candidate_count, table_10_2_domain_ladder
from hypergraph_experiment.streamlit_common import (
    batch_cell,
    build_ch10_column_name_map,
    render_batch_per_run_tables,
    render_parameters_table,
    render_sidebar_warehouse,
    render_table_with_copy_csv,
    run_batch_per_run_rows,
    streamlit_progress_callback,
)

st.set_page_config(page_title="§10.2 配置域", layout="wide")
render_sidebar_warehouse()

st.title("§10.2　配置域與域型約束收縮")
render_section_glossary(st, "10.2")
st.markdown(
    r"""
    在固定候選母集 $\mathcal{C}_{cand}$ 上對嵌套 **$\Lambda_{\mathrm{dom}}^{(i)}$** 逐層篩選，觀察
    **$|\mathrm{Cfg}_{\Lambda^{(i)}}|$**（合法配置數）與保留比例（對應論文表 10-2）。
    符號對照：$\mathcal{C}_{cand}$=候選母集、$\Lambda_{\mathrm{dom}}$=域型約束、$\mathrm{Cfg}_\Lambda$=合法配置域。
    **建議**（§10.2.6）：$k_{\min}=2$，$k_{\max}=3$，$N_{\mathrm{cand}}\approx 5000$；
    $n\in\{6,8,10,12\}$，掃描弱／中／強域型約束。
    """
)

with st.expander("enum｜枚舉規模預估（sample_limit=0 等價近似）"):
    st.caption(f"{L10_CAND_ENUM}；枚舉總數 $\\sum_r \\binom{{M}}{{r}}$（§10.2 候選母集）")
    en = st.number_input(L10_N, 2, 16, 8, key="t2est_n")
    ek_min = st.number_input(L10_K_MIN, 2, 6, 2, key="t2est_kmin")
    ek = st.number_input(L10_K_MAX, 2, 6, 3, key="t2est_k")
    em = st.number_input(L10_M_MAX, 1, 20, 10, key="t2est_m")
    if st.button("估算候選總數", key="t2est_btn"):
        tot, m_e = exhaustive_candidate_count(int(en), int(ek), int(em), int(ek_min))
        st.info(f"M={m_e}，枚舉候選配置總數 ≈ {tot:,}")

st.subheader("建議固定參數（可覆寫）")
c1, c2, c3 = st.columns(3)
with c1:
    n = st.number_input(L10_N, 2, 12, 8, key="t2n1")
    k_min = st.number_input(L10_K_MIN, 2, 6, 2, key="t2kmin1")
    k_max = st.number_input(L10_K_MAX, 2, 6, 3, key="t2k1")
with c2:
    m_max = st.number_input(L10_M_MAX, 1, 20, 10, key="t2m1")
    sl = st.number_input(L10_SAMPLE_LIMIT, 0, 50_000, 5000, key="t2sl1")
with c3:
    sd = st.number_input(L10_SEED, 0, 2_000_000_000, 20, key="t2sd1")
    n_seed = st.number_input(L10_N_SEED_BATCHES, 1, 100, 20, key="t2nseed1")

st.divider()
st.subheader("建議變量掃描參數（可覆寫）")
st.markdown(r"**三層域型約束 $\Lambda_{\mathrm{dom}}^{(i)}$（由弱至強；§10.2.4）**")
level_count = st.number_input("層數｜$L$ — 域型梯子層數（建議 3；§10.2）", 2, 4, 3, key="t2_levels")
lv_rows: list[dict] = []
for idx in range(int(level_count)):
    defaults = [
        ("Λ_dom^(1) 弱", False, 8, False),
        ("Λ_dom^(2) 中", True, 5, False),
        ("Λ_dom^(3) 強", True, 4, True),
        ("Λ_dom^(4) 極強", True, 3, True),
    ][idx]
    a, b, d, e = st.columns([2, 1, 1, 1])
    with a:
        lb = st.text_input("層標籤", value=defaults[0], key=f"l2{idx}a", help="對應論文表中約束層名稱")
    with b:
        co = st.checkbox(
            "2-section 連通｜域型條件",
            value=defaults[1],
            key=f"l2{idx}b",
            help=L10_CONNECTED,
        )
    with d:
        md = st.number_input(L10_D_MAX, 1, 20, defaults[2], key=f"l2{idx}c")
    with e:
        ft = st.checkbox(
            "禁二元△｜forbidden motif",
            value=defaults[3],
            key=f"l2{idx}d",
            help=L10_FORBID_TRI,
        )
    lv_rows.append(
        {
            "label": lb,
            "connected": co,
            "max_degree": int(md),
            "forbid_pair_triangles": ft,
        }
    )

if st.button("執行表 10-2", key="run_t2_once"):
    bar = st.progress(0)
    stt = st.empty()

    def _p(done: int, total: int, msg: str) -> None:
        streamlit_progress_callback(bar, stt)(done, total, msg)

    rows, n_cand = table_10_2_domain_ladder(
        n=int(n),
        min_edge_size=int(k_min),
        max_edge_size=int(k_max),
        max_edges=int(m_max),
        sample_limit=int(sl),
        seed=int(sd),
        num_seeds=int(n_seed),
        levels=lv_rows,
        progress=_p,
    )
    st.session_state["t2_last"] = (rows, n_cand)
    st.session_state["t2_last_params"] = {
        "n": int(n),
        "k_min": int(k_min),
        "k_max": int(k_max),
        "m_max": int(m_max),
        "sample_limit": int(sl),
        "seed": int(sd),
        "n_seed": int(n_seed),
        "levels": lv_rows,
    }
    st.success(f"完成；候選母集大小 $|\\mathcal{{C}}_{{cand}}|={n_cand}$")

if "t2_last" in st.session_state:
    rr, nc = st.session_state["t2_last"]
    st.caption(rf"$|\mathcal{{C}}_{{cand}}|={nc}$；表中合法配置數即 $|\mathrm{{Cfg}}_{{\Lambda^{{(i)}}}}|$")
    df_t2_once = pd.DataFrame(rr)
    render_table_with_copy_csv(
        df_t2_once,
        key_prefix="t2_once",
        csv_filename="table_10_2_single.csv",
        hide_index=True,
        column_name_map=build_ch10_column_name_map(df_t2_once.columns),
    )
    render_parameters_table(
        st.session_state.get("t2_last_params"),
        key_prefix="t2_once",
        csv_filename="table_10_2_single_params.csv",
    )
    st.caption(
        "違規主因欄位：依與篩選程式相同之順序，記錄每筆被排除配置**第一個**未通過之域型條件；"
        "集中度為各主因筆數中最大占比（愈接近 1 表示愈集中於單一違規類型）。"
        "「本層新增排除」僅計入上一層合法、本層始遭排除之配置。"
    )

st.subheader("批次參數表")
st.caption(
    "欄位含固定與掃描參數：節點數、最小與最大超邊大小、最大邊數、候選採樣上限、偽隨機基底種子、種子批次數（論文 $N_{seed}$）、層數，"
    "以及各層域型條件。**偽隨機基底種子**與 **種子批次數（論文 $N_{seed}$）** 不同，見摺疊符號對照。"
)
stop_err = st.checkbox("遇錯即停", value=False)
if st.button("載入論文建議批次模板（§10.2）", key="load_t2_template"):
    st.session_state["t2_batch_template_df"] = pd.DataFrame(
        [
            {
                "節點數": 6,
                "最小超邊階數": 2,
                "最大超邊階數": 3,
                "最大超邊數": 10,
                "候選採樣上限": 5000,
                "偽隨機基底種子": 20,
                "種子批次數（論文N_seed）": 20,
                "約束層數": 3,
                "層1連通": False,
                "層1度上限": 8,
                "層1禁二元三角": False,
                "層2連通": True,
                "層2度上限": 5,
                "層2禁二元三角": False,
                "層3連通": True,
                "層3度上限": 4,
                "層3禁二元三角": True,
                "層4連通": True,
                "層4度上限": 3,
                "層4禁二元三角": True,
            },
            {
                "節點數": 8,
                "最小超邊階數": 2,
                "最大超邊階數": 3,
                "最大超邊數": 10,
                "候選採樣上限": 5000,
                "偽隨機基底種子": 21,
                "種子批次數（論文N_seed）": 20,
                "約束層數": 3,
                "層1連通": False,
                "層1度上限": 8,
                "層1禁二元三角": False,
                "層2連通": True,
                "層2度上限": 5,
                "層2禁二元三角": False,
                "層3連通": True,
                "層3度上限": 4,
                "層3禁二元三角": True,
                "層4連通": True,
                "層4度上限": 3,
                "層4禁二元三角": True,
            },
        ]
    )
default_batch = st.session_state.get(
    "t2_batch_template_df",
    pd.DataFrame(
        [
            {
                "節點數": 8,
                "最小超邊階數": 2,
                "最大超邊階數": 3,
                "最大超邊數": 10,
                "候選採樣上限": 5000,
                "偽隨機基底種子": 20,
                "種子批次數（論文N_seed）": 20,
                "約束層數": 3,
                "層1連通": False,
                "層1度上限": 8,
                "層1禁二元三角": False,
                "層2連通": True,
                "層2度上限": 5,
                "層2禁二元三角": False,
                "層3連通": True,
                "層3度上限": 4,
                "層3禁二元三角": True,
                "層4連通": True,
                "層4度上限": 3,
                "層4禁二元三角": True,
            }
        ]
    ),
)
df_ed = st.data_editor(default_batch, num_rows="dynamic", key="t2_batch")

if st.button("批次執行表 10-2", key="t2_batch_go"):

    def _row(r: pd.Series, prog):
        level_defaults = [
            ("弱", False, 8, False),
            ("中", True, 5, False),
            ("強", True, 4, True),
            ("極強", True, 3, True),
        ]
        lv_count = max(2, min(4, int(batch_cell(r, "約束層數", "level_count", 3))))
        row_levels = []
        for i in range(1, lv_count + 1):
            lb, co, dm, ft = level_defaults[i - 1]
            row_levels.append(
                {
                    "label": lb,
                    "connected": bool(batch_cell(r, f"層{i}連通", f"l{i}_connected", co)),
                    "max_degree": int(batch_cell(r, f"層{i}度上限", f"l{i}_d_max", dm)),
                    "forbid_pair_triangles": bool(batch_cell(r, f"層{i}禁二元三角", f"l{i}_forbid_tri", ft)),
                }
            )
        rows_b, nc = table_10_2_domain_ladder(
            n=int(batch_cell(r, "節點數", "n")),
            min_edge_size=int(batch_cell(r, "最小超邊階數", "k_min", 2)),
            max_edge_size=int(batch_cell(r, "最大超邊階數", "k_max")),
            max_edges=int(batch_cell(r, "最大超邊數", "m_max")),
            sample_limit=int(batch_cell(r, "候選採樣上限", "sample_limit")),
            seed=int(batch_cell(r, "偽隨機基底種子", "seed")),
            num_seeds=int(batch_cell(r, "種子批次數（論文N_seed）", "n_seed", 20)),
            levels=row_levels,
            progress=prog,
        )
        # 與單次「執行表 10-2」相同：多列梯子表（中文欄位），不依批次扁平為寬表
        return {
            "table_10_2_ladder_rows": rows_b,
            "n_candidates": nc,
        }

    st.session_state["t2_batch_runs"] = run_batch_per_run_rows(
        df_ed,
        _row,
        stop_on_error=stop_err,
        use_progress=True,
    )


def _t2_batch_display_parts(res: object, run_idx: int):
    """組出與單次 pd.DataFrame(rr) 同形之表 10-2。"""
    if not isinstance(res, dict):
        return []
    rows = res.get("table_10_2_ladder_rows") or []
    df = pd.DataFrame(rows)
    nc = res.get("n_candidates")
    title = f"批次執行表 10-2（候選母集 |𝒞_cand|={nc}）"
    return [(title, df, f"table_10_2_batch_run_{run_idx}.csv", True)]


if "t2_batch_runs" in st.session_state:
    render_batch_per_run_tables(
        st.session_state["t2_batch_runs"],
        _t2_batch_display_parts,
        key_prefix="t2_batch",
        column_name_map_for_df=lambda d: build_ch10_column_name_map(d.columns),
    )
