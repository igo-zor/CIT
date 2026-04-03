"""
§10.6 局部解析族與拼合障礙（循環 parity 視窗；對齊 §10.6.3–10.6.6）。
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

from hypergraph_experiment.ch10_symbol_glossary import (
    L10_ETA_CTX,
    L10_M_CTX,
    L10_MODE_CTX,
    L10_N_CTX,
    L10_N_SEARCH_106,
    L10_NODES_BIT,
    L10_SEED,
    L10_T_LOC,
    L10_W_CTX,
    render_section_glossary,
)
from hypergraph_experiment.experiments.exp_10_6_contexts import (
    run_canonical_demo_10_6,
    run_experiment_10_6,
    section_10_6_output_parameters_df,
)
from hypergraph_experiment.ch10_paper_presets import CH10_6_BASELINE
from hypergraph_experiment.streamlit_common import (
    batch_cell,
    build_ch10_column_name_map,
    flatten_result_row,
    render_batch_per_run_tables,
    render_parameters_table,
    render_sidebar_warehouse,
    render_table_with_copy_csv,
    run_batch_per_run_rows,
)

st.set_page_config(page_title="§10.6 拼合障礙", layout="wide")
render_sidebar_warehouse()

st.title("§10.6　拼合障礙與不可全域延拓")
render_section_glossary(st, "10.6")
st.markdown(
    r"""
    以 **GF(2) XOR 約束** 操作化 §10.6.6 之「循環交疊視窗」：參數 **$M,w,\eta,n,N_{\mathrm{ctx}}$** 與論文 §10.6.4／建議表一致；
    **obstruction** 模式在併入方程後翻轉奇性以得到不可全域延拓樣本；**satisfiable** 於 $w=2$ 時閉合環上 XOR 使主線可滿足。
    輸出對齊 §10.6.5：$\chi_{\mathrm{glue}}$、$\rho_{\mathrm{glue}}$、$r_{\min}$（以「上下文」子族計）、$N_{\mathrm{val}}$、$\rho_{\mathrm{val}}$。
    """
)

if st.button("顯示論文最小主例（canonical）"):
    _demo66 = run_canonical_demo_10_6()
    st.json(_demo66)
    st.subheader("§10.6.5 輸出參數（論文主表，canonical）")
    st.caption("與《約束世界論 30》§10.6.5 欄位對齊之三列簡表。")
    render_table_with_copy_csv(
        section_10_6_output_parameters_df(_demo66),
        key_prefix="t6_canonical_main",
        csv_filename="table_10_6_canonical_section_165_main.csv",
        hide_index=True,
    )

with st.form("f66"):
    st.subheader("§10.6.3 固定／建議表（可覆寫）")
    st.caption(L10_N_SEARCH_106 + "；parity 並查集實際極快，此欄主要與論文表對齊。")
    x0, x1 = st.columns(2)
    with x0:
        sd = st.number_input(L10_SEED, 0, 2_000_000_000, int(CH10_6_BASELINE["seed"]))
        n_search = st.number_input(
            "N_search｜全域搜尋步數上限", 100, 200_000, int(CH10_6_BASELINE["n_search"])
        )
    with x1:
        t_loc = st.number_input(
            L10_T_LOC + "（2–4；現僅 parity，數值小幅影響位元摺疊）",
            2,
            4,
            int(CH10_6_BASELINE["T_loc"]),
        )

    st.divider()
    st.subheader("§10.6.4 變量與執行參數")
    a1, a2, a3 = st.columns(3)
    with a1:
        n_ctx = st.number_input(L10_N_CTX, 10, 5000, int(CH10_6_BASELINE["n_ctx"]))
        n_nodes = st.number_input(L10_NODES_BIT, 4, 24, int(CH10_6_BASELINE["n_nodes"]))
    with a2:
        m_ctx = st.number_input(L10_M_CTX + "（論文建議 4–8）", 2, 16, int(CH10_6_BASELINE["M"]))
        w_ctx = st.number_input(L10_W_CTX + "（2–3）", 2, 6, int(CH10_6_BASELINE["w_ctx"]))
    with a3:
        eta_ctx = st.number_input(L10_ETA_CTX + "（1–2 為主）", 1, 5, int(CH10_6_BASELINE["eta_ctx"]))
        mode = st.selectbox(L10_MODE_CTX, ["obstruction", "satisfiable"])
    sub = st.form_submit_button("執行批次")

if sub:
    if int(w_ctx) <= int(eta_ctx):
        st.error("請確保 w > η（§10.6.4）。")
    else:
        st.session_state["out66"] = run_experiment_10_6(
            n_ctx=int(n_ctx),
            n_nodes=int(n_nodes),
            mode=str(mode),
            seed=int(sd),
            M=int(m_ctx),
            w=int(w_ctx),
            eta=int(eta_ctx),
            t_loc=int(t_loc),
            n_search=int(n_search),
        )
        st.session_state["out66_params"] = {
            "n_ctx": int(n_ctx),
            "n_nodes": int(n_nodes),
            "M": int(m_ctx),
            "w_ctx": int(w_ctx),
            "eta_ctx": int(eta_ctx),
            "T_loc": int(t_loc),
            "n_search": int(n_search),
            "mode": str(mode),
            "seed": int(sd),
        }

out66 = st.session_state.get("out66")
if out66:
    st.subheader("§10.6.5 輸出參數（主表，可複製／下載 CSV）")
    st.caption(
        "與《約束世界論 30》§10.6.5（一）–（四）對齊；下方「扁平欄位」為程式匯出寬表（欄名帶前綴 e6s_）。"
    )
    if (out66.get("metrics") or {}).get("error"):
        st.error(str(out66["metrics"]["error"]))
    df_t6_main = section_10_6_output_parameters_df(out66)
    render_table_with_copy_csv(
        df_t6_main,
        key_prefix="t6_single_main",
        csv_filename="table_10_6_single_section_165_main.csv",
        hide_index=True,
    )
    st.subheader("單次執行指標（扁平欄位，可複製／下載 CSV）")
    df_t6_single = pd.DataFrame([flatten_result_row("e6s", out66)])
    render_table_with_copy_csv(
        df_t6_single,
        key_prefix="t6_single",
        csv_filename="table_10_6_single.csv",
        column_name_map=build_ch10_column_name_map(df_t6_single.columns),
    )
    render_parameters_table(
        st.session_state.get("out66_params"),
        key_prefix="t6_single",
        csv_filename="table_10_6_single_params.csv",
    )
    with st.expander("完整 JSON（除錯）"):
        st.json(out66)

st.subheader("參數掃描批次")
st.caption(
    "欄位：上下文樣本數、節點數 n、上下文數 M、視窗 w、交疊 η、局部型別 T_loc、"
    "全域搜尋上限、執行模式、偽隨機基底種子。"
)
if st.button("載入論文建議批次模板（§10.6）", key="load_t6_template"):
    st.session_state["t6_batch_template_df"] = pd.DataFrame(
        [
            {
                "上下文樣本數": int(CH10_6_BASELINE["n_ctx"]),
                "位元節點數": int(CH10_6_BASELINE["n_nodes"]),
                "上下文數 M": int(CH10_6_BASELINE["M"]),
                "視窗大小 w": int(CH10_6_BASELINE["w_ctx"]),
                "視窗交疊 η": int(CH10_6_BASELINE["eta_ctx"]),
                "局部型別數 T_loc": int(CH10_6_BASELINE["T_loc"]),
                "全域搜尋上限": int(CH10_6_BASELINE["n_search"]),
                "執行模式": "obstruction",
                "偽隨機基底種子": int(CH10_6_BASELINE["seed"]),
            },
            {
                "上下文樣本數": int(CH10_6_BASELINE["n_ctx"]),
                "位元節點數": int(CH10_6_BASELINE["n_nodes"]),
                "上下文數 M": int(CH10_6_BASELINE["M"]),
                "視窗大小 w": int(CH10_6_BASELINE["w_ctx"]),
                "視窗交疊 η": int(CH10_6_BASELINE["eta_ctx"]),
                "局部型別數 T_loc": int(CH10_6_BASELINE["T_loc"]),
                "全域搜尋上限": int(CH10_6_BASELINE["n_search"]),
                "執行模式": "satisfiable",
                "偽隨機基底種子": int(CH10_6_BASELINE["seed"]) + 1,
            },
        ]
    )
df = st.data_editor(
    st.session_state.get(
        "t6_batch_template_df",
        pd.DataFrame(
            [
                {
                    "上下文樣本數": int(CH10_6_BASELINE["n_ctx"]),
                    "位元節點數": int(CH10_6_BASELINE["n_nodes"]),
                    "上下文數 M": int(CH10_6_BASELINE["M"]),
                    "視窗大小 w": int(CH10_6_BASELINE["w_ctx"]),
                    "視窗交疊 η": int(CH10_6_BASELINE["eta_ctx"]),
                    "局部型別數 T_loc": int(CH10_6_BASELINE["T_loc"]),
                    "全域搜尋上限": int(CH10_6_BASELINE["n_search"]),
                    "執行模式": "obstruction",
                    "偽隨機基底種子": int(CH10_6_BASELINE["seed"]),
                }
            ]
        ),
    ),
    num_rows="dynamic",
    key="b66",
)


def _t6_batch_display_parts(res: object, run_idx: int):
    if not isinstance(res, dict):
        return []
    df_main = section_10_6_output_parameters_df(res)
    df = pd.DataFrame([flatten_result_row("e6s", res)])
    return [
        (
            "§10.6.5 輸出參數（主表）",
            df_main,
            f"table_10_6_batch_run_{run_idx}_section_165_main.csv",
            True,
        ),
        (
            "單次執行指標（與本頁單次區塊同形）",
            df,
            f"table_10_6_batch_run_{run_idx}.csv",
        ),
    ]


if st.button("批次 §10.6"):

    def _r(r: pd.Series, prog):
        return run_experiment_10_6(
            n_ctx=int(batch_cell(r, "上下文樣本數", "n_ctx")),
            n_nodes=int(batch_cell(r, "位元節點數", "n_nodes")),
            M=int(batch_cell(r, "上下文數 M", "M", CH10_6_BASELINE["M"])),
            w=int(batch_cell(r, "視窗大小 w", "w_ctx", CH10_6_BASELINE["w_ctx"])),
            eta=int(batch_cell(r, "視窗交疊 η", "eta_ctx", CH10_6_BASELINE["eta_ctx"])),
            t_loc=int(batch_cell(r, "局部型別數 T_loc", "T_loc", CH10_6_BASELINE["T_loc"])),
            n_search=int(batch_cell(r, "全域搜尋上限", "n_search", CH10_6_BASELINE["n_search"])),
            mode=str(batch_cell(r, "執行模式", "mode")),
            seed=int(batch_cell(r, "偽隨機基底種子", "seed")),
            progress=prog,
        )

    st.session_state["t6_batch_runs"] = run_batch_per_run_rows(
        df, _r, stop_on_error=False, use_progress=True
    )

if "t6_batch_runs" in st.session_state:
    render_batch_per_run_tables(
        st.session_state["t6_batch_runs"],
        _t6_batch_display_parts,
        key_prefix="t6_batch",
        column_name_map_for_df=lambda d: build_ch10_column_name_map(d.columns),
    )
