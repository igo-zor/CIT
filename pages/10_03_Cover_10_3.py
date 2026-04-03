"""
§10.3 相容覆蓋、穩定化與解析分割（靜態模式 + 表 10-3）。
"""

from __future__ import annotations

import statistics
from typing import Any, Dict, List

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from hypergraph_experiment.ch10_symbol_glossary import (
    L10_CFG_ADMISSIBLE,
    L10_CONNECTED,
    L10_D_MAX,
    L10_DELTA,
    L10_FORBID_TRI,
    L10_K_MAX,
    L10_M_MAX,
    L10_N,
    L10_N_CFG,
    L10_N_REP,
    L10_OVERLAP,
    L10_SAMPLE_LIMIT,
    L10_SEED,
    L10_S_LAMBDA,
    render_section_glossary,
)
from hypergraph_experiment.core import (
    SIGNATURES,
    analyze_static,
    run_full_experiment,
    sample_candidates_and_filter,
    subsample_obs_configs,
)
from hypergraph_experiment.ch10_paper_presets import CH10_3_BASELINE
from hypergraph_experiment.streamlit_common import (
    batch_cell,
    build_ch10_column_name_map,
    render_batch_per_run_tables,
    render_parameters_table,
    render_sidebar_warehouse,
    render_table_with_copy_csv,
    run_batch_per_run_rows,
    streamlit_progress_callback,
    zip_download_hypergraph_run,
)

st.set_page_config(page_title="§10.3 解析覆蓋", layout="wide")
render_sidebar_warehouse()

st.title("§10.3　相容覆蓋與解析分割")
render_section_glossary(st, "10.3")
st.markdown(
    r"""
    在 **$\mathrm{Cfg}_\Lambda$** 上建 **$\approx_\delta$** 鄰域，穩定化得 **$\sim_\delta$** 與 **$S_{\Lambda,\delta}$**；
    輸出重疊率、傳遞違反率、熵等。
    符號對照：$\approx_\delta$=簽名距離不超過 $\delta$ 的相容關係、$\sim_\delta$=穩定化等價、$S_{\Lambda,\delta}$=解析單元空間。
    **建議**：$N_{\mathrm{cfg}}$ 固定先比較 weak/medium/strong 與 $\delta$ 掃描；$\delta\in\{0,1,2,\ldots\}$（整數簽名距離）。
    本實作之 **δ** 與 ``signature_distance`` 一致，為**整數**離散尺度；若與正文以小數標度敘述並觀，僅為表述尺度差異，操作上以本頁數值為準。
    """
)

with st.form("f10_3"):
    st.subheader("建議固定參數（可覆寫）")
    a1, a2 = st.columns(2)
    with a1:
        n = st.number_input(L10_N, 2, 12, 8)
        k_max = st.number_input(L10_K_MAX, 2, 6, 3)
        m_max = st.number_input(L10_M_MAX, 1, 20, 10)
        d_max = st.number_input(L10_D_MAX, 1, 20, 4)
    with a2:
        sl = st.number_input(L10_SAMPLE_LIMIT, 0, 50_000, 5000)
        n_cfg_in = st.number_input(
            L10_N_CFG,
            1,
            50_000,
            int(CH10_3_BASELINE["n_cfg"]),
            help="與候選上限分工：先由 N_cand 生成可採用域，再抽 N_cfg 筆作解析觀測集。",
        )
        n_rep_in = st.number_input(
            L10_N_REP,
            1,
            200,
            20,
            help="同一組參數下重抽觀測集次數；1 表示不重抽。",
        )
        sd = st.number_input(L10_SEED, 0, 2_000_000_000, 20)
        conn = st.checkbox("2-section 連通｜域型條件", True, help=L10_CONNECTED)
        ft = st.checkbox("禁二元△｜forbidden motif", False, help=L10_FORBID_TRI)
    st.divider()
    st.subheader("建議變量掃描參數（可覆寫）")
    st.caption("本頁固定同時掃描 weak / medium / strong 三種解析簽名版本（對應論文 Sig_Λ）。")
    b1, b2 = st.columns(2)
    with b1:
        delta_min = st.number_input("解析閾值整數最小值 δ_min", 0, 50, 0)
        delta_max = st.number_input("解析閾值整數最大值 δ_max（含）", 0, 50, 2)
    with b2:
        s_min = st.number_input(
            "重疊率鄰域最小支持 s_min（僅 |T(c)|≥s_min 參與配對；0 表示不過濾）",
            0,
            50,
            int(CH10_3_BASELINE["s_min"]),
        )
    sub = st.form_submit_button("執行靜態實驗")

if sub:
    bar = st.progress(0.0)
    txt = st.empty()

    signatures: List[str] = sorted(SIGNATURES.keys())
    rows: List[Dict[str, Any]] = []

    if int(delta_max) < int(delta_min):
        st.error("請確保 δ_max ≥ δ_min。")
    else:
        try:
            deltas: List[int] = list(range(int(delta_min), int(delta_max) + 1))
            n_rep_eff = int(max(1, n_rep_in))
            total = len(signatures) * len(deltas) * n_rep_eff
            step = 0

            txt.write("建立候選與可採用域（僅一次）...")
            _, admissible_configs = sample_candidates_and_filter(
                n=int(n),
                max_edge_size=int(k_max),
                max_edges=int(m_max),
                sample_limit=int(sl),
                seed=int(sd),
                connected=conn,
                max_degree=int(d_max),
                forbid_pair_triangles=ft,
                progress=None,
            )
            if not admissible_configs:
                st.error("可採用配置集合為空，請放寬域型條件或加大候選採樣上限。")
                rows = []
                raise ValueError("empty admissible set")

            for sig_name in signatures:
                for d_val in deltas:
                    overlap_vals: List[float] = []
                    trans_vals: List[float] = []
                    u_vals: List[float] = []
                    iso_vals: List[float] = []
                    entropy_vals: List[float] = []
                    class_vals: List[float] = []
                    obs_vals: List[int] = []
                    for rep_idx in range(n_rep_eff):
                        step += 1
                        rep_seed = int(sd) + rep_idx * 1009
                        txt.write(
                            f"簽名 {sig_name}，δ={d_val}，重抽 {rep_idx + 1}/{n_rep_eff}（{step}/{total}）"
                        )
                        obs, _n_req, n_obs_actual, _notice = subsample_obs_configs(
                            admissible_configs, int(n_cfg_in), seed=rep_seed
                        )
                        a = analyze_static(obs, sig_name, int(d_val), s_min=int(s_min))
                        overlap_vals.append(float(a.get("overlap_rate", 0.0)))
                        trans_vals.append(float(a.get("transitivity_violation_rate", 0.0)))
                        u_vals.append(float(a.get("compression_ratio_U", 0.0)))
                        iso_vals.append(float(a.get("isol_rate_compat_graph", 0.0)))
                        entropy_vals.append(float(a.get("entropy_bits", 0.0)))
                        class_vals.append(float(a.get("num_equivalence_classes", 0.0)))
                        obs_vals.append(int(n_obs_actual))
                        bar.progress(step / float(total))
                    rows.append(
                        {
                            "解析簽名": sig_name,
                            "解析閾值整數": int(d_val),
                            "觀測配置數": int(round(statistics.mean(obs_vals))),
                            "解析單元數": float(statistics.mean(class_vals)),
                            "平均單元大小": (
                                (float(statistics.mean(obs_vals)) / float(statistics.mean(class_vals)))
                                if statistics.mean(class_vals) > 0
                                else 0.0
                            ),
                            "解析壓縮比 U_Λ": float(statistics.mean(u_vals)),
                            "重疊率 R_overlap": float(statistics.mean(overlap_vals)),
                            "相容孤立率 R_iso": float(statistics.mean(iso_vals)),
                            "傳遞違反率 R_trans_viol": float(statistics.mean(trans_vals)),
                            "解析熵位元": float(statistics.mean(entropy_vals)),
                        }
                    )
            st.success("§10.3 掃描完成。")
        except Exception as e:
            st.error(f"執行 §10.3 掃描時發生錯誤：{e}")
            rows = []

    if rows:
        st.subheader("§10.3 主表：解析簽名 × 解析閾值掃描")
        df_main = pd.DataFrame(rows).sort_values(["解析簽名", "解析閾值整數"])
        render_table_with_copy_csv(
            df_main,
            key_prefix="t3_main",
            csv_filename="table_10_3_main.csv",
            column_name_map=build_ch10_column_name_map(df_main.columns),
            hide_index=True,
        )
        render_parameters_table(
            {
                "n": int(n),
                "k_max": int(k_max),
                "m_max": int(m_max),
                "d_max": int(d_max),
                "connected": bool(conn),
                "forbid_pair_triangles": bool(ft),
                "sample_limit": int(sl),
                "seed": int(sd),
                "n_cfg": int(n_cfg),
                "n_rep": int(n_rep),
                "s_min": int(s_min),
                "signatures": list(signatures),
                "delta_values": list(deltas),
            },
            key_prefix="t3_main",
            csv_filename="table_10_3_main_params.csv",
        )
        st.caption(
            f"指標對照：候選勢 **|𝒞_cand|**、{L10_CFG_ADMISSIBLE}、{L10_OVERLAP}、{L10_S_LAMBDA}，"
            "此處主表僅顯示平均指標。"
        )

        def _plot_metric(metric_col: str, title: str, ylabel: str) -> None:
            fig, ax = plt.subplots(figsize=(6, 3))
            for sig_name in signatures:
                sub = df_main[df_main["解析簽名"] == sig_name].sort_values("解析閾值整數")
                ax.plot(
                    sub["解析閾值整數"],
                    sub[metric_col],
                    marker="o",
                    label=sig_name,
                )
            ax.set_title(title)
            ax.set_xlabel("解析閾值 δ")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.3)
            ax.legend()
            st.pyplot(fig)
            plt.close(fig)

        st.subheader("輔圖：指標隨解析閾值 δ 之變化")
        _plot_metric("重疊率 R_overlap", "R_overlap vs δ", "R_overlap")
        _plot_metric("傳遞違反率 R_trans_viol", "R_trans-viol vs δ", "R_trans-viol")
        _plot_metric("解析壓縮比 U_Λ", "U_Λ vs δ", "U_Λ = |S|/N_cfg")
        _plot_metric("相容孤立率 R_iso", "R_iso vs δ", "R_iso")

        st.subheader("簽名比較表：固定 δ 下之 weak / medium / strong")
        delta_choices = sorted(df_main["解析閾值整數"].unique())
        delta_choice = st.selectbox(
            "選擇展示用解析閾值 δ（對應論文簽名比較表）",
            delta_choices,
            index=delta_choices.index(0) if 0 in delta_choices else 0,
        )
        df_delta = df_main[df_main["解析閾值整數"] == int(delta_choice)].copy()
        if not df_delta.empty:
            df_sig = []
            for sig_name in signatures:
                row_sig = df_delta[df_delta["解析簽名"] == sig_name]
                if row_sig.empty:
                    continue
                r = row_sig.iloc[0]
                df_sig.append(
                    {
                        "解析簽名": sig_name,
                        "解析單元數": float(r["解析單元數"]),
                        "平均單元大小": float(r["平均單元大小"]),
                        "解析壓縮比 U_Λ": float(r["解析壓縮比 U_Λ"]),
                        "重疊率 R_overlap": float(r["重疊率 R_overlap"]),
                        "相容孤立率 R_iso": float(r["相容孤立率 R_iso"]),
                        "傳遞違反率 R_trans-viol": float(r["傳遞違反率 R_trans_viol"]),
                        "解析熵位元": float(r["解析熵位元"]),
                    }
                )
            if df_sig:
                df_sig_table = pd.DataFrame(df_sig)
                render_table_with_copy_csv(
                    df_sig_table,
                    key_prefix="t3_sig_compare",
                    csv_filename=f"table_10_3_signature_compare_delta_{delta_choice}.csv",
                    hide_index=True,
                    column_name_map=build_ch10_column_name_map(df_sig_table.columns),
                )

st.subheader("批次（多組 n、候選上限、δ、偽隨機基底種子）")
st.caption(
    "欄位包含節點數、超邊大小上限、邊數上限、度數上限、候選採樣上限、解析閾值、偽隨機基底種子與解析簽名層級；"
    "本批次列預設不啟用連通與禁三角域型，僅供快速掃描。"
)
if st.button("載入論文建議批次模板（§10.3）", key="load_t3_template"):
    st.session_state["t3_batch_template_df"] = pd.DataFrame(
        [
            {
                "節點數": 8,
                "最大超邊階數": 3,
                "最大超邊數": 10,
                "頂點度上限": 4,
                "二部圖連通": True,
                "禁止二元三角": False,
                "候選採樣上限": 5000,
                "輸入配置數": 300,
                "重複次數": 20,
                "解析閾值整數": 0,
                "偽隨機基底種子": 20,
                "解析簽名": "weak",
                "重疊率鄰域最小支持": 2,
            },
            {
                "節點數": 8,
                "最大超邊階數": 3,
                "最大超邊數": 10,
                "頂點度上限": 4,
                "二部圖連通": True,
                "禁止二元三角": False,
                "候選採樣上限": 5000,
                "輸入配置數": 300,
                "重複次數": 20,
                "解析閾值整數": 0,
                "偽隨機基底種子": 20,
                "解析簽名": "medium",
                "重疊率鄰域最小支持": 2,
            },
            {
                "節點數": 8,
                "最大超邊階數": 3,
                "最大超邊數": 10,
                "頂點度上限": 4,
                "二部圖連通": True,
                "禁止二元三角": False,
                "候選採樣上限": 5000,
                "輸入配置數": 300,
                "重複次數": 20,
                "解析閾值整數": 0,
                "偽隨機基底種子": 20,
                "解析簽名": "strong",
                "重疊率鄰域最小支持": 2,
            },
        ]
    )
df_b = st.data_editor(
    st.session_state.get(
        "t3_batch_template_df",
        pd.DataFrame(
            [
                {
                    "節點數": 8,
                    "最大超邊階數": 3,
                    "最大超邊數": 10,
                    "頂點度上限": 4,
                    "二部圖連通": True,
                    "禁止二元三角": False,
                    "候選採樣上限": 5000,
                    "輸入配置數": 300,
                    "重複次數": 20,
                    "解析閾值整數": 0,
                    "偽隨機基底種子": 20,
                    "解析簽名": "medium",
                    "重疊率鄰域最小支持": 2,
                }
            ]
        ),
    ),
    num_rows="dynamic",
    key="b10_3_ed",
)


def _row_static(r: pd.Series, prog):
    return run_full_experiment(
        mode="static",
        n=int(batch_cell(r, "節點數", "n")),
        max_edge_size=int(batch_cell(r, "最大超邊階數", "k_max")),
        max_edges=int(batch_cell(r, "最大超邊數", "m_max")),
        max_degree=int(batch_cell(r, "頂點度上限", "d_max")),
        connected=bool(batch_cell(r, "二部圖連通", "connected", True)),
        forbid_pair_triangles=bool(batch_cell(r, "禁止二元三角", "forbid_pair_triangles", False)),
        sample_limit=int(batch_cell(r, "候選採樣上限", "sample_limit")),
        n_cfg=int(batch_cell(r, "輸入配置數", "n_cfg", 300)),
        n_rep=int(batch_cell(r, "重複次數", "n_rep", 20)),
        signature=str(batch_cell(r, "解析簽名", "signature", "medium")),
        delta=int(batch_cell(r, "解析閾值整數", "delta")),
        s_min=int(batch_cell(r, "重疊率鄰域最小支持", "s_min", 2)),
        seed=int(batch_cell(r, "偽隨機基底種子", "seed")),
        show_sample_configs=0,
        progress=prog,
    )


def _t3_batch_display_parts(res: object, run_idx: int):
    if not isinstance(res, dict):
        return []
    df = pd.DataFrame([flatten_result_row("r3s", res)])
    return [
        (
            "單次執行扁平指標表（與本頁單次區塊同形）",
            df,
            f"table_10_3_batch_run_{run_idx}.csv",
        )
    ]


if st.button("批次靜態實驗"):
    st.session_state["t3_batch_runs"] = run_batch_per_run_rows(
        df_b, _row_static, stop_on_error=False, use_progress=True
    )

if "t3_batch_runs" in st.session_state:
    render_batch_per_run_tables(
        st.session_state["t3_batch_runs"],
        _t3_batch_display_parts,
        key_prefix="t3_batch",
        column_name_map_for_df=lambda d: build_ch10_column_name_map(d.columns),
    )
