"""
§10.4 解析細化、纖維與推前一致性。
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

from hypergraph_experiment.ch10_symbol_glossary import (
    L10_DELTA_COARSE,
    L10_DELTA_FINE,
    L10_D_MAX,
    L10_H_COARSE,
    L10_H_FINE,
    L10_JS_TERM,
    L10_K_MAX,
    L10_KERNEL,
    L10_M_MAX,
    L10_N,
    L10_PUSH_ERR,
    L10_SAMPLE_LIMIT,
    L10_SEED,
    L10_SIG_COARSE,
    L10_SIG_FINE,
    render_section_glossary,
)
from hypergraph_experiment.core import run_full_experiment, sample_candidates_and_filter
from hypergraph_experiment.refinement import compare_ordered_refinement_paths, preset_layer_to_substep_spec
from hypergraph_experiment.ch10_paper_presets import CH10_4_BASELINE
from hypergraph_experiment.streamlit_common import (
    batch_cell,
    build_ch10_column_name_map,
    flatten_result_row,
    render_batch_per_run_tables,
    render_hypergraph_preview,
    render_parameters_table,
    render_sidebar_warehouse,
    render_table_with_copy_csv,
    run_batch_per_run_rows,
    streamlit_progress_callback,
    zip_download_hypergraph_run,
)

# 實驗 B 子步驟種類（對應 refinement.RefinementSubstepSpec.kind，勿與主表 weak/medium/strong 標籤混淆）
SUBSTEP_KIND_KEYS = ("edge_scale_motif", "degree_split", "adjacency_motif_fine")
SUBSTEP_KIND_LABELS_ZH = {
    "edge_scale_motif": "超邊尺度（邊階多重集）— 動機粗分",
    "degree_split": "度序列細分（邊階+度數）",
    "adjacency_motif_fine": "鄰接與動機計數（精細）",
}


REFINEMENT_PRESETS = {
    "weak→medium（主表組合一）": {
        "coarse_sig": "weak",
        "coarse_delta": 3,
        "fine_sig": "medium",
        "fine_delta": 1,
    },
    "medium→strong（主表組合二）": {
        "coarse_sig": "medium",
        "coarse_delta": 2,
        "fine_sig": "strong",
        "fine_delta": 0,
    },
}
st.set_page_config(page_title="§10.4 細化", layout="wide")
render_sidebar_warehouse()

st.title("§10.4　解析細化與可交換性觀察（$\\Lambda\\to\\Lambda'$）")
render_section_glossary(st, "10.4")
st.markdown(
    r"""
    檢查 **$\pi_{\Lambda'\to\Lambda}$**（粗細投影）、纖維條件核、推前誤差與雙路徑終端差異（**JS_term**）。
    本頁將主表實驗（10.4.5/10.4.6）與 A/B 實驗（10.4.7）拆為兩組獨立輸入與執行。
    """
)

st.subheader("共用基礎參數")
b1, b2 = st.columns(2)
with b1:
    n = st.number_input(L10_N, 2, 12, 8)
    k_max = st.number_input(L10_K_MAX, 2, 6, 3)
    m_max = st.number_input(L10_M_MAX, 1, 20, 10)
    d_max = st.number_input(L10_D_MAX, 1, 20, 4)
with b2:
    sl = st.number_input(L10_SAMPLE_LIMIT, 100, 20_000, 5000)
    sd = st.number_input(L10_SEED, 0, 2_000_000_000, 20)
    conn = st.checkbox("2-section 連通｜域型條件", True)
    ft = st.checkbox("禁二元△｜forbidden motif", False)
st.caption("下列兩組實驗共用上述候選生成與域型條件。")

st.divider()
st.subheader("實驗 A：主表（10.4.5 / 10.4.6）")
with st.form("f10_4_main"):
    preset_name = st.selectbox("細化鏈 preset（粗/細層與閾值綁定）", list(REFINEMENT_PRESETS.keys()))
    preset = REFINEMENT_PRESETS[preset_name]
    r_cs = str(preset["coarse_sig"])
    r_cd = int(preset["coarse_delta"])
    r_fs = str(preset["fine_sig"])
    r_fd = int(preset["fine_delta"])
    st.caption(
        f"{L10_SIG_COARSE}={r_cs}，{L10_DELTA_COARSE}={r_cd}；"
        f"{L10_SIG_FINE}={r_fs}，{L10_DELTA_FINE}={r_fd}"
    )
    c1, c2, c3 = st.columns(3)
    with c1:
        km_main = st.selectbox(L10_KERNEL, ["uniform", "proportional"], key="km_main")
        n_cfg_coarse = st.number_input(
            "粗層樣本數 N_cfg_coarse", 1, 50_000, int(CH10_4_BASELINE["coarse_sample_size"])
        )
    with c2:
        n_cfg_fine = st.number_input(
            "細層樣本數 N_cfg_fine", 1, 50_000, int(CH10_4_BASELINE["fine_sample_size"])
        )
        fiber_sample_size = st.number_input("纖維樣本數上限（每粗單元）", 1, 5000, 200)
    with c3:
        eps_push_threshold = st.number_input(
            "推前誤差閾值 ε_push", 0.0, 1.0, float(CH10_4_BASELINE["eps_push_threshold"]), step=0.01
        )
        js_threshold = st.number_input(
            "JS 差異閾值", 0.0, 1.0, float(CH10_4_BASELINE["js_threshold"]), step=0.01
        )
    run_main = st.form_submit_button("執行主表實驗（10.4.5 / 10.4.6）")

if run_main:
    bar = st.progress(0)
    pr = streamlit_progress_callback(bar, st.empty())
    try:
        res_main = run_full_experiment(
            mode="static",
            n=int(n),
            max_edge_size=int(k_max),
            max_edges=int(m_max),
            max_degree=int(d_max),
            connected=bool(conn),
            forbid_pair_triangles=bool(ft),
            sample_limit=int(sl),
            signature=str(r_fs),
            delta=int(r_fd),
            seed=int(sd),
            refinement_enabled=True,
            refine_coarse_signature=str(r_cs),
            refine_coarse_delta=int(r_cd),
            refine_fine_signature=str(r_fs),
            refine_fine_delta=int(r_fd),
            refine_kernel=str(km_main),
            refine_compare_chains=False,
            refine_coarse_sample_size=int(n_cfg_coarse),
            refine_fine_sample_size=int(min(n_cfg_fine, n_cfg_coarse)),
            refine_fiber_sample_size=int(fiber_sample_size),
            progress=pr,
        )
        res_main["thresholds_10_4"] = {
            "epsilon_push_threshold": float(eps_push_threshold),
            "js_threshold": float(js_threshold),
        }
        st.session_state["res_10_4_main"] = res_main
        st.session_state["res_10_4_main_params"] = {
            "n": int(n),
            "max_edge_size": int(k_max),
            "max_edges": int(m_max),
            "max_degree": int(d_max),
            "connected": bool(conn),
            "forbid_pair_triangles": bool(ft),
            "sample_limit": int(sl),
            "seed": int(sd),
            "refine_coarse_signature": str(r_cs),
            "refine_coarse_delta": int(r_cd),
            "refine_fine_signature": str(r_fs),
            "refine_fine_delta": int(r_fd),
            "refine_kernel": str(km_main),
            "refine_coarse_sample_size": int(n_cfg_coarse),
            "refine_fine_sample_size": int(min(n_cfg_fine, n_cfg_coarse)),
            "refine_fiber_sample_size": int(fiber_sample_size),
            "epsilon_push_threshold": float(eps_push_threshold),
            "js_threshold": float(js_threshold),
        }
        st.success("主表實驗完成")
    except Exception as e:
        st.error(str(e))

st.divider()
st.subheader("實驗 B：A/B 順序可交換性（10.4.7）")
with st.form("f10_4_ab"):
    preset_name_ab = st.selectbox(
        "A/B 預設來源 preset（僅填入預設子步驟種類與 δ）", list(REFINEMENT_PRESETS.keys())
    )
    preset_ab = REFINEMENT_PRESETS[preset_name_ab]
    step_a_def = preset_layer_to_substep_spec(str(preset_ab["coarse_sig"]), int(preset_ab["coarse_delta"]))
    step_b_def = preset_layer_to_substep_spec(str(preset_ab["fine_sig"]), int(preset_ab["fine_delta"]))
    st.caption("A 與 B 參數固定拆分為兩組獨立輸入；預設值由 preset 帶入，可直接調整。")
    a1, a2 = st.columns(2)
    with a1:
        st.markdown("**細化 A 參數（step_a）**")
        kind_a = st.selectbox(
            "A 子步驟映射 R_A",
            list(SUBSTEP_KIND_KEYS),
            index=list(SUBSTEP_KIND_KEYS).index(step_a_def["kind"]),
            format_func=lambda k: SUBSTEP_KIND_LABELS_ZH[str(k)],
        )
        delta_a = st.number_input("A 相容閾值 δ（整數）", 0, 50, int(step_a_def["delta"]))
    with a2:
        st.markdown("**細化 B 參數（step_b）**")
        kind_b = st.selectbox(
            "B 子步驟映射 R_B",
            list(SUBSTEP_KIND_KEYS),
            index=list(SUBSTEP_KIND_KEYS).index(step_b_def["kind"]),
            format_func=lambda k: SUBSTEP_KIND_LABELS_ZH[str(k)],
        )
        delta_b = st.number_input("B 相容閾值 δ（整數）", 0, 50, int(step_b_def["delta"]))
    step_a = {"kind": str(kind_a), "delta": int(delta_a)}
    step_b = {"kind": str(kind_b), "delta": int(delta_b)}
    k_ab1, k_ab2 = st.columns(2)
    with k_ab1:
        km_ab = st.selectbox("A/B 細化核模式", ["uniform", "proportional"], key="km_ab")
        n_cfg_coarse_ab = st.number_input(
            "A/B 粗層樣本數 N_cfg_coarse",
            1,
            50_000,
            int(CH10_4_BASELINE["coarse_sample_size"]),
            key="n_cfg_coarse_ab",
        )
    with k_ab2:
        n_cfg_fine_ab = st.number_input(
            "A/B 細層樣本數 N_cfg_fine",
            1,
            50_000,
            int(CH10_4_BASELINE["fine_sample_size"]),
            key="n_cfg_fine_ab",
        )
        fiber_sample_size_ab = st.number_input("A/B 纖維樣本數上限", 1, 5000, 200, key="fiber_sample_size_ab")
    run_ab = st.form_submit_button("執行 A/B 實驗（10.4.7）")

if run_ab:
    try:
        _, cfgs = sample_candidates_and_filter(
            n=int(n),
            max_edge_size=int(k_max),
            max_edges=int(m_max),
            sample_limit=int(sl),
            seed=int(sd),
            connected=bool(conn),
            max_degree=int(d_max),
            forbid_pair_triangles=bool(ft),
        )
        cmp_out = compare_ordered_refinement_paths(
            cfgs,
            step_a=step_a,  # type: ignore[arg-type]
            step_b=step_b,  # type: ignore[arg-type]
            kernel_mode=str(km_ab),  # type: ignore[arg-type]
            coarse_sample_size=int(n_cfg_coarse_ab),
            fine_sample_size=int(min(n_cfg_fine_ab, n_cfg_coarse_ab)),
            sample_seed=int(sd),
            max_fiber_size=int(fiber_sample_size_ab),
        )
        st.session_state["res_10_4_ab"] = cmp_out
        st.session_state["res_10_4_ab_params"] = {
            "n": int(n),
            "max_edge_size": int(k_max),
            "max_edges": int(m_max),
            "max_degree": int(d_max),
            "connected": bool(conn),
            "forbid_pair_triangles": bool(ft),
            "sample_limit": int(sl),
            "seed": int(sd),
            "step_a": dict(step_a),
            "step_b": dict(step_b),
            "kernel_mode": str(km_ab),
            "coarse_sample_size": int(n_cfg_coarse_ab),
            "fine_sample_size": int(min(n_cfg_fine_ab, n_cfg_coarse_ab)),
            "max_fiber_size": int(fiber_sample_size_ab),
        }
        st.success("A/B 實驗完成")
    except Exception as e:
        st.error(str(e))

res_main = st.session_state.get("res_10_4_main")
res_ab = st.session_state.get("res_10_4_ab")
thr_global = (res_main or {}).get("thresholds_10_4", {}) if isinstance(res_main, dict) else {}
ep_th_global = float(thr_global.get("epsilon_push_threshold", float(CH10_4_BASELINE["eps_push_threshold"])))
js_th_global = float(thr_global.get("js_threshold", float(CH10_4_BASELINE["js_threshold"])))
if res_main:
    st.subheader("單次執行扁平指標表（研究記錄）")
    df_t4_single = pd.DataFrame([flatten_result_row("r4s", res_main)])
    render_table_with_copy_csv(
        df_t4_single,
        key_prefix="t4_single",
        csv_filename="table_10_4_single.csv",
        column_name_map=build_ch10_column_name_map(df_t4_single.columns),
    )
    render_parameters_table(
        st.session_state.get("res_10_4_main_params"),
        key_prefix="t4_single",
        csv_filename="table_10_4_single_params.csv",
    )
    ref = res_main.get("refinement_10_4") or {}
    single = ref.get("single_step_Lambda_to_Lambda_prime") or {}
    tw_main = res_ab if isinstance(res_ab, dict) else {}
    thr = res_main.get("thresholds_10_4") or {}
    ep_th = float(thr.get("epsilon_push_threshold", float(CH10_4_BASELINE["eps_push_threshold"])))
    js_th = float(thr.get("js_threshold", float(CH10_4_BASELINE["js_threshold"])))

    # 論文 10.4.5 對齊主表（H 粗/細、ε_push、JS_term、ΔH_term、可交換性判定）
    main_row = {
        "粗層熵 H(p_Λ)": single.get("H_p_Lambda_bits"),
        "細層熵 H(p_Λ')": single.get("H_p_Lambda_prime_bits"),
        "推前誤差 ε_push": single.get("pushforward_max_error"),
        "終端 JS 差異 JS_term": tw_main.get("js_divergence_bits_terminal_ab_ba"),
        "終端熵差 ΔH_term": tw_main.get("entropy_abs_diff_terminal_ab_ba"),
    }
    ep_val = float(main_row["推前誤差 ε_push"] or 0.0)
    js_val_any = main_row["終端 JS 差異 JS_term"]
    js_val = float(js_val_any) if js_val_any is not None else None
    if js_val is None:
        ex_label = "無法判定（終端 JS 維度不一致或未計算）"
    else:
        ex_label = "近似可交換" if (ep_val <= ep_th and js_val <= js_th) else "具順序依賴"
    main_row["可交換性判定"] = ex_label
    st.subheader("10.4.5 對齊主表（輸出參數）")
    df_t45 = pd.DataFrame([main_row]).round(3)
    render_table_with_copy_csv(
        df_t45,
        key_prefix="t4_out_main",
        csv_filename="table_10_4_outputs_main.csv",
        hide_index=True,
        column_name_map=build_ch10_column_name_map(df_t45.columns),
    )
    st.caption(
        f"10.4.5 判準門檻：ε_push*={ep_th:.3f}，ε_JS*={js_th:.3f}；"
        "可交換性判定以兩者同時達標為準。"
    )

    st.caption(f"10.4.7 觀察一/二：單步細化 {L10_PUSH_ERR}；熵指標 {L10_H_COARSE}、{L10_H_FINE}。")
    if "error" not in single:
        ep = float(single.get("pushforward_max_error", 0.0) or 0.0)
        st.metric("ε_push 推前最大誤差", f"{ep:.3f}", help=f"判準閾值：{ep_th:.3f}")
        st.caption("判準：" + ("推前一致（within ε）" if ep <= ep_th else "推前偏差較大（超過 ε）"))
        x1, x2, x3 = st.columns(3)
        x1.metric("H(p_Λ) 粗層熵 (bit)", f'{single.get("H_p_Lambda_bits", 0):.4f}')
        x2.metric("H(p_Λ′) 細層熵 (bit)", f'{single.get("H_p_Lambda_prime_bits", 0):.4f}')
        x3.metric("細化成立（推前一致）", "是" if single.get("refinement_valid") else "否")
        det = single.get("detail") or {}
        ks = det.get("kernel_stability") or {}
        if ks and int(ks.get("kernel_n_fibers", 0) or 0) > 0:
            y1, y2, y3 = st.columns(3)
            y1.metric("纖維核穩定度：平均非零比例", f'{float(ks.get("kernel_mean_nnz_ratio", 0)):.4f}')
            y2.metric("纖維核穩定度：平均熵 bit", f'{float(ks.get("kernel_mean_entropy_bits", 0)):.4f}')
            y3.metric("纖維樣本上限", det.get("fiber_sample_size"))
        st.caption(
            f"樣本數：粗層={det.get('coarse_sample_size_used')}，細層={det.get('fine_sample_size_used')}"
        )
    if tw_main and "error" not in tw_main:
        js_v = tw_main.get("js_divergence_bits_terminal_ab_ba")
        st.write(
            f"JS_term={tw_main.get('js_divergence_bits_terminal_ab_ba')}；"
            f"ΔH_term={tw_main.get('entropy_abs_diff_terminal_ab_ba')}"
        )
        if js_v is not None:
            st.caption("10.4.7 判讀：" + ("雙路徑近似可交換（JS_term 低）" if float(js_v) <= js_th else "雙路徑差異可觀（JS_term 偏高）"))
if isinstance(res_ab, dict):
    st.subheader("10.4.7 A/B 進階比較表")
    cmp_out = res_ab
    rows_cmp = []
    for key in ("A_to_B", "B_to_A"):
        d = cmp_out.get(key) or {}
        ks = (d.get("kernel_stability") or {}) if isinstance(d, dict) else {}
        rows_cmp.append(
            {
                "鏈": d.get("path_key"),
                "推前誤差": d.get("pushforward_max_error"),
                "終端熵 H_term (bit)": d.get("entropy_fine_bits"),
                "纖維核穩定度：平均非零比例": ks.get("kernel_mean_nnz_ratio"),
                "纖維核穩定度：平均熵 bit": ks.get("kernel_mean_entropy_bits"),
                "終端維度": d.get("terminal_dim"),
            }
        )
    df_cmp = pd.DataFrame(rows_cmp).round(3)
    render_table_with_copy_csv(
        df_cmp,
        key_prefix="t4_ab_cmp",
        csv_filename="table_10_4_ab_compare.csv",
        hide_index=True,
        column_name_map=build_ch10_column_name_map(df_cmp.columns),
    )
    render_parameters_table(
        st.session_state.get("res_10_4_ab_params"),
        key_prefix="t4_ab_cmp",
        csv_filename="table_10_4_ab_params.csv",
    )
    st.caption(
        f"10.4.7 跨鏈摘要：JS_term={cmp_out.get('js_divergence_bits_terminal_ab_ba')}；"
        f"ΔH_term={cmp_out.get('entropy_abs_diff_terminal_ab_ba')}"
    )
    js_cross = cmp_out.get("js_divergence_bits_terminal_ab_ba")
    if js_cross is None:
        st.caption("跨鏈可交換性：無法判定（終端維度不一致）。")
    else:
        ep_ab = float((rows_cmp[0].get("推前誤差") or 0.0))
        ep_ba = float((rows_cmp[1].get("推前誤差") or 0.0))
        ok = ep_ab <= ep_th_global and ep_ba <= ep_th_global and float(js_cross) <= js_th_global
        st.caption("10.4.7 跨鏈可交換性判定：" + ("近似可交換" if ok else "具順序依賴"))
if res_main:
    with st.expander("refinement JSON（除錯）"):
        st.json(ref)
    render_hypergraph_preview(res_main.get("sample_configs") or [], key_prefix="p04")
    st.download_button(
        "下載 ZIP",
        zip_download_hypergraph_run(res_main, "x.zip"),
        "refine_run.zip",
        "application/zip",
    )

st.subheader("批次細化參數表")
st.caption(
    "欄位含固定與掃描參數：節點數、超邊與邊數上限、度數上限、連通與禁三角條件、候選採樣上限、偽隨機基底種子、粗細解析簽名與閾值、細化核。"
)
if st.button("載入論文建議批次模板（§10.4）", key="load_t4_template"):
    st.session_state["t4_batch_template_df"] = pd.DataFrame(
        [
            {
                "節點數": 8,
                "最大超邊階數": 3,
                "最大超邊數": 10,
                "頂點度上限": 4,
                "二部圖連通": True,
                "禁止二元三角": False,
                "候選採樣上限": 5000,
                "偽隨機基底種子": 20,
                "細化鏈 preset": "weak→medium（主表組合一）",
                "細化核模式": "uniform",
                "粗層樣本數": 2000,
                "細層樣本數": 2000,
                "纖維樣本數上限": 200,
                "推前誤差閾值": 0.01,
                "JS 差異閾值": 0.01,
            },
            {
                "節點數": 8,
                "最大超邊階數": 3,
                "最大超邊數": 10,
                "頂點度上限": 4,
                "二部圖連通": True,
                "禁止二元三角": False,
                "候選採樣上限": 5000,
                "偽隨機基底種子": 21,
                "細化鏈 preset": "medium→strong（主表組合二）",
                "細化核模式": "proportional",
                "粗層樣本數": 2000,
                "細層樣本數": 2000,
                "纖維樣本數上限": 200,
                "推前誤差閾值": 0.01,
                "JS 差異閾值": 0.01,
            },
        ]
    )
df_r4 = st.data_editor(
    st.session_state.get(
        "t4_batch_template_df",
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
                    "偽隨機基底種子": 20,
                    "細化鏈 preset": "weak→medium（主表組合一）",
                    "細化核模式": "uniform",
                    "粗層樣本數": 2000,
                    "細層樣本數": 2000,
                    "纖維樣本數上限": 200,
                    "推前誤差閾值": 0.01,
                    "JS 差異閾值": 0.01,
                }
            ]
        ),
    ),
    num_rows="dynamic",
    key="b10_4_ed",
)
if st.button("批次 §10.4", key="run_batch_t4"):

    def _row_refine(row: pd.Series, prog):
        return run_full_experiment(
            mode="static",
            n=int(batch_cell(row, "節點數", "n")),
            max_edge_size=int(batch_cell(row, "最大超邊階數", "k_max")),
            max_edges=int(batch_cell(row, "最大超邊數", "m_max")),
            max_degree=int(batch_cell(row, "頂點度上限", "d_max")),
            connected=bool(batch_cell(row, "二部圖連通", "connected", True)),
            forbid_pair_triangles=bool(batch_cell(row, "禁止二元三角", "forbid_pair_triangles", False)),
            sample_limit=int(batch_cell(row, "候選採樣上限", "sample_limit")),
            signature=str(
                REFINEMENT_PRESETS[str(batch_cell(row, "細化鏈 preset", "preset", "weak→medium（主表組合一）"))][
                    "fine_sig"
                ]
            ),
            delta=int(
                REFINEMENT_PRESETS[str(batch_cell(row, "細化鏈 preset", "preset", "weak→medium（主表組合一）"))][
                    "fine_delta"
                ]
            ),
            seed=int(batch_cell(row, "偽隨機基底種子", "seed")),
            refinement_enabled=True,
            refine_coarse_signature=str(
                REFINEMENT_PRESETS[str(batch_cell(row, "細化鏈 preset", "preset", "weak→medium（主表組合一）"))][
                    "coarse_sig"
                ]
            ),
            refine_coarse_delta=int(
                REFINEMENT_PRESETS[str(batch_cell(row, "細化鏈 preset", "preset", "weak→medium（主表組合一）"))][
                    "coarse_delta"
                ]
            ),
            refine_fine_signature=str(
                REFINEMENT_PRESETS[str(batch_cell(row, "細化鏈 preset", "preset", "weak→medium（主表組合一）"))][
                    "fine_sig"
                ]
            ),
            refine_fine_delta=int(
                REFINEMENT_PRESETS[str(batch_cell(row, "細化鏈 preset", "preset", "weak→medium（主表組合一）"))][
                    "fine_delta"
                ]
            ),
            refine_kernel=str(batch_cell(row, "細化核模式", "kernel")),
            refine_compare_chains=True,
            refine_coarse_sample_size=int(batch_cell(row, "粗層樣本數", "refine_coarse_sample_size", 2000)),
            refine_fine_sample_size=int(batch_cell(row, "細層樣本數", "refine_fine_sample_size", 2000)),
            refine_fiber_sample_size=int(batch_cell(row, "纖維樣本數上限", "refine_fiber_sample_size", 200)),
            progress=prog,
        )

    st.session_state["t4_batch_runs"] = run_batch_per_run_rows(
        df_r4, _row_refine, stop_on_error=False, use_progress=True
    )


def _t4_batch_display_parts(res: object, run_idx: int):
    if not isinstance(res, dict):
        return []
    df = pd.DataFrame([flatten_result_row("r4s", res)])
    return [
        (
            "單次執行扁平指標表（與本頁單次區塊同形）",
            df,
            f"table_10_4_batch_run_{run_idx}.csv",
        )
    ]


if "t4_batch_runs" in st.session_state:
    render_batch_per_run_tables(
        st.session_state["t4_batch_runs"],
        _t4_batch_display_parts,
        key_prefix="t4_batch",
        column_name_map_for_df=lambda d: build_ch10_column_name_map(d.columns),
    )
