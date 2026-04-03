"""
§10.7 超圖動力學：實驗 A（單路徑）與實驗 B（雙解析路徑終端比較）。
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

from hypergraph_experiment.ch10_paper_presets import (
    CH10_7_BASELINE,
    CH10_7_N0_CHOICES,
    CH10_7_T_CHOICES,
)
from hypergraph_experiment.ch10_section_10_7_tables import (
    build_experiment_10_7_b_output_table,
    build_section_10_7_5_output_table,
)
from hypergraph_experiment.ch10_symbol_glossary import (
    L10_DELTA,
    L10_DELTA_COARSE,
    L10_DELTA_FINE,
    L10_K_MAX,
    L10_M_MAX,
    L10_N,
    L10_N0_RUNS,
    L10_OBS_SIG,
    L10_SAMPLE_LIMIT,
    L10_SEED,
    L10_SIG_COARSE,
    L10_SIG_FINE,
    L10_T_STEPS,
    render_section_glossary,
)
from hypergraph_experiment.core import SIGNATURES, run_full_experiment
from hypergraph_experiment.experiments.exp_10_7_paths import run_experiment_10_7
from hypergraph_experiment.streamlit_common import (
    batch_cell,
    build_ch10_column_name_map,
    flatten_result_row,
    render_batch_per_run_tables,
    render_parameters_table,
    render_result_charts,
    render_sidebar_warehouse,
    render_table_with_copy_csv,
    run_batch_per_run_rows,
    streamlit_progress_callback,
    zip_download_hypergraph_run,
)

# 論文 §10.7.6 主線預設（單一真相來源）
_B = CH10_7_BASELINE
_SIG_KEYS: list[str] = sorted(SIGNATURES.keys())
_DEFAULT_SIG = str(_B["signature"])
_SIG_INDEX = _SIG_KEYS.index(_DEFAULT_SIG) if _DEFAULT_SIG in _SIG_KEYS else 0


def _idx_in_tuple(t: tuple[int, ...], val: int, fallback: int = 0) -> int:
    """回傳 val 在有序元組 t 中的索引；不存在則回傳 fallback。"""
    try:
        return t.index(int(val))
    except ValueError:
        return fallback


st.set_page_config(page_title="§10.7 動力學", layout="wide")
render_sidebar_warehouse()

st.title("§10.7　動力學與解析熵")
render_section_glossary(st, "10.7")
st.markdown(
    r"""
    **實驗 A**（§10.7 主線）：單一解析觀測層 $\Lambda_{\mathrm{obs}}$ 下之軌道與熵序列
    $H_\Lambda^{(\ell)}$，對齊論文 §10.7.5 輸出表。
    **實驗 B**（擴充）：同一局部更新規則下，兩組 $(\mathrm{Sig}_\Lambda,\delta)$ 之**終端**類別分布比較（JS 等），
    非 §10.7.5 主表，另附對照表。
    """
)

# ---------------------------------------------------------------------------
# 實驗 A：單路徑動力學
# ---------------------------------------------------------------------------
st.header("實驗 A：單路徑動力學（§10.7 主線）")

with st.form("form_10_7_experiment_a"):
    st.subheader("域型與採樣（程式實作；論文固定合法配置域）")
    a1, a2 = st.columns(2)
    with a1:
        n_a = st.number_input(L10_N, 2, 12, int(_B["n"]), key="a7_n")
        k_max_a = st.number_input(L10_K_MAX, 2, 6, int(_B["k_max"]), key="a7_k")
        m_max_a = st.number_input(L10_M_MAX, 1, 20, int(_B["m_max"]), key="a7_m")
        d_max_a = st.number_input(
            "$d_{max}$ — 頂點最大允許度數（域型；§10.7）",
            1,
            20,
            int(_B["d_max"]),
            key="a7_d",
        )
    with a2:
        sl_a = st.number_input(
            L10_SAMPLE_LIMIT, 200, 50_000, int(_B["sample_limit"]), key="a7_sl"
        )
        conn_a = st.checkbox("2-section 連通｜域型條件", True, key="a7_conn")
        ft_a = st.checkbox("禁二元△｜forbidden motif", False, key="a7_ft")

    st.subheader("論文變量（§10.7.4）與建議主線（§10.7.6；部分僅記錄於參數表）")
    a3, a4, a5 = st.columns(3)
    with a3:
        use_paper_t = st.checkbox(
            "$T$ 使用論文建議離散集合", True, key="a7_use_paper_t"
        )
        if use_paper_t:
            steps_a = st.selectbox(
                L10_T_STEPS,
                list(CH10_7_T_CHOICES),
                index=_idx_in_tuple(CH10_7_T_CHOICES, int(_B["steps"]), 2),
                key="a7_steps_sel",
            )
        else:
            steps_a = st.number_input(
                L10_T_STEPS, 5, 2000, int(_B["steps"]), key="a7_steps_num"
            )
        use_paper_n0 = st.checkbox(
            "$N_0$ 使用論文建議離散集合", True, key="a7_use_paper_n0"
        )
        if use_paper_n0:
            runs_a = st.selectbox(
                L10_N0_RUNS,
                list(CH10_7_N0_CHOICES),
                index=_idx_in_tuple(CH10_7_N0_CHOICES, int(_B["runs"]), 1),
                key="a7_runs_sel",
            )
        else:
            runs_a = st.number_input(
                L10_N0_RUNS, 1, 300, int(_B["runs"]), key="a7_runs_num"
            )
    with a4:
        sig_a = st.selectbox(
            L10_OBS_SIG, _SIG_KEYS, index=_SIG_INDEX, key="a7_sig"
        )
        delta_a = st.number_input(
            L10_DELTA, 0, 50, int(_B["delta"]), key="a7_delta"
        )
        sd_a = st.number_input(
            L10_SEED, 0, 2_000_000_000, int(_B["seed"]), key="a7_seed"
        )
    with a5:
        eps_plat_a = st.number_input(
            r"$\varepsilon_{\mathrm{plat}}$ — 熵平台判定（|ΔH|≤ε；§10.7.3）",
            0.0,
            1.0,
            float(_B["eps_plat"]),
            step=0.005,
            key="a7_eps",
        )
        st.caption("下列欄位對齊 §10.7.6 建議表；已接線至動力學模擬與彙總輸出。")
        m_trial_a = st.number_input(
            "$M_{\\mathrm{trial}}$ — 每步候選更新數（已接線）",
            1,
            100,
            int(_B["m_trial"]),
            key="a7_mt",
        )
        w_h_a = st.number_input(
            "$w_H$ — 熵滑動窗口", 1, 100, int(_B["w_h"]), key="a7_wh"
        )
        w_a_a = st.number_input(
            "$w_A$ — 吸引子／週期檢測窗口", 1, 100, int(_B["w_a"]), key="a7_wa"
        )
        p_max_a = st.number_input(
            "$P_{\\max}$ — 最大週期長度", 1, 100, int(_B["p_max"]), key="a7_pm"
        )
        n_seed_107_a = st.number_input(
            "$N_{\\mathrm{seed}}$ — §10.7.6 重跑次數（已接線）",
            1,
            200,
            int(_B["n_seed_107"]),
            key="a7_ns",
        )

    submitted_a = st.form_submit_button("執行實驗 A（單路徑動力學）")

if submitted_a:
    bar = st.progress(0)
    pr = streamlit_progress_callback(bar, st.empty())
    try:
        res = run_full_experiment(
            mode="dynamics",
            n=int(n_a),
            max_edge_size=int(k_max_a),
            max_edges=int(m_max_a),
            max_degree=int(d_max_a),
            connected=bool(conn_a),
            forbid_pair_triangles=bool(ft_a),
            sample_limit=int(sl_a),
            signature=sig_a,
            delta=int(delta_a),
            runs=int(runs_a),
            steps=int(steps_a),
            m_trial=int(m_trial_a),
            w_h=int(w_h_a),
            w_a=int(w_a_a),
            p_max=int(p_max_a),
            n_seed_107=int(n_seed_107_a),
            seed=int(sd_a),
            epsilon_plat=float(eps_plat_a),
            refinement_enabled=False,
            progress=pr,
        )
        st.session_state["res_dy"] = res
        st.session_state["res_dy_params"] = {
            "mode": "dynamics",
            "n": int(n_a),
            "max_edge_size": int(k_max_a),
            "max_edges": int(m_max_a),
            "max_degree": int(d_max_a),
            "connected": bool(conn_a),
            "forbid_pair_triangles": bool(ft_a),
            "sample_limit": int(sl_a),
            "signature": str(sig_a),
            "delta": int(delta_a),
            "runs": int(runs_a),
            "steps": int(steps_a),
            "seed": int(sd_a),
            "epsilon_plat": float(eps_plat_a),
            "m_trial": int(m_trial_a),
            "w_h": int(w_h_a),
            "w_a": int(w_a_a),
            "p_max": int(p_max_a),
            "n_seed_107": int(n_seed_107_a),
        }
        st.success("實驗 A 完成")
    except Exception as e:
        st.error(str(e))

res = st.session_state.get("res_dy")
if res:
    _an7 = res.get("analysis") or {}
    if "error" not in _an7:
        m1, m2, m3, m4 = st.columns(4)
        m1.metric(
            "合法更新步占比（軌跡平均）",
            f'{float(_an7.get("legal_update_step_fraction_mean", 0.0)):.4f}',
        )
        es = _an7.get("entropy_summary") or {}
        m2.metric(
            "熵平台最長段長",
            es.get("plateau_max_length", "—"),
        )
        m3.metric(
            "熵平台段數",
            es.get("plateau_num_segments", "—"),
        )
        m4.metric(
            r"$N_{\mathrm{reach}}$（軌跡平均）",
            f'{_an7.get("n_reach_mean", "—")}',
        )
    with st.expander("接線自檢（顯示實際採用之參數與重跑彙總）", expanded=False):
        st.markdown(
            f"- analysis.m_trial = `{_an7.get('m_trial', '—')}`\n"
            f"- analysis.w_h = `{_an7.get('w_h', '—')}`\n"
            f"- analysis.w_a = `{_an7.get('w_a', '—')}`\n"
            f"- analysis.p_max = `{_an7.get('p_max', '—')}`\n"
            f"- result.analysis_seed_summary exists = `{bool(res.get('analysis_seed_summary'))}`\n"
            f"- result.analysis_seeds count = `{len(res.get('analysis_seeds') or [])}`\n"
        )
        if res.get("analysis_seed_summary"):
            st.json(res.get("analysis_seed_summary"))
    st.subheader("§10.7.5 論文輸出參數對照表")
    df_paper_a = build_section_10_7_5_output_table(_an7)
    render_table_with_copy_csv(
        df_paper_a,
        key_prefix="t7_paper_10_7_5",
        csv_filename="table_10_7_section_10_7_5_outputs.csv",
        hide_index=True,
    )
    st.subheader("單次動力學指標（扁平欄位，可複製／下載 CSV）")
    df_t7_dyn_single = pd.DataFrame([flatten_result_row("d7s", res)])
    render_table_with_copy_csv(
        df_t7_dyn_single,
        key_prefix="t7_dyn_single",
        csv_filename="table_10_7_dynamics_single.csv",
        column_name_map=build_ch10_column_name_map(df_t7_dyn_single.columns),
    )
    render_parameters_table(
        st.session_state.get("res_dy_params"),
        key_prefix="t7_dyn_single",
        csv_filename="table_10_7_dynamics_single_params.csv",
    )
    render_result_charts(res, key_prefix="d7")
    st.download_button(
        "下載 ZIP",
        zip_download_hypergraph_run(res, "d.zip"),
        "dynamics_run.zip",
        "application/zip",
    )

# ---------------------------------------------------------------------------
# 實驗 B：雙解析路徑
# ---------------------------------------------------------------------------
st.header("實驗 B：雙解析路徑終端比較（§10.7 擴充）")
st.caption(
    "與實驗 A 共用同一類型之域型與採樣設定；終端分布差異見下表（非 §10.7.5 主線）。"
)

with st.form("form_10_7_experiment_b"):
    b1, b2 = st.columns(2)
    with b1:
        n_b = st.number_input(L10_N, 2, 12, int(_B["n"]), key="b7_n")
        k_max_b = st.number_input(L10_K_MAX, 2, 6, int(_B["k_max"]), key="b7_k")
        m_max_b = st.number_input(L10_M_MAX, 1, 20, int(_B["m_max"]), key="b7_m")
        d_max_b = st.number_input(
            "$d_{max}$ — 頂點最大允許度數（域型；§10.7）",
            1,
            20,
            int(_B["d_max"]),
            key="b7_d",
        )
    with b2:
        sl_b = st.number_input(
            L10_SAMPLE_LIMIT, 200, 50_000, int(_B["sample_limit"]), key="b7_sl"
        )
        conn_b = st.checkbox("2-section 連通｜域型條件", True, key="b7_conn")
        ft_b = st.checkbox("禁二元△｜forbidden motif", False, key="b7_ft")

    b3, b4, b5 = st.columns(3)
    with b3:
        use_paper_t_b = st.checkbox(
            "$T$ 使用論文建議離散集合", True, key="b7_use_paper_t"
        )
        if use_paper_t_b:
            st_b = st.selectbox(
                L10_T_STEPS,
                list(CH10_7_T_CHOICES),
                index=_idx_in_tuple(CH10_7_T_CHOICES, int(_B["steps"]), 2),
                key="b7_steps_sel",
            )
        else:
            st_b = st.number_input(
                L10_T_STEPS, 5, 2000, int(_B["steps"]), key="b7_steps_num"
            )
        use_paper_n0_b = st.checkbox(
            "$N_0$ 使用論文建議離散集合", True, key="b7_use_paper_n0"
        )
        if use_paper_n0_b:
            runs_b = st.selectbox(
                L10_N0_RUNS,
                list(CH10_7_N0_CHOICES),
                index=_idx_in_tuple(CH10_7_N0_CHOICES, int(_B["runs"]), 1),
                key="b7_runs_sel",
            )
        else:
            runs_b = st.number_input(
                L10_N0_RUNS, 1, 300, int(_B["runs"]), key="b7_runs_num"
            )
    with b4:
        sd_b = st.number_input(
            L10_SEED, 0, 2_000_000_000, int(_B["seed"]), key="b7_seed"
        )
    with b5:
        sa = st.selectbox(
            "路徑 A｜" + L10_SIG_COARSE, _SIG_KEYS, index=0, key="b7_sa"
        )
        da = st.number_input(
            "路徑 A｜" + L10_DELTA_COARSE, 0, 50, 3, key="b7_da"
        )
        sb = st.selectbox(
            "路徑 B｜" + L10_SIG_FINE, _SIG_KEYS, index=min(2, len(_SIG_KEYS) - 1), key="b7_sb"
        )
        db = st.number_input(
            "路徑 B｜" + L10_DELTA_FINE, 0, 50, 0, key="b7_db"
        )

    submitted_b = st.form_submit_button("執行實驗 B（雙路徑終端比較）")

if submitted_b:
    bar = st.progress(0)
    pr = streamlit_progress_callback(bar, st.empty())
    try:
        st.session_state["cmp7"] = run_experiment_10_7(
            n=int(n_b),
            max_edge_size=int(k_max_b),
            max_edges=int(m_max_b),
            sample_limit=int(sl_b),
            seed=int(sd_b),
            runs=int(runs_b),
            steps=int(st_b),
            sig_path_a=sa,
            delta_a=int(da),
            sig_path_b=sb,
            delta_b=int(db),
            connected=bool(conn_b),
            max_degree=int(d_max_b),
            forbid_pair_triangles=bool(ft_b),
            progress=pr,
        )
        st.session_state["cmp7_params"] = {
            "n": int(n_b),
            "max_edge_size": int(k_max_b),
            "max_edges": int(m_max_b),
            "max_degree": int(d_max_b),
            "connected": bool(conn_b),
            "forbid_pair_triangles": bool(ft_b),
            "sample_limit": int(sl_b),
            "seed": int(sd_b),
            "runs": int(runs_b),
            "steps": int(st_b),
            "sig_path_a": str(sa),
            "delta_a": int(da),
            "sig_path_b": str(sb),
            "delta_b": int(db),
        }
        st.success("實驗 B 完成")
    except Exception as e:
        st.error(str(e))

cmp7 = st.session_state.get("cmp7")
if cmp7:
    st.subheader("實驗 B：擴充輸出對照表（非 §10.7.5 主表）")
    df_paper_b = build_experiment_10_7_b_output_table(cmp7)
    render_table_with_copy_csv(
        df_paper_b,
        key_prefix="t7_paper_b",
        csv_filename="table_10_7_experiment_b_outputs.csv",
        hide_index=True,
    )
    st.subheader("雙路徑比較指標（扁平欄位，可複製／下載 CSV）")
    df_t7_path_single = pd.DataFrame([flatten_result_row("p7", cmp7)])
    render_table_with_copy_csv(
        df_t7_path_single,
        key_prefix="t7_path_cmp_single",
        csv_filename="table_10_7_paths_single.csv",
        column_name_map=build_ch10_column_name_map(df_t7_path_single.columns),
    )
    render_parameters_table(
        st.session_state.get("cmp7_params"),
        key_prefix="t7_path_cmp_single",
        csv_filename="table_10_7_paths_single_params.csv",
    )
    with st.expander("完整 JSON（除錯）"):
        st.json(cmp7)

# ---------------------------------------------------------------------------
# 實驗 A：批次掃描
# ---------------------------------------------------------------------------
st.header("實驗 A：參數掃描批次")
st.caption(
    "欄位對應實驗 A 之已接線參數：節點數、超邊與邊數上限、候選採樣上限、解析簽名、"
    "解析閾值、軌道數、步數、偽隨機基底種子、域型條件、熵平台閾值。"
    "（$M_{\\mathrm{trial}}$、$w_H$ 等僅見單次表單，批次未接線者不列入。）"
)
if st.button("載入論文建議批次模板（§10.7）", key="load_t7_template"):
    st.session_state["t7_batch_template_df"] = pd.DataFrame(
        [
            {
                "節點數": int(_B["n"]),
                "最大超邊階數": int(_B["k_max"]),
                "最大超邊數": int(_B["m_max"]),
                "候選採樣上限": int(_B["sample_limit"]),
                "解析簽名": str(_B["signature"]),
                "解析閾值整數": int(_B["delta"]),
                "軌道數": int(_B["runs"]),
                "演化步數": int(_B["steps"]),
                "偽隨機基底種子": int(_B["seed"]),
                "頂點度上限": int(_B["d_max"]),
                "二部圖連通": True,
                "禁止二元三角": False,
                "熵平台判定閾值": float(_B["eps_plat"]),
            },
            {
                "節點數": int(_B["n"]),
                "最大超邊階數": int(_B["k_max"]),
                "最大超邊數": int(_B["m_max"]),
                "候選採樣上限": int(_B["sample_limit"]),
                "解析簽名": str(_B["signature"]),
                "解析閾值整數": int(_B["delta"]),
                "軌道數": int(_B["runs"]),
                "演化步數": 500,
                "偽隨機基底種子": int(_B["seed"]) + 1,
                "頂點度上限": int(_B["d_max"]),
                "二部圖連通": True,
                "禁止二元三角": False,
                "熵平台判定閾值": float(_B["eps_plat"]),
            },
        ]
    )
dfd = st.data_editor(
    st.session_state.get(
        "t7_batch_template_df",
        pd.DataFrame(
            [
                {
                    "節點數": int(_B["n"]),
                    "最大超邊階數": int(_B["k_max"]),
                    "最大超邊數": int(_B["m_max"]),
                    "候選採樣上限": int(_B["sample_limit"]),
                    "解析簽名": str(_B["signature"]),
                    "解析閾值整數": int(_B["delta"]),
                    "軌道數": int(_B["runs"]),
                    "演化步數": int(_B["steps"]),
                    "偽隨機基底種子": int(_B["seed"]),
                    "頂點度上限": int(_B["d_max"]),
                    "二部圖連通": True,
                    "禁止二元三角": False,
                    "熵平台判定閾值": float(_B["eps_plat"]),
                }
            ]
        ),
    ),
    num_rows="dynamic",
    key="bd7",
)


def _t7_dyn_batch_display_parts(res: object, run_idx: int):
    if not isinstance(res, dict):
        return []
    df = pd.DataFrame([flatten_result_row("d7s", res)])
    return [
        (
            "單次動力學指標（與實驗 A 單次區塊同形）",
            df,
            f"table_10_7_dynamics_batch_run_{run_idx}.csv",
        )
    ]


if st.button("批次動力學（實驗 A）", key="batch_dyn_7"):

    def _r(r: pd.Series, prog):
        return run_full_experiment(
            mode="dynamics",
            n=int(batch_cell(r, "節點數", "n")),
            max_edge_size=int(batch_cell(r, "最大超邊階數", "k_max")),
            max_edges=int(batch_cell(r, "最大超邊數", "m_max")),
            max_degree=int(batch_cell(r, "頂點度上限", "max_degree", int(_B["d_max"]))),
            connected=bool(batch_cell(r, "二部圖連通", "connected", True)),
            forbid_pair_triangles=bool(
                batch_cell(r, "禁止二元三角", "forbid_pair_triangles", False)
            ),
            sample_limit=int(batch_cell(r, "候選採樣上限", "sample_limit")),
            signature=str(batch_cell(r, "解析簽名", "signature")),
            delta=int(batch_cell(r, "解析閾值整數", "delta")),
            runs=int(batch_cell(r, "軌道數", "runs")),
            steps=int(batch_cell(r, "演化步數", "steps")),
            seed=int(batch_cell(r, "偽隨機基底種子", "seed")),
            epsilon_plat=float(
                batch_cell(r, "熵平台判定閾值", "epsilon_plat", float(_B["eps_plat"]))
            ),
            progress=prog,
        )

    st.session_state["t7_dyn_batch_runs"] = run_batch_per_run_rows(
        dfd, _r, stop_on_error=False, use_progress=True
    )

if "t7_dyn_batch_runs" in st.session_state:
    render_batch_per_run_tables(
        st.session_state["t7_dyn_batch_runs"],
        _t7_dyn_batch_display_parts,
        key_prefix="t7_batch",
        column_name_map_for_df=lambda d: build_ch10_column_name_map(d.columns),
    )
