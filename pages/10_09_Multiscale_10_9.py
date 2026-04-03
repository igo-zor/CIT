"""
§10.9 解析視窗與多尺度宏觀量（論文參數與 §10.9.5 輸出表對齊）。
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

from hypergraph_experiment.ch10_paper_presets import CH10_9_BASELINE
from hypergraph_experiment.ch10_symbol_glossary import (
    L10_CONNECTED,
    L10_D_MAX,
    L10_DELTA_T,
    L10_FORBID_TRI,
    L10_K_MAX,
    L10_M_MAX,
    L10_M_TRIAL,
    L10_N,
    L10_N_HIST,
    L10_N_SEED_RUNS_109,
    L10_OBS_SIG,
    L10_R_DEPTH,
    L10_SAMPLE_LIMIT,
    L10_SEED,
    L10_T_STEPS,
    L10_W_LIST,
    render_section_glossary,
)
from hypergraph_experiment.core import SIGNATURES
from hypergraph_experiment.experiments.exp_10_9_multiscale import (
    run_experiment_10_9,
    section_10_9_output_parameters_df,
)
from hypergraph_experiment.streamlit_common import (
    batch_cell,
    build_ch10_column_name_map,
    flatten_result_row,
    per_window_metrics_dataframe_zh,
    render_batch_per_run_tables,
    render_parameters_table,
    render_sidebar_warehouse,
    render_table_with_copy_csv,
    run_batch_per_run_rows,
    streamlit_progress_callback,
)

_B9: dict = CH10_9_BASELINE  # 論文 §10.9.6 單一真相來源
_SIG_KEYS: list[str] = sorted(SIGNATURES.keys())
_DEFAULT_SIG = str(_B9["sig_obs"])
_SIG_INDEX = _SIG_KEYS.index(_DEFAULT_SIG) if _DEFAULT_SIG in _SIG_KEYS else 1

st.set_page_config(page_title="§10.9 多尺度", layout="wide")
render_sidebar_warehouse()

st.title("§10.9　解析窗口、多尺度觀察與宏觀穩定殘影")
render_section_glossary(st, "10.9")
st.markdown(
    r"""
    固定**同一合法域**上之一條或多條微觀演化歷史，對不同時間聚合寬度 **$w$** 與步進 **$\Delta t$** 比較宏觀量；
    **$\sigma_{\mathrm{obs}}$** 決定宏觀標籤所用解析簽名（論文 §10.9.4）。
    **$C_{\mathrm{eff}}$** 之本實作固定為「窗口末態超邊數 $|E|$」代理（§10.9.2（五）允許多種操作化之一）。
    """
)

st.subheader("域型與採樣（§10.9.3）")
d1, d2, d3 = st.columns(3)
with d1:
    n = st.number_input(L10_N, 3, 16, int(_B9["n"]))
    k_max = st.number_input(L10_K_MAX, 2, 8, int(_B9["k_max"]))
    m_max = st.number_input(L10_M_MAX, 2, 30, int(_B9["m_max"]))
with d2:
    sample_limit = st.number_input(L10_SAMPLE_LIMIT, 50, 100_000, int(_B9["sample_limit"]))
    d_max = st.number_input(L10_D_MAX, 1, 20, int(_B9["d_max"]))
    conn = st.checkbox(L10_CONNECTED, bool(_B9["connected"]))
with d3:
    forbid_tri = st.checkbox(L10_FORBID_TRI, False)
    sd = st.number_input(L10_SEED, 0, 2_000_000_000, int(_B9["seed"]))

st.divider()
st.subheader("動力學與論文主變量（§10.9.4–10.9.6）")
c1, c2, c3 = st.columns(3)
with c1:
    steps = st.number_input(L10_T_STEPS, 20, 2000, int(_B9["steps"]))
    m_trial = st.number_input(L10_M_TRIAL, 1, 50, int(_B9["m_trial"]))
    ws_in = st.text_input(f"{L10_W_LIST}", str(_B9["window_list"]))
with c2:
    r = st.number_input(L10_R_DEPTH, 0, 3, int(_B9["r"]))
    sig_obs_9 = st.selectbox(L10_OBS_SIG, _SIG_KEYS, index=_SIG_INDEX)
    delta_t = st.number_input(L10_DELTA_T, 1, 20, int(_B9["delta_t"]))
with c3:
    eps_plat = st.slider(
        r"$\varepsilon_{\mathrm{plat}}^{(w)}$ — 平台閾（§10.9；施於 $C_{\mathrm{eff}}$ 代理）",
        0.001,
        0.1,
        float(_B9["epsilon_plat"]),
        step=0.001,
        help="完整語義見摺疊「論文符號對照」。",
    )
    n_hist = st.number_input(L10_N_HIST, 1, 200, int(_B9["n_hist"]))
    n_seed_runs = st.number_input(
        r"$N_{\mathrm{seed}}$ — 穩健性重跑（§10.9.6；程式鍵 n_seed_runs）",
        1,
        100,
        int(_B9["n_seed_runs"]),
    )

if st.button("執行 §10.9"):
    # 解析逗號分隔之視窗寬度，並以 Streamlit 進度條承接實驗內 progress 回呼
    ws = tuple(int(x.strip()) for x in ws_in.split(",") if x.strip())
    bar9 = st.progress(0)
    st_9 = st.empty()
    pr9 = streamlit_progress_callback(bar9, st_9)
    st.session_state["out69"] = run_experiment_10_9(
        n=int(n),
        max_edge_size=int(k_max),
        max_edges=int(m_max),
        sample_limit=int(sample_limit),
        seed=int(sd),
        connected=bool(conn),
        max_degree=int(d_max),
        forbid_pair_triangles=bool(forbid_tri),
        steps=int(steps),
        window_sizes=ws,
        r=int(r),
        sig_obs=str(sig_obs_9),
        m_trial=int(m_trial),
        delta_t=int(delta_t),
        epsilon_plat=float(eps_plat),
        n_hist=int(n_hist),
        n_seed=int(n_seed_runs),
        progress=pr9,
    )
    st.session_state["out69_params"] = {
        "n": int(n),
        "max_edge_size": int(k_max),
        "max_edges": int(m_max),
        "sample_limit": int(sample_limit),
        "seed": int(sd),
        "connected": bool(conn),
        "max_degree": int(d_max),
        "forbid_pair_triangles": bool(forbid_tri),
        "steps": int(steps),
        "window_sizes": list(ws),
        "r": int(r),
        "sig_obs": str(sig_obs_9),
        "m_trial": int(m_trial),
        "delta_t": int(delta_t),
        "epsilon_plat": float(eps_plat),
        "n_hist": int(n_hist),
        "n_seed_runs": int(n_seed_runs),
    }

out69 = st.session_state.get("out69")
if out69:
    st.subheader("§10.9.5 論文輸出參數對照表")
    df_t9_paper = section_10_9_output_parameters_df(out69)
    render_table_with_copy_csv(
        df_t9_paper,
        key_prefix="t9_paper_outputs",
        csv_filename="table_10_9_section_10_9_5_outputs.csv",
        hide_index=True,
    )
    pw = out69.get("per_window")
    if pw:
        st.subheader("各視窗宏觀指標（對照論文 §10.9.5；欄名為中文說明）")
        # 與其他表同形：一鍵複製 Markdown + 下載 CSV（欄名已為中文，不重複映射）
        df_pw_zh = per_window_metrics_dataframe_zh(pd.DataFrame(pw))
        render_table_with_copy_csv(
            df_pw_zh,
            key_prefix="t9_per_window",
            csv_filename="table_10_9_per_window.csv",
            hide_index=True,
        )
    st.subheader("單次執行指標（扁平欄位，可複製／下載 CSV）")
    df_t9_single = pd.DataFrame([flatten_result_row("e9s", out69)])
    render_table_with_copy_csv(
        df_t9_single,
        key_prefix="t9_single",
        csv_filename="table_10_9_single.csv",
        column_name_map=build_ch10_column_name_map(df_t9_single.columns),
    )
    # 優先顯示本次執行之表單參數；若無則回退至結果內 parameters
    _params = st.session_state.get("out69_params") or out69.get("parameters")
    render_parameters_table(
        _params,
        key_prefix="t9_single",
        csv_filename="table_10_9_single_params.csv",
    )
    with st.expander("完整 JSON（除錯）"):
        st.json(out69)

st.caption(
    "論文 **N_seed=20**（穩健性）見 §10.9.6；若與 **N_hist=30** 同開可能較耗時，請視需求調低。"
)

st.subheader("批次 §10.9")
st.caption(
    "欄位與單次表單對齊：域型、演化步數、視窗清單、解析簽名、Δt、ε_plat、N_hist、N_seed_runs 等。"
)


def _default_batch_row() -> dict:
    return {
        "節點數": int(_B9["n"]),
        "最大超邊階數": int(_B9["k_max"]),
        "最大超邊數": int(_B9["m_max"]),
        "候選採樣上限": int(_B9["sample_limit"]),
        "頂點度上限": int(_B9["d_max"]),
        "二部圖連通": bool(_B9["connected"]),
        "禁止二元三角": False,
        "演化步數": int(_B9["steps"]),
        "每步候選更新數": int(_B9["m_trial"]),
        "視窗寬度清單": str(_B9["window_list"]),
        "鄰域深度": int(_B9["r"]),
        "整體觀測簽名": str(_B9["sig_obs"]),
        "聚合步長 Δt": int(_B9["delta_t"]),
        "宏觀平台閾 ε_plat": float(_B9["epsilon_plat"]),
        "微觀歷史條數": int(_B9["n_hist"]),
        "穩健性重跑次數": int(_B9["n_seed_runs"]),
        "偽隨機基底種子": int(_B9["seed"]),
    }


if st.button("載入論文建議批次模板（§10.9）", key="load_t9_template"):
    st.session_state["t9_batch_template_df"] = pd.DataFrame(
        [
            _default_batch_row(),
            {
                **_default_batch_row(),
                "偽隨機基底種子": int(_B9["seed"]) + 1,
                "演化步數": 500,
                "視窗寬度清單": "1,4,8",
            },
        ]
    )
dfb = st.data_editor(
    st.session_state.get(
        "t9_batch_template_df",
        pd.DataFrame([_default_batch_row()]),
    ),
    num_rows="dynamic",
    key="b69",
)


def _t9_batch_display_parts(res: object, run_idx: int):
    """與單次區塊同序：論文 §10.9.5 表、per_window、扁平欄位表。"""
    if not isinstance(res, dict):
        return []
    parts: list = []
    df_p = section_10_9_output_parameters_df(res)
    parts.append(
        (
            "§10.9.5 論文輸出參數對照表",
            df_p,
            f"table_10_9_paper_outputs_batch_run_{run_idx}.csv",
            True,
        )
    )
    pw = res.get("per_window")
    if pw:
        dfw = per_window_metrics_dataframe_zh(pd.DataFrame(pw))
        parts.append(
            (
                "各視窗宏觀指標（中文欄名）",
                dfw,
                f"table_10_9_per_window_batch_run_{run_idx}.csv",
                True,
            )
        )
    dff = pd.DataFrame([flatten_result_row("e9s", res)])
    parts.append(
        (
            "單次執行扁平指標（與本頁單次區塊同形）",
            dff,
            f"table_10_9_flat_batch_run_{run_idx}.csv",
        )
    )
    return parts


if st.button("批次 §10.9"):

    def _r(row: pd.Series, prog):
        ws_raw = batch_cell(row, "視窗寬度清單", "window_sizes", str(_B9["window_list"]))
        ws = tuple(int(x.strip()) for x in str(ws_raw).split(",") if x.strip())
        return run_experiment_10_9(
            n=int(batch_cell(row, "節點數", "n")),
            max_edge_size=int(batch_cell(row, "最大超邊階數", "k_max", _B9["k_max"])),
            max_edges=int(batch_cell(row, "最大超邊數", "m_max", _B9["m_max"])),
            sample_limit=int(batch_cell(row, "候選採樣上限", "sample_limit", _B9["sample_limit"])),
            seed=int(batch_cell(row, "偽隨機基底種子", "seed")),
            connected=bool(batch_cell(row, "二部圖連通", "connected", _B9["connected"])),
            max_degree=int(batch_cell(row, "頂點度上限", "d_max", _B9["d_max"])),
            forbid_pair_triangles=bool(batch_cell(row, "禁止二元三角", "forbid_pair_triangles", False)),
            steps=int(batch_cell(row, "演化步數", "steps")),
            window_sizes=ws,
            r=int(batch_cell(row, "鄰域深度", "r", _B9["r"])),
            sig_obs=str(batch_cell(row, "整體觀測簽名", "sig_obs", _B9["sig_obs"])),
            m_trial=int(batch_cell(row, "每步候選更新數", "m_trial", _B9["m_trial"])),
            delta_t=int(batch_cell(row, "聚合步長 Δt", "delta_t", _B9["delta_t"])),
            epsilon_plat=float(batch_cell(row, "宏觀平台閾 ε_plat", "epsilon_plat", _B9["epsilon_plat"])),
            n_hist=int(batch_cell(row, "微觀歷史條數", "n_hist", _B9["n_hist"])),
            n_seed=int(batch_cell(row, "穩健性重跑次數", "n_seed_runs", _B9["n_seed_runs"])),
            progress=prog,
        )

    st.session_state["t9_batch_runs"] = run_batch_per_run_rows(
        dfb, _r, stop_on_error=False, use_progress=True
    )

if "t9_batch_runs" in st.session_state:
    render_batch_per_run_tables(
        st.session_state["t9_batch_runs"],
        _t9_batch_display_parts,
        key_prefix="t9_batch",
        column_name_map_for_df=lambda d: build_ch10_column_name_map(d.columns),
    )
