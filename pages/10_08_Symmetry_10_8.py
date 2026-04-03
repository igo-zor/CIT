"""
§10.8 對稱初態與局部破缺。
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

from hypergraph_experiment.ch10_symbol_glossary import (
    L10_ETA,
    L10_INIT_FAMILY,
    L10_M_EDGES,
    L10_N,
    L10_N_SAMPLES,
    L10_R_DEPTH,
    L10_SEED,
    L10_T_SB,
    render_section_glossary,
)
from hypergraph_experiment.core import SIGNATURES
from hypergraph_experiment.experiments.exp_10_8_symmetry import (
    run_experiment_10_8,
    run_experiment_10_8_three_arm,
    section_10_8_output_parameters_df,
)
from hypergraph_experiment.ch10_paper_presets import CH10_8_BASELINE
from hypergraph_experiment.streamlit_common import (
    batch_cell,
    build_ch10_column_name_map,
    flatten_result_row,
    render_batch_per_run_tables,
    render_parameters_table,
    render_sidebar_warehouse,
    render_table_with_copy_csv,
    run_batch_per_run_rows,
    streamlit_progress_callback,
)

st.set_page_config(page_title="§10.8 對稱破缺", layout="wide")
render_sidebar_warehouse()

st.title("§10.8　對稱初態與可比性擴張")
render_section_glossary(st, "10.8")
st.markdown(
    r"""
    比較 **sym / pert / rand** 三類初態之局部型別數與型別熵（**$N_{\mathrm{type}}$** 相關）。
    符號對照：$\eta$=微擾強度、$T_{\mathrm{sb}}$=對稱破缺後演化步數、$r$=局部鄰域深度。
    **建議**：$n=12$，$m=18$，$r=2$，掃描 **$\eta$** 與 **$T_{\mathrm{sb}}$**（破缺步數可對應本頁動力學步數）。
    """
)

st.subheader("建議固定參數（可覆寫）")
c1, c2 = st.columns(2)
with c1:
    n = st.number_input(L10_N, 4, 16, int(CH10_8_BASELINE["n"]))
    m = st.number_input(L10_M_EDGES, 4, 30, int(CH10_8_BASELINE["m"]))
    r = st.number_input(L10_R_DEPTH, 1, 3, int(CH10_8_BASELINE["r"]))
with c2:
    mode_108 = st.radio(
        "執行模式",
        ["三臂並列（sym／pert／rand）", "單一族別"],
        index=0,
        horizontal=True,
        help="三臂：同參數連跑三類初態族並給出 ρ_type（pert∥sym）與 rand 輔助對照。",
    )
    fam = "sym"
    if mode_108 == "單一族別":
        fam = st.selectbox(L10_INIT_FAMILY, ["sym", "pert", "rand"])
    ns = st.number_input(L10_N_SAMPLES, 1, 80, int(CH10_8_BASELINE["n_samples"]))
    sd8 = st.number_input(L10_SEED, 0, 2_000_000_000, int(CH10_8_BASELINE["seed"]), key="t8sd")
st.divider()
st.subheader("建議變量掃描參數（可覆寫）")
v1, v2 = st.columns(2)
with v1:
    sig_obs_8 = st.selectbox(
        "整體觀測簽名 sig_obs",
        list(sorted(SIGNATURES.keys())),
        index=1,
    )
    eta = st.slider(L10_ETA, 0.0, 0.5, float(CH10_8_BASELINE["eta"]))
with v2:
    Tsb = st.number_input(L10_T_SB, 0, 200, int(CH10_8_BASELINE["T_sb"]))

if st.button("執行 §10.8"):
    # 單次執行與 §10.7 一致：建立 Streamlit 進度條並傳入實驗回呼（三臂時總刻度為 3×n_samples）
    bar_108 = st.progress(0)
    status_108 = st.empty()
    pr_108 = streamlit_progress_callback(bar_108, status_108)
    if mode_108 == "三臂並列（sym／pert／rand）":
        out = run_experiment_10_8_three_arm(
            n=int(n),
            m=int(m),
            eta=float(eta),
            T_sb=int(Tsb),
            r=int(r),
            sig_obs=str(sig_obs_8),
            n_samples=int(ns),
            dynamics_steps=int(Tsb),
            seed=int(sd8),
            progress=pr_108,
        )
    else:
        out = run_experiment_10_8(
            n=int(n),
            m=int(m),
            init_family=fam,
            eta=float(eta),
            T_sb=int(Tsb),
            r=int(r),
            sig_obs=str(sig_obs_8),
            n_samples=int(ns),
            dynamics_steps=int(Tsb),
            seed=int(sd8),
            progress=pr_108,
        )
    st.session_state["o68"] = out
    st.session_state["o68_params"] = {
        "執行模式": str(mode_108),
        "n": int(n),
        "m": int(m),
        "init_family": str(fam) if mode_108 == "單一族別" else "sym+pert+rand",
        "eta": float(eta),
        "T_sb": int(Tsb),
        "r": int(r),
        "sig_obs": str(sig_obs_8),
        "n_samples": int(ns),
        "dynamics_steps": int(Tsb),
        "seed": int(sd8),
    }

o68 = st.session_state.get("o68")
if o68:
    st.subheader("§10.8.5 論文輸出參數對照表")
    df_t8_paper = section_10_8_output_parameters_df(o68)
    render_table_with_copy_csv(
        df_t8_paper,
        key_prefix="t8_paper_outputs",
        csv_filename="table_10_8_section_10_8_5_outputs.csv",
        hide_index=True,
    )
    st.subheader("單次執行指標（扁平欄位，可複製／下載 CSV）")
    df_t8_single = pd.DataFrame([flatten_result_row("e8s", o68)])
    render_table_with_copy_csv(
        df_t8_single,
        key_prefix="t8_single",
        csv_filename="table_10_8_single.csv",
        column_name_map=build_ch10_column_name_map(df_t8_single.columns),
    )
    render_parameters_table(
        st.session_state.get("o68_params"),
        key_prefix="t8_single",
        csv_filename="table_10_8_single_params.csv",
    )
    if o68.get("experiment") == "10.8_three_arm":
        comp = o68.get("comparison") or {}
        st.subheader("三臂並列摘要（comparison）")
        st.dataframe(pd.DataFrame([comp]), use_container_width=True)
        z0, z1, z2, z3, z4 = st.columns(5)
        z0.metric("N_type 均值 sym", comp.get("N_type_mean_sym", "—"))
        z1.metric("N_type 均值 pert", comp.get("N_type_mean_pert", "—"))
        z2.metric("N_type 均值 rand", comp.get("N_type_mean_rand", "—"))
        z3.metric("ρ_type（pert∥sym）", comp.get("rho_type_pert_over_sym", "—"))
        z4.metric("輔助 ρ（rand∥sym）", comp.get("rho_type_rand_over_sym", "—"))
        st.caption(
            "論文 **ρ_type** 僅對應 pert 與 sym；rand 為第三對照臂。"
            " 各臂種子見結果 JSON 內 parameters.seeds_by_family。"
        )
        with st.expander("各臂 per_sample（除錯；資料量大）"):
            for arm in ("sym", "pert", "rand"):
                bf = (o68.get("by_family") or {}).get(arm)
                if isinstance(bf, dict):
                    st.markdown(f"**{arm}**")
                    st.json(bf.get("per_sample", []))
    else:
        _m8 = o68.get("metrics") or {}
        z1, z2, z3 = st.columns(3)
        z1.metric("N_iso 平均（代理）", _m8.get("N_iso_mean", "—"))
        z2.metric("A_reach 可達對平均距離", _m8.get("A_reach_mean", "—"))
        z3.metric("可達無序對比例平均", _m8.get("A_reach_pair_fraction_mean", "—"))
        st.caption(
            "單一族別時，ρ_type 請以 sym 與 pert 各跑一次後以 N_type_mean 代入 "
            "``rho_type_expansion``，或改用上方「三臂並列」模式一次取得。"
        )
    with st.expander("完整 JSON（除錯）"):
        st.json(o68)

st.subheader("三族對照相批次（建議各跑一列後自行彙總 $N_{type}$／熵）")
st.caption(
    "欄位含固定與掃描參數：初態族別、節點數、邊數、鄰域深度、微擾強度、破缺步數、樣本數、演化步數與偽隨機基底種子。"
)
if st.button("載入論文建議批次模板（§10.8）", key="load_t8_template"):
    st.session_state["t8_batch_template_df"] = pd.DataFrame(
        [
            {
                "初態族別": "sym",
                "節點數": 12,
                "超邊數": 18,
                "鄰域深度": 2,
                "微擾強度": 0.10,
                "對稱破缺步數占位": 20,
                "重複樣本數": 30,
                "動力學步數": 20,
                "整體觀測簽名": "medium",
                "偽隨機基底種子": 20,
            },
            {
                "初態族別": "pert",
                "節點數": 12,
                "超邊數": 18,
                "鄰域深度": 2,
                "微擾強度": 0.10,
                "對稱破缺步數占位": 20,
                "重複樣本數": 30,
                "動力學步數": 20,
                "整體觀測簽名": "medium",
                "偽隨機基底種子": 21,
            },
            {
                "初態族別": "rand",
                "節點數": 12,
                "超邊數": 18,
                "鄰域深度": 2,
                "微擾強度": 0.20,
                "對稱破缺步數占位": 50,
                "重複樣本數": 30,
                "動力學步數": 50,
                "整體觀測簽名": "medium",
                "偽隨機基底種子": 22,
            },
        ]
    )
df = st.data_editor(
    st.session_state.get(
        "t8_batch_template_df",
        pd.DataFrame(
            [
                {
                    "初態族別": "sym",
                    "節點數": 12,
                    "超邊數": 18,
                    "鄰域深度": 2,
                    "微擾強度": 0.10,
                    "對稱破缺步數占位": 20,
                    "重複樣本數": 30,
                    "動力學步數": 20,
                    "整體觀測簽名": "medium",
                    "偽隨機基底種子": 20,
                }
            ]
        ),
    ),
    num_rows="dynamic",
    key="b68",
)

def _t8_batch_display_parts(res: object, run_idx: int):
    if not isinstance(res, dict):
        return []
    df = pd.DataFrame([flatten_result_row("e8s", res)])
    return [
        (
            "單次執行指標（與本頁單次區塊同形）",
            df,
            f"table_10_8_batch_run_{run_idx}.csv",
        )
    ]


if st.button("批次 §10.8"):

    def _r(row: pd.Series, prog):
        tsb = int(batch_cell(row, "對稱破缺步數占位", "T_sb", 20))
        return run_experiment_10_8(
            n=int(batch_cell(row, "節點數", "n")),
            m=int(batch_cell(row, "超邊數", "m")),
            init_family=str(batch_cell(row, "初態族別", "init_family")),
            eta=float(batch_cell(row, "微擾強度", "eta", 0.10)),
            T_sb=tsb,
            r=int(batch_cell(row, "鄰域深度", "r", 2)),
            sig_obs=str(batch_cell(row, "整體觀測簽名", "sig_obs", "medium")),
            n_samples=int(batch_cell(row, "重複樣本數", "n_samples")),
            dynamics_steps=int(batch_cell(row, "動力學步數", "dynamics_steps", tsb)),
            seed=int(batch_cell(row, "偽隨機基底種子", "seed", 0)),
            progress=prog,
        )

    st.session_state["t8_batch_runs"] = run_batch_per_run_rows(
        df, _r, stop_on_error=False, use_progress=True
    )

if "t8_batch_runs" in st.session_state:
    render_batch_per_run_tables(
        st.session_state["t8_batch_runs"],
        _t8_batch_display_parts,
        key_prefix="t8_batch",
        column_name_map_for_df=lambda d: build_ch10_column_name_map(d.columns),
    )
