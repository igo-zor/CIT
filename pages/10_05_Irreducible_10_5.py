"""
§10.5 解析不可分解與非可分結構。
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from hypergraph_experiment.viz import matplotlib_figure_to_png_bytes

from hypergraph_experiment.ch10_symbol_glossary import (
    L10_ALPHA_CROSS,
    L10_DELTA_ENT,
    L10_K_MAX,
    L10_K_MIN,
    L10_M_EDGES,
    L10_N_CFG_10_5,
    L10_SAMPLE_LIMIT,
    L10_SEED,
    render_section_glossary,
)
from hypergraph_experiment.experiments.exp_10_5_bipartite import (
    run_experiment_10_5,
    section_10_5_batch_sweep_dataframe,
    section_10_5_output_parameters_df,
    section_10_5_sweep_axis_candidates,
)
from hypergraph_experiment.ch10_paper_presets import CH10_5_BASELINE
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

st.set_page_config(page_title="§10.5 不可分解", layout="wide")
render_sidebar_warehouse()

st.title("§10.5　解析不可分解性")
render_section_glossary(st, "10.5")
st.markdown(
    r"""
    對應《約束世界論 30》**§10.5.1**：在固定二分切分 $V=V_A\sqcup V_B$ 下，檢查是否存在穩定之**解析不可分解單元**，
    且其比例是否與**非可分**解析層統計（$D_{\mathrm{sep}}$、$I(A;B)$）同向變化——非量子態模擬，而是 9.4 之有限離散操作化。
    控制 **$\alpha_{\mathrm{cross}}$**、估 **$\rho_{\mathrm{irred}}$**、**$\rho_{\mathrm{cross}}$**、**$D_{\mathrm{sep}}$**、**$I(A;B)$**。
    **建議**：$n_A=n_B\in\{4,6,8\}$，$m\in\{12,16,20,24\}$，$k_{\min}=2$、$k_{\max}=3$，$N_{\mathrm{cand}}\sim 10^3$–$10^4$、$N_{\mathrm{cfg}}$ 依 §10.5.10。
    """
)

with st.form("f55"):
    st.subheader("建議固定參數（可覆寫）")
    a, b, c = st.columns(3)
    with a:
        na = st.number_input(
            "n_A｜$|\mathcal{V}_A|$ — A 側節點數（§10.5.4（二））", 4, 12, int(CH10_5_BASELINE["n_a"])
        )
        nb = st.number_input(
            "n_B｜$|\mathcal{V}_B|$ — B 側節點數（§10.5.4（二））", 4, 12, int(CH10_5_BASELINE["n_b"])
        )
        m_e = st.number_input(L10_M_EDGES, 2, 40, int(CH10_5_BASELINE["m_edges"]))
    with b:
        k_lo = st.number_input(L10_K_MIN, 2, 6, int(CH10_5_BASELINE["k_min"]))
        k_hi = st.number_input(L10_K_MAX, 2, 8, int(CH10_5_BASELINE["k_max"]))
    with c:
        d_e = st.number_input(L10_DELTA_ENT, 0, 20, int(CH10_5_BASELINE["delta_ent"]))
        sd = st.number_input(L10_SEED, 0, 2_000_000_000, int(CH10_5_BASELINE["seed"]), step=1)
    st.divider()
    st.subheader("觀測集與掃描（可覆寫）")
    r1, r2 = st.columns(2)
    with r1:
        sl = st.number_input(L10_SAMPLE_LIMIT, 50, 20_000, int(CH10_5_BASELINE["sample_limit"]))
        use_all_cfg = st.checkbox(
            "以全部合法配置作為解析觀測集（略過 $N_{cfg}$ 上限；論文穩健性以外可選）",
            value=False,
        )
    with r2:
        n_cfg_in = st.number_input(
            L10_N_CFG_10_5,
            1,
            50_000,
            int(CH10_5_BASELINE["n_cfg"]),
            help="若勾選上方「以全部合法配置…」，執行時將傳入 n_cfg=None。",
        )
        ax = st.slider(L10_ALPHA_CROSS, 0.0, 1.0, float(CH10_5_BASELINE["alpha_cross"]))
    sub = st.form_submit_button("執行")

if sub:
    if int(k_hi) < int(k_lo):
        st.error("請確保 k_max ≥ k_min。")
    else:
        bar = st.progress(0)
        pr = streamlit_progress_callback(bar, st.empty())
        n_cfg_val = None if use_all_cfg else int(n_cfg_in)
        out = run_experiment_10_5(
            n_a=int(na),
            n_b=int(nb),
            m_edges=int(m_e),
            max_edge_size=int(k_hi),
            k_min=int(k_lo),
            alpha_cross=float(ax),
            sample_limit=int(sl),
            n_cfg=n_cfg_val,
            seed=int(sd),
            delta_ent=int(d_e),
            progress=pr,
        )
        st.session_state["out55"] = out
        st.session_state["out55_params"] = {
            "n_a": int(na),
            "n_b": int(nb),
            "m_edges": int(m_e),
            "k_min": int(k_lo),
            "k_max": int(k_hi),
            "alpha_cross": float(ax),
            "sample_limit": int(sl),
            "n_cfg": n_cfg_val,
            "seed": int(sd),
            "delta_ent": int(d_e),
        }

out55 = st.session_state.get("out55")
if out55:
    st.subheader("§10.5.5 輸出參數（主表，可複製／下載 CSV）")
    st.caption(
        "列舉論文 §10.5.5（一）–（四）之記號與數值；與下方扁平化寬表為同一跑一次實驗之不同呈現。"
    )
    df_t5_main = section_10_5_output_parameters_df(out55)
    render_table_with_copy_csv(
        df_t5_main,
        key_prefix="t5_single_main",
        csv_filename="table_10_5_single_section_155_main.csv",
        hide_index=True,
    )
    if out55.get("n_cfg_notice"):
        st.info(str(out55["n_cfg_notice"]))
    st.subheader("單次執行指標（扁平欄位，可複製／下載 CSV）")
    df_t5_single = pd.DataFrame([flatten_result_row("e5s", out55)])
    render_table_with_copy_csv(
        df_t5_single,
        key_prefix="t5_single",
        csv_filename="table_10_5_single.csv",
        column_name_map=build_ch10_column_name_map(df_t5_single.columns),
    )
    render_parameters_table(
        st.session_state.get("out55_params"),
        key_prefix="t5_single",
        csv_filename="table_10_5_single_params.csv",
    )
    # 與 §10.3 相同：matplotlib 輔圖 + 可下載 PNG
    _m55 = out55.get("metrics")
    if isinstance(_m55, dict) and not _m55.get("error"):
        st.subheader("輔圖：§10.5.5 主要指標（單次）")
        st.caption("以長條圖摘要本次執行之 ρ_irred、ρ_cross、D_sep、I(A;B)；可下載 PNG。")
        _bar_specs: list[tuple[str, str]] = [
            ("ρ_irred", "rho_irred"),
            ("ρ_cross", "rho_cross_mean"),
            ("D_sep（TV）", "D_sep_total_variation"),
            ("D_sep（JS，bits）", "D_sep_JS_bits"),
            ("I(A;B)（bits）", "I_A_B_bits"),
        ]
        _blabels: list[str] = []
        _bvals: list[float] = []
        for _zh, _mk in _bar_specs:
            _v = _m55.get(_mk)
            if _v is None:
                continue
            try:
                _bvals.append(float(_v))
                _blabels.append(_zh)
            except (TypeError, ValueError):
                continue
        if _blabels:
            _fig0, _ax0 = plt.subplots(figsize=(6, 3.2))
            _ax0.barh(_blabels, _bvals, color="#4C72B0")
            _ax0.set_title("§10.5.5 主要指標（單次執行）")
            _ax0.grid(axis="x", alpha=0.3)
            _ax0.set_xlabel("數值")
            _png0 = matplotlib_figure_to_png_bytes(_fig0)
            st.pyplot(_fig0)
            plt.close(_fig0)
            st.download_button(
                label="下載 PNG：fig_10_5_single_metrics.png",
                data=_png0,
                file_name="fig_10_5_single_metrics.png",
                mime="image/png",
                key="t5_single_metrics_png_dl",
            )
    with st.expander("完整 JSON（除錯）"):
        st.json(out55)

st.subheader("批次（α_cross 等掃描列）")
st.caption(
    "欄位：甲乙區節點數、邊數、最小／最大超邊階數、跨區傾向、候選上限、"
    "合法樣本數（論文 $N_{cfg}$；填 **0** 表示以全部合法配置為觀測集）、熵差閾值、偽隨機基底種子。"
)

# 論文 §10.5.4（一）：離散掃描軸範例
_ALPHA_SCAN_10_5 = (0.00, 0.10, 0.20, 0.30, 0.50, 0.80)


def _t5_default_row(alpha_cross: float) -> dict[str, object]:
    return {
        "甲區節點數": int(CH10_5_BASELINE["n_a"]),
        "乙區節點數": int(CH10_5_BASELINE["n_b"]),
        "超邊數": int(CH10_5_BASELINE["m_edges"]),
        "最小超邊階數": int(CH10_5_BASELINE["k_min"]),
        "最大超邊階數": int(CH10_5_BASELINE["k_max"]),
        "跨區傾向": float(alpha_cross),
        "候選採樣上限": int(CH10_5_BASELINE["sample_limit"]),
        "合法樣本數（論文N_cfg）": int(CH10_5_BASELINE["n_cfg"]),
        "熵差閾值": int(CH10_5_BASELINE["delta_ent"]),
        "偽隨機基底種子": int(CH10_5_BASELINE["seed"]),
    }


if st.button("載入論文 α_cross 掃描列（§10.5.4（一））", key="load_t5_alpha_scan"):
    st.session_state["t5_batch_template_df"] = pd.DataFrame([_t5_default_row(a) for a in _ALPHA_SCAN_10_5])
if st.button("載入論文建議批次模板（§10.5.10 基準列兩筆）", key="load_t5_template"):
    st.session_state["t5_batch_template_df"] = pd.DataFrame(
        [_t5_default_row(0.30), _t5_default_row(0.60)]
    )
df = st.data_editor(
    st.session_state.get(
        "t5_batch_template_df",
        pd.DataFrame([_t5_default_row(float(CH10_5_BASELINE["alpha_cross"]))]),
    ),
    num_rows="dynamic",
    key="b55",
)


def _t5_batch_display_parts(res: object, run_idx: int):
    if not isinstance(res, dict):
        return []
    df_main = section_10_5_output_parameters_df(res)
    df_flat = pd.DataFrame([flatten_result_row("e5s", res)])
    return [
        (
            "§10.5.5 輸出參數（主表）",
            df_main,
            f"table_10_5_batch_run_{run_idx}_section_155_main.csv",
        ),
        (
            "單次執行指標（與本頁單次區塊同形）",
            df_flat,
            f"table_10_5_batch_run_{run_idx}.csv",
        ),
    ]


if st.button("批次 §10.5"):

    def _r(r: pd.Series, prog):
        n_cfg_raw = batch_cell(
            r, "合法樣本數（論文N_cfg）", "n_cfg", int(CH10_5_BASELINE["n_cfg"])
        )
        try:
            n_cfg_int = int(n_cfg_raw)
        except (TypeError, ValueError):
            n_cfg_int = int(CH10_5_BASELINE["n_cfg"])
        n_cfg_opt = None if n_cfg_int <= 0 else n_cfg_int
        return run_experiment_10_5(
            n_a=int(batch_cell(r, "甲區節點數", "n_a")),
            n_b=int(batch_cell(r, "乙區節點數", "n_b")),
            m_edges=int(batch_cell(r, "超邊數", "m_edges")),
            k_min=int(batch_cell(r, "最小超邊階數", "k_min", CH10_5_BASELINE["k_min"])),
            max_edge_size=int(batch_cell(r, "最大超邊階數", "k_max", CH10_5_BASELINE["k_max"])),
            alpha_cross=float(batch_cell(r, "跨區傾向", "alpha_cross")),
            sample_limit=int(batch_cell(r, "候選採樣上限", "sample_limit")),
            n_cfg=n_cfg_opt,
            seed=int(batch_cell(r, "偽隨機基底種子", "seed")),
            delta_ent=int(batch_cell(r, "熵差閾值", "delta_ent", 0)),
            progress=prog,
        )

    st.session_state["t5_batch_runs"] = run_batch_per_run_rows(
        df, _r, stop_on_error=False, use_progress=True
    )

if "t5_batch_runs" in st.session_state:
    render_batch_per_run_tables(
        st.session_state["t5_batch_runs"],
        _t5_batch_display_parts,
        key_prefix="t5_batch",
        column_name_map_for_df=lambda d: build_ch10_column_name_map(d.columns),
    )

    # 論文 §10.5.1／§10.5.5：批次掃描列之觀察重點（與 §10.3 相同以 matplotlib 繪製並可下載 PNG）
    st.subheader("批次掃描曲線（論文觀察重點，PNG）")
    st.caption(
        "橫軸為參數表中實際變動之數值欄（預設選 **跨區傾向** ＝ "
        r"$\alpha_{\mathrm{cross}}$ 掃描）；僅繪製無錯誤之成功列，並依橫軸排序。"
        " 圖表以 matplotlib 產生（同 §10.3 輔圖）；每張圖下方附**同資料之表格**，"
        "含一鍵複製 Markdown 與 CSV 下載，以及對應 PNG。"
    )
    _t5_sweep_df = section_10_5_batch_sweep_dataframe(st.session_state["t5_batch_runs"])
    if _t5_sweep_df.empty:
        st.info("目前無可繪圖之成功批次列（請確認各列實驗未回傳錯誤）。")
    else:
        _t5_x_opts = section_10_5_sweep_axis_candidates(_t5_sweep_df)
        if not _t5_x_opts:
            st.warning(
                "無法偵測掃描軸：請在批次表中讓至少一個**數值參數欄**（如跨區傾向、超邊數）在不同列之間有變化。"
            )
        else:
            _pref = "跨區傾向"
            _ix = _t5_x_opts.index(_pref) if _pref in _t5_x_opts else 0
            _x_col = st.selectbox(
                "拆線圖橫軸（掃描參數）",
                options=_t5_x_opts,
                index=_ix,
                key="t5_batch_sweep_x",
            )
            # （論文小節、縱軸欄位、圖表說明）
            _t5_chart_specs: list[tuple[str, str, str]] = [
                (
                    "10.5.5（一）",
                    "rho_irred",
                    r"$\rho_{\mathrm{irred}}$：解析不可分解單元占比",
                ),
                (
                    "10.5.5（二）",
                    "rho_cross_mean",
                    r"$\rho_{\mathrm{cross}}$：跨區塊耦合密度（代表態平均）",
                ),
                (
                    "10.5.5（三）",
                    "D_sep_TV",
                    r"$D_{\mathrm{sep}}$（TV）：非可分偏離，正文建議主距離",
                ),
                (
                    "10.5.5（三）補充",
                    "D_sep_JS_bits",
                    r"$D_{\mathrm{sep}}$（JS，bits）：程式補充對照",
                ),
                (
                    "10.5.5（四）",
                    "I_A_B_bits",
                    r"$I(A;B)$：區塊互資訊（bits）",
                ),
            ]
            def _t5_axis_slug(name: str) -> str:
                """產生檔名用簡短後綴（避免路徑非法字元）。"""
                s = re.sub(r"[^\w\-]+", "_", str(name), flags=re.ASCII)
                return (s.strip("_")[:24] or "axis").lower()

            # 掃描曲線縱軸欄位 → 表格表頭用中文說明（與 §10.5.5 小節對齊）
            _t5_sweep_y_col_zh: dict[str, str] = {
                "rho_irred": "ρ_irred（解析不可分解單元占比；§10.5.5（一））",
                "rho_cross_mean": "ρ_cross 平均（跨區塊耦合密度；§10.5.5（二））",
                "D_sep_TV": "D_sep（TV；非可分偏離主距離；§10.5.5（三））",
                "D_sep_JS_bits": "D_sep（JS，bits；程式補充；§10.5.5（三））",
                "I_A_B_bits": "I(A;B)（bits；區塊互資訊；§10.5.5（四））",
            }

            _x_slug = _t5_axis_slug(_x_col)
            for _j, (_sec, _y_col, _cap) in enumerate(_t5_chart_specs):
                st.markdown(f"**§{_sec}**　{_cap}")
                _sub = _t5_sweep_df[[_x_col, _y_col]].copy()
                _sub[_x_col] = pd.to_numeric(_sub[_x_col], errors="coerce")
                _sub[_y_col] = pd.to_numeric(_sub[_y_col], errors="coerce")
                _sub = _sub.dropna(subset=[_x_col, _y_col]).sort_values(_x_col)
                if _sub.empty:
                    st.caption("此指標無有效數值可繪製（可能全為空或無法轉為數值）。")
                else:
                    _fig, _ax = plt.subplots(figsize=(6, 3))
                    _ax.plot(
                        _sub[_x_col],
                        _sub[_y_col],
                        marker="o",
                        color="#C44E52",
                    )
                    _ax.set_title(f"§{_sec}　{_cap}")
                    _ax.set_xlabel(str(_x_col))
                    _ax.set_ylabel(_y_col)
                    _ax.grid(alpha=0.3)
                    _png_b = matplotlib_figure_to_png_bytes(_fig)
                    st.pyplot(_fig)
                    plt.close(_fig)
                    _fn = f"fig_10_5_batch_{_y_col}_{_x_slug}.png"
                    st.download_button(
                        label=f"下載 PNG：{_fn}",
                        data=_png_b,
                        file_name=_fn,
                        mime="image/png",
                        key=f"t5_batch_png_dl_{_j}_{_x_slug}",
                    )
                    # 與圖同資料：表格 + 一鍵複製 Markdown（與全站 render_table_with_copy_csv 同形）
                    _df_tbl = _sub[[_x_col, _y_col]].copy()
                    _x_hdr = f"掃描參數（橫軸）｜{_x_col}"
                    _y_hdr = _t5_sweep_y_col_zh.get(
                        _y_col, f"§{_sec}｜{_y_col}"
                    )
                    st.markdown("**本圖對應數據表**（列序與橫軸排序一致）")
                    render_table_with_copy_csv(
                        _df_tbl,
                        key_prefix=f"t5_batch_sweep_{_j}_{_x_slug}",
                        csv_filename=f"table_10_5_batch_sweep_{_y_col}_{_x_slug}.csv",
                        hide_index=True,
                        column_name_map={_x_col: _x_hdr, _y_col: _y_hdr},
                    )
