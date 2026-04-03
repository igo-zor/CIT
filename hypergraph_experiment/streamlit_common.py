"""
Streamlit 共用：進度條包裝、批次執行、圖表與結果扁平化（不依賴單一頁面腳本）。
"""

from __future__ import annotations

import io
import json
import re
import zipfile
from typing import Any, Callable, Dict, List, Literal, Optional, Sequence, Tuple, Union

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from hypergraph_experiment.ch10_symbol_glossary import (
    L10_ALPHA_CROSS,
    L10_CONNECTED,
    L10_D_MAX,
    L10_DELTA,
    L10_DELTA_COARSE,
    L10_DELTA_ENT,
    L10_DELTA_FINE,
    L10_DELTA_T,
    L10_ETA,
    L10_EPS_PLAT_109,
    L10_FORBID_TRI,
    L10_INIT_FAMILY,
    L10_K_MAX,
    L10_K_MIN,
    L10_KERNEL,
    L10_M_EDGES,
    L10_M_MAX,
    L10_M_TRIAL,
    L10_M_CTX,
    L10_MODE_CTX,
    L10_N,
    L10_N0_RUNS,
    L10_N_CFG,
    L10_N_CFG_10_5,
    L10_N_CTX,
    L10_N_HIST,
    L10_N_REP,
    L10_N_SAMPLES,
    L10_N_SEED_BATCHES,
    L10_N_SEED_RUNS_109,
    L10_ETA_CTX,
    L10_N_SEARCH_106,
    L10_N_VAL_106,
    L10_NODES_BIT,
    L10_T_LOC,
    L10_OBS_SIG,
    L10_P_MAX,
    L10_N_SEED_107,
    L10_R_DEPTH,
    L10_SAMPLE_LIMIT,
    L10_SEED,
    L10_SIG_COARSE,
    L10_SIG_FINE,
    L10_SIG_WEAK,
    L10_T_SB,
    L10_T_STEPS,
    L10_W_A,
    L10_W_H,
    L10_W_CTX,
    L10_W_LIST,
    render_full_ch10_sidebar_note,
)
from hypergraph_experiment.output_rounding import round_floats_for_output
from hypergraph_experiment.progress_types import ProgressCallback
from hypergraph_experiment.storage import default_store_root, save_run
from hypergraph_experiment.viz import (
    plot_entropy_time_series,
    plot_hypergraph,
    plot_largest_class_bars,
)


# 新版中文欄名 → 舊版 session 曾使用之同義欄名（向後相容）
_BATCH_ZH_LEGACY: dict[str, tuple[str, ...]] = {
    "偽隨機基底種子": ("隨機種子",),
    "種子批次數（論文N_seed）": ("種子批次數",),
}


def batch_cell(row: pd.Series, zh: str, en: str, default: Optional[Any] = None) -> Any:
    """
    讀取批次列：優先使用中文欄名，否則回退英文欄名（舊版 session 範本相容）。

    Args:
        row: 單列參數。
        zh: 中文欄名（無符號、與論文語義對齊）。
        en: 程式內部英文鍵名。
        default: 兩者皆缺時之預設值。
    """
    order = (zh,) + _BATCH_ZH_LEGACY.get(zh, ())
    for k in order:
        if k in row.index:
            return row[k]
    if en in row.index:
        return row[en]
    return default


# 單次／批次扁平化欄名前綴（與 flatten_result_row、run_batch_dataframe 一致）
_CH10_FLAT_PREFIXES: tuple[str, ...] = (
    "r3s_",
    "r4s_",
    "e5s_",
    "e6s_",
    "d7s_",
    "e8s_",
    "e9s_",
    "p7_",
    "dy_",
    "r4_",
    "r_",
    "e5_",
    "e6_",
    "e8_",
    "e9_",
    "t2_",
)


def _strip_ch10_flat_column_prefix(col: str) -> str:
    """移除扁平化單次／批次列前綴，得到與論文表結構對齊之鍵幹。"""
    for prefix in _CH10_FLAT_PREFIXES:
        if col.startswith(prefix):
            return col[len(prefix) :]
    return col


# §10.9 per_window 寬表：程式內部鍵 → 介面／CSV 用中文欄名（與論文 §10.9.5 語義對齊）
CH10_9_PER_WINDOW_COL_ZH: dict[str, str] = {
    "w": "時間聚合視窗寬度 w（步）",
    "H_macro_bits": "宏觀型別熵（bits；宏觀標籤分布）",
    "H_macro_bits_std": "宏觀型別熵·標準差（多條歷史彙總時）",
    "Var_C_eff_edges": "Var(C_eff)（代理：末態超邊數 |E| 之序列變異）",
    "Var_C_eff_edges_std": "Var(C_eff)·標準差（多條歷史彙總時）",
    "n_window_positions": "可放置視窗位置數（沿時間滑動之列數）",
    "JS_vs_w1_bits": "JS_w（相對最細 w=1；bits）",
    "JS_vs_w1_bits_std": "JS_w·標準差（多條歷史彙總時）",
    "R_edge_bar": "微觀邊周轉率整段平均",
    "R_edge_bar_std": "邊周轉率平均·標準差（多條歷史彙總時）",
    "tau_unit_max": "宏觀標籤最長連續段長度 τ（最大）",
    "tau_unit_max_std": "τ_max·標準差（多條歷史彙總時）",
    "tau_unit_mean": "宏觀標籤平均連續段長度 τ（平均）",
    "tau_unit_mean_std": "τ_mean·標準差（多條歷史彙總時）",
    "L_plat_max_len_Ceff": "C_eff 序列最長平台長度（ε_plat 判定）",
    "L_plat_max_len_Ceff_std": "最長平台長度·標準差",
    "L_plat_num_seg_Ceff": "C_eff 序列平台段數",
    "L_plat_num_seg_Ceff_std": "平台段數·標準差",
}


def per_window_metrics_dataframe_zh(df: pd.DataFrame) -> pd.DataFrame:
    """
    將 §10.9 ``per_window`` 寬表之英文鍵換成中文欄名（顯示／匯出 CSV 皆可）。

    未列於對照表之欄位保留原名，避免未預期欄位被靜默丟失。
    """
    if df.empty:
        return df
    rename = {c: CH10_9_PER_WINDOW_COL_ZH.get(str(c), str(c)) for c in df.columns}
    return df.rename(columns=rename)


def build_ch10_column_name_map(columns: list[str] | pd.Index) -> dict[str, str]:
    """建立第十章表格欄名映射（中文、無符號）。"""
    token_map: dict[str, str] = {
        "n": "節點數",
        "k": "超邊大小",
        "kmin": "最小超邊大小",
        "kmax": "最大超邊大小",
        "m": "邊數",
        "mmax": "最大邊數",
        "sample": "樣本",
        "limit": "上限",
        "seed": "偽隨機基底種子",
        "level": "層",
        "count": "數量",
        "label": "標籤",
        "cfg": "合法配置數",
        "cand": "候選",
        "candidate": "候選",
        "ratio": "比例",
        "shrink": "收縮率",
        "prev": "前層",
        "forbid": "排除",
        "subset": "子集",
        "chain": "鏈式",
        "ok": "成立",
        "connected": "連通",
        "triangles": "二元三角",
        "tri": "二元三角",
        "degree": "度數",
        "max": "最大",
        "min": "最小",
        "analysis": "分析",
        "metrics": "指標",
        "entropy": "熵",
        "bits": "位元",
        "overlap": "重疊",
        "transitivity": "傳遞",
        "violation": "違反",
        "compression": "壓縮",
        "classes": "類別",
        "equivalence": "等價",
        "path": "路徑",
        "terminal": "終端",
        "window": "視窗",
        "macro": "宏觀",
        "runs": "軌道數",
        "steps": "步數",
        "delta": "解析閾值",
        "signature": "解析簽名",
        "mode": "模式",
        "experiment": "實驗編號",
        "run": "執行",
        "index": "序號",
        "error": "錯誤訊息",
        "js": "分布差異",
        "mean": "平均",
        "len": "長度",
        "first": "首值",
        "last": "末值",
        "plat": "平台",
        "plateau": "平台",
        "fraction": "占比",
        "turnover": "周轉",
        "isolated": "孤立",
        "compat": "相容",
        "graph": "圖",
        "isol": "孤立",
        "legal": "合法",
        "update": "更新",
        "step": "步",
        "epsilon": "閾值",
        "tau": "持續",
        "reach": "可達",
        "pair": "成對",
        "smin": "鄰域下限",
        "kernel": "核",
        "nnz": "非零占比",
        "stability": "穩定度",
        "fibers": "纖維數",
        "edge": "邊",
        "bar": "平均欄",
        "unit": "單元",
    }
    exact_map: dict[str, str] = {
        "run_index": "執行序號",
        "error": "錯誤訊息",
        "n_candidates": "候選母集數",
        "候選母集數": "候選母集數",
        "num_obs_configs": (
            "觀測配置數（解析統計實際觀測筆數；§10.3／§10.5；"
            "若參數設定 N_cfg 則為抽樣後筆數）"
        ),
        "n_cfg_requested": "輸入配置數請求值",
        "n_cfg_notice": "觀測集抽樣提示",
        "n_cfg": "輸入配置數",
        "parameters_n_cfg": "參數 輸入配置數",
        "n_rep": "重複次數",
        "parameters_n_rep": "參數 重複次數",
        "analysis_rep_summary": "重抽摘要",
        "num_admissible_filtered": "域型過濾後可採用配置數（§10.5）",
        "n_val": L10_N_VAL_106,
    }
    level_tail_map = {
        "label": "約束層級",
        "cfg": "合法配置數",
        "cand_ratio": "佔候選比例",
        "cand_shrink_ratio": "相對候選收縮率",
        "prev_ratio": "相對前層保留率",
        "prev_shrink_ratio": "相對前層收縮率",
        "N_forbid_i": "累計排除數",
        "N_forbid_delta": "本層新增排除數",
        "subset_prev": "是否為前層子集",
        "chain_subset_ok": "鏈式子集成立",
        "rho_i_to_prev": "相對前層保留率",
        "viol_all_too_many_edges": "排除主因超邊數超限筆數",
        "viol_all_edge_size_bad": "排除主因超邊大小不符筆數",
        "viol_all_degree_excess": "排除主因度數超限筆數",
        "viol_all_disconnected": "排除主因連通未滿足筆數",
        "viol_all_pair_triangle": "排除主因禁二元三角筆數",
        "viol_new_too_many_edges": "新增排除主因超邊數超限筆數",
        "viol_new_edge_size_bad": "新增排除主因超邊大小不符筆數",
        "viol_new_degree_excess": "新增排除主因度數超限筆數",
        "viol_new_disconnected": "新增排除主因連通未滿足筆數",
        "viol_new_pair_triangle": "新增排除主因禁二元三角筆數",
        "viol_concentration": "違規主因集中度",
        "viol_new_concentration": "新增排除違規主因集中度",
    }
    win_tail_zh = CH10_9_PER_WINDOW_COL_ZH

    out: dict[str, str] = {}
    for col in [str(c) for c in columns]:
        if col in exact_map:
            out[col] = exact_map[col]
            continue
        stripped = _strip_ch10_flat_column_prefix(col)
        if stripped in exact_map:
            out[col] = exact_map[stripped]
            continue
        m_level = re.match(r"^level_(\d+)_(.+)$", stripped)
        if m_level:
            idx = int(m_level.group(1))
            tail = m_level.group(2)
            out[col] = f"第{idx}層{level_tail_map.get(tail, tail.replace('_', ' '))}"
            continue
        m_win = re.match(r"^per_window_(\d+)_(.+)$", stripped)
        if m_win:
            idx = int(m_win.group(1))
            tail_raw = m_win.group(2)
            tail_zh = win_tail_zh.get(tail_raw, tail_raw.replace("_", " "))
            # idx 為扁平化時之列序（第幾個 w），與該列「w」欄數值一致時可互相對照
            out[col] = f"第{idx + 1}列（per_window）·{tail_zh}"
            continue
        parts = [p for p in stripped.split("_") if p]
        zh_parts: list[str] = []
        for p in parts:
            key = p.lower()
            zh_parts.append(token_map.get(key, p))
        out[col] = " ".join(zh_parts)
    return out


def streamlit_progress_callback(bar: Any, status: Any | None = None) -> ProgressCallback:
    """
    建立適用於 ``run_full_experiment`` 等函式之回呼，更新 Streamlit 進度條。

    Args:
        bar: ``st.progress`` 回傳之物件。
        status: 可選 ``st.status`` 或 ``st.empty()`` 容器供文字訊息。

    Returns:
        ``(done,total,message) -> None`` 之回呼；total 為 0 時略過更新。
    """

    def _cb(done: int, total: int, message: str) -> None:
        if total <= 0:
            return
        bar.progress(min(1.0, float(done) / float(total)))
        if status is not None:
            status.write(message)

    return _cb


def render_sidebar_warehouse() -> None:
    """顯示實驗倉儲根目錄與論文符號提示。"""
    store = default_store_root()
    st.sidebar.markdown("**倉儲目錄**")
    st.sidebar.code(str(store), language="text")
    render_full_ch10_sidebar_note(st)


def render_result_charts(result: dict[str, Any], key_prefix: str = "") -> None:
    """依模式顯示靜態柱狀圖或動力學熵曲線。"""
    mode = (result.get("parameters") or {}).get("mode", "static")
    analysis = result.get("analysis") or {}
    if not isinstance(analysis, dict) or "error" in analysis:
        st.warning(analysis.get("error", "無分析結果可繪圖。"))
        return

    if mode == "dynamics":
        w_h = analysis.get("w_h")
        series_w = analysis.get("entropy_time_series_wH") or []
        series = analysis.get("entropy_time_series") or []
        pe = (result.get("parameters") or {}).get("epsilon_plat")
        plat_e = float(pe) if pe is not None else None
        if isinstance(w_h, int) and w_h > 1 and series_w:
            st.image(
                io.BytesIO(plot_entropy_time_series(series_w, plateau_epsilon=plat_e)),
                caption=(
                    "解析熵時間序列（窗口熵 $H_{\\Lambda,w_H}^{(\\ell)}$；"
                    f"w_H={w_h}；§10.7.2（六）、§10.7.5）"
                ),
            )
        elif series:
            st.image(
                io.BytesIO(plot_entropy_time_series(series, plateau_epsilon=plat_e)),
                caption="解析熵時間序列（對應 $H_\\Lambda^{(\\ell)}$ 操作化；§10.7.5）",
            )
        else:
            st.info("無熵時間序列。")
    else:
        st.image(
            io.BytesIO(plot_largest_class_bars(analysis)),
            caption="前 10 大解析單元（等價類）大小（§10.3、$S_{\\Lambda,\\delta}$）",
        )


def render_hypergraph_preview(sample_configs: list[dict[str, Any]], key_prefix: str) -> None:
    """自範例配置選一筆顯示超圖。"""
    if not sample_configs:
        st.info("無範例超圖可顯示。")
        return
    style = st.selectbox(
        "超圖繪製樣式（配置 $c=(V,E)$ 之視覺化；§10.1.1）",
        ["spatial", "incidence"],
        index=0,
        format_func=lambda s: "空間嵌入（2-section + 凸包）"
        if s == "spatial"
        else "二部圖（頂點–超邊，incidence）",
        key=f"{key_prefix}_hyper_style",
    )
    labels = [
        f"範例 {i + 1}（|E|={len(c.get('hyperedges') or [])}）"
        for i, c in enumerate(sample_configs)
    ]
    idx = st.selectbox(
        "選擇要視覺化之範例超圖",
        range(len(labels)),
        format_func=lambda i: labels[i],
        key=f"{key_prefix}_hyper_idx",
    )
    png = plot_hypergraph(sample_configs[idx], style=style)
    st.image(io.BytesIO(png), caption="超圖預覽")


# by_family：§10.8 三臂結果內含三份完整 per_sample，扁平化會過大且重複。
_SKIP_FLATTEN_KEYS = frozenset({"per_sample", "sample_configs", "by_family"})


def _safe_flat_key(s: str) -> str:
    """扁平化欄位名：避免空白與過長符號造成 CSV 困擾。"""
    return str(s).replace(" ", "_").replace("|", "").replace(".", "_")


def _deep_flatten_into_row(
    key: str,
    v: Any,
    row: dict[str, Any],
    *,
    max_list_items: int = 48,
    max_json_len: int = 12_000,
) -> None:
    """遞迴將任意巢狀結構壓入單列字典（純量、dict、list 可辨識型別）。"""
    if v is None or isinstance(v, (bool, int, float)):
        row[key] = v
        return
    if isinstance(v, str):
        row[key] = v if len(v) <= max_json_len else v[: max_json_len - 1] + "…"
        return
    if isinstance(v, dict):
        if not v:
            row[key] = ""
            return
        for sk, sv in v.items():
            _deep_flatten_into_row(
                f"{key}_{_safe_flat_key(sk)}",
                sv,
                row,
                max_list_items=max_list_items,
                max_json_len=max_json_len,
            )
        return
    if isinstance(v, (list, tuple)):
        lst = list(v)
        if not lst:
            row[key] = ""
            return
        if all(isinstance(x, dict) for x in lst):
            for i, d in enumerate(lst[:max_list_items]):
                if isinstance(d, dict):
                    for dk, dv in d.items():
                        _deep_flatten_into_row(
                            f"{key}_{i}_{_safe_flat_key(dk)}",
                            dv,
                            row,
                            max_list_items=max_list_items,
                            max_json_len=max_json_len,
                        )
            if len(lst) > max_list_items:
                row[f"{key}_truncated_remaining"] = len(lst) - max_list_items
            return
        if all(isinstance(x, (int, float)) for x in lst):
            row[f"{key}_len"] = len(lst)
            row[f"{key}_min"] = min(lst)
            row[f"{key}_max"] = max(lst)
            row[f"{key}_mean"] = sum(lst) / len(lst)
            row[f"{key}_first"] = lst[0]
            row[f"{key}_last"] = lst[-1]
            return
        try:
            s = json.dumps(lst, ensure_ascii=False)
        except TypeError:
            s = str(lst)
        row[key] = s if len(s) <= max_json_len else s[: max_json_len - 1] + "…"
        return
    row[key] = str(v)[:max_json_len]


def flatten_result_row(prefix: str, data: dict[str, Any]) -> dict[str, Any]:
    """將單次實驗結果深度壓平為一列（鍵加前綴；略過超大欄位）。"""
    row: dict[str, Any] = {}
    for k, v in data.items():
        if k in _SKIP_FLATTEN_KEYS:
            continue
        _deep_flatten_into_row(f"{prefix}_{_safe_flat_key(k)}", v, row)
    return round_floats_for_output(row)


def run_batch_dataframe(
    df_params: pd.DataFrame,
    row_runner: Callable[[pd.Series, ProgressCallback | None], dict[str, Any]],
    *,
    flatten_prefix: str,
    stop_on_error: bool,
    use_progress: bool,
) -> pd.DataFrame:
    """
    對參數表逐列執行 ``row_runner``，累積扁平化列為 DataFrame。

    Args:
        df_params: 每列一組參數。
        row_runner: ``(row, progress) -> result_dict``。
        flatten_prefix: 扁平化欄位前綴。
        stop_on_error: 若 True 則遇例外即中斷。
        use_progress: 是否建立進度條。

    Returns:
        合併後之結果表（含 ``run_index``、``error`` 欄位若適用）。
    """
    rows_out: List[dict[str, Any]] = []
    prog_holder = st.progress(0) if use_progress else None
    stat_holder = st.empty() if use_progress else None
    prog_cb = (
        streamlit_progress_callback(prog_holder, stat_holder) if use_progress else None
    )
    n = len(df_params)
    for i in range(n):
        if prog_cb and n > 0:
            prog_cb(i, max(1, n), f"批次 {i + 1}/{n}")
        ser = df_params.iloc[i]
        try:
            res = row_runner(ser, prog_cb)
            flat = flatten_result_row(flatten_prefix, res)
            flat["run_index"] = i
            rows_out.append(flat)
        except Exception as e:
            rows_out.append({"run_index": i, "error": str(e)})
            if stop_on_error:
                break
        if prog_cb and n > 0:
            prog_cb(i + 1, n, f"完成 {i + 1}/{n}")
    return pd.DataFrame(rows_out)


# 批次逐組顯示：(標題, DataFrame, CSV 檔名)；可選第四個元素為 hide_index（與單次表 10-2 一致）
BatchDisplayPart = Union[
    Tuple[str, pd.DataFrame, str],
    Tuple[str, pd.DataFrame, str, bool],
]


def _batch_display_hide_index(part: BatchDisplayPart) -> bool:
    return bool(part[3]) if len(part) > 3 else False


def _param_row_caption(ser: pd.Series, run_index: int, max_fields: int = 10) -> str:
    """產生參數列簡要說明（對照 data_editor 第幾列）。"""
    parts: list[str] = [f"參數表列索引={run_index}"]
    n = 0
    for k, v in ser.items():
        if n >= max_fields:
            parts.append("…")
            break
        parts.append(f"{k}={v}")
        n += 1
    return "；".join(parts)


def run_batch_per_run_rows(
    df_params: pd.DataFrame,
    row_runner: Callable[[pd.Series, ProgressCallback | None], Any],
    *,
    stop_on_error: bool,
    use_progress: bool,
) -> List[Dict[str, Any]]:
    """
    對參數表逐列執行實驗，回傳可供 ``render_batch_per_run_tables`` 使用之結構（不扁平化）。

    每筆元素含 ``run_index``、``param_caption``、``param_row``，
    以及成功時 ``result`` 或失敗時 ``error``。
    """
    prog_holder = st.progress(0) if use_progress else None
    stat_holder = st.empty() if use_progress else None
    prog_cb = (
        streamlit_progress_callback(prog_holder, stat_holder) if use_progress else None
    )
    n = len(df_params)
    out: List[Dict[str, Any]] = []
    for i in range(n):
        if prog_cb and n > 0:
            prog_cb(i, max(1, n), f"批次 {i + 1}/{n}")
        ser = df_params.iloc[i]
        cap = _param_row_caption(ser, i)
        param_row = {str(k): ser[k] for k in ser.index}
        try:
            res = row_runner(ser, prog_cb)
            out.append(
                {
                    "run_index": i,
                    "error": None,
                    "result": res,
                    "param_caption": cap,
                    "param_row": param_row,
                }
            )
        except Exception as e:
            out.append(
                {
                    "run_index": i,
                    "error": str(e),
                    "result": None,
                    "param_caption": cap,
                    "param_row": param_row,
                }
            )
            if stop_on_error:
                if prog_cb and n > 0:
                    prog_cb(i + 1, n, f"完成 {i + 1}/{n}")
                return out
        if prog_cb and n > 0:
            prog_cb(i + 1, n, f"完成 {i + 1}/{n}")
    return out


def render_batch_per_run_tables(
    batch_results: Sequence[Dict[str, Any]],
    display_parts: Callable[[Any, int], Sequence[BatchDisplayPart]],
    *,
    key_prefix: str,
    column_name_map_for_df: Optional[Callable[[pd.DataFrame], Dict[str, str]]] = None,
    group_title_template: str = "批次第 {n} 組",
) -> None:
    """
    依 ``run_batch_per_run_rows`` 之輸出，為每組參數渲染與單次實驗同結構之一張或多張表。

    Args:
        batch_results: 內含 ``run_index``、``error``、``result``、``param_caption``、``param_row``。
        display_parts: ``(result, run_index) ->`` 各表之標題、``DataFrame``、下載檔名（及可選 ``hide_index``）。
        key_prefix: Streamlit 元件鍵前綴（會加上組別與表序避免衝突）。
        column_name_map_for_df: 若提供，則依各表 ``DataFrame`` 分別建立中文欄名映射。
        group_title_template: 每組標題模板，``{n}`` 為 1-based 序號。
    """
    for item in batch_results:
        i = int(item["run_index"])
        st.subheader(group_title_template.format(n=i + 1))
        st.caption(str(item.get("param_caption", "")))
        render_parameters_table(
            item.get("param_row"),
            key_prefix=f"{key_prefix}_params_{i}",
            csv_filename=f"{key_prefix}_params_run_{i}.csv",
            title="本組輸入實驗參數表（固定參數＋變數）",
        )
        if item.get("error"):
            st.error(str(item["error"]))
            continue
        res = item.get("result")
        if res is None:
            st.warning("無結果可顯示。")
            continue
        parts = display_parts(res, i)
        for j, part in enumerate(parts):
            title, df, csv_fn = part[0], part[1], part[2]
            hide_ix = _batch_display_hide_index(part)
            if df is None or getattr(df, "empty", True):
                continue
            st.markdown(f"**{title}**")
            cmap = column_name_map_for_df(df) if column_name_map_for_df else None
            render_table_with_copy_csv(
                df,
                key_prefix=f"{key_prefix}_{i}_{j}",
                csv_filename=csv_fn,
                hide_index=hide_ix,
                column_name_map=cmap,
            )


def zip_download_hypergraph_run(result: dict[str, Any], file_name: str) -> bytes:
    """打包 result 為 ZIP（JSON + 主要 PNG）。"""
    an = result.get("analysis") or {}
    samples = result.get("sample_configs") or []
    zbuf = io.BytesIO()
    tmp_manifest = {
        "parameters": result.get("parameters"),
        **{k: v for k, v in result.items() if k != "parameters"},
    }
    with zipfile.ZipFile(zbuf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("result.json", json.dumps(tmp_manifest, ensure_ascii=False, indent=2))
        if isinstance(an, dict) and "error" not in an:
            if (result.get("parameters") or {}).get("mode") == "dynamics" and an.get(
                "entropy_time_series"
            ):
                _pe = (result.get("parameters") or {}).get("epsilon_plat")
                _plat = float(_pe) if _pe is not None else None
                zf.writestr(
                    "figures/entropy_time_series.png",
                    plot_entropy_time_series(
                        an["entropy_time_series"], plateau_epsilon=_plat
                    ),
                )
            else:
                zf.writestr("figures/largest_classes.png", plot_largest_class_bars(an))
        if samples:
            zf.writestr(
                "figures/sample_hypergraph.png",
                plot_hypergraph(samples[0], style="spatial"),
            )
    zbuf.seek(0)
    return zbuf.getvalue()


def save_run_with_figures(result: dict[str, Any], notes: str, hg_style: str) -> str:
    """寫入實驗庫並附帶圖檔。"""
    an = result.get("analysis") or {}
    samples = result.get("sample_configs") or []
    extra: dict[str, bytes] = {}
    if isinstance(an, dict) and "error" not in an:
        if (result.get("parameters") or {}).get("mode") == "dynamics" and an.get(
            "entropy_time_series"
        ):
            _pe = (result.get("parameters") or {}).get("epsilon_plat")
            _plat = float(_pe) if _pe is not None else None
            extra["figures/entropy_time_series.png"] = plot_entropy_time_series(
                an["entropy_time_series"], plateau_epsilon=_plat
            )
        else:
            extra["figures/largest_classes.png"] = plot_largest_class_bars(an)
    if samples:
        extra["figures/sample_hypergraph.png"] = plot_hypergraph(samples[0], style=hg_style)
    return save_run(result, notes=notes, extra_files=extra)


def _df_to_markdown_table(df: pd.DataFrame) -> str:
    """將 DataFrame 轉成可複製的 Markdown 表格字串。"""

    def _safe_header(col: str) -> str:
        if col == "|Cfg|":
            return r"$\lvert Cfg \rvert$"
        return col.replace("|", r"\|")

    cols = [str(c) for c in df.columns]
    head = "| " + " | ".join(_safe_header(c) for c in cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows: list[str] = []
    for _, row in df.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            txt = "" if pd.isna(v) else str(v)
            vals.append(txt.replace("|", r"\|").replace("\n", " "))
        rows.append("| " + " | ".join(vals) + " |")
    return "\n".join([head, sep, *rows])


def _render_copy_markdown_button(md_text: str, key: str) -> None:
    """渲染一鍵複製 Markdown 的前端按鈕。"""
    js_text = md_text.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")
    components.html(
        f"""
<button id="{key}" style="padding:0.35rem 0.8rem;border:1px solid #666;border-radius:0.4rem;cursor:pointer;">
  一鍵複製表格 Markdown
</button>
<script>
  const btn = document.getElementById("{key}");
  if (btn) {{
    btn.onclick = async () => {{
      try {{
        await navigator.clipboard.writeText(`{js_text}`);
        btn.innerText = "已複製 Markdown";
      }} catch (e) {{
        btn.innerText = "複製失敗（請手動）";
      }}
      setTimeout(() => btn.innerText = "一鍵複製表格 Markdown", 1500);
    }};
  }}
</script>
""",
        height=42,
    )


def _split_symbol_and_semantics(label: str) -> tuple[str, str]:
    """由 glossary 標籤拆出「論文記號」與「語義」。"""
    # 註：glossary 標籤格式多為「記號 — 語義」，若無分隔符則整段視為語義。
    if "—" in label:
        left, right = label.split("—", 1)
        return left.strip(), right.strip()
    return "", label.strip()


_PARAM_KEY_TO_GLOSSARY_LABEL: dict[str, str] = {
    # 通用與 §10.2/10.3
    "n": L10_N,
    "k_min": L10_K_MIN,
    "min_edge_size": L10_K_MIN,
    "k_max": L10_K_MAX,
    "max_edge_size": L10_K_MAX,
    "m_max": L10_M_MAX,
    "max_edges": L10_M_MAX,
    "m": L10_M_EDGES,
    "sample_limit": L10_SAMPLE_LIMIT,
    "n_cfg": L10_N_CFG,
    "num_seeds": L10_N_SEED_BATCHES,
    "n_seed": L10_N_SEED_BATCHES,
    "n_rep": L10_N_REP,
    "seed": L10_SEED,
    "d_max": L10_D_MAX,
    "max_degree": L10_D_MAX,
    "connected": L10_CONNECTED,
    "forbid_pair_triangles": L10_FORBID_TRI,
    "delta": L10_DELTA,
    "signature": L10_OBS_SIG,
    "s_min": "s_min — 重疊率統計最小支持門檻（§10.3）",
    # §10.4
    "refine_coarse_signature": r"$\mathrm{Sig}_{\Lambda}$ — 粗層解析簽名（§10.4）",
    "refine_coarse_delta": r"$\delta$ — 粗層解析閾值（§10.4）",
    "refine_fine_signature": r"$\mathrm{Sig}_{\Lambda'}$ — 細層解析簽名（§10.4）",
    "refine_fine_delta": r"$\delta'$ — 細層解析閾值（§10.4）",
    "refine_kernel": L10_KERNEL,
    "refine_coarse_sample_size": r"$N_{cfg}^{\Lambda}$ — 粗層觀測樣本數（§10.4 實作）",
    "refine_fine_sample_size": r"$N_{cfg}^{\Lambda'}$ — 細層觀測樣本數（§10.4 實作）",
    "refine_fiber_sample_size": r"$N_{fiber}$ — 每粗單元纖維樣本上限（§10.4 實作）",
    "epsilon_push_threshold": r"$\varepsilon_{push}^{*}$ — 推前誤差判定閾值（§10.4.5/10.4.7）",
    "js_threshold": r"$\varepsilon_{JS}^{*}$ — 終端 JS 判定閾值（§10.4.5/10.4.7）",
    "refinement_enabled": "refinement_enabled — 是否啟用 §10.4 細化／纖維分析（程式開關）",
    "refine_compare_chains": "refine_compare_chains — 是否計算雙路徑 A→B／B→A（§10.4.7；程式開關）",
    "kernel_mode": L10_KERNEL,
    "coarse_sample_size": r"$N_{cfg}^{\Lambda}$ — 粗層觀測樣本數（§10.4.7；compare_ordered_refinement_paths）",
    "fine_sample_size": r"$N_{cfg}^{\Lambda'}$ — 細層觀測樣本數（§10.4.7；compare_ordered_refinement_paths）",
    "max_fiber_size": r"$N_{fiber}$ — 纖維取樣上限（§10.4.7；compare_ordered_refinement_paths）",
    "step_a": r"$R_A$ — 子步驟 A 之映射規格（§9.6-C、§10.4.7；kind+δ）",
    "step_b": r"$R_B$ — 子步驟 B 之映射規格（§9.6-C、§10.4.7；kind+δ）",
    "signatures": r"$\{\sigma_{res}\}$ — 掃描用解析簽名版本集合（weak／medium／strong；§10.3）",
    "delta_values": r"$\{\delta\}$ — 掃描用相容閾值集合（§10.3）",
    "levels": r"$\{\Lambda_{\mathrm{dom}}^{(i)}\}$ — 域型約束梯子各層設定（§10.2）",
    "show_sample_configs": "show_sample_configs — 範例配置輸出筆數（程式；非論文參數）",
    # §10.5（《約束世界論 30》小節編號）
    "n_a": r"$n_A$ — 二分圖甲區 $|V_A|$（§10.5.4（二））",
    "n_b": r"$n_B$ — 二分圖乙區 $|V_B|$（§10.5.4（二））",
    "m_edges": L10_M_EDGES,
    "alpha_cross": L10_ALPHA_CROSS,
    "delta_ent": L10_DELTA_ENT,
    # §10.6
    "n_ctx": L10_N_CTX,
    "n_nodes": L10_NODES_BIT,
    "M": L10_M_CTX,
    "w_ctx": L10_W_CTX,
    "eta_ctx": L10_ETA_CTX,
    "T_loc": L10_T_LOC,
    "n_search": L10_N_SEARCH_106,
    "n_val": L10_N_VAL_106,
    "mode": (
        "mode — 執行模式：§10.6 為 obstruction／satisfiable；"
        "run_full_experiment 為 static（靜態解析）／dynamics（動力學；§10.3／§10.7）"
    ),
    # §10.7
    "runs": L10_N0_RUNS,
    "steps": L10_T_STEPS,
    "epsilon_plat": r"$\varepsilon_{\mathrm{plat}}$ — 熵平台判定閾值（§10.7.3、§10.7.6）",
    "m_trial": L10_M_TRIAL,
    "w_h": L10_W_H,
    "w_a": L10_W_A,
    "p_max": L10_P_MAX,
    "n_seed_107": L10_N_SEED_107,
    "sig_path_a": L10_SIG_COARSE,
    "delta_a": L10_DELTA_COARSE,
    "sig_path_b": L10_SIG_FINE,
    "delta_b": L10_DELTA_FINE,
    # §10.8
    "init_family": L10_INIT_FAMILY,
    "eta": L10_ETA,
    "T_sb": L10_T_SB,
    "r": L10_R_DEPTH,
    "n_samples": L10_N_SAMPLES,
    "dynamics_steps": L10_T_STEPS,
    "sig_obs": L10_OBS_SIG,
    # §10.9
    "window_sizes": L10_W_LIST,
    "delta_t": L10_DELTA_T,
    "n_hist": L10_N_HIST,
    "n_seed_runs": L10_N_SEED_RUNS_109,
    "epsilon_plat_109": L10_EPS_PLAT_109,
}

# 批次表 data_editor 常見繁中欄名（與程式內英文鍵同義）
_ZH_PARAM_KEY_TO_GLOSSARY_LABEL: dict[str, str] = {
    "節點數": L10_N,
    "最小超邊階數": L10_K_MIN,
    "最大超邊階數": L10_K_MAX,
    "最大超邊數": L10_M_MAX,
    "候選採樣上限": L10_SAMPLE_LIMIT,
    "偽隨機基底種子": L10_SEED,
    "種子批次數（論文N_seed）": L10_N_SEED_BATCHES,
    "約束層數": "約束層數 — 域型約束梯子層數（§10.2）",
    "頂點度上限": L10_D_MAX,
    "二部圖連通": L10_CONNECTED,
    "禁止二元三角": L10_FORBID_TRI,
    "輸入配置數": L10_N_CFG,
    "重複次數": L10_N_REP,
    "解析簽名": L10_SIG_WEAK,
    "解析閾值整數": L10_DELTA,
    "重疊率鄰域最小支持": "s_min — 重疊率統計最小支持門檻（§10.3）",
    "甲區節點數": r"$n_A$ — 二分圖甲區節點數（§10.5.4（二））",
    "乙區節點數": r"$n_B$ — 二分圖乙區節點數（§10.5.4（二））",
    "合法樣本數（論文N_cfg）": L10_N_CFG_10_5,
    "超邊數": L10_M_EDGES,
    "跨區傾向": L10_ALPHA_CROSS,
    "熵差閾值": L10_DELTA_ENT,
    "上下文樣本數": L10_N_CTX,
    "位元節點數": L10_NODES_BIT,
    "執行模式": L10_MODE_CTX,
    "上下文數 M": L10_M_CTX,
    "視窗大小 w": L10_W_CTX,
    "視窗交疊 η": L10_ETA_CTX,
    "局部型別數 T_loc": L10_T_LOC,
    "全域搜尋上限": L10_N_SEARCH_106,
    "軌道數": L10_N0_RUNS,
    "演化步數": L10_T_STEPS,
    "熵平台判定閾值": "ε_plat — 熵平台判定閾值（§10.7）",
    "初態族別": L10_INIT_FAMILY,
    "鄰域深度": L10_R_DEPTH,
    "微擾強度": L10_ETA,
    "對稱破缺步數占位": L10_T_SB,
    "重複樣本數": L10_N_SAMPLES,
    "動力學步數": L10_T_STEPS,
    "整體觀測簽名": L10_OBS_SIG,
    "視窗寬度清單": L10_W_LIST,
    "聚合步長 Δt": L10_DELTA_T,
    "宏觀平台閾 ε_plat": L10_EPS_PLAT_109,
    "微觀歷史條數": L10_N_HIST,
    "穩健性重跑次數": L10_N_SEED_RUNS_109,
    "每步候選更新數": L10_M_TRIAL,
    "細化鏈 preset": "細化鏈 preset — 論文 10.4.6 固定粗細鏈之介面標籤（§10.4）",
    "細化核模式": L10_KERNEL,
    "粗層樣本數": r"$N_{cfg}^{\Lambda}$ — 粗層觀測樣本數（§10.4 實作）",
    "細層樣本數": r"$N_{cfg}^{\Lambda'}$ — 細層觀測樣本數（§10.4 實作）",
    "纖維樣本數上限": r"$N_{fiber}$ — 每粗單元纖維樣本上限（§10.4 實作）",
    "推前誤差閾值": r"$\varepsilon_{push}^{*}$ — 推前誤差判定閾值（§10.4.5/10.4.7）",
    "JS 差異閾值": r"$\varepsilon_{JS}^{*}$ — 終端 JS 判定閾值（§10.4.5/10.4.7）",
}

_LEVEL_ZH_RE = re.compile(r"^層(\d+)(連通|度上限|禁二元三角)$")


def _resolve_param_glossary_label(key: str) -> str:
    """依參數鍵（英文或批次中文）取得 glossary 對照字串；無則回傳空字串。"""
    k = str(key)
    if k in _PARAM_KEY_TO_GLOSSARY_LABEL:
        return _PARAM_KEY_TO_GLOSSARY_LABEL[k]
    if k in _ZH_PARAM_KEY_TO_GLOSSARY_LABEL:
        return _ZH_PARAM_KEY_TO_GLOSSARY_LABEL[k]
    m = _LEVEL_ZH_RE.match(k)
    if m:
        idx, suffix = m.group(1), m.group(2)
        if suffix == "連通":
            return f"第{idx}層連通 — 域型梯子第 {idx} 層之 2-section 連通條件（§10.2）"
        if suffix == "度上限":
            return f"第{idx}層度上限 — 域型梯子第 {idx} 層之頂點最大度數 $d_{{max}}$（§10.2）"
        if suffix == "禁二元三角":
            return f"第{idx}層禁二元三角 — 域型梯子第 {idx} 層之禁二元三角形 motif（§10.2）"
    return ""


def parameters_to_df(parameters: dict[str, Any]) -> pd.DataFrame:
    """將輸入參數字典展平成四欄資料表（鍵、論文記號、語義、值）。"""
    rows: list[dict[str, str]] = []
    for k, v in parameters.items():
        if isinstance(v, (dict, list, tuple)):
            try:
                v_txt = json.dumps(v, ensure_ascii=False, sort_keys=True)
            except TypeError:
                v_txt = str(v)
        else:
            v_txt = str(v)
        label = _resolve_param_glossary_label(str(k))
        if label:
            symbol, semantics = _split_symbol_and_semantics(label)
        else:
            symbol = "—"
            semantics = (
                f"介面欄位「{k}」：請對照該節頁面摺疊「論文符號對照」或 README_CH10_FIELDS.md；"
                "若為自訂欄位可逕以參數鍵與參數值記錄。"
            )
        rows.append(
            {
                "參數鍵": str(k),
                "論文記號": symbol,
                "論文語義": semantics,
                "參數值": v_txt,
            }
        )
    return pd.DataFrame(rows)


def render_parameters_table(
    parameters: Any,
    *,
    key_prefix: str,
    csv_filename: str,
    title: str = "輸入實驗參數表（固定參數＋變數）",
    column_name_map: dict[str, str] | None = None,
) -> None:
    """渲染輸入參數表，沿用既有一鍵複製 Markdown 與 CSV 下載。"""
    st.markdown(f"**{title}**")
    if not isinstance(parameters, dict):
        st.warning("未提供可用的輸入參數（parameters）。")
        return
    df_params = parameters_to_df(parameters)
    if df_params.empty:
        st.warning("輸入參數為空，無可顯示資料。")
        return
    render_table_with_copy_csv(
        df_params,
        key_prefix=f"{key_prefix}_params",
        csv_filename=csv_filename,
        hide_index=True,
        column_name_map=column_name_map,
    )


def render_table_with_copy_csv(
    df: pd.DataFrame,
    *,
    key_prefix: str,
    csv_filename: str,
    hide_index: bool = False,
    width: Literal["stretch", "content"] = "stretch",
    column_name_map: dict[str, str] | None = None,
) -> None:
    """渲染資料表，並附一鍵複製 Markdown 與 CSV 下載。

    Args:
        width: Streamlit 新版寬度語意；``stretch`` 對應已棄用之 ``use_container_width=True``，
            ``content`` 對應 ``False``。
    """
    df_view = df.rename(columns=column_name_map) if column_name_map else df
    st.dataframe(df_view, width=width, hide_index=hide_index)
    c_l, c_r = st.columns([1, 1])
    with c_l:
        _render_copy_markdown_button(_df_to_markdown_table(df_view), f"copy_md_{key_prefix}")
    with c_r:
        st.download_button(
            "下載結果 CSV（含表頭）",
            df_view.to_csv(index=False).encode("utf-8-sig"),
            csv_filename,
            "text/csv",
            key=f"dl_csv_{key_prefix}",
        )

