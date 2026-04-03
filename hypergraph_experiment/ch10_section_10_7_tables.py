"""
§10.7 論文輸出參數對照表（與《約束世界論》§10.7.5 及雙路徑擴充對齊）。

將程式內分析字典轉為寬表，供 Streamlit 與 CSV 匯出；與扁平化欄位並存。
"""

from __future__ import annotations

from typing import Any

import pandas as pd


def build_section_10_7_5_output_table(analysis: dict[str, Any]) -> pd.DataFrame:
    """
    由 ``analyze_dynamics`` 回傳之 ``analysis`` 字典建立 §10.7.5 輸出參數對照列。

    Args:
        analysis: 動力學分析結果；若含 ``error`` 鍵則仍回傳說明列。

    Returns:
        四欄 DataFrame：論文小節、符號、數值或摘要、備註（程式鍵／接線狀態）。
    """
    if not isinstance(analysis, dict):
        return _single_row_df("—", "—", "—", "分析結果型別異常。")

    if analysis.get("error"):
        return _single_row_df("—", "—", "—", str(analysis.get("error")))

    es = analysis.get("entropy_summary") or {}
    es_wH = analysis.get("entropy_summary_wH") or {}
    series = analysis.get("entropy_time_series")
    series_wH = analysis.get("entropy_time_series_wH")
    n = len(series) if isinstance(series, list) else 0
    n_wH = len(series_wH) if isinstance(series_wH, list) else 0
    w_h = analysis.get("w_h")

    def _fmt(x: Any) -> str:
        if x is None:
            return "—"
        if isinstance(x, float):
            return f"{x:.8g}"
        return str(x)

    rows: list[dict[str, str]] = []

    # （一）熵序列
    if n > 0 and isinstance(series, list):
        s0 = series[0]
        s1 = series[-1]
        sm = es.get("mean")
        summ = _fmt(sm) if sm is not None else "—"
        val = (
            f"長度={n}（時步 0…T 之聚合分布熵）；"
            f"H(首)={_fmt(s0)}；H(末)={_fmt(s1)}；時間平均={summ}"
        )
    else:
        val = "—"
    rows.append(
        {
            "論文小節": "§10.7.5（一）",
            "符號": r"$H_\Lambda^{(\ell)}$（逐步；程式鍵 entropy_time_series）",
            "數值或摘要": val,
            "備註": "完整序列見圖表；扁平 CSV 含 first／last／mean 等摘要欄。",
        }
    )

    # （一-b）窗口熵序列（若 w_H>1）
    if isinstance(w_h, int) and w_h > 1:
        if n_wH > 0 and isinstance(series_wH, list):
            s0 = series_wH[0]
            s1 = series_wH[-1]
            sm = es_wH.get("mean")
            summ = _fmt(sm) if sm is not None else "—"
            val_w = (
                f"長度={n_wH}（w_H={w_h}，窗口 pooled 分布熵）；"
                f"H(首)={_fmt(s0)}；H(末)={_fmt(s1)}；時間平均={summ}"
            )
        else:
            val_w = "—"
        rows.append(
            {
                "論文小節": "§10.7.5（一）",
                "符號": r"$H_{\Lambda,w_H}^{(\ell)}$（程式鍵 entropy_time_series_wH）",
                "數值或摘要": val_w,
                "備註": "w_H=1 時與逐步熵一致；w_H>1 可降低波動（§10.7.2（六））。",
            }
        )

    # （二）平台區長度
    plat_len = (
        es_wH.get("plateau_max_length")
        if isinstance(w_h, int) and w_h > 1
        else es.get("plateau_max_length")
    )
    rows.append(
        {
            "論文小節": "§10.7.5（二）",
            "符號": r"$L_{\mathrm{plat}}$",
            "數值或摘要": _fmt(plat_len),
            "備註": (
                "以 |ΔH|≤ε_plat 連續段之操作化；"
                + ("使用窗口熵 entropy_summary_wH。" if isinstance(w_h, int) and w_h > 1 else "使用逐步熵 entropy_summary。")
            ),
        }
    )

    # （三）週期長度
    p_sum = analysis.get("p_cycle_summary") or {}
    rows.append(
        {
            "論文小節": "§10.7.5（三）",
            "符號": r"$P_{\mathrm{cycle}}$",
            "數值或摘要": (
                f"平均={_fmt(p_sum.get('mean'))}；最小={_fmt(p_sum.get('min'))}；最大={_fmt(p_sum.get('max'))}"
                if isinstance(p_sum, dict)
                else "—"
            ),
            "備註": "以解析類別標籤序列檢測；窗口 w_A 與上限 P_max 見 analysis.w_a／analysis.p_max。",
        }
    )

    # （四）吸引子進入時刻
    ell_sum = analysis.get("ell_a_summary") or {}
    rows.append(
        {
            "論文小節": "§10.7.5（四）",
            "符號": r"$\ell_A$",
            "數值或摘要": (
                f"平均={_fmt(ell_sum.get('mean'))}；最小={_fmt(ell_sum.get('min'))}；最大={_fmt(ell_sum.get('max'))}"
                if isinstance(ell_sum, dict)
                else "—"
            ),
            "備註": "取最早進入「平台（熵差分）」或「週期（解析標籤）」者；平台候選 ell_a_plat 與週期候選 ell_a_cycle_per_run 可供除錯。",
        }
    )

    # （五）可達狀態數
    nr_mean = analysis.get("n_reach_mean")
    nr_min = analysis.get("n_reach_min")
    nr_max = analysis.get("n_reach_max")
    if nr_mean is not None:
        nrv = (
            f"軌跡平均={_fmt(nr_mean)}；單軌最小={_fmt(nr_min)}；單軌最大={_fmt(nr_max)}"
        )
        note = "n_reach_mean／min／max：各軌跡相異解析單元數（類別索引）之彙總。"
    else:
        nrv = "—"
        note = "預期鍵 n_reach_mean（請更新 hypergraph_experiment.core.analyze_dynamics）。"
    rows.append(
        {
            "論文小節": "§10.7.5（五）",
            "符號": r"$N_{\mathrm{reach}}$",
            "數值或摘要": nrv,
            "備註": note,
        }
    )

    # （六）平均合法更新率
    r_adm = analysis.get("r_adm_mean")
    m_trial = analysis.get("m_trial")
    proxy = analysis.get("legal_update_step_fraction_mean")
    rows.append(
        {
            "論文小節": "§10.7.5（六）",
            "符號": r"$\bar r_{\mathrm{adm}}$",
            "數值或摘要": (
                f"{_fmt(r_adm)}（r_adm_mean；M_trial={_fmt(m_trial)}）" if r_adm is not None else "—"
            ),
            "備註": (
                "論文定義為每步 M_adm/M_trial；本版已接線為 r_adm_mean。"
                + (
                    f" 另保留代理量 legal_update_step_fraction_mean={_fmt(proxy)}（相鄰步是否變更）供向後相容。"
                    if proxy is not None
                    else ""
                )
            ),
        }
    )

    # （七）熵波動幅度
    use_es = es_wH if isinstance(w_h, int) and w_h > 1 else es
    hmax = use_es.get("max")
    hmin = use_es.get("min")
    if hmax is not None and hmin is not None:
        d_osc = float(hmax) - float(hmin)
        d_txt = f"{d_osc:.8g}"
    else:
        d_txt = "—"
    rows.append(
        {
            "論文小節": "§10.7.5（七）",
            "符號": r"$\Delta H_{\mathrm{osc}}$",
            "數值或摘要": d_txt,
            "備註": (
                "max(H)−min(H)，取自 "
                + ("entropy_summary_wH（窗口熵）。" if isinstance(w_h, int) and w_h > 1 else "entropy_summary（逐步熵）。")
            ),
        }
    )

    return pd.DataFrame(rows)


def build_experiment_10_7_b_output_table(cmp: dict[str, Any]) -> pd.DataFrame:
    """
    雙解析路徑終端比較之輸出對照（論文擴充，非 §10.7.5 主表）。

    Args:
        cmp: ``run_experiment_10_7`` 回傳字典。

    Returns:
        四欄 DataFrame，與 :func:`build_section_10_7_5_output_table` 同形。
    """
    if not isinstance(cmp, dict):
        return _single_row_df("—", "—", "—", "結果型別異常。")

    if cmp.get("error"):
        return _single_row_df("—", "—", "—", str(cmp.get("error")))

    m = cmp.get("metrics") or {}

    def _fmt(x: Any) -> str:
        if x is None:
            return "—"
        if isinstance(x, float):
            return f"{x:.8g}"
        return str(x)

    rows = [
        {
            "論文小節": "擴充（終端）",
            "符號": r"$\mathrm{JS}$（bit）",
            "數值或摘要": _fmt(m.get("js_terminal_bits")),
            "備註": "metrics.js_terminal_bits；兩路徑終端類別分布之 Jensen–Shannon。",
        },
        {
            "論文小節": "擴充（終端）",
            "符號": r"$H$（路徑 A，bit）",
            "數值或摘要": _fmt(m.get("H_terminal_a_bits")),
            "備註": "終端分布熵（路徑 A 之 Sig、δ）。",
        },
        {
            "論文小節": "擴充（終端）",
            "符號": r"$H$（路徑 B，bit）",
            "數值或摘要": _fmt(m.get("H_terminal_b_bits")),
            "備註": "終端分布熵（路徑 B 之 Sig、δ）。",
        },
        {
            "論文小節": "擴充（終端）",
            "符號": r"$|H_A-H_B|$（bit）",
            "數值或摘要": _fmt(m.get("entropy_abs_diff_terminal")),
            "備註": "metrics.entropy_abs_diff_terminal。",
        },
        {
            "論文小節": "擴充（終端）",
            "符號": "終端類別索引一致率",
            "數值或摘要": _fmt(m.get("terminal_class_agree_rate")),
            "備註": "兩觀測層下終點是否落在同一解析類別索引（離散一致）。",
        },
        {
            "論文小節": "—",
            "符號": "—",
            "數值或摘要": "—",
            "備註": "此表非 §10.7.5 主線輸出；主線見實驗 A 之 §10.7.5 對照表。",
        },
    ]
    return pd.DataFrame(rows)


def _single_row_df(a: str, b: str, c: str, d: str) -> pd.DataFrame:
    """建立單列說明用 DataFrame。"""
    return pd.DataFrame([{"論文小節": a, "符號": b, "數值或摘要": c, "備註": d}])
