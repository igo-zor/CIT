"""
實驗庫：列出已封存 run、表格式比對、熵曲線重疊、CSV／ZIP 匯出。
"""

from __future__ import annotations

import io
import zipfile
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from hypergraph_experiment.streamlit_common import render_parameters_table, render_table_with_copy_csv
from hypergraph_experiment.storage import default_store_root, list_all_runs, load_manifest


def _flatten_manifest_for_csv(m: dict[str, Any]) -> dict[str, Any]:
    """將 manifest 壓成單列字典，供 CSV 匯出。"""
    params = m.get("parameters") or {}
    analysis = m.get("analysis") if isinstance(m.get("analysis"), dict) else {}
    row: dict[str, Any] = {
        "run_id": m.get("run_id", ""),
        "created_at": m.get("created_at", ""),
        "notes": m.get("notes", ""),
        "num_candidates": m.get("num_candidates"),
        "num_admissible_configs": m.get("num_admissible_configs"),
        "num_obs_configs": m.get("num_obs_configs"),
    }
    for k, v in params.items():
        row[f"param_{k}"] = v
    if "error" in analysis:
        row["analysis_error"] = analysis.get("error")
    else:
        for key in (
            "num_equivalence_classes",
            "entropy_bits",
            "compression_ratio_U",
            "overlap_rate",
            "transitivity_violation_rate",
            "avg_branching_factor",
            "avg_return_time",
        ):
            if key in analysis:
                row[f"analysis_{key}"] = analysis.get(key)
        summ = analysis.get("entropy_summary") or {}
        for sk in ("start", "end", "max", "min", "mean"):
            if sk in summ:
                row[f"analysis_entropy_{sk}"] = summ.get(sk)
    return row


def _params_diff(m1: dict[str, Any], m2: dict[str, Any]) -> pd.DataFrame:
    """兩份 manifest 之 parameters 差異表。"""
    p1 = m1.get("parameters") or {}
    p2 = m2.get("parameters") or {}
    keys = sorted(set(p1) | set(p2))
    rows = []
    for k in keys:
        v1, v2 = p1.get(k), p2.get(k)
        if v1 != v2:
            rows.append({"欄位": k, "run_A": v1, "run_B": v2})
    return pd.DataFrame(rows)


def main() -> None:
    st.title("實驗庫與比對")
    root = default_store_root()
    st.caption(f"倉儲根目錄：`{root}`")

    runs = list_all_runs(root)
    if not runs:
        st.warning("尚無封存實驗。請於主頁執行後按「封存到實驗庫」。")
        return

    df_summary = pd.DataFrame(
        [
            {
                "run_id": r.run_id,
                "created_at": r.created_at,
                "mode": r.mode,
                "n": r.n,
                "delta": r.delta,
                "|Cfg|": r.num_admissible,
                "熵(可比)": r.entropy_or_none,
                "notes": r.notes,
                "path": str(r.path),
            }
            for r in runs
        ]
    )
    render_table_with_copy_csv(
        df_summary,
        key_prefix="exp_lib_summary",
        csv_filename="experiment_summary.csv",
        hide_index=True,
    )

    options = {r.run_id: r.path for r in runs}
    selected_ids = st.multiselect("選取要比對的實驗（可多選）", list(options.keys()), default=[])

    if len(selected_ids) < 1:
        st.info("請至少選取一筆以查看詳情與匯出。")
        return

    manifests: list[dict[str, Any]] = []
    for rid in selected_ids:
        path = options[rid]
        try:
            manifests.append(load_manifest(Path(path)))
        except Exception as e:
            st.error(f"讀取 {rid} 失敗：{e}")

    if not manifests:
        return

    compare_df = pd.DataFrame([_flatten_manifest_for_csv(m) for m in manifests])
    st.subheader("並排指標（表格）")
    st.caption(
        "靜態分析之 **analysis_compression_ratio_U** 為論文 **U_Λ = |S| / N_cfg**（解析單元數／觀測集大小）；"
        "**N_cfg** 即封存檔之 **num_obs_configs**（未設定 n_cfg 抽樣時等於 |Cfg_Λ|）。"
    )
    render_table_with_copy_csv(
        compare_df,
        key_prefix="exp_lib_compare",
        csv_filename="experiment_compare.csv",
        hide_index=True,
    )
    st.subheader("各實驗輸入參數表")
    for idx, m in enumerate(manifests):
        run_id = str(m.get("run_id", f"run_{idx}"))
        st.markdown(f"**{run_id}**")
        render_parameters_table(
            m.get("parameters"),
            key_prefix=f"exp_lib_params_{idx}",
            csv_filename=f"experiment_params_{run_id}.csv",
        )

    st.subheader("動力學：熵曲線重疊")
    dyn_manifests = [m for m in manifests if (m.get("parameters") or {}).get("mode") == "dynamics"]
    series_list: list[tuple[str, list[float]]] = []
    for m in dyn_manifests:
        an = m.get("analysis") or {}
        s = an.get("entropy_time_series") if isinstance(an, dict) else None
        if isinstance(s, list) and s:
            series_list.append((str(m.get("run_id")), [float(x) for x in s]))
    if len(series_list) >= 1:
        fig, ax = plt.subplots(figsize=(8, 4))
        for label, seq in series_list:
            ax.plot(range(len(seq)), seq, marker="o", markersize=2, linewidth=1, label=label[:24])
        ax.set_xlabel("t")
        ax.set_ylabel("H (bits)")
        ax.set_title("解析熵時間序列（多 run 重疊）")
        ax.legend(fontsize=7, loc="best")
        ax.grid(True, alpha=0.3)
        st.pyplot(fig)
        plt.close(fig)
    else:
        st.info("選取之 run 中無動力學熵序列可重疊。")

    st.subheader("兩兩參數差異")
    if len(manifests) >= 2:
        a, b = manifests[0], manifests[1]
        st.markdown(f"比較 **{a.get('run_id')}** 與 **{b.get('run_id')}**（其餘組合請自行選取順序於前兩項）")
        diff_df = _params_diff(a, b)
        if diff_df.empty:
            st.success("兩者 parameters 完全相同。")
        else:
            render_table_with_copy_csv(
                diff_df,
                key_prefix="exp_lib_params_diff",
                csv_filename="experiment_params_diff.csv",
                hide_index=True,
            )
    else:
        st.info("選取至少兩筆可顯示參數差異表。")

    st.subheader("單筆 ZIP 下載")
    for m in manifests:
        rid = m.get("run_id", "run")
        p = Path(options.get(rid, ""))
        if not p.is_dir():
            continue
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for f in sorted(p.rglob("*")):
                if f.is_file():
                    zf.write(f, f.relative_to(p))
        buf.seek(0)
        st.download_button(
            f"下載 {rid}.zip",
            data=buf.getvalue(),
            file_name=f"{rid}.zip",
            mime="application/zip",
            key=f"zip_{rid}",
        )


main()
