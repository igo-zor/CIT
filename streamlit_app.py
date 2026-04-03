"""
約束世界論 · 第十章實驗平台入口。

啟動：於專案根目錄執行 ``streamlit run streamlit_app.py``（建議使用虛擬環境）。
左側導覽可切換 §10.2–§10.9 各實驗頁與實驗庫。
"""

from __future__ import annotations

import subprocess
import sys

import streamlit as st

from hypergraph_experiment.ch10_symbol_glossary import MD_10_1_1_CORE
from hypergraph_experiment.storage import default_store_root


def main() -> None:
    st.set_page_config(page_title="約束世界論 · 第十章實驗", layout="wide")
    st.title("第十章　有限超圖實驗平台")
    st.markdown(
        """
        本平台依論文章節分頁：**§10.2** 配置域梯子、**§10.3** 靜態解析覆蓋、**§10.4** 細化、
        **§10.5** 不可分解、**§10.6** 拼合障礙、**§10.7** 動力學、**§10.8** 對稱破缺、**§10.9** 多尺度視窗。
        請自左欄選取頁面；各頁頂部附實驗簡述與建議參數，並支援 **批次參數表** 與計算 **進度條**。
        """
    )
    root = default_store_root()
    st.info(f"實驗倉儲目錄：**{root}**")
    with st.expander("論文符號總表（§10.1.1　統一符號與操作型定義）", expanded=False):
        st.markdown(MD_10_1_1_CORE)
    st.markdown(
        """
        - **實驗庫**：`pages/2_Experiment_library.py`（封存 run 比對）
        - **輕量 CLI**：`python experiment.py --help`
        """
    )


if __name__ == "__main__":
    try:
        from streamlit.runtime.scriptrunner import get_script_run_ctx

        if get_script_run_ctx() is not None:
            main()
        else:
            rc = subprocess.call([sys.executable, "-m", "streamlit", "run", __file__, *sys.argv[1:]])
            raise SystemExit(rc)
    except ImportError:
        rc = subprocess.call([sys.executable, "-m", "streamlit", "run", __file__, *sys.argv[1:]])
        raise SystemExit(rc)
