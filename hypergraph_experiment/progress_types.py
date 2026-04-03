"""
進度回呼型別：供 Streamlit 進度條或 CLI 靜默模式共用。

callback(done, total, message)：done 為已完成步數（含當前步），total 為總步數（≥1）。
"""

from __future__ import annotations

from typing import Callable, Optional

# 已完成步數、總步數、人類可讀說明
ProgressCallback = Optional[Callable[[int, int, str], None]]
