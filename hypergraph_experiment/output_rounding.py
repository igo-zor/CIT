"""
實驗輸出浮點數位數統一處理。

供核心分析、第十章實驗模組、Streamlit 批次扁平化與 JSON 序列化共用，
避免介面、CSV 與論文表格式不一致。
"""

from __future__ import annotations

import math
import numbers
from typing import Any

# 全專案實驗結果中，有限浮點數統一保留之小數位數
OUTPUT_FLOAT_DECIMALS: int = 3


def round_floats_for_output(obj: Any, *, ndigits: int | None = None) -> Any:
    """
    遞迴將字典／列表／元組中的有限浮點數四捨五入至指定位數。

    Args:
        obj: 任意巢狀結構；通常為實驗回傳之 ``dict``。
        ndigits: 小數位數；``None`` 時使用 ``OUTPUT_FLOAT_DECIMALS``。

    Returns:
        結構形狀不變之新物件。``None``、``str``、``bool``、整數不改變；
        ``inf``／``nan`` 不改變；無法辨識之型別（如 ``HypergraphConfig``）原樣回傳。

    Note:
        布林須先於 ``numbers.Integral`` 判斷，否則 ``True`` 會被當成 ``1``。
    """
    if ndigits is None:
        ndigits = OUTPUT_FLOAT_DECIMALS
    if obj is None or isinstance(obj, str):
        return obj
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, numbers.Integral):
        return int(obj)
    if isinstance(obj, numbers.Real):
        x = float(obj)
        if not math.isfinite(x):
            return x
        return round(x, ndigits)
    if isinstance(obj, dict):
        return {k: round_floats_for_output(v, ndigits=ndigits) for k, v in obj.items()}
    if isinstance(obj, list):
        return [round_floats_for_output(v, ndigits=ndigits) for v in obj]
    if isinstance(obj, tuple):
        return tuple(round_floats_for_output(v, ndigits=ndigits) for v in obj)
    return obj
