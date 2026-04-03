"""
時間序列摘要：平台區長度與段數（§10.7、§10.9 共用）。

本模組將連續或離散序列之「近平坦區段」與標籤游程長度抽成可 JSON 之純量，
供動力學熵曲線與多尺度宏觀邊數序列共用。
"""

from __future__ import annotations

from typing import Any, Hashable, List, Sequence


def plateau_length_max_abs_diff(
    series: Sequence[float],
    *,
    epsilon: float,
) -> dict[str, float | int | None]:
    """
    依相鄰差分絕對值判定平台：僅當 |x[t+1]-x[t]|<=epsilon 時視為同一段。

    Args:
        series: 實數時間序列（例如熵或有效邊數）。
        epsilon: 平台閾值；超過則切段。

    Returns:
        max_plateau_length：最長段之點數（含端點）。
        num_plateau_segments：長度至少 2 之段數。
        mean_plateau_length：上述段長之算術平均（無段時為 None）。
    """
    if not series:
        return {
            "max_plateau_length": None,
            "num_plateau_segments": 0,
            "mean_plateau_length": None,
        }
    if len(series) == 1:
        return {
            "max_plateau_length": 1,
            "num_plateau_segments": 1,
            "mean_plateau_length": 1.0,
        }

    lengths: List[int] = []
    run_start = 0
    for t in range(len(series) - 1):
        if abs(float(series[t + 1]) - float(series[t])) > float(epsilon):
            seg_len = t + 1 - run_start + 1
            if seg_len >= 2:
                lengths.append(seg_len)
            run_start = t + 1
    last_len = len(series) - run_start
    if last_len >= 2:
        lengths.append(last_len)

    if not lengths:
        return {
            "max_plateau_length": 1,
            "num_plateau_segments": 0,
            "mean_plateau_length": None,
        }

    return {
        "max_plateau_length": int(max(lengths)),
        "num_plateau_segments": len(lengths),
        "mean_plateau_length": round(sum(lengths) / len(lengths), 6),
    }


def state_change_fraction_traj(traj: Sequence[Hashable]) -> float:
    """
    軌跡上「狀態有變」之步數占比；在僅允許合法後繼之動力學下，作為合法更新率之操作化。
    """
    if len(traj) < 2:
        return 0.0
    n = sum(1 for i in range(len(traj) - 1) if traj[i] != traj[i + 1])
    return n / (len(traj) - 1)


def run_length_stats_for_labels(labels: Sequence[Hashable]) -> dict[str, Any]:
    """
    離散標籤序列之游程統計：連續相同標籤視為一個宏觀「單元」。

    Returns:
        tau_unit_max／tau_unit_mean／n_runs 等，對應 §10.9 單元持續時間之簡化代理。
    """
    if not labels:
        return {
            "tau_unit_max": None,
            "tau_unit_mean": None,
            "n_runs": 0,
        }
    runs: List[int] = []
    i = 0
    while i < len(labels):
        j = i + 1
        while j < len(labels) and labels[j] == labels[i]:
            j += 1
        runs.append(j - i)
        i = j
    return {
        "tau_unit_max": int(max(runs)),
        "tau_unit_mean": round(sum(runs) / len(runs), 6) if runs else None,
        "n_runs": len(runs),
    }


def mean_edge_turnover_rate(traj: Sequence[Any]) -> float:
    """
    微觀邊周轉率：相鄰兩態超邊集合對稱差 |E'△E''|，除以 |E'|（至少為 1）後對時間取平均。

    Args:
        traj: 具 ``hyperedges`` 屬性之 ``HypergraphConfig`` 序列。
    """
    if len(traj) < 2:
        return 0.0
    acc = 0.0
    for i in range(len(traj) - 1):
        e0 = getattr(traj[i], "hyperedges", frozenset())
        e1 = getattr(traj[i + 1], "hyperedges", frozenset())
        sym = len(e0 ^ e1)
        denom = max(1, len(e0))
        acc += sym / denom
    return acc / (len(traj) - 1)
