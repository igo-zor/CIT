"""
第十章論文建議參數（baseline + 掃描集合）。

僅管理頁面預設值與批次模板，不影響核心計算。
"""

from __future__ import annotations

from typing import Final

# §10.2
CH10_2_BASELINE: Final[dict[str, int]] = {
    "n": 8,
    "k_min": 2,
    "k_max": 3,
    "m_max": 10,
    "sample_limit": 5000,
    "n_seed": 20,
}

# §10.3
CH10_3_BASELINE: Final[dict[str, int]] = {
    "n": 8,
    "k_max": 3,
    "m_max": 10,
    "d_max": 4,
    "sample_limit": 5000,
    "n_cfg": 300,
    "n_rep": 20,
    "seed": 20,
    "delta_min": 0,
    "delta_max": 2,
    "s_min": 2,
}

# §10.4
CH10_4_BASELINE: Final[dict[str, float | int]] = {
    "n": 8,
    "k_max": 3,
    "m_max": 10,
    "d_max": 4,
    "sample_limit": 5000,
    "seed": 20,
    "coarse_sample_size": 2000,
    "fine_sample_size": 2000,
    "fiber_sample_size": 200,
    "eps_push_threshold": 0.01,
    "js_threshold": 0.01,
}

# §10.5（《約束世界論 30》§10.5.10 建議參數表）
CH10_5_BASELINE: Final[dict[str, float | int]] = {
    "n_a": 6,
    "n_b": 6,
    "m_edges": 16,
    "k_min": 2,
    "k_max": 3,
    "alpha_cross": 0.30,
    "sample_limit": 5000,
    "n_cfg": 2000,
    "delta_ent": 0,
    "seed": 20,
}

# §10.6（《約束世界論 30》§10.6.6 建議參數表基準）
CH10_6_BASELINE: Final[dict[str, int | str]] = {
    "n_ctx": 1000,
    "n_nodes": 8,
    "M": 4,
    "w_ctx": 2,
    "eta_ctx": 1,
    "T_loc": 2,
    "n_search": 5000,
    "mode": "obstruction",
    "seed": 20,
}

# §10.7（《約束世界論 30》§10.7.6 建議主線固定組合 + 域型預設）
CH10_7_BASELINE: Final[dict[str, float | int | str]] = {
    "n": 8,
    "k_max": 3,
    "m_max": 10,
    "sample_limit": 5000,
    "signature": "medium",
    "delta": 1,
    "runs": 50,
    "steps": 200,
    "seed": 30,
    "eps_plat": 0.01,
    "d_max": 4,
    # 論文主線：N_0=50, T=200, M_trial=10, Λ_obs=medium, w_H=10, w_A=20, N_seed=30, P_max=20
    "m_trial": 10,
    "w_h": 10,
    "w_a": 20,
    "p_max": 20,
    "n_seed_107": 30,
}

# §10.7.4 建議離散掃描集合（介面選單用）
CH10_7_T_CHOICES: Final[tuple[int, ...]] = (50, 100, 200, 500, 1000)
CH10_7_N0_CHOICES: Final[tuple[int, ...]] = (20, 50, 100, 200)

# §10.8
CH10_8_BASELINE: Final[dict[str, float | int | str]] = {
    "n": 12,
    "m": 18,
    "r": 2,
    "eta": 0.10,
    "T_sb": 20,
    "n_samples": 30,
    "sig_obs": "medium",
    "seed": 20,
}

# §10.9（《約束世界論 30》§10.9.6 建議主組合 + 域型預設）
CH10_9_BASELINE: Final[dict[str, int | float | str | bool]] = {
    "n": 6,
    "k_max": 3,
    "m_max": 10,
    "sample_limit": 5000,
    "d_max": 5,
    "connected": False,
    "r": 2,
    "steps": 300,
    "window_list": "1,2,4,8,16",
    "seed": 20,
    "m_trial": 10,
    "sig_obs": "medium",
    "delta_t": 1,
    "epsilon_plat": 0.01,
    "n_hist": 30,
    # 論文 §10.9.6 建議 N_seed=20（穩健性）；鍵名 n_seed_runs 避免與 §10.2 種子批次數混淆
    "n_seed_runs": 1,
}
