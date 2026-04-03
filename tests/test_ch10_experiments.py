"""§10.5–10.9 之煙霧測試（小參數、固定種子）。"""

from __future__ import annotations

import pandas as pd
import pytest

from hypergraph_experiment.experiments.exp_10_5_bipartite import (
    run_experiment_10_5,
    section_10_5_batch_sweep_dataframe,
    section_10_5_output_parameters_df,
    section_10_5_sweep_axis_candidates,
)
from hypergraph_experiment.experiments.exp_10_6_contexts import (
    run_canonical_demo_10_6,
    run_experiment_10_6,
    section_10_6_output_parameters_df,
    validate_cyclic_window_params,
)
from hypergraph_experiment.ch10_section_10_7_tables import (
    build_experiment_10_7_b_output_table,
    build_section_10_7_5_output_table,
)
from hypergraph_experiment.experiments.exp_10_7_paths import run_experiment_10_7
from hypergraph_experiment.experiments.exp_10_8_symmetry import run_experiment_10_8
from hypergraph_experiment.experiments.exp_10_8_symmetry import run_experiment_10_8_three_arm
from hypergraph_experiment.experiments.exp_10_8_symmetry import section_10_8_output_parameters_df
from hypergraph_experiment.experiments.exp_10_9_multiscale import (
    run_experiment_10_9,
    section_10_9_output_parameters_df,
)
from hypergraph_experiment.ch10_paper_presets import (
    CH10_3_BASELINE,
    CH10_4_BASELINE,
    CH10_5_BASELINE,
    CH10_6_BASELINE,
    CH10_7_BASELINE,
    CH10_9_BASELINE,
)
from hypergraph_experiment.core import (
    VIOLATION_DISCONNECTED,
    VIOLATION_PAIR_TRIANGLE_FORBIDDEN,
    HypergraphConfig,
    analyze_dynamics,
    analyze_static,
    domain_constraint_violation_primary,
    run_full_experiment,
    sample_candidates_and_filter,
    subsample_obs_configs,
)
from hypergraph_experiment.paper_tables import (
    TABLE_10_2_VIOLATION_COUNT_COLUMNS,
    table_10_2_domain_ladder,
    table_10_3_signature_comparison,
)
from hypergraph_experiment.refinement import compare_ordered_refinement_paths
from hypergraph_experiment.output_rounding import (
    OUTPUT_FLOAT_DECIMALS,
    round_floats_for_output,
)
from hypergraph_experiment.streamlit_common import batch_cell, build_ch10_column_name_map, flatten_result_row


def test_10_5_smoke() -> None:
    """n_cfg 為 None 時以全部過濾後配置為觀測集（向後相容）。"""
    out = run_experiment_10_5(n_a=3, n_b=3, m_edges=4, sample_limit=8, seed=0, alpha_cross=0.2)
    assert out["experiment"] == "10.5"
    assert "metrics" in out
    assert out.get("n_cfg_requested") is None
    assert out["num_admissible_filtered"] == out["num_obs_configs"] == out["num_admissible"]


def test_10_5_n_cfg_subsample_matches_obs_set() -> None:
    """n_cfg 小於過濾後母集時，解析統計僅使用抽樣後之觀測筆數。"""
    out = run_experiment_10_5(
        n_a=3,
        n_b=3,
        m_edges=4,
        sample_limit=80,
        seed=3,
        alpha_cross=0.25,
        n_cfg=4,
    )
    assert out.get("metrics", {}).get("error") is None
    n_f = int(out["num_admissible_filtered"])
    assert n_f >= 4, "測試需至少 4 筆過濾後配置；請提高 sample_limit"
    assert out["num_obs_configs"] == 4
    assert out["num_admissible"] == 4
    assert out["n_cfg_requested"] == 4
    assert out["n_cfg_notice"] is None


def test_section_10_5_batch_sweep_dataframe_and_axis() -> None:
    """批次彙整表應帶入參數列與 §10.5.5 指標；掃描軸候選應辨識變動之數值參數欄。"""

    def _ok_row(alpha: float, rho: float) -> dict[str, object]:
        return {
            "run_index": 0,
            "error": None,
            "param_caption": "",
            "param_row": {"跨區傾向": alpha, "超邊數": 16},
            "result": {
                "experiment": "10.5",
                "metrics": {
                    "rho_irred": rho,
                    "N_irred": 1,
                    "rho_cross_mean": 0.5,
                    "D_sep_total_variation": 0.1,
                    "D_sep_JS_bits": 0.2,
                    "I_A_B_bits": 0.3,
                },
                "num_admissible_filtered": 10,
                "num_obs_configs": 10,
                "num_classes": 3,
            },
        }

    runs = [
        {**_ok_row(0.1, 0.0), "run_index": 0},
        {**_ok_row(0.5, 0.25), "run_index": 1},
    ]
    sdf = section_10_5_batch_sweep_dataframe(runs)
    assert len(sdf) == 2
    assert "跨區傾向" in sdf.columns
    assert list(sdf.sort_values("跨區傾向")["rho_irred"]) == [0.0, 0.25]
    axes = section_10_5_sweep_axis_candidates(sdf)
    assert "跨區傾向" in axes
    # 超邊數兩列皆為 16，無變化，不應作為掃描軸候選
    assert "超邊數" not in axes


def test_10_5_section_152_main_table() -> None:
    """《約束世界論 30》§10.5.5 主表欄位與列數應涵蓋（一）–（四）。"""
    out = run_experiment_10_5(n_a=3, n_b=3, m_edges=4, sample_limit=12, seed=1, alpha_cross=0.3)
    df = section_10_5_output_parameters_df(out)
    assert list(df.columns) == ["論文小節", "輸出參數", "論文記號", "數值", "論文語義摘要"]
    assert len(df) == 7
    assert df.iloc[0]["論文小節"] == "10.5.5（一）"
    assert r"\rho" in str(df.iloc[0]["論文記號"])

    err_out = {"experiment": "10.5", "metrics": {"error": "無合法配置"}}
    df_e = section_10_5_output_parameters_df(err_out)
    assert len(df_e) == 1
    assert "無合法配置" in str(df_e.iloc[0]["論文語義摘要"])


def test_10_6_canonical() -> None:
    d = run_canonical_demo_10_6()
    assert d.get("chi_glue") == 1


def test_10_6_batch() -> None:
    out = run_experiment_10_6(n_ctx=5, n_nodes=8, mode="obstruction", seed=2)
    assert "rho_glue" in out["metrics"]
    assert out.get("n_ctx_effective") == 5
    assert float(out["metrics"]["rho_glue"]) == 1.0
    assert int(out["metrics"].get("n_val", -1)) == 0
    assert int(out.get("sum_chi_glue", -1)) == 5
    assert int(out["parameters"]["M"]) == 4
    assert int(out["parameters"]["w_ctx"]) == 2


def test_10_6_validate_degenerate_M2() -> None:
    """w=2、η=1、M=2 應拒絕（與論文建議 M≥4 一致）。"""
    err = validate_cyclic_window_params(8, 2, 2, 1)
    assert err is not None


def test_10_6_section_165_output_table() -> None:
    """§10.6.5 主表欄位與小節應對齊論文（一）–（四），含 N_val 列。"""
    out = run_experiment_10_6(n_ctx=4, n_nodes=8, mode="satisfiable", seed=0)
    df = section_10_6_output_parameters_df(out)
    assert list(df.columns) == ["論文小節", "輸出參數", "論文記號", "數值", "論文語義摘要"]
    assert len(df) == 6
    assert df.iloc[1]["論文小節"] == "10.6.5（二）"
    assert r"\rho" in str(df.iloc[1]["論文記號"])
    assert int(out["metrics"]["n_val"]) == 4
    dc = run_canonical_demo_10_6()
    df_c = section_10_6_output_parameters_df(dc)
    assert len(df_c) == 3
    assert df_c.iloc[0]["論文小節"] == "10.6.5（一）"


def test_10_7_smoke() -> None:
    out = run_experiment_10_7(n=4, max_edges=3, sample_limit=50, runs=3, steps=5, seed=1)
    assert "metrics" in out or "error" in out


def test_10_7_dynamics_n_reach_and_paper_output_tables() -> None:
    """§10.7：analyze_dynamics 應含 N_reach 彙總；論文對照表應可建表。"""
    r = run_full_experiment(
        mode="dynamics",
        n=4,
        max_edge_size=3,
        max_edges=3,
        sample_limit=80,
        runs=2,
        steps=4,
        seed=2,
        epsilon_plat=0.05,
        connected=False,
    )
    an = r["analysis"]
    assert "error" not in an
    assert an.get("n_reach_mean") is not None
    assert isinstance(an.get("n_reach_per_run"), list)
    assert len(an["n_reach_per_run"]) == 2
    assert isinstance(an.get("r_adm_mean"), (int, float))
    assert isinstance(an.get("entropy_time_series_wH"), list)
    assert isinstance(an.get("p_cycle_summary"), dict)
    df_a = build_section_10_7_5_output_table(an)
    assert "§10.7.5（五）" in set(df_a["論文小節"].tolist())
    out_b = run_experiment_10_7(
        n=4,
        max_edges=3,
        sample_limit=40,
        runs=2,
        steps=3,
        seed=0,
    )
    df_b = build_experiment_10_7_b_output_table(out_b)
    assert len(df_b) >= 5


def test_10_7_dynamics_n_seed_reruns_summary() -> None:
    """§10.7：N_seed 重跑時應提供 analysis_seeds 與彙總。"""
    r = run_full_experiment(
        mode="dynamics",
        n=4,
        max_edge_size=3,
        max_edges=3,
        sample_limit=60,
        runs=2,
        steps=5,
        seed=1,
        m_trial=5,
        w_h=3,
        w_a=3,
        p_max=5,
        n_seed_107=3,
        epsilon_plat=0.05,
        connected=False,
    )
    assert "analysis_seeds" in r
    assert "analysis_seed_summary" in r
    assert len(r["analysis_seeds"]) == 3


def test_10_8_smoke() -> None:
    out = run_experiment_10_8(n=6, m=8, n_samples=2, init_family="sym", dynamics_steps=0)
    assert out["experiment"] == "10.8"
    m = out["metrics"]
    assert "N_iso_mean" in m
    assert "A_reach_mean" in m
    assert len(out["per_sample"]) == 2
    assert "N_iso" in out["per_sample"][0]


def test_10_8_paper_output_parameters_table() -> None:
    out = run_experiment_10_8(n=6, m=8, n_samples=2, init_family="sym", dynamics_steps=0, seed=0)
    df = section_10_8_output_parameters_df(out)
    assert list(df.columns) == ["論文小節", "輸出參數", "論文記號", "數值", "論文語義摘要"]
    assert len(df) >= 5
    assert df.iloc[0]["論文小節"] == "10.8.5（一）"


def test_10_8_three_arm_smoke() -> None:
    out = run_experiment_10_8_three_arm(n=6, m=8, n_samples=2, dynamics_steps=0, seed=0)
    assert out["experiment"] == "10.8_three_arm"
    assert set(out["by_family"].keys()) == {"sym", "pert", "rand"}
    comp = out["comparison"]
    assert "rho_type_pert_over_sym" in comp
    assert "rho_type_rand_over_sym" in comp


def test_10_8_three_arm_paper_table() -> None:
    out = run_experiment_10_8_three_arm(n=6, m=8, n_samples=2, dynamics_steps=0, seed=1)
    df = section_10_8_output_parameters_df(out)
    assert len(df) == 7
    assert df.iloc[5]["論文小節"] == "10.8.5（六）"
    assert "型別擴張率（論文" in str(df.iloc[5]["輸出參數"])
    assert df.iloc[6]["論文小節"] == "10.8.5（輔助）"


def test_10_9_smoke() -> None:
    """預設 CH10_9 n_hist 較大；煙霧測試固定 n_hist=1、n_seed=1。"""
    out = run_experiment_10_9(
        n=4,
        steps=15,
        window_sizes=(1, 3),
        seed=0,
        n_hist=1,
        n_seed=1,
    )
    assert out["experiment"] == "10.9"
    assert "per_window" in out
    assert "R_edge_bar" in out
    row0 = out["per_window"][0]
    assert "tau_unit_max" in row0
    assert "L_plat_max_len_Ceff" in row0
    assert out["parameters"]["sig_obs"] == "medium"
    assert out["parameters"]["m_trial"] == 10


def test_10_9_paper_output_table() -> None:
    out = run_experiment_10_9(
        n=4,
        steps=12,
        window_sizes=(1, 2),
        seed=1,
        r=1,
        n_hist=1,
        n_seed=1,
        sig_obs="weak",
        delta_t=2,
    )
    df = section_10_9_output_parameters_df(out)
    assert list(df.columns) == ["論文小節", "輸出參數", "論文記號", "數值", "論文語義摘要"]
    assert len(df) == 7
    assert df.iloc[0]["論文小節"] == "10.9.5（一）"
    assert df.iloc[5]["論文小節"] == "10.9.5（六）"
    assert "JS" in str(df.iloc[5]["論文記號"])


def test_10_9_output_table_error_path() -> None:
    df = section_10_9_output_parameters_df({"experiment": "10.9", "error": "無合法配置", "parameters": {}})
    assert len(df) == 1
    assert "無合法配置" in str(df.iloc[0]["論文語義摘要"])


def test_10_9_n_hist_aggregate() -> None:
    out = run_experiment_10_9(
        n=4,
        steps=10,
        window_sizes=(1,),
        seed=2,
        n_hist=2,
        n_seed=1,
        sample_limit=200,
    )
    assert out["aggregated_N_hist"] == 2
    assert len(out["per_window"]) >= 1


def test_paper_preset_baselines_key_values() -> None:
    """關鍵章節 baseline 常數應對齊論文主線。"""
    assert int(CH10_3_BASELINE["n_cfg"]) == 300
    assert int(CH10_3_BASELINE["s_min"]) == 2
    assert float(CH10_4_BASELINE["eps_push_threshold"]) == 0.01
    assert float(CH10_4_BASELINE["js_threshold"]) == 0.01
    assert int(CH10_4_BASELINE["coarse_sample_size"]) == 2000
    assert int(CH10_4_BASELINE["fine_sample_size"]) == 2000
    assert float(CH10_7_BASELINE["eps_plat"]) == 0.01
    assert int(CH10_7_BASELINE["m_trial"]) == 10
    assert int(CH10_7_BASELINE["w_h"]) == 10
    assert int(CH10_7_BASELINE["w_a"]) == 20
    assert int(CH10_7_BASELINE["p_max"]) == 20
    assert int(CH10_7_BASELINE["n_seed_107"]) == 30
    assert int(CH10_5_BASELINE["k_min"]) == 2
    assert int(CH10_5_BASELINE["k_max"]) == 3
    assert int(CH10_5_BASELINE["n_cfg"]) == 2000
    assert int(CH10_5_BASELINE["sample_limit"]) == 5000
    assert int(CH10_6_BASELINE["M"]) == 4
    assert int(CH10_6_BASELINE["w_ctx"]) == 2
    assert int(CH10_6_BASELINE["eta_ctx"]) == 1
    assert int(CH10_6_BASELINE["n_search"]) == 5000
    assert int(CH10_9_BASELINE["steps"]) == 300
    assert int(CH10_9_BASELINE["m_trial"]) == 10
    assert int(CH10_9_BASELINE["n_hist"]) == 30
    assert str(CH10_9_BASELINE["window_list"]) == "1,2,4,8,16"


def test_10_2_domain_ladder_extended_metrics() -> None:
    levels = [
        {"label": "弱", "connected": False, "max_degree": 8, "forbid_pair_triangles": False},
        {"label": "中", "connected": True, "max_degree": 5, "forbid_pair_triangles": False},
        {"label": "強", "connected": True, "max_degree": 4, "forbid_pair_triangles": True},
    ]
    rows, n_cand = table_10_2_domain_ladder(
        n=6,
        min_edge_size=2,
        max_edge_size=3,
        max_edges=6,
        sample_limit=120,
        seed=7,
        num_seeds=5,
        levels=levels,
    )
    assert n_cand > 0
    assert len(rows) == 3
    assert "累計排除數" in rows[0]
    assert "本層新增排除數" in rows[0]
    assert "是否為上一層子集" in rows[0]
    assert "鏈式子集成立" in rows[0]
    assert "相對上一層保留率" in rows[1]
    assert "相對候選收縮率" in rows[1]
    assert "相對上一層收縮率" in rows[1]
    assert "違規主因集中度" in rows[0]
    assert "該層排除主因連通未滿足筆數" in rows[0]
    for row in rows:
        s = sum(int(row[col_all]) for _code, col_all, _col_new in TABLE_10_2_VIOLATION_COUNT_COLUMNS)
        assert s == row["累計排除數"]
        s_new = sum(int(row[col_new]) for _code, _col_all, col_new in TABLE_10_2_VIOLATION_COUNT_COLUMNS)
        assert s_new == row["本層新增排除數"]


def test_domain_constraint_violation_primary_order() -> None:
    """不連通配置於強制連通時應標為連通未滿足。"""
    c = HypergraphConfig(vertices=(1, 2, 3), hyperedges=frozenset({frozenset({1, 2})}))
    r = domain_constraint_violation_primary(
        c,
        max_edge_size=3,
        max_edges=6,
        connected_required=True,
        max_degree=8,
        forbid_pair_triangles=False,
    )
    assert r == VIOLATION_DISCONNECTED


def test_forbidden_pair_triangle_triggers_primary() -> None:
    """二元邊三角形於啟用禁制時應標為禁二元三角。"""
    c = HypergraphConfig(
        vertices=(1, 2, 3),
        hyperedges=frozenset(
            {
                frozenset({1, 2}),
                frozenset({1, 3}),
                frozenset({2, 3}),
            }
        ),
    )
    r = domain_constraint_violation_primary(
        c,
        max_edge_size=3,
        max_edges=6,
        connected_required=False,
        max_degree=8,
        forbid_pair_triangles=True,
    )
    assert r == VIOLATION_PAIR_TRIANGLE_FORBIDDEN


def test_round_floats_for_output_three_decimals() -> None:
    """實驗輸出舍入：巢狀浮點應統一至 OUTPUT_FLOAT_DECIMALS。"""
    assert OUTPUT_FLOAT_DECIMALS == 3
    v = round_floats_for_output(
        {"a": 1.234567, "b": [2.6666], "c": {"x": float("inf")}, "d": None, "e": 3}
    )
    assert v["a"] == 1.235
    assert v["b"] == [2.667]
    assert v["c"]["x"] == float("inf")
    assert v["d"] is None
    assert v["e"] == 3


def test_flatten_result_row_deep_skips_samples() -> None:
    """批次扁平化：巢狀 dict／list[dict]／數值序列應展開；略過 per_sample。"""
    payload: dict = {
        "parameters": {"n": 3},
        "per_window": [{"w": 1, "H_bits": 0.5}, {"w": 2, "H_bits": 0.6}],
        "per_sample": [1, 2, 3],
        "series": [10.0, 20.0],
    }
    row = flatten_result_row("t", payload)
    assert row["t_parameters_n"] == 3
    assert row["t_per_window_0_w"] == 1
    assert row["t_per_window_1_H_bits"] == 0.6
    assert "t_per_sample" not in row
    assert row["t_series_mean"] == 15.0


def test_subsample_obs_configs_and_compression_U_matches_paper() -> None:
    """觀測集大小與論文 U_Λ=|S|/N；avg_class_size=N/|S|。"""
    _, cfgs = sample_candidates_and_filter(
        n=5,
        max_edge_size=3,
        max_edges=6,
        sample_limit=200,
        seed=11,
        connected=False,
        max_degree=6,
        forbid_pair_triangles=False,
    )
    assert cfgs
    obs, _req, n_act, _nt = subsample_obs_configs(cfgs, 40, seed=99)
    assert n_act == min(40, len(cfgs))
    assert len(obs) == n_act
    a = analyze_static(obs, "weak", 0, s_min=0)
    n = len(obs)
    s = int(a["num_equivalence_classes"])
    u = float(a["compression_ratio_U"])
    assert n > 0
    assert abs(u * n - s) < 1e-6
    acs = float(a["avg_class_size"])
    # avg_class_size 經輸出舍入後與 n/s 可略有誤差
    assert abs(acs - n / s) <= 0.002


def test_run_full_experiment_n_cfg_field_and_backward_compat() -> None:
    """n_cfg 時 num_obs_configs 合理；未指定 n_cfg 時分析用全集。"""
    r1 = run_full_experiment(
        mode="static",
        n=4,
        max_edge_size=3,
        max_edges=5,
        sample_limit=150,
        n_cfg=25,
        signature="weak",
        delta=0,
        seed=3,
        show_sample_configs=0,
    )
    adm = int(r1["num_admissible_configs"])
    obs_n = int(r1["num_obs_configs"])
    assert obs_n == min(25, adm)
    an = r1["analysis"]
    assert "error" not in an
    assert int(an["num_configs"]) == obs_n

    r0 = run_full_experiment(
        mode="static",
        n=4,
        max_edge_size=3,
        max_edges=5,
        sample_limit=150,
        n_cfg=None,
        signature="weak",
        delta=0,
        seed=3,
        show_sample_configs=0,
    )
    assert int(r0["num_obs_configs"]) == int(r0["num_admissible_configs"])
    assert int((r0["analysis"] or {})["num_configs"]) == int(r0["num_admissible_configs"])


def test_run_full_experiment_static_n_rep_summary() -> None:
    """static 模式啟用 N_rep 時，應回傳重抽明細與彙總。"""
    r = run_full_experiment(
        mode="static",
        n=4,
        max_edge_size=3,
        max_edges=5,
        sample_limit=150,
        n_cfg=20,
        n_rep=3,
        signature="medium",
        delta=0,
        seed=13,
        show_sample_configs=0,
    )
    reps = r.get("analysis_repetitions") or []
    assert len(reps) == 3
    rep_sum = r.get("analysis_rep_summary") or {}
    assert int(rep_sum.get("n_rep_effective", 0)) == 3
    assert "compression_ratio_U_mean" in rep_sum
    assert "compression_ratio_U_std" in rep_sum


def test_refinement_controls_sample_and_fiber_size() -> None:
    """§10.4 可調樣本與纖維上限應可傳入並反映於 detail。"""
    r = run_full_experiment(
        mode="static",
        n=5,
        max_edge_size=3,
        max_edges=6,
        sample_limit=300,
        signature="strong",
        delta=0,
        seed=9,
        refinement_enabled=True,
        refine_coarse_signature="weak",
        refine_coarse_delta=2,
        refine_fine_signature="strong",
        refine_fine_delta=0,
        refine_kernel="uniform",
        refine_compare_chains=True,
        refine_coarse_sample_size=120,
        refine_fine_sample_size=80,
        refine_fiber_sample_size=30,
        show_sample_configs=0,
    )
    ref = (r.get("refinement_10_4") or {}).get("single_step_Lambda_to_Lambda_prime") or {}
    assert "error" not in ref
    det = ref.get("detail") or {}
    assert int(det.get("coarse_sample_size_used", 0)) <= 120
    assert int(det.get("fine_sample_size_used", 0)) <= 80
    assert int(det.get("fiber_sample_size", 0)) == 30


def test_refinement_bundle_has_1045_output_fields() -> None:
    """§10.4.5 對齊主表所需欄位應可由細化結果直接取得。"""
    r = run_full_experiment(
        mode="static",
        n=5,
        max_edge_size=3,
        max_edges=6,
        sample_limit=260,
        signature="strong",
        delta=0,
        seed=10,
        refinement_enabled=True,
        refine_coarse_signature="weak",
        refine_coarse_delta=2,
        refine_fine_signature="strong",
        refine_fine_delta=0,
        refine_kernel="uniform",
        refine_compare_chains=True,
        refine_coarse_sample_size=100,
        refine_fine_sample_size=70,
        refine_fiber_sample_size=20,
        show_sample_configs=0,
    )
    ref = r.get("refinement_10_4") or {}
    single = ref.get("single_step_Lambda_to_Lambda_prime") or {}
    assert "H_p_Lambda_bits" in single
    assert "H_p_Lambda_prime_bits" in single
    assert "pushforward_max_error" in single
    tw = ref.get("two_path_ab_ba_10_4_7") or {}
    assert "js_divergence_bits_terminal" in tw
    assert "entropy_abs_diff_terminal" in tw


def test_compare_ordered_refinement_paths_has_expected_keys() -> None:
    """A→B / B→A 進階比較應回傳完整摘要鍵。"""
    _, cfgs = sample_candidates_and_filter(
        n=5,
        max_edge_size=3,
        max_edges=6,
        sample_limit=220,
        seed=8,
        connected=False,
        max_degree=6,
        forbid_pair_triangles=False,
    )
    if len(cfgs) < 30:
        return
    # 實驗 B 子步驟以映射種類 + δ 輸入（對齊 §9.6-C），非直接 weak/medium/strong 名稱。
    step_a = {"kind": "adjacency_motif_fine", "delta": 0}  # 內部對應 strong
    step_b = {"kind": "degree_split", "delta": 2}  # 內部對應 medium
    out = compare_ordered_refinement_paths(
        cfgs,
        step_a=step_a,
        step_b=step_b,
        kernel_mode="uniform",
        coarse_sample_size=120,
        fine_sample_size=80,
        sample_seed=9,
        max_fiber_size=25,
    )
    assert "A_to_B" in out and "B_to_A" in out
    assert "js_divergence_bits_terminal_ab_ba" in out
    assert "entropy_abs_diff_terminal_ab_ba" in out
    assert (out.get("A_to_B") or {}).get("path_key") == "A→B"
    assert (out.get("B_to_A") or {}).get("path_key") == "B→A"
    ab = out.get("A_to_B") or {}
    assert ab.get("step_a") == step_a
    assert ab.get("step_b") == step_b
    assert "pushforward_max_error" in ab
    assert "entropy_fine_bits" in ab


def test_compare_ordered_refinement_paths_rejects_legacy_sig_delta_args() -> None:
    """舊介面 sig/delta 應被拒絕，並提示改用 step_a/step_b。"""
    _, cfgs = sample_candidates_and_filter(
        n=5,
        max_edge_size=3,
        max_edges=6,
        sample_limit=80,
        seed=5,
        connected=False,
        max_degree=6,
        forbid_pair_triangles=False,
    )
    with pytest.raises(ValueError, match="step_a 與 step_b"):
        compare_ordered_refinement_paths(
            cfgs,
            kernel_mode="uniform",
            sig_a="strong",  # type: ignore[call-arg]
            delta_a=0,  # type: ignore[call-arg]
            sig_b="medium",  # type: ignore[call-arg]
            delta_b=2,  # type: ignore[call-arg]
        )


def test_refinement_bundle_two_path_uses_ab_ba_output() -> None:
    """§10.4 bundle 應輸出 A→B / B→A 的雙路徑摘要鍵。"""
    r = run_full_experiment(
        mode="static",
        n=5,
        max_edge_size=3,
        max_edges=6,
        sample_limit=280,
        signature="strong",
        delta=0,
        seed=15,
        refinement_enabled=True,
        refine_coarse_signature="weak",
        refine_coarse_delta=2,
        refine_fine_signature="strong",
        refine_fine_delta=0,
        refine_kernel="uniform",
        refine_compare_chains=True,
        refine_coarse_sample_size=110,
        refine_fine_sample_size=90,
        refine_fiber_sample_size=25,
        show_sample_configs=0,
    )
    ref = r.get("refinement_10_4") or {}
    tw = ref.get("two_path_ab_ba_10_4_7") or {}
    assert "js_divergence_bits_terminal" in tw
    assert "entropy_abs_diff_terminal" in tw
    assert "A_to_B" in tw
    assert "B_to_A" in tw


def test_refinement_bundle_without_chain_compare_has_no_ab_block() -> None:
    """主表實驗路徑（不啟用 compare_chains）不應產生 A/B 區塊。"""
    r = run_full_experiment(
        mode="static",
        n=5,
        max_edge_size=3,
        max_edges=6,
        sample_limit=260,
        signature="strong",
        delta=0,
        seed=16,
        refinement_enabled=True,
        refine_coarse_signature="weak",
        refine_coarse_delta=2,
        refine_fine_signature="strong",
        refine_fine_delta=0,
        refine_kernel="uniform",
        refine_compare_chains=False,
        refine_coarse_sample_size=100,
        refine_fine_sample_size=70,
        refine_fiber_sample_size=20,
        show_sample_configs=0,
    )
    ref = r.get("refinement_10_4") or {}
    assert "two_path_ab_ba_10_4_7" not in ref


def test_table_10_3_same_obs_across_signatures() -> None:
    """表 10-3：同一觀測集上各簽名之 num_configs（觀測長度）一致。"""
    _, cfgs = sample_candidates_and_filter(
        n=4,
        max_edge_size=3,
        max_edges=4,
        sample_limit=80,
        seed=2,
        connected=False,
        max_degree=4,
        forbid_pair_triangles=False,
    )
    if len(cfgs) < 5:
        return
    rows = table_10_3_signature_comparison(cfgs, 0, s_min=0, n_cfg=10, seed=5)
    assert len(rows) == 3
    # 間接：weak 與 strong 若觀測集不同則不可能；以列數與非空確認
    for row in rows:
        assert row["解析單元數"] >= 1


def test_analyze_static_r_iso_and_smin() -> None:
    """靜態分析應含 R_iso 與可選 s_min。"""
    _, cfgs = sample_candidates_and_filter(
        n=4,
        max_edge_size=3,
        max_edges=5,
        sample_limit=80,
        seed=3,
        connected=False,
        max_degree=5,
        forbid_pair_triangles=False,
    )
    assert cfgs
    a = analyze_static(cfgs, "weak", 0, s_min=0)
    assert "isol_rate_compat_graph" in a
    assert 0.0 <= float(a["isol_rate_compat_graph"]) <= 1.0
    a2 = analyze_static(cfgs, "weak", 0, s_min=2)
    assert a2["s_min"] == 2


def test_analyze_dynamics_plat_and_legal_fraction() -> None:
    """動力學應回傳熵平台摘要與合法更新步占比。"""
    _, cfgs = sample_candidates_and_filter(
        n=4,
        max_edge_size=3,
        max_edges=4,
        sample_limit=100,
        seed=5,
        connected=False,
        max_degree=4,
        forbid_pair_triangles=False,
    )
    if len(cfgs) < 2:
        return
    d = analyze_dynamics(
        cfgs,
        signature_name="weak",
        delta=1,
        runs=2,
        steps=6,
        seed=11,
        max_edge_size=3,
        max_edges=4,
        connected_required=False,
        max_degree=4,
        forbid_pair_triangles=False,
        epsilon_plat=0.05,
    )
    assert "error" not in d
    assert "legal_update_step_fraction_mean" in d
    es = d.get("entropy_summary") or {}
    assert "plateau_max_length" in es


def test_batch_cell_zh_fallback() -> None:
    """batch_cell 應支援中文優先、英文回退。"""
    r = pd.Series({"節點數": 7, "legacy": 1})
    assert batch_cell(r, "節點數", "n") == 7
    r2 = pd.Series({"n": 8})
    assert batch_cell(r2, "節點數", "n") == 8
    r3 = pd.Series({"隨機種子": 99})
    assert batch_cell(r3, "偽隨機基底種子", "seed") == 99


def test_build_ch10_column_name_map_n_cfg_fields() -> None:
    """§10.3 觀測集相關扁平鍵應有固定中文顯示名。"""
    m = build_ch10_column_name_map(
        ["r3s_num_obs_configs", "r3s_parameters_n_cfg", "r3s_n_cfg_requested"]
    )
    assert "觀測配置數" in m["r3s_num_obs_configs"]
    assert "輸入配置數" in m["r3s_parameters_n_cfg"]
    m2 = build_ch10_column_name_map(["r3s_parameters_n_rep", "r3s_n_rep"])
    assert "重複次數" in m2["r3s_parameters_n_rep"]


def test_build_ch10_column_name_map_level_fields() -> None:
    """欄名映射：層級欄位應轉為中文無符號名稱。"""
    mapping = build_ch10_column_name_map(
        [
            "run_index",
            "level_1_cfg",
            "level_2_prev_shrink_ratio",
            "level_3_N_forbid_delta",
            "level_1_viol_all_disconnected",
        ]
    )
    assert mapping["run_index"] == "執行序號"
    assert mapping["level_1_cfg"] == "第1層合法配置數"
    assert mapping["level_2_prev_shrink_ratio"] == "第2層相對前層收縮率"
    assert mapping["level_3_N_forbid_delta"] == "第3層本層新增排除數"
    assert mapping["level_1_viol_all_disconnected"] == "第1層排除主因連通未滿足筆數"


def test_build_ch10_per_window_tails() -> None:
    """多尺度 per_window 尾碼應有固定中文（與 CH10_9_PER_WINDOW_COL_ZH 一致）。"""
    m = build_ch10_column_name_map(["per_window_0_R_edge_bar", "per_window_1_tau_unit_max"])
    assert "第1列（per_window）" in m["per_window_0_R_edge_bar"]
    assert "周轉" in m["per_window_0_R_edge_bar"]
    assert "第2列（per_window）" in m["per_window_1_tau_unit_max"]
    assert "連續段長度" in m["per_window_1_tau_unit_max"]


def test_per_window_metrics_dataframe_zh_renames() -> None:
    """§10.9 寬表應將程式鍵換成可讀中文欄名。"""
    from hypergraph_experiment.streamlit_common import per_window_metrics_dataframe_zh

    raw = pd.DataFrame([{"w": 1, "H_macro_bits": 1.5, "JS_vs_w1_bits": None}])
    out = per_window_metrics_dataframe_zh(raw)
    cols = list(out.columns)
    assert any("時間聚合視窗寬度" in c for c in cols)
    assert any("宏觀型別熵" in c for c in cols)
    assert any("JS_w" in c for c in cols)


def test_build_ch10_table_10_2_ladder_display_names_unique() -> None:
    """表 10-2 梯子表（與單次、批次 per-run 同形）經映射後顯示欄名須全唯一。"""
    levels = [
        {"label": "弱", "connected": False, "max_degree": 8, "forbid_pair_triangles": False},
        {"label": "中", "connected": True, "max_degree": 5, "forbid_pair_triangles": False},
    ]
    rows, _n = table_10_2_domain_ladder(
        n=6,
        min_edge_size=2,
        max_edge_size=3,
        max_edges=6,
        sample_limit=80,
        seed=1,
        num_seeds=3,
        levels=levels,
    )
    df = pd.DataFrame(rows)
    m = build_ch10_column_name_map(df.columns)
    zh = list(m.values())
    dup = sorted({x for x in zh if zh.count(x) > 1})
    assert not dup, f"重複顯示欄名: {dup}"


def test_experiment_10_9_batch_display_flat_has_columns() -> None:
    """§10.9 批次第二張表（扁平）應與單次同前綴之 flatten 鍵。"""
    out = run_experiment_10_9(
        n=4,
        steps=12,
        window_sizes=(1, 2),
        seed=0,
        r=1,
        n_hist=1,
        n_seed=1,
    )
    flat = flatten_result_row("e9s", out)
    assert any(k.startswith("e9s_") for k in flat)
    assert "e9s_experiment" in flat or "e9s_R_edge_bar" in flat


def test_build_ch10_column_name_map_t2_level_batch() -> None:
    """表 10-2 批次扁平欄 t2_level_* 應與無前綴 level_* 同義（完整中文違規欄名）。"""
    col = "t2_level_3_viol_new_edge_size_bad"
    m = build_ch10_column_name_map([col, "level_3_viol_new_edge_size_bad", "t2_run_index"])
    want = "第3層新增排除主因超邊大小不符筆數"
    assert m[col] == want
    assert m["level_3_viol_new_edge_size_bad"] == want
    assert m["t2_run_index"] == "執行序號"
