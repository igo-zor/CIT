"""
論文 §10.5–§10.9 之獨立實驗入口（與 core 之通用候選／動力學互補）。
"""

from __future__ import annotations

from hypergraph_experiment.experiments.exp_10_5_bipartite import run_experiment_10_5
from hypergraph_experiment.experiments.exp_10_6_contexts import (
    run_canonical_demo_10_6,
    run_experiment_10_6,
)
from hypergraph_experiment.experiments.exp_10_7_paths import run_experiment_10_7
from hypergraph_experiment.experiments.exp_10_8_symmetry import (
    rho_type_expansion,
    run_experiment_10_8,
    run_experiment_10_8_three_arm,
)
from hypergraph_experiment.experiments.exp_10_9_multiscale import (
    run_experiment_10_9,
    section_10_9_output_parameters_df,
)

__all__ = [
    "run_experiment_10_5",
    "run_experiment_10_6",
    "run_canonical_demo_10_6",
    "run_experiment_10_7",
    "run_experiment_10_8",
    "run_experiment_10_8_three_arm",
    "run_experiment_10_9",
    "section_10_9_output_parameters_df",
    "rho_type_expansion",
]
