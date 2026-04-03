"""
第十章超圖實驗之可程式化核心：域型約束篩選、解析簽名、等價類與動力學示範。

本套件供 CLI（experiment.py）與 Streamlit GUI 共用。
"""

from __future__ import annotations

from hypergraph_experiment.core import (
    HypergraphConfig,
    SIGNATURES,
    analyze_dynamics,
    analyze_static,
    filter_configs,
    generate_candidate_configs,
    run_full_experiment,
    satisfies_domain_constraints,
)
from hypergraph_experiment.refinement import (
    analyze_refinement_pair,
    analyze_section_10_4_bundle,
    apply_refinement_step,
    partition_analytic_units,
)

__all__ = [
    "HypergraphConfig",
    "SIGNATURES",
    "analyze_dynamics",
    "analyze_static",
    "analyze_refinement_pair",
    "analyze_section_10_4_bundle",
    "apply_refinement_step",
    "filter_configs",
    "generate_candidate_configs",
    "partition_analytic_units",
    "run_full_experiment",
    "satisfies_domain_constraints",
]
