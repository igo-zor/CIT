#!/usr/bin/env python3
"""
第十章超圖 toy model：命令列入口（計算邏輯見 hypergraph_experiment.core）。

Usage examples:
  python experiment.py --n 5 --max-edge-size 3 --max-edges 4 --delta 0 --signature weak
  python experiment.py --n 6 --sample-limit 5000 --connected --runs 50 --steps 40
  python experiment.py --mode dynamics --n 5 --steps 30 --runs 100 --signature medium
"""

from __future__ import annotations

import argparse
import json
import random

from hypergraph_experiment.core import SIGNATURES, run_full_experiment


def build_parser() -> argparse.ArgumentParser:
    """建構與先前版本相同參數之命令列解析器。"""
    p = argparse.ArgumentParser(description="Minimal hypergraph experiment for Chapter 10")
    p.add_argument("--mode", choices=["static", "dynamics"], default="static")
    p.add_argument("--n", type=int, default=5, help="Number of vertices")
    p.add_argument("--max-edge-size", type=int, default=3)
    p.add_argument("--max-edges", type=int, default=4)
    p.add_argument("--max-degree", type=int, default=4)
    p.add_argument("--connected", action="store_true", help="Require connected 2-section graph")
    p.add_argument("--forbid-pair-triangles", action="store_true", help="Forbid 2-edge triangles")
    p.add_argument(
        "--sample-limit",
        type=int,
        default=2000,
        help="Random sample size; use 0 for exhaustive enumeration when feasible",
    )
    p.add_argument(
        "--n-cfg",
        type=int,
        default=None,
        help="§10.3：自可採用域抽取之觀測配置數 N_cfg；省略則以全部可採用配置分析",
    )
    p.add_argument(
        "--n-rep",
        type=int,
        default=1,
        help="§10.3：觀測集重抽次數 N_rep（僅 static；1 表示不重抽）",
    )
    p.add_argument("--signature", choices=sorted(SIGNATURES.keys()), default="medium")
    p.add_argument("--delta", type=int, default=0)
    p.add_argument(
        "--s-min",
        type=int,
        default=0,
        help="§10.3 重疊率：僅 |T(c)|>=s_min 之配置參與平均；0 表示不過濾",
    )
    p.add_argument(
        "--epsilon-plat",
        type=float,
        default=0.02,
        help="§10.7 熵序列平台閾值 |ΔH|≤ε",
    )
    p.add_argument("--runs", type=int, default=50)
    p.add_argument("--steps", type=int, default=30)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--show-sample-configs", type=int, default=3)
    p.add_argument(
        "--refinement",
        action="store_true",
        help="Append §10.4 refinement: π, fibers, K_s, pushforward check (static mode only)",
    )
    p.add_argument("--refine-coarse-sig", choices=sorted(SIGNATURES.keys()), default="weak")
    p.add_argument("--refine-coarse-delta", type=int, default=3)
    p.add_argument("--refine-fine-sig", choices=sorted(SIGNATURES.keys()), default="medium")
    p.add_argument("--refine-fine-delta", type=int, default=0)
    p.add_argument("--refine-kernel", choices=["uniform", "proportional"], default="uniform")
    p.add_argument(
        "--no-refine-compare-chains",
        action="store_true",
        help="Skip §10.4.1 two-path JS / entropy comparison",
    )
    return p


def main() -> None:
    args = build_parser().parse_args()
    random.seed(args.seed)
    result = run_full_experiment(
        mode=args.mode,
        n=args.n,
        max_edge_size=args.max_edge_size,
        max_edges=args.max_edges,
        max_degree=args.max_degree,
        connected=args.connected,
        forbid_pair_triangles=args.forbid_pair_triangles,
        sample_limit=args.sample_limit,
        n_cfg=args.n_cfg,
        n_rep=args.n_rep,
        signature=args.signature,
        delta=args.delta,
        s_min=args.s_min,
        epsilon_plat=args.epsilon_plat,
        runs=args.runs,
        steps=args.steps,
        seed=args.seed,
        show_sample_configs=args.show_sample_configs,
        refinement_enabled=args.refinement,
        refine_coarse_signature=args.refine_coarse_sig,
        refine_coarse_delta=args.refine_coarse_delta,
        refine_fine_signature=args.refine_fine_sig,
        refine_fine_delta=args.refine_fine_delta,
        refine_kernel=args.refine_kernel,
        refine_compare_chains=not args.no_refine_compare_chains,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
