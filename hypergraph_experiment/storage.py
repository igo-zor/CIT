"""
實驗封存：每筆 run 獨立目錄、manifest.json、選用圖檔與範例超圖 JSON。

倉儲根路徑預設為專案下 ``experiments_data``，可藉環境變數 ``CIT_EXPERIMENTS_ROOT`` 覆寫。
"""

from __future__ import annotations

import json
import os
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from hypergraph_experiment.output_rounding import round_floats_for_output


def default_store_root() -> Path:
    """回傳實驗倉儲根目錄（不存在則建立）。"""
    env = os.environ.get("CIT_EXPERIMENTS_ROOT", "").strip()
    if env:
        root = Path(env).resolve()
    else:
        # 預設與套件同層之專案根
        root = Path(__file__).resolve().parent.parent / "experiments_data"
    root.mkdir(parents=True, exist_ok=True)
    return root


def new_run_id() -> str:
    """產生用於目錄名稱之 run 識別碼（UTC 時間 + 短 UUID）。"""
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    short = uuid.uuid4().hex[:8]
    return f"{ts}_{short}"


@dataclass
class RunRecord:
    """單筆實驗在磁碟上之摘要（供列表與比對 UI 使用）。"""

    run_id: str
    path: Path
    created_at: str
    mode: str
    n: int
    delta: int
    num_admissible: int
    entropy_or_none: float | None
    notes: str


def _safe_entropy_from_analysis(analysis: dict[str, Any]) -> float | None:
    """自 analysis 字典抽出單一可比較之熵欄位（靜態用 entropy_bits，動態用摘要 mean）。"""
    if "error" in analysis:
        return None
    if "entropy_bits" in analysis:
        v = analysis["entropy_bits"]
        return float(v) if v is not None else None
    summ = analysis.get("entropy_summary") or {}
    m = summ.get("mean")
    return float(m) if m is not None else None


def build_manifest(
    result: dict[str, Any],
    *,
    run_id: str,
    notes: str = "",
) -> dict[str, Any]:
    """
    合併實驗結果與封存中繼資料，形成寫入 manifest.json 之完整字典。

    Args:
        result: ``run_full_experiment`` 回傳值。
        run_id: 目錄名稱所使用之識別碼。
        notes: 使用者備註（可為空）。

    Returns:
        可 JSON 序列化之 manifest。
    """
    params = result.get("parameters") or {}
    created = datetime.now(timezone.utc).isoformat()
    manifest: dict[str, Any] = {
        "run_id": run_id,
        "created_at": created,
        "notes": notes,
        "parameters": params,
        "num_candidates": result.get("num_candidates"),
        "num_admissible_configs": result.get("num_admissible_configs"),
        "sample_configs": result.get("sample_configs"),
        "analysis": result.get("analysis"),
    }
    return manifest


def save_run(
    result: dict[str, Any],
    *,
    store_root: Path | None = None,
    notes: str = "",
    extra_files: dict[str, bytes] | None = None,
) -> Path:
    """
    建立新 run 目錄並寫入 manifest、範例超圖與選用二進位附件（如 PNG）。

    Args:
        result: 完整實驗結果。
        store_root: 倉儲根；預設 ``default_store_root()``。
        notes: 備註。
        extra_files: 相對路徑 -> 位元組（例如 ``{\"figures/entropy.png\": data}``）。

    Returns:
        該次 run 之目錄路徑。
    """
    root = store_root or default_store_root()
    rid = new_run_id()
    run_dir = root / rid
    run_dir.mkdir(parents=True, exist_ok=False)

    manifest = round_floats_for_output(build_manifest(result, run_id=rid, notes=notes))
    manifest_path = run_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    sample = result.get("sample_configs") or []
    sample_path = run_dir / "sample_hypergraphs.json"
    sample_path.write_text(json.dumps(sample, ensure_ascii=False, indent=2), encoding="utf-8")

    if extra_files:
        for rel, data in extra_files.items():
            target = run_dir / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)

    return run_dir


def load_manifest(run_dir: Path) -> dict[str, Any]:
    """讀取單一 run 目錄內之 manifest.json。"""
    p = run_dir / "manifest.json"
    if not p.is_file():
        raise FileNotFoundError(f"找不到 manifest：{p}")
    return json.loads(p.read_text(encoding="utf-8"))


def iter_run_directories(store_root: Path | None = None) -> Iterator[Path]:
    """依目錄名排序迭代倉儲內各 run 資料夾。"""
    root = store_root or default_store_root()
    if not root.is_dir():
        return
    for child in sorted(root.iterdir(), key=lambda p: p.name, reverse=True):
        if child.is_dir() and (child / "manifest.json").is_file():
            yield child


def summarize_run(run_dir: Path) -> RunRecord | None:
    """由 manifest 建立 ``RunRecord``；若檔案損毀則回傳 None。"""
    try:
        m = load_manifest(run_dir)
    except (OSError, json.JSONDecodeError):
        return None
    params = m.get("parameters") or {}
    analysis = m.get("analysis") or {}
    return RunRecord(
        run_id=m.get("run_id", run_dir.name),
        path=run_dir,
        created_at=str(m.get("created_at", "")),
        mode=str(params.get("mode", "")),
        n=int(params.get("n", 0)),
        delta=int(params.get("delta", 0)),
        num_admissible=int(m.get("num_admissible_configs") or 0),
        entropy_or_none=_safe_entropy_from_analysis(analysis if isinstance(analysis, dict) else {}),
        notes=str(m.get("notes", "")),
    )


def list_all_runs(store_root: Path | None = None) -> list[RunRecord]:
    """列出倉儲內所有有效 run 之摘要。"""
    out: list[RunRecord] = []
    for d in iter_run_directories(store_root):
        rec = summarize_run(d)
        if rec is not None:
            out.append(rec)
    return out
