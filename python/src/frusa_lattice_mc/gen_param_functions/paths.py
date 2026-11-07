"""
Path construction for frusa_mc simulation runs.

Single source of truth for the directory layout:

    root/
    ├── input/{series}/model_params_{run}[_stage_{stage}].json
    │                 /mc_params_{run}[_stage_{stage}].json
    └── data/{series}/run_{run}[_stage_{stage}]/
                                ├── average_e/
                                ├── e_records/
                                └── structures/
"""

from pathlib import Path


def c_path(p: Path) -> str:
    """Format a Path as the C++ side expects: resolved, with a trailing slash."""
    return str(p.resolve()) + "/"



# ── directory builders ──────────────────────────────────────────────


def input_dir(root: Path, run_name: str) -> Path:
    return root / "input" / run_name


def _stage_name(stage_index: int | None, n_stages: int | None) -> str:
    """Return a zero-padded stage suffix, or nothing for a one-stage run."""
    if stage_index is None and n_stages is None:
        return ""
    if stage_index is None or n_stages is None:
        raise ValueError("stage_index and n_stages must be provided together")
    if n_stages < 1:
        raise ValueError("n_stages must be positive")
    if not 0 <= stage_index < n_stages:
        raise ValueError(
            f"stage_index must be in [0, {n_stages}), got {stage_index}"
        )

    width = len(str(n_stages - 1))
    return f"_stage_{stage_index:0{width}d}"


def input_file_paths(
    root: Path,
    series_name: str,
    run_index: int,
    stage_index: int | None = None,
    n_stages: int | None = None,
) -> tuple[Path, Path]:
    d = input_dir(root, series_name)
    stage_name = _stage_name(stage_index, n_stages)
    return (
        d / f"model_params_{run_index}{stage_name}.json",
        d / f"mc_params_{run_index}{stage_name}.json",
    )


def data_dir(
    root: Path,
    series_name: str,
    run_index: int,
    stage_index: int | None = None,
    n_stages: int | None = None,
) -> Path:
    return (
        root
        / "data"
        / series_name
        / f"run_{run_index}{_stage_name(stage_index, n_stages)}"
    ).resolve()


def structures_dir(data: Path) -> Path:
    return data / "structures"


def energy_av_dir(data: Path) -> Path:
    return data / "average_e"


def energy_records_dir(data: Path) -> Path:
    return data / "e_records"
