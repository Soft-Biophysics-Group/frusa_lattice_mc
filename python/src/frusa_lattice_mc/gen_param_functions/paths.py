"""
Path construction for frusa_mc simulation runs.

Single source of truth for the directory layout:

    root/
    ├── input/{prefix}/{slug}/model_params_{i}.json
    │                        /mc_params_{i}.json
    └── data/{prefix}/{slug}/run_{i}/
                                    ├── average_e/
                                    ├── e_records/
                                    └── structures/
"""

from pathlib import Path


def c_path(p: Path) -> str:
    """Format a Path as the C++ side expects: resolved, with a trailing slash."""
    return str(p.resolve()) + "/"


def run_slug(series_index: int, label: str | None = None) -> str:
    if label is not None:
        return f"series_{series_index:02}_{label}"
    return f"series_{series_index:02}"


# ── directory builders ──────────────────────────────────────────────


def input_dir(root: Path, prefix: str, slug: str) -> Path:
    return root / "input" / prefix / slug


def input_file_paths(
    root: Path, prefix: str, slug: str, run_index: int
) -> tuple[Path, Path]:
    d = input_dir(root, prefix, slug)
    return d / f"model_params_{run_index}.json", d / f"mc_params_{run_index}.json"


def data_dir(root: Path, prefix: str, slug: str, run_index: int) -> Path:
    return (root / "data" / prefix / slug / f"run_{run_index}").resolve()


def structures_dir(data: Path) -> Path:
    return data / "structures"


def energy_av_dir(data: Path) -> Path:
    return data / "average_e"


def energy_records_dir(data: Path) -> Path:
    return data / "e_records"
