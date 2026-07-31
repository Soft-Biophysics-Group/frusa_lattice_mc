import json
from pathlib import Path
from typing import NamedTuple

class Stage(NamedTuple):
    model_file: Path
    mc_file: Path


def write_stages(path: Path, jobs: list[list[Stage]]) -> None:
    """Write a series of jobs, each made of multiple stages (individual simulations).
    Scope: writing run files for multi-stage simulations"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for stages in jobs:
            f.write(
                " ".join(f"{model.resolve()} {mc.resolve()}" for model, mc in stages)
                + "\n"
            )

def write_single_stage(path: Path, jobs: list[Stage]) -> None:
    """Write a series of single-stage jobs, "mc model". Used for simple simulation
    campaings and data analysis."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for model_file, mc_file in jobs:
            f.write(f"{mc_file.resolve()} {model_file.resolve()}\n")

def read_manifest(path:Path) -> list[list[Stage]]:
    """Read the manifest located at path, whether it is in several or single stage format"""
    jobs = []
    for line in path.read_text().splitlines():
        fields = [Path(f) for f in line.split()]
        if not fields:
            continue
        if len(fields) % 2:
            raise ValueError(f"Odd column count in {path}: {line!r}")
        jobs.append(
            [_identify_pair(fields[i], fields[i + 1]) for i in range(0, len(fields), 2)]
        )
    return jobs


def _identify_pair(first: Path, second: Path) -> Stage:
    """Sort a column pair into (model, mc) by what the JSON holds. Agnostic to the order in
    which files were written."""
    is_mc = ["mcs_eq" in json.loads(path.read_text()) for path in (first, second)]
    match is_mc:
        case [False, True]:
            return Stage(first, second)
        case [True, False]:
            return Stage(second, first)
        case _:
            kind = "MC" if is_mc[0] else "model"
            raise ValueError(f"{first} and {second} are both {kind} parameter files")
