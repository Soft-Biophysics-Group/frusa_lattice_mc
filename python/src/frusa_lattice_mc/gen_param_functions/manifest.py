"""The manifest: one line per SLURM array task, "model mc" per stage in order.

Older manifests put mc first. read_manifest identifies columns by their JSON
contents rather than by position, so those files still read correctly.
"""

import json
import os
from pathlib import Path
from typing import NamedTuple, Sequence


class Stage(NamedTuple):
    model_file: Path
    mc_file: Path


def run_dir(stage: Stage) -> Path:
    """The directory this stage writes into: the parent of its structures/ folder."""
    mc = json.loads(stage.mc_file.read_text())
    address = mc.get("final_structure_address")
    if address is None:
        # MCParams.to_dict drops None fields, so a run configured without a
        # final structure has no such key at all. Name the file that is missing
        # it: a bare KeyError says nothing about which of hundreds it was.
        raise ValueError(f"{stage.mc_file} has no final_structure_address")
    return Path(address).parent

def write_manifest(path: Path, jobs: Sequence[Sequence[Stage]]) -> None:
    """Write one line per job, a job being the stages of one array task.

    The inverse of read_manifest. A run with nothing chained after it is a
    one-stage job: [[stage], [stage], ...].
    """
    n_stages = len(jobs[0])
    for index, stages in enumerate(jobs):
        # generate_array_script_by_stages sizes the script from the first line,
        # so a ragged manifest silently drops the extra stages of every other job.
        if len(stages) != n_stages:
            raise ValueError(
                f"Job {index} of {path} has {len(stages)} stages, job 0 has {n_stages}"
            )

    text = "".join(
        " ".join(
            f"{stage.model_file.resolve()} {stage.mc_file.resolve()}" for stage in stages
        )
        + "\n"
        for stages in jobs
    )

    # Written whole, then moved into place: a failure here leaves the previous
    # manifest intact rather than truncating it.
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


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
