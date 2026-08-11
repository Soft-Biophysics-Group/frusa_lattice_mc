"""Vincent Ouazan-Reboul, 2026-07-31. Co-written by Claude Opus 5.

Running the stages of one manifest line, wherever the job happens to run.

Shared by the local runner (every line, in a process pool) and the SLURM task
runner (one line per array task). Each stage decides for itself whether to be
skipped, continued or run from scratch, which is what lets a single script both
start a campaign and pick it up after an interruption.
"""

import json
import subprocess
from pathlib import Path
from typing import NamedTuple

from ..gen_param_functions import continuation
from ..gen_param_functions.manifest import Stage

CONTINUE_DIR_NAME = "_continued"


class TaskResult(NamedTuple):
    name: str
    n_ran: int
    n_skipped: int
    n_resumed: int
    failed_stage: int | None
    returncode: int


def task_name(stages: list[Stage]) -> str:
    mc_file = stages[0].mc_file
    return f"{mc_file.parent.name}_{mc_file.stem.replace('mc_params_', 'run_')}"


def missing_initial_state(stages: list[Stage]) -> Path | None:
    """The state the first stage resumes from, if named but absent."""
    model = json.loads(stages[0].model_file.read_text())
    if model.get("initialize_option") != "from_file":
        return None
    state = Path(model["state_input"])
    return None if state.is_file() else state


def independent_stages(stages: list[Stage]) -> bool:
    """
    Whether a manifest line's stages can run in parallel rather than in order.

    A line is a chain only because a later stage initialises from an earlier
    one's final structure. If no stage reads a structure from disk at all, the
    ordering carries no meaning: these are independent runs that were grouped
    onto one line, typically to keep a SLURM array small.

    Deliberately conservative — any stage reading from a file keeps the whole
    line sequential, since it may be reading a sibling's output.
    """
    return not any(
        json.loads(stage.model_file.read_text()).get("initialize_option") == "from_file"
        for stage in stages
    )


def split_independent(jobs: list[list[Stage]]) -> list[list[Stage]]:
    """One job per stage for every line whose stages do not depend on each other."""
    split: list[list[Stage]] = []
    for stages in jobs:
        if len(stages) > 1 and independent_stages(stages):
            split.extend([stage] for stage in stages)
        else:
            split.append(stages)
    return split


def newest_inputs(stage: Stage, attempt_dir: Path, continue_root: Path) -> Stage:
    """This stage's latest param files: newest earlier attempt, else the original."""
    rel = stage.mc_file.relative_to(stage.mc_file.parents[2])
    for previous in sorted(continue_root.glob("*"), reverse=True):
        if previous != attempt_dir and (previous / rel).is_file():
            mc_file = previous / rel
            model_file = json.loads(mc_file.read_text())["model_params_file"]
            return Stage(Path(model_file), mc_file)
    return stage


def resume_stage(stage: Stage, attempt_dir: Path, continue_root: Path) -> Stage | None:
    """Temp param files resuming from the last checkpoint, or None if there is none."""
    source = newest_inputs(stage, attempt_dir, continue_root)
    return continuation.continue_stage(source, attempt_dir)


def discard_if_empty(attempt_dir: Path | None) -> None:
    """Drop an attempt directory nothing was written to, so reruns stay tidy."""
    if attempt_dir is not None and attempt_dir.is_dir() and not any(attempt_dir.iterdir()):
        attempt_dir.rmdir()


def run_stages(
    stages: list[Stage],
    log_dir: Path,
    attempt_dir: Path | None,
    continue_root: Path,
    force: bool,
    executable: Path,
) -> TaskResult:
    """Run every stage of one task in order, stopping at the first failure."""
    name = task_name(stages)
    n_ran = n_skipped = n_resumed = 0

    with (log_dir / f"{name}.log").open("a") as log:
        for stage_index, stage in enumerate(stages, start=1):
            if not force and continuation.is_finished(stage.mc_file):
                n_skipped += 1
                continue

            to_run = stage
            if attempt_dir is not None and not force:
                resumed = resume_stage(stage, attempt_dir, continue_root)
                if resumed is not None:
                    to_run = resumed
                    n_resumed += 1

            log.write(f"=== stage {stage_index}/{len(stages)}: {to_run.mc_file}\n")
            log.flush()
            result = subprocess.run(
                [
                    str(executable),
                    "-m", str(to_run.model_file),
                    "-M", str(to_run.mc_file),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
            )
            if result.returncode != 0:
                return TaskResult(
                    name, n_ran, n_skipped, n_resumed, stage_index, result.returncode
                )

            n_ran += 1

    return TaskResult(name, n_ran, n_skipped, n_resumed, None, 0)
