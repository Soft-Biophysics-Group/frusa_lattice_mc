"""
Inspect and continue unfinished frusa_mc jobs.

Design: query functions (find_unfinished, checkpoint_progress, get_t_schedule)
are pure inspections that never write.  write_continuation_inputs is the only
function that touches the filesystem; continue_stage and prepare_continuation
wrap it for whole stages and manifests, and are what the local and SLURM
runners share so both produce identical inputs.

A continuation writes into the run's existing output directories: checkpoints
carry a structure_index_offset so the series stays a single globally-indexed
set of files, and energies were already keyed by temperature. One run, one
dataset, however many times it was interrupted.
"""


import re
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .manifest import Stage, read_manifest, write_stages
from .params import write_json, load_json


# ── Data returned by queries ────────────────────────────────────────


@dataclass
class UnfinishedJob:
    """All the info needed to decide whether / how to continue a job."""

    mc_file: Path
    mc_file_relative: Path  # relative to input_root
    last_completed_index: int  # into this file's schedule, not the whole run
    total_steps: int
    structure_index_offset: int = 0  # steps this file's schedule starts after

    @property
    def global_last_completed(self) -> int:
        """Last completed step counted over the whole run, continuations included."""
        return self.structure_index_offset + self.last_completed_index

    @property
    def global_total_steps(self) -> int:
        return self.structure_index_offset + self.total_steps

    @property
    def fraction_done(self) -> float:
        return self.global_last_completed / max(self.global_total_steps - 1, 1)


# ── Query functions (read-only) ─────────────────────────────────────

_STRUCT_RE = re.compile(r"structure_(\d+)\.dat")


def checkpoint_progress(checkpoint_dir: Path, offset: int = 0) -> int:
    """
    Highest structure index found, counted from `offset`, or -1 if none reach it.

    Files are named by their global index; `offset` is the run's
    structure_index_offset, so the result indexes that run's own schedule.
    """
    max_idx = -1
    if not checkpoint_dir.is_dir():
        return max_idx
    for p in checkpoint_dir.iterdir():
        m = _STRUCT_RE.match(p.name)
        if m:
            max_idx = max(max_idx, int(m.group(1)))
    return max_idx - offset if max_idx >= offset else -1


def get_t_schedule(mc_params: dict) -> tuple[np.ndarray, np.ndarray]:
    """
    Reconstruct the temperature profile from MC params.

    Returns (physical_temperatures, schedule_values) where schedule_values
    are what the C++ code steps through (e.g. inverse temperatures).
    """
    schedule = mc_params["cooling_schedule"]

    if schedule == "arbitrary":
        # The temperatures are given explicitly, so they are also what the
        # C++ code steps through.
        t_arr = np.asarray(mc_params["T_array"], dtype=float)
        return t_arr, t_arr

    t_i = mc_params["Ti"]
    t_f = mc_params["Tf"]
    n_t = mc_params["Nt"]

    t_sched = np.linspace(t_i, t_f, n_t)

    if schedule == "inverse":
        return 1.0 / t_sched, t_sched
    elif schedule == "exponential":
        return 10.0**t_sched, t_sched
    elif schedule == "linear":
        return t_sched, t_sched
    else:
        raise ValueError(f"Unknown cooling schedule: {schedule!r}")

def is_finished(mc_file: Path) -> bool:
    mc = load_json(mc_file)
    return (Path(mc["final_structure_address"]) / "final_structure.dat").is_file()

def find_unfinished(
    input_root: Path,
    mc_glob: str = "mc_params_*.json",
) -> list[UnfinishedJob]:
    """
    Scan input_root for MC param files whose final_structure.dat is missing.

    Returns structured UnfinishedJob objects — no side effects.
    """
    input_root = Path(input_root)
    jobs: list[UnfinishedJob] = []

    for mc_file in sorted(input_root.glob(mc_glob)):
        mc = load_json(mc_file)
        final = Path(mc["final_structure_address"]) / "final_structure.dat"

        if not final.is_file():
            offset = mc.get("structure_index_offset", 0)
            jobs.append(
                UnfinishedJob(
                    mc_file=mc_file,
                    mc_file_relative=mc_file.relative_to(input_root),
                    last_completed_index=checkpoint_progress(
                        Path(mc["checkpoint_address"]), offset
                    ),
                    total_steps=mc["Nt"],
                    structure_index_offset=offset,
                )
            )

    return jobs


def print_unfinished_summary(jobs: list[UnfinishedJob]) -> None:
    """Pretty-print a list of unfinished jobs. Replaces the old do_nothing=True path."""
    if not jobs:
        print("All jobs are finished.")
        return

    print(f"Found {len(jobs)} unfinished job(s):\n")
    for job in jobs:
        print(
            f"  {job.mc_file}\n"
            f"    last checkpoint: {job.global_last_completed}"
            f" / {job.global_total_steps - 1}"
            f"  ({job.fraction_done:.0%})\n"
        )


# ── Action function (writes files) ─────────────────────────────────


def write_continuation_inputs(
    jobs: list[UnfinishedJob],
    old_input_root: Path,
    new_input_root: Path,
) -> list[Path]:
    """
    For each unfinished job, write new MC and model param files that
    resume from the last checkpoint.

    Returns the list of newly created MC param file paths.
    """
    old_input_root = Path(old_input_root).resolve()
    new_input_root = Path(new_input_root).resolve()
    created: list[Path] = []

    for job in jobs:
        # Skip any files that are already inside the continuation tree
        if job.mc_file.resolve().is_relative_to(new_input_root):
            continue

        mc = load_json(job.mc_file)
        _, t_sched = get_t_schedule(mc)

        resume_index = job.last_completed_index + 1
        if mc["cooling_schedule"] == "arbitrary":
            # Shrinking Nt is not enough when the temperatures are explicit:
            # drop the steps already done and re-derive Ti/Tf/Nt from what is left.
            remaining = [float(t) for t in t_sched[resume_index:]]
            mc["T_array"] = remaining
            mc["Nt"] = len(remaining)
            mc["Ti"] = remaining[0]
            mc["Tf"] = remaining[-1]
        else:
            mc["Ti"] = float(t_sched[resume_index])
            mc["Nt"] = mc["Nt"] - resume_index

        # Structures keep accumulating in the run's one checkpoint directory: the
        # continuation numbers its own from where the interrupted run stopped, so
        # the series stays a single, gapless, globally-indexed dataset. Only
        # never-completed indices are ever written.
        offset = mc.get("structure_index_offset", 0)
        last_global_index = offset + job.last_completed_index
        mc["structure_index_offset"] = offset + resume_index

        # Point initialization at the last checkpoint
        old_model_file = Path(mc["model_params_file"])
        model = load_json(old_model_file)
        model["initialize_option"] = "from_file"
        model["state_input"] = str(
            Path(mc["checkpoint_address"]) / f"structure_{last_global_index}.dat"
        )

        # Write new files mirroring the relative structure under new_input_root
        new_mc_file = new_input_root / job.mc_file_relative
        model_rel = old_model_file.relative_to(old_input_root)
        new_model_file = new_input_root / model_rel

        mc["model_params_file"] = str(new_model_file.resolve())

        write_json(model, new_model_file)
        write_json(mc, new_mc_file)
        created.append(new_mc_file)

    return created


# ── Manifest-level continuation (shared by the local and SLURM runners) ──


def new_attempt_dir(continue_root: Path, label: str | None = None) -> Path:
    """
    Create and return a fresh directory for one resume attempt's param files.

    Named by timestamp, plus `label` to keep concurrent array tasks apart, and
    suffixed on collision. Attempts must never share a directory, since each one
    reads the previous attempt's inputs to work out where to restart — so the
    directory is claimed with mkdir itself, which is atomic, rather than by
    testing for existence first.
    """
    continue_root = Path(continue_root)
    stamp = time.strftime("%y%m%d_%H%M%S")
    if label:
        stamp = f"{stamp}_{label}"
    attempt, n = continue_root / stamp, 1
    while True:
        try:
            attempt.mkdir(parents=True)
            return attempt
        except FileExistsError:
            attempt = continue_root / f"{stamp}_{n:02}"  # zero-padded so names sort
            n += 1


def continue_stage(stage: Stage, attempt_dir: Path) -> Stage | None:
    """Continuation param files for one stage, or None if it cannot be resumed."""
    if is_finished(stage.mc_file):
        return None

    mc = load_json(stage.mc_file)
    offset = mc.get("structure_index_offset", 0)
    progress = checkpoint_progress(Path(mc["checkpoint_address"]), offset)
    if progress < 0:
        return None

    # input/<prefix>/<slug>/mc_params_i.json — mirror from the input root down, so
    # two stages of one job (same slug, different prefix) never share a temp path.
    input_root = stage.mc_file.parents[2]
    job = UnfinishedJob(
        mc_file=stage.mc_file,
        mc_file_relative=stage.mc_file.relative_to(input_root),
        last_completed_index=progress,
        total_steps=mc["Nt"],
        structure_index_offset=offset,
    )
    created = write_continuation_inputs([job], input_root, attempt_dir)
    if not created:
        return None
    return Stage(Path(load_json(created[0])["model_params_file"]), created[0])


def prepare_continuation(manifest_path: Path, attempt_dir: Path) -> Path | None:
    """
    Write continuation inputs for every resumable stage of a manifest.

    Stages that already finished, and those that never wrote a checkpoint, keep
    their original param files. Returns a new manifest pointing at the mix, or
    None if nothing needs continuing.
    """
    manifest_path = Path(manifest_path)
    attempt_dir = Path(attempt_dir)

    new_jobs: list[list[Stage]] = []
    n_continued = 0
    for stages in read_manifest(manifest_path):
        new_stages = []
        for stage in stages:
            continued = continue_stage(stage, attempt_dir)
            n_continued += continued is not None
            new_stages.append(continued if continued is not None else stage)
        new_jobs.append(new_stages)

    if not n_continued:
        return None

    new_manifest = attempt_dir / manifest_path.name
    write_stages(new_manifest, new_jobs)
    return new_manifest
