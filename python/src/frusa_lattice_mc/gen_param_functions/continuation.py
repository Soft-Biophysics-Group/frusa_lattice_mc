"""
Inspect and continue unfinished frusa_mc jobs.

Design: query functions (find_unfinished, get_progress, get_t_range) are
pure inspections that never write.  write_continuation_inputs is the
only function that touches the filesystem.
"""


import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .params import write_json, load_json
from .paths import c_path


# ── Data returned by queries ────────────────────────────────────────


@dataclass
class UnfinishedJob:
    """All the info needed to decide whether / how to continue a job."""

    mc_file: Path
    mc_file_relative: Path  # relative to input_root
    last_completed_index: int
    total_steps: int

    @property
    def fraction_done(self) -> float:
        return self.last_completed_index / max(self.total_steps - 1, 1)


# ── Query functions (read-only) ─────────────────────────────────────

_STRUCT_RE = re.compile(r"structure_(\d+)\.dat")


def _checkpoint_progress(checkpoint_dir: Path) -> int:
    """Return the highest structure index found, or -1 if none exist."""
    max_idx = -1
    if not checkpoint_dir.is_dir():
        return max_idx
    for p in checkpoint_dir.iterdir():
        m = _STRUCT_RE.match(p.name)
        if m:
            max_idx = max(max_idx, int(m.group(1)))
    return max_idx


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
            progress = _checkpoint_progress(Path(mc["checkpoint_address"]))
            jobs.append(
                UnfinishedJob(
                    mc_file=mc_file,
                    mc_file_relative=mc_file.relative_to(input_root),
                    last_completed_index=progress,
                    total_steps=mc["Nt"],
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
            f"    last checkpoint: {job.last_completed_index} / {job.total_steps - 1}"
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

        # Save continued structures in a sibling directory to avoid overwriting
        old_checkpoint = Path(mc["checkpoint_address"])
        continued_structs = old_checkpoint.parent / "structures_continued"
        continued_structs.mkdir(parents=True, exist_ok=True)
        mc["checkpoint_address"] = c_path(continued_structs)
        mc["final_structure_address"] = c_path(continued_structs)

        # Point initialization at the last checkpoint
        old_model_file = Path(mc["model_params_file"])
        model = load_json(old_model_file)
        model["initialize_option"] = "from_file"
        model["state_input"] = str(
            old_checkpoint / f"structure_{job.last_completed_index}.dat"
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
