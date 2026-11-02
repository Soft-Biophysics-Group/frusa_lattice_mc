"""Gzip the structure checkpoints of a campaign's runs, one file at a time.

Same shape as `frusa-mc-task`: line N of the manifest is array task N, so it can
follow the simulation in the same SLURM task. Without an index, every line is done.
Run from the project root, where the paths inside the parameter files resolve.

    frusa-mc-compress input/04_Q2_M90_very_short_times/model_mc_files.txt "$SLURM_ARRAY_TASK_ID"
    frusa-mc-compress input/04_Q2_M90_very_short_times/model_mc_files.txt

Safe on runs in any state, and to repeat: an unfinished run keeps its newest
checkpoint plain so it can still be resumed, and a file any stage of the line
starts from is never touched. See `frusa_lattice_mc.checkpoints`.
"""

import argparse
from pathlib import Path

from ..checkpoints import compress_checkpoints, structure_path
from ..gen_param_functions.continuation import checkpoint_progress, is_finished
from ..gen_param_functions.manifest import Stage, read_manifest
from ..json_dump import load_json


def checkpoint_dir(stage: Stage) -> Path | None:
    address = load_json(stage.mc_file).get("checkpoint_address")
    return None if address is None else Path(address)


def files_to_keep(stages: list[Stage]) -> set[Path]:
    """Structures the C++ may still read: resume points and state inputs."""
    keep = set()
    for stage in stages:
        model = load_json(stage.model_file)
        if model.get("initialize_option") == "from_file":
            keep.add(Path(model["state_input"]))
        directory = checkpoint_dir(stage)
        if directory is not None and not is_finished(stage.mc_file):
            newest = checkpoint_progress(directory)
            if newest >= 0:
                keep.add(structure_path(directory, newest))
    return keep


def compress_line(stages: list[Stage], dry_run: bool) -> int:
    """Compress every checkpoint directory of one manifest line. Returns files done."""
    keep = files_to_keep(stages)
    directories = {d for d in map(checkpoint_dir, stages) if d is not None}
    n_done = 0
    for directory in sorted(directories):
        done = compress_checkpoints(directory, keep, dry_run)
        verb = "would compress" if dry_run else "compressed"
        print(f"{verb} {len(done)} structure(s) in {directory}")
        n_done += len(done)
    return n_done


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="Manifest naming the runs.")
    parser.add_argument(
        "index", type=int, nargs="?",
        help="1-based manifest line, e.g. $SLURM_ARRAY_TASK_ID. All lines if omitted.",
    )
    parser.add_argument(
        "-n", "--dry-run", action="store_true",
        help="Say what would be compressed without touching anything.",
    )
    args = parser.parse_args()

    jobs = read_manifest(args.manifest)
    if args.index is None:
        lines = jobs
    elif 1 <= args.index <= len(jobs):
        lines = [jobs[args.index - 1]]
    else:
        parser.error(f"index {args.index} outside 1..{len(jobs)}")

    for stages in lines:
        compress_line(stages, args.dry_run)


if __name__ == "__main__":
    main()
