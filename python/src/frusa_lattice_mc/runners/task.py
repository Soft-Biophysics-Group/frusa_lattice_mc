"""Vincent Ouazan-Reboul, 2026-07-31. Co-written by Claude Opus 5.

Run one line of a manifest, continuing it if it was interrupted.

One array task per manifest line. Each stage is skipped if it already wrote its
final structure, continued from its last checkpoint if it has one, and run from
scratch otherwise — so the same SLURM script starts a campaign and picks it up
after a timeout, with nothing to edit in between. Resubmit it as often as you
like; once everything is finished it becomes a no-op.

    frusa-mc-task input/03_long_simus_T_1/manifest.txt "$SLURM_ARRAY_TASK_ID"
"""

import argparse
from pathlib import Path

from ..config import find_executable
from ..gen_param_functions import manifest
from ..gen_param_functions.continuation import new_attempt_dir
from .stages import (
    CONTINUE_DIR_NAME,
    discard_if_empty,
    missing_initial_state,
    run_stages,
    task_name,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument(
        "index", type=int, help="1-based manifest line, e.g. $SLURM_ARRAY_TASK_ID"
    )
    parser.add_argument(
        "--no-resume", action="store_true",
        help="Run interrupted stages from the start instead of continuing them.",
    )
    parser.add_argument(
        "--force", action="store_true", help="Redo stages that already completed."
    )
    parser.add_argument("--log-dir", type=Path, default=None)
    parser.add_argument("--continue-root", type=Path, default=None)
    parser.add_argument("--executable", type=Path, default=None)
    args = parser.parse_args()

    jobs = manifest.read_manifest(args.manifest)
    if not 1 <= args.index <= len(jobs):
        raise SystemExit(
            f"Line {args.index} is out of range: {args.manifest} has {len(jobs)} line(s)."
        )
    stages = jobs[args.index - 1]

    missing = missing_initial_state(stages)
    if missing is not None:
        raise SystemExit(
            f"This run starts from a structure that does not exist: {missing}\n"
            "Its preceding stage has to run first."
        )

    log_dir = args.log_dir or args.manifest.parent / "logs_stages"
    log_dir.mkdir(parents=True, exist_ok=True)
    continue_root = args.continue_root or args.manifest.parent / CONTINUE_DIR_NAME
    attempt_dir = (
        None
        if args.no_resume
        else new_attempt_dir(continue_root, label=f"task{args.index:03}")
    )

    result = run_stages(
        stages,
        log_dir,
        attempt_dir,
        continue_root,
        args.force,
        find_executable(args.executable),
    )
    discard_if_empty(attempt_dir)

    print(
        f"{task_name(stages)}: {result.n_ran} run,"
        f" {result.n_skipped} skipped, {result.n_resumed} continued"
    )
    if result.returncode != 0:
        print(f"FAILED at stage {result.failed_stage} (rc={result.returncode})")
        if result.note:
            print(result.note)
    print(f"Log: {log_dir / f'{task_name(stages)}.log'}")
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
