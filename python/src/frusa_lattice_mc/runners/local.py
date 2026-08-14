"""Vincent Ouazan-Reboul, 2026-07-30. Co-written by Claude Opus 5.

Run a manifest on this machine instead of on the cluster.

Local counterpart of the SLURM array job: one worker per manifest line, running
that line's stages in order via runners.stages. Completed stages are skipped;
with --resume, an interrupted stage continues from its last checkpoint using
temporary param files under <manifest dir>/_continued/. Those inputs come from
continuation.continue_stage, the same call the cluster path makes, so a run
resumed here and one resumed there produce the same dataset.

    frusa-mc-local input/03_long_simus_T_1/manifest.txt -j 6 --resume
"""

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from ..config import find_executable
from ..gen_param_functions import continuation, manifest
from .stages import (
    CONTINUE_DIR_NAME,
    TaskResult,
    discard_if_empty,
    missing_initial_state,
    run_stages,
    split_independent,
    task_name,
)


def run_manifest(
    manifest_path: Path,
    *,
    workers: int | None = None,
    resume: bool = False,
    force: bool = False,
    split: bool = True,
    log_dir: Path | None = None,
    continue_root: Path | None = None,
    executable: Path | None = None,
    verbose: bool = True,
) -> list[TaskResult]:
    """Run every job of a manifest, `workers` at a time. One result per job."""
    manifest_path = Path(manifest_path)
    workers = workers or max(1, (os.cpu_count() or 2) // 2)
    log_dir = Path(log_dir) if log_dir else manifest_path.parent / "logs_local"
    continue_root = (
        Path(continue_root)
        if continue_root
        else manifest_path.parent / CONTINUE_DIR_NAME
    )
    executable = find_executable(executable)

    lines = manifest.read_manifest(manifest_path)
    if not lines:
        raise ValueError(f"{manifest_path} is empty")

    # Stages on one line run in order, so a line is the unit of parallelism.
    # Lines whose stages do not actually depend on each other get split, or a
    # grouping chosen to keep a SLURM array small would cap the pool here too.
    jobs = split_independent(lines) if split else lines

    orphans = [
        (task_name(stages), missing)
        for stages in jobs
        if (missing := missing_initial_state(stages)) is not None
    ]
    if orphans:
        listing = "\n".join(f"  {name}: {path}" for name, path in orphans[:5])
        raise FileNotFoundError(
            f"{len(orphans)} run(s) resume from a structure that does not exist:\n"
            f"{listing}\nRun the preceding stage first."
        )

    attempt_dir = continuation.new_attempt_dir(continue_root) if resume else None

    log_dir.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    results: list[TaskResult] = []

    if verbose:
        n_stages = max(len(stages) for stages in jobs)
        if len(jobs) != len(lines):
            print(
                f"{len(lines)} manifest line(s) split into {len(jobs)} independent run(s)"
                " — their stages do not read each other's output"
            )
        print(f"{len(jobs)} run(s), up to {n_stages} stage(s) each, {workers} at a time")
        if workers > len(jobs):
            print(
                f"Note: {workers} workers requested but only {len(jobs)} job(s) to run;"
                " parallelism is per manifest line."
            )
        if attempt_dir is not None:
            print(f"Resuming; continuation inputs in {attempt_dir}")

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(
                run_stages,
                stages,
                log_dir,
                attempt_dir,
                continue_root,
                force,
                executable,
            )
            for stages in jobs
        ]
        for n_done, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            if verbose:
                status = (
                    f"ok ({result.n_ran} run, {result.n_skipped} skipped,"
                    f" {result.n_resumed} resumed)"
                    if result.returncode == 0
                    else f"FAILED at stage {result.failed_stage}"
                    f" (rc={result.returncode})"
                )
                if result.note:
                    status = f"{status}\n    {result.note}"
                print(
                    f"[{n_done}/{len(jobs)}] {result.name}: {status}"
                    f"  ({(time.monotonic() - start) / 60:.1f} min elapsed)"
                )

    discard_if_empty(attempt_dir)

    if verbose:
        n_failed = sum(r.returncode != 0 for r in results)
        print(
            f"Done in {(time.monotonic() - start) / 60:.1f} min, {n_failed} failure(s)"
        )
        print(f"Logs under {log_dir}")

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument(
        "-j", "--workers", type=int, default=None,
        help="Concurrent runs. Defaults to half the logical cores.",
    )
    parser.add_argument("--log-dir", type=Path, default=None)
    parser.add_argument("--continue-root", type=Path, default=None)
    parser.add_argument(
        "-n", "--dry-run", action="store_true", help="List what would run, then exit."
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Restart interrupted stages from their last checkpoint.",
    )
    parser.add_argument(
        "--force", action="store_true", help="Redo stages that already completed."
    )
    parser.add_argument(
        "--no-split", action="store_true",
        help="Keep every manifest line as one job, even when its stages are"
             " independent and could run in parallel.",
    )
    args = parser.parse_args()

    if args.dry_run:
        lines = manifest.read_manifest(args.manifest)
        jobs = lines if args.no_split else split_independent(lines)
        if len(jobs) != len(lines):
            print(
                f"{len(lines)} manifest line(s) split into {len(jobs)} independent run(s)"
            )
        for stages in jobs:
            todo = [
                f"stage {i}"
                for i, stage in enumerate(stages, start=1)
                if args.force or not continuation.is_finished(stage.mc_file)
            ]
            print(f"  {task_name(stages)}: {', '.join(todo) if todo else 'all done'}")
        print(f"{len(jobs)} job(s)")
        return

    results = run_manifest(
        args.manifest,
        workers=args.workers,
        resume=args.resume,
        force=args.force,
        split=not args.no_split,
        log_dir=args.log_dir,
        continue_root=args.continue_root,
    )
    raise SystemExit(1 if any(r.returncode != 0 for r in results) else 0)


if __name__ == "__main__":
    main()
