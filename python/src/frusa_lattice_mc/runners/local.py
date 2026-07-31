"""Vincent Ouazan-Reboul, 2026-07-30. Co-written by Claude Opus 5.

Run a manifest on this machine instead of on the cluster.

Local counterpart of slurm.generate_array_script_by_stages: one worker per
manifest line, running that line's stages in order. Completed stages are
skipped; with --resume, an interrupted stage restarts from its last checkpoint
via temporary param files under <manifest dir>/_local_continued/.

    frusa-mc-local input/03_long_simus_T_1/manifest.txt -j 6 --resume
"""

import argparse
import json
import os
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import NamedTuple

from ..config import find_executable
from ..gen_param_functions import continuation, manifest
from ..gen_param_functions.manifest import Stage

CONTINUE_DIR_NAME = "_local_continued"


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


def newest_inputs(stage: Stage, attempt_dir: Path, continue_root: Path) -> Stage:
    """This stage's latest param files: newest earlier attempt, else the original."""
    rel = stage.mc_file.relative_to(stage.mc_file.parent.parent)
    for previous in sorted(continue_root.glob("*"), reverse=True):
        if previous != attempt_dir and (previous / rel).is_file():
            mc_file = previous / rel
            model_file = json.loads(mc_file.read_text())["model_params_file"]
            return Stage(Path(model_file), mc_file)
    return stage


def resume_stage(stage: Stage, attempt_dir: Path, continue_root: Path) -> Stage | None:
    """Temp param files resuming from the last checkpoint, or None if there is none."""
    source = newest_inputs(stage, attempt_dir, continue_root)
    mc = json.loads(source.mc_file.read_text())
    progress = continuation.checkpoint_progress(Path(mc["checkpoint_address"]))
    if progress < 0:
        return None

    input_root = source.mc_file.parent.parent
    job = continuation.UnfinishedJob(
        mc_file=source.mc_file,
        mc_file_relative=source.mc_file.relative_to(input_root),
        last_completed_index=progress,
        total_steps=mc["Nt"],
    )
    created = continuation.write_continuation_inputs([job], input_root, attempt_dir)
    if not created:
        return None
    model_file = json.loads(created[0].read_text())["model_params_file"]
    return Stage(Path(model_file), created[0])


def publish_final_structure(original: Stage, resumed: Stage) -> None:
    """Copy the resumed run's final structure where the original stage promised it."""
    source = Path(json.loads(resumed.mc_file.read_text())["final_structure_address"])
    target = Path(json.loads(original.mc_file.read_text())["final_structure_address"])
    source, target = source / "final_structure.dat", target / "final_structure.dat"
    if source.is_file() and source != target:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


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

            if to_run is not stage:
                publish_final_structure(stage, to_run)
            n_ran += 1

    return TaskResult(name, n_ran, n_skipped, n_resumed, None, 0)


def run_manifest(
    manifest_path: Path,
    *,
    workers: int | None = None,
    resume: bool = False,
    force: bool = False,
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

    jobs = manifest.read_manifest(manifest_path)
    if not jobs:
        raise ValueError(f"{manifest_path} is empty")

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

    attempt_dir = None
    if resume:
        attempt_dir = continue_root / time.strftime("%y%m%d_%H%M%S")
        attempt_dir.mkdir(parents=True)

    log_dir.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    results: list[TaskResult] = []

    if verbose:
        n_stages = max(len(stages) for stages in jobs)
        print(f"{len(jobs)} run(s), up to {n_stages} stage(s) each, {workers} at a time")
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
                print(
                    f"[{n_done}/{len(jobs)}] {result.name}: {status}"
                    f"  ({(time.monotonic() - start) / 60:.1f} min elapsed)"
                )

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
    args = parser.parse_args()

    if args.dry_run:
        for stages in manifest.read(args.manifest):
            todo = [
                f"stage {i}"
                for i, stage in enumerate(stages, start=1)
                if args.force or not continuation.is_finished(stage.mc_file)
            ]
            print(f"  {task_name(stages)}: {', '.join(todo) if todo else 'all done'}")
        return

    results = run_manifest(
        args.manifest,
        workers=args.workers,
        resume=args.resume,
        force=args.force,
        log_dir=args.log_dir,
        continue_root=args.continue_root,
    )
    raise SystemExit(1 if any(r.returncode != 0 for r in results) else 0)


if __name__ == "__main__":
    main()
