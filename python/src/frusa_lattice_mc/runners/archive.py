"""Vincent Ouazan-Reboul, 2026-08-11. Co-written by Claude Opus 5.

Compress a campaign's run directories into .tar.zst archives — one archive per
run, holding all its stage directories — and unpack them again when the data is
needed on disk.

One run per SLURM array task, the same shape as `frusa-mc-task`: line N of the
manifest is array task N.

    frusa-mc-archive input/02_long_simus_more_particles/model_mc_files.txt "$SLURM_ARRAY_TASK_ID"
    frusa-mc-archive input/02_long_simus_more_particles/model_mc_files.txt "$SLURM_ARRAY_TASK_ID" --decompress
    frusa-mc-archive input/02_long_simus_more_particles/model_mc_files.txt --todo

Resubmitting is safe in both directions. A run already in the state being asked
for is skipped, a run with nothing to work from is skipped, and an interrupted
task leaves neither a half-written archive nor a half-written run directory
behind. The source directory is only ever deleted with --delete-source, and then
only once the archive has been read back and checked; likewise the archive is
only ever deleted with --delete-archive, after a successful extraction.

`--todo` prints the SLURM array specification of the runs still waiting, so a
partially finished campaign is resubmitted for exactly those:

    sbatch --array="$(frusa-mc-archive MANIFEST --todo)" scripts/archive_run.slurm

Nothing here knows about any particular campaign: the manifest names the runs,
and the run directory of each is read from its own MC parameter file.
"""

import argparse
import shutil
import sys
from pathlib import Path

from ..archive import (
    DEFAULT_LEVEL,
    RunDirs,
    archive_for,
    compress_run,
    dir_size,
    extract_run,
    source_files,
    verify_archive,
)
from ..gen_param_functions.manifest import read_manifest, run_dir


def run_directories(manifest_file: Path) -> list[RunDirs]:
    """The directories of every run of the campaign, one entry per manifest line.

    Input files are regenerated on each machine, so the addresses inside the MC
    params are the right ones for wherever this is running.
    """
    return [[run_dir(stage) for stage in job] for job in read_manifest(manifest_file)]


def array_spec(indices: list[int]) -> str:
    """Condense 1-based task numbers into a SLURM --array specification."""
    if not indices:
        return ""
    runs: list[list[int]] = [[indices[0], indices[0]]]
    for i in indices[1:]:
        if i == runs[-1][1] + 1:
            runs[-1][1] = i
        else:
            runs.append([i, i])
    return ",".join(str(a) if a == b else f"{a}-{b}" for a, b in runs)


def compress_one(
    run_dirs: RunDirs, level: int, delete_source: bool, dry_run: bool
) -> int:
    """Compress one run, every stage of it, into one archive. Returns a process exit code."""
    archive = archive_for(run_dirs[0])

    if archive.is_file():
        print(f"already compressed, nothing to do: {archive}")
        return 0
    present = [d for d in run_dirs if d.is_dir()]
    if not present:
        # Most runs of a campaign do not exist on any one machine, and some may
        # never have run. That is not a failure of this job.
        print(f"no such run directory, skipping: {run_dirs[0]}")
        return 0
    # An unfinished run may still be running or be resumed: its archive would go
    # stale, and --delete-source would pull its checkpoints out from under it.
    # final_structure.dat marks a finished stage, as in continuation.is_finished.
    unfinished = [d for d in run_dirs if not (d / "structures/final_structure.dat").is_file()]
    if unfinished:
        print(
            f"REFUSED: {len(unfinished)} of {len(run_dirs)} stage(s) unfinished,"
            f" first {unfinished[0]}",
            file=sys.stderr,
        )
        return 1

    files = source_files(present)
    if not files:
        print(f"run directory is empty, skipping: {run_dirs[0]}")
        return 0

    before = dir_size(files)
    summary = f"{len(files)} files, {before / 1e6:.1f} MB"
    if dry_run:
        print(f"would compress into {archive} ({summary})")
        if delete_source:
            print(f"would then delete the source directories of {archive}")
        return 0

    print(f"compressing into {archive} ({summary})")
    archive = compress_run(present, level)

    problems = verify_archive(archive, files)
    if problems:
        print(f"ARCHIVE FAILED VERIFICATION: {archive}", file=sys.stderr)
        for problem in problems[:10]:
            print(f"  {problem}", file=sys.stderr)
        print("leaving the source directory untouched", file=sys.stderr)
        return 1

    after = archive.stat().st_size
    print(
        f"wrote {archive} ({after / 1e6:.1f} MB, {before / max(after, 1):.1f}x smaller)"
    )

    if delete_source:
        for directory in present:
            shutil.rmtree(directory)
        print(f"deleted the source directories of {archive}")
    return 0


def decompress_one(
    run_dirs: RunDirs, delete_archive: bool, dry_run: bool, dest: Path | None = None
) -> int:
    """Restore one run, every stage of it, from its archive. Returns a process exit code.

    With `dest`, the run lands in `dest/<series>/` and the canonical location is
    left alone — the way to get at archived data without disturbing data/. The
    series directory is kept in the path because run names repeat across
    series, and two `run_5`s would otherwise collide in one destination.
    """
    archive = archive_for(run_dirs[0])
    series_dir = run_dirs[0].parent
    target = series_dir if dest is None else dest / series_dir.name

    if all((target / d.name).is_dir() for d in run_dirs):
        print(f"already on disk, nothing to do: {target / run_dirs[0].name}")
        return 0
    if not archive.is_file():
        # Same story as compression: most runs live on some other machine.
        print(f"no such archive, skipping: {archive}")
        return 0

    if dry_run:
        size = archive.stat().st_size
        print(f"would extract {archive} ({size / 1e6:.1f} MB) into {target}")
        if delete_archive:
            print(f"would then delete {archive}")
        return 0

    print(f"extracting {archive} ({archive.stat().st_size / 1e6:.1f} MB)")
    try:
        files = extract_run(archive, target)
    except Exception as exc:
        # A truncated or corrupt archive trips the zstd frame checksum or the
        # tar reader; either way nothing has landed at the target.
        print(f"EXTRACTION FAILED: {archive}: {exc}", file=sys.stderr)
        print("leaving the archive untouched", file=sys.stderr)
        return 1

    print(f"wrote {target} ({len(files)} files, {dir_size(files) / 1e6:.1f} MB)")

    if delete_archive:
        archive.unlink()
        print(f"deleted {archive}")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="Manifest naming the runs.")
    parser.add_argument(
        "index",
        type=int,
        nargs="?",
        help="1-based manifest line, e.g. $SLURM_ARRAY_TASK_ID",
    )
    parser.add_argument(
        "-d",
        "--decompress",
        action="store_true",
        help="Unpack the run back out of its archive instead of compressing it.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="Print every run this would cover, then exit.",
    )
    parser.add_argument(
        "--todo",
        action="store_true",
        help="Print the --array spec of the runs still to do, then exit.",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Say what would happen without writing or deleting anything.",
    )
    parser.add_argument(
        "--level",
        type=int,
        default=DEFAULT_LEVEL,
        help=f"zstd compression level (default {DEFAULT_LEVEL}).",
    )
    parser.add_argument(
        "--delete-source",
        action="store_true",
        help="Delete the run directory once its archive verifies. Off by default.",
    )
    parser.add_argument(
        "--delete-archive",
        action="store_true",
        help="With --decompress: delete the archive once extracted. Off by default.",
    )
    parser.add_argument(
        "--dest",
        type=Path,
        help="With --decompress: unpack under DEST instead of back into data/.",
    )
    args = parser.parse_args()

    if args.delete_source and args.decompress:
        parser.error("--delete-source belongs to compression, not --decompress")
    if args.delete_archive and not args.decompress:
        parser.error("--delete-archive only makes sense with --decompress")
    if args.dest is not None and not args.decompress:
        parser.error("--dest only makes sense with --decompress")
    if args.dest is not None and args.delete_archive:
        # Deleting the archive after unpacking it somewhere throwaway would
        # leave the campaign with no copy of the run at all.
        parser.error("--delete-archive with --dest would destroy the only copy")

    run_dirs = run_directories(args.manifest)
    # A run's first directory stands for the whole run in the overviews below:
    # its stages are archived, and deleted, together.
    runs = [dirs[0] for dirs in run_dirs]

    if args.list:
        for i, run in enumerate(runs, start=1):
            on_disk, compressed = run.is_dir(), archive_for(run).is_file()
            match on_disk, compressed:
                case True, True:
                    state = "both"
                case True, False:
                    state = "present"
                case False, True:
                    state = "compressed"
                case _:
                    state = "-"
            print(f"{i:4}  {state:10}  {run}")
        present = sum(d.is_dir() for d in runs)
        done = sum(archive_for(d).is_file() for d in runs)
        print(f"\n{len(runs)} run(s): {present} present, {done} compressed")
        return

    if args.todo:
        # A run needs compressing while it has a directory and no archive, and
        # needs restoring in the mirror case. Either way, a run that is absent
        # in both forms is nobody's work to do.
        if args.decompress:
            todo = [
                i
                for i, run in enumerate(runs, start=1)
                if archive_for(run).is_file() and not run.is_dir()
            ]
        else:
            todo = [
                i
                for i, run in enumerate(runs, start=1)
                if run.is_dir() and not archive_for(run).is_file()
            ]
        print(array_spec(todo))
        # An empty array spec is not something sbatch will accept, and there is
        # nothing to submit anyway.
        raise SystemExit(0 if todo else 1)

    if args.index is None:
        parser.error("give a manifest line number, or --list / --todo")
    if not 1 <= args.index <= len(runs):
        raise SystemExit(
            f"Line {args.index} is out of range: {args.manifest} has {len(runs)} line(s)."
        )

    dirs = run_dirs[args.index - 1]
    if args.decompress:
        raise SystemExit(
            decompress_one(dirs, args.delete_archive, args.dry_run, args.dest)
        )
    raise SystemExit(compress_one(dirs, args.level, args.delete_source, args.dry_run))


if __name__ == "__main__":
    main()
