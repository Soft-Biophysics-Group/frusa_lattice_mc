"""Vincent Ouazan-Reboul, 2026-08-03. Co-written by Claude Opus 5.

Reading and writing run data whether it sits on disk or inside a .tar.zst archive.

Runs are archived one per archive, `<run>.tar.zst`, next to the run directory.
A single-stage run's archive holds one top-level `<run>/` entry; a multi-stage
run's holds its stage directories `<run>_stage_NNN/` side by side, so a whole
run is one file however many stages it has. Member names are always relative to
the series directory, exactly as on disk. Analysis scripts should not care which
of the two forms a run happens to be in on the machine they run on:

    with run_data(run_dir) as data_dir:
        if data_dir is None:
            continue                      # neither directory nor archive here
        state = LatticeState(data_dir / "structures/final_structure.dat", model)

`data_dir` is the run directory itself when it exists (nothing is copied), and a
temporary extraction otherwise, deleted when the block ends.

This module owns the archive format — how a run is packed, checked and unpacked
— and nothing else. It imports the standard library and the zstd binding only,
so it stays usable on nodes where the geometry stack (and its scipy extension
modules) cannot be imported. The CLI that drives it over a campaign lives in
`frusa_lattice_mc.runners.archive`.
"""

import os
import re
import shutil
import tarfile
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from backports.zstd import ZstdFile

ARCHIVE_SUFFIX = ".tar.zst"
# @ AGENT: What motivated this value? Is it the standard for zstd?
DEFAULT_LEVEL = 19
# Mirrors gen_param_functions.paths._stage_name; the two must change together.
_STAGE_SUFFIX = re.compile(r"_stage_\d+$")

# The directories of one run, all in one series directory: its stage directories
# in stage order, or just the run directory of a single-stage run.
RunDirs = list[Path]


def archive_for(run_dir: Path) -> Path:
    """The archive that holds (or would hold) this run: <run>.tar.zst.

    `run_dir` may be a stage directory: every stage of a run maps to one archive.
    """
    run_name = _STAGE_SUFFIX.sub("", run_dir.name)
    return run_dir.with_name(run_name + ARCHIVE_SUFFIX)


def extract_into(archive: Path, dest: Path) -> list[Path]:
    """Unpack `archive` under `dest`. Returns the top-level directories it wrote.

    Streams the tar through the decompressor rather than holding it in memory,
    the mirror of how the archive was written.
    """
    dest.mkdir(parents=True, exist_ok=True)
    tops = set()
    with (
        ZstdFile(archive, "rb") as compressed,
        tarfile.open(fileobj=compressed, mode="r|") as tar,
    ):
        for member in tar:
            # filter="data" refuses members with absolute or ../ paths, so a
            # damaged archive cannot write outside dest.
            tar.extract(member, dest, filter="data")
            tops.add(Path(member.name).parts[0])
    return sorted(dest / top for top in tops)


def read_member(archive: Path, relpath: str) -> bytes | None:
    """Stream one file out of the archive without unpacking the rest.

    `relpath` is the member's full name, starting with its run or stage
    directory: `run_5/structures/final_structure.dat`. Compared whole rather
    than by suffix, so the same file name in two stages is never confused.
    """
    with (
        ZstdFile(archive, "rb") as compressed,
        # @ AGENT What does r| mean?
        tarfile.open(fileobj=compressed, mode="r|") as tar,
    ):
        for member in tar:
            if member.isfile() and member.name == str(relpath):
                extracted = tar.extractfile(member)
                return None if extracted is None else extracted.read()
    return None


@contextmanager
def run_data(run_dir: Path, tmp_root: Path | None = None) -> Iterator[Path | None]:
    """Yield a directory holding this run's files, or None if the run is absent.

    The uncompressed directory wins when both forms are present: it costs
    nothing to read and is the one a running simulation writes to. Extraction
    goes to a temporary directory that is removed on the way out, so an analysis
    pass over a compressed campaign never needs the disk space for more than one
    run at a time.

    `tmp_root` picks where that temporary directory lives — worth setting to
    node-local scratch on the cluster, where $TMPDIR may be small or shared.
    """
    if run_dir.is_dir():
        yield run_dir
        return

    archive = archive_for(run_dir)
    if not archive.is_file():
        yield None
        return

    with tempfile.TemporaryDirectory(dir=tmp_root, prefix="run_data_") as tmp:
        extract_into(archive, Path(tmp))
        yield Path(tmp) / run_dir.name


def source_files(run_dirs: RunDirs) -> list[Path]:
    """Every file under the `run_dirs`, in a stable order."""
    return sorted(p for d in run_dirs for p in d.rglob("*") if p.is_file())


def dir_size(paths: list[Path]) -> int:
    return sum(p.stat().st_size for p in paths)


def compress_run(run_dirs: RunDirs, level: int = DEFAULT_LEVEL) -> Path:
    """Write one run's directories into <run>.tar.zst. Returns the archive path.

    Streams through the compressor rather than building the tar in memory, and
    lands on the final name with a rename, so a cancelled task never leaves
    something a later pass would mistake for a finished archive.
    """
    archive = archive_for(run_dirs[0])
    partial = archive.with_name(archive.name + ".tmp")
    partial.unlink(missing_ok=True)

    with (
        ZstdFile(partial, "wb", level=level) as compressed,
        tarfile.open(fileobj=compressed, mode="w|") as tar,
    ):
        for run_dir in run_dirs:
            tar.add(run_dir, arcname=run_dir.name)

    os.replace(partial, archive)
    return archive


def verify_archive(archive: Path, expected: list[Path]) -> list[str]:
    """Read the archive back. Returns a list of problems, empty when sound."""
    # The archive sits in the series directory, which members are relative to.
    want = {str(p.relative_to(archive.parent)): p.stat().st_size for p in expected}
    found: dict[str, int] = {}
    with (
        ZstdFile(archive, "rb") as compressed,
        tarfile.open(fileobj=compressed, mode="r|") as tar,
    ):
        for member in tar:
            if member.isfile():
                found[member.name] = member.size

    problems = []
    for name, size in want.items():
        if name not in found:
            problems.append(f"missing from archive: {name}")
        elif found[name] != size:
            problems.append(
                f"size differs for {name}: {size} on disk, {found[name]} in archive"
            )
    return problems


def extract_run(archive: Path, dest_dir: Path) -> list[Path]:
    """Unpack `archive`'s run or stage directories into `dest_dir`, usually the
    series directory. Returns the files written.

    Unpacks into a staging directory next to the destination and renames each
    directory into place, so an interrupted task never leaves a partial run or
    stage directory that a later pass — or an analysis script — would take for
    the real thing.
    """
    staging = dest_dir / (archive.name + ".extracting")
    shutil.rmtree(staging, ignore_errors=True)

    try:
        tops = extract_into(archive, staging)
        targets = [dest_dir / top.name for top in tops]
        for top, target in zip(tops, targets):
            # Never overwrite: a directory already there is either the live copy
            # or was moved in whole by an earlier, interrupted pass.
            if not target.exists():
                os.replace(top, target)
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    return source_files(targets)


def relocate(path: Path, run_dir: Path, data_dir: Path) -> Path:
    """Rewrite a path recorded under `run_dir` to point inside `data_dir`.

    Paths in the MC parameter files are absolute and name the canonical run
    directory; when the data was extracted to a temporary folder they have to be
    re-anchored there.
    """
    return data_dir / Path(path).relative_to(run_dir)
