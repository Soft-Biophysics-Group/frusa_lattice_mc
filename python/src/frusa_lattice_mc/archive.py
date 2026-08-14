"""Vincent Ouazan-Reboul, 2026-08-03. Co-written by Claude Opus 5.

Reading and writing run data whether it sits on disk or inside a .tar.zst archive.

Runs are archived one directory per archive, `<run_dir>.tar.zst`, holding a
single top-level `<run_dir.name>/` entry. Analysis scripts should not care which
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
import shutil
import tarfile
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from backports.zstd import ZstdFile

ARCHIVE_SUFFIX = ".tar.zst"
DEFAULT_LEVEL = 19


def archive_for(run_dir: Path) -> Path:
    """The archive that holds (or would hold) this run: <run_dir>.tar.zst."""
    return run_dir.with_name(run_dir.name + ARCHIVE_SUFFIX)


def extract_into(archive: Path, dest: Path) -> Path:
    """Unpack `archive` under `dest`. Returns the run directory it wrote.

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

    if len(tops) != 1:
        raise FileNotFoundError(
            f"{archive} should hold exactly one top-level directory, found {sorted(tops)}"
        )
    return dest / tops.pop()


def read_member(archive: Path, relpath: str) -> bytes | None:
    """Stream one file out of the archive without unpacking the rest.
    `relpath` is the path to the wanted file relative to the archive root.
    """
    with (
        ZstdFile(archive, "rb") as compressed,
        tarfile.open(fileobj=compressed, mode="r|") as tar,
    ):
        for member in tar:
            # Members carry the <run_dir>/ prefix; callers think run-relative.
            # Compared whole rather than by suffix: a run holding a nested
            # structures/ would otherwise hand back whichever copy the tar
            # happens to list first.
            _, _, rel = member.name.partition("/")
            if member.isfile() and rel == str(relpath):
                extracted = tar.extractfile(member)
                return None if extracted is None else extracted.read()
    return None

def read_final(archive: Path) -> bytes | None:
    """The run's final structure, straight out of the archive."""
    return read_member(archive, "structures/final_structure.dat")


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
        yield extract_into(archive, Path(tmp))


def source_files(run_dir: Path) -> list[Path]:
    """Every file under `run_dir`, in a stable order."""
    return sorted(p for p in run_dir.rglob("*") if p.is_file())


def dir_size(paths: list[Path]) -> int:
    return sum(p.stat().st_size for p in paths)


def compress_run(run_dir: Path, level: int = DEFAULT_LEVEL) -> Path:
    """Write <run_dir>.tar.zst. Returns the archive path.

    Streams through the compressor rather than building the tar in memory, and
    lands on the final name with a rename, so a cancelled task never leaves
    something a later pass would mistake for a finished archive.
    """
    archive = archive_for(run_dir)
    partial = archive.with_name(archive.name + ".tmp")
    partial.unlink(missing_ok=True)

    with (
        ZstdFile(partial, "wb", level=level) as compressed,
        tarfile.open(fileobj=compressed, mode="w|") as tar,
    ):
        tar.add(run_dir, arcname=run_dir.name)

    os.replace(partial, archive)
    return archive


def verify_archive(archive: Path, run_dir: Path, expected: list[Path]) -> list[str]:
    """Read the archive back. Returns a list of problems, empty when sound."""
    want = {
        str(Path(run_dir.name) / p.relative_to(run_dir)): p.stat().st_size
        for p in expected
    }
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


def extract_run(archive: Path, run_dir: Path) -> list[Path]:
    """Unpack `archive` into `run_dir`. Returns the files written.

    Unpacks into a staging directory next to the destination and renames it into
    place, so an interrupted task never leaves a partial run directory that a
    later pass — or an analysis script — would take for the real thing.
    """
    staging = run_dir.with_name(run_dir.name + ".extracting")
    shutil.rmtree(staging, ignore_errors=True)

    try:
        run_dir.parent.mkdir(parents=True, exist_ok=True)
        os.replace(extract_into(archive, staging), run_dir)
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    return source_files(run_dir)


def relocate(path: Path, run_dir: Path, data_dir: Path) -> Path:
    """Rewrite a path recorded under `run_dir` to point inside `data_dir`.

    Paths in the MC parameter files are absolute and name the canonical run
    directory; when the data was extracted to a temporary folder they have to be
    re-anchored there.
    """
    return data_dir / Path(path).relative_to(run_dir)
