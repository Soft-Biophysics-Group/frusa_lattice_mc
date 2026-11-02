"""Per-file gzip compression of a run's structure checkpoints.

Each structure_N.dat becomes structure_N.dat.gz in place. Structures are dense text
dumps of a mostly empty lattice and shrink ~250x this way. Unlike one archive per
run, every frame stays a separate file, readable without decompressing the rest:
np.loadtxt opens .gz transparently, and the resume logic's structure_(\\d+)\\.dat
pattern still counts compressed checkpoints.

The C++ only reads plain text, so what it may read back stays plain:
final_structure.dat (whose presence also marks a run finished), the newest
checkpoint of an unfinished run (what a resume starts from) and any file a stage
names as its state_input.
"""

import gzip
import os
import re
import shutil
from collections.abc import Collection
from pathlib import Path

GZ_SUFFIX = ".gz"
_CHECKPOINT_RE = re.compile(r"structure_\d+\.dat")


def structure_path(struct_dir: str | Path, index: int | None = None) -> Path:
    """The file holding structure `index` (the final structure if None).

    Plain when present, else its .gz; the plain name when neither exists, so
    errors name the file a reader expects.
    """
    name = "final_structure.dat" if index is None else f"structure_{index}.dat"
    plain = Path(struct_dir) / name
    compressed = plain.with_name(name + GZ_SUFFIX)
    if not plain.is_file() and compressed.is_file():
        return compressed
    return plain


def gzip_in_place(path: Path) -> Path:
    """Replace `path` by `path`.gz.

    The .gz is written under a temporary name and moved into place before the
    original is removed, so an interruption leaves at worst a stray .tmp or both
    copies, never a truncated .gz and no original.
    """
    target = path.with_name(path.name + GZ_SUFFIX)
    tmp = path.with_name(target.name + ".tmp")
    with path.open("rb") as src, gzip.open(tmp, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.replace(tmp, target)
    path.unlink()
    return target


def compress_checkpoints(
    struct_dir: Path, keep: Collection[Path] = (), dry_run: bool = False
) -> list[Path]:
    """Gzip every plain structure_N.dat in `struct_dir` not in `keep`.

    Returns the plain files compressed (or that would be, with `dry_run`).
    final_structure.dat never matches. Leftover .tmp files from an interrupted
    pass are removed: the original they were made from is still on disk.
    """
    keep = {Path(p).resolve() for p in keep}
    if not struct_dir.is_dir():
        return []
    todo = sorted(
        p for p in struct_dir.iterdir()
        if _CHECKPOINT_RE.fullmatch(p.name) and p.resolve() not in keep
    )
    if dry_run:
        return todo
    for stale in struct_dir.glob(f"structure_*.dat{GZ_SUFFIX}.tmp"):
        stale.unlink()
    for path in todo:
        gzip_in_place(path)
    return todo
