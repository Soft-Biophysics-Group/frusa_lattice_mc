"""Vincent Ouazan-Reboul, 2026-07-31. Co-written by Claude Opus 5.

Merge legacy `structures_continued/` directories back into `structures/`.

Continuations used to write into a sibling directory, renumbering from 0, so an
interrupted run ended up split across two folders with two overlapping index
series. They now share one directory and one global numbering
(structure_index_offset). This migrates data written under the old scheme:
continued structure_i becomes structures/structure_{boundary + i}.dat, where
`boundary` is the global index the continuation resumed at.

Dry run by default — nothing moves until you pass --apply.

    python scripts/migrate_continued_structures.py ../../data
    python scripts/migrate_continued_structures.py ../../data --apply

`boundary` is read from the run's own parameter files when they can be found
(original Nt minus continuation Nt) and cross-checked against the highest index
in structures/. The two must agree, or the run is refused.

If any run looks wrong, nothing is migrated — moving data is hard to undo, so a
surprise should be understood before the rest proceeds. To migrate the healthy
runs meanwhile, point the script at a narrower directory.
"""

import argparse
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

STRUCT_RE = re.compile(r"structure_(\d+)\.dat")
CONTINUED_DIR = "structures_continued"
FINAL_NAME = "final_structure.dat"


def indices(directory: Path) -> list[int]:
    """Every structure index in `directory`, sorted."""
    return sorted(
        int(m.group(1))
        for p in directory.iterdir()
        if (m := STRUCT_RE.fullmatch(p.name))
    )


def gaps(idx: list[int]) -> list[int]:
    """Indices missing from an otherwise contiguous 0..max series."""
    return sorted(set(range(idx[-1] + 1)) - set(idx)) if idx else []


@dataclass
class Migration:
    """One run's planned merge, plus whatever makes it unsafe."""

    run_dir: Path
    structures: Path
    continued: Path
    boundary: int = 0
    boundary_from_params: int | None = None
    moves: list[tuple[Path, Path]] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.problems


# ── Locating the parameter files that describe a run ────────────────


def mc_files_by_checkpoint(input_roots: list[Path]) -> dict[Path, list[dict]]:
    """Index every mc param file under `input_roots` by its checkpoint directory."""
    by_dir: dict[Path, list[dict]] = {}
    for root in input_roots:
        for mc_file in root.rglob("mc_params_*.json"):
            try:
                mc = json.loads(mc_file.read_text())
                key = Path(mc["checkpoint_address"]).resolve()
            except (json.JSONDecodeError, KeyError, OSError):
                continue
            by_dir.setdefault(key, []).append(mc)
    return by_dir


def boundary_from_params(m: Migration, by_dir: dict[Path, list[dict]]) -> int | None:
    """
    The index the continuation resumed at, from Nt shrinkage.

    The old continuation kept Tf and dropped the steps already done, so the
    number of steps it lost is exactly where it picked up.
    """
    originals = by_dir.get(m.structures.resolve(), [])
    continuations = by_dir.get(m.continued.resolve(), [])
    if not originals or not continuations:
        return None

    candidates = {
        o["Nt"] - c["Nt"]
        for o in originals
        for c in continuations
        if "Nt" in o and "Nt" in c
    }
    if not candidates:
        return None
    if len(candidates) != 1:
        m.problems.append(
            f"parameter files disagree on where the continuation resumed: {sorted(candidates)}"
        )
        return None
    return candidates.pop()


# ── Planning ────────────────────────────────────────────────────────


def plan(run_dir: Path, by_dir: dict[Path, list[dict]]) -> Migration:
    """Work out one run's merge without touching anything."""
    m = Migration(run_dir, run_dir / "structures", run_dir / CONTINUED_DIR)

    if not m.structures.is_dir():
        m.problems.append(f"no structures/ beside {CONTINUED_DIR}/")
        return m

    original, continued = indices(m.structures), indices(m.continued)
    if not continued:
        m.notes.append(f"{CONTINUED_DIR}/ holds no structures")
    if not original:
        m.problems.append("structures/ is empty, so the boundary cannot be derived")
        return m

    m.boundary = original[-1] + 1
    m.boundary_from_params = boundary_from_params(m, by_dir)
    if m.boundary_from_params is None:
        m.notes.append("no parameter files found; boundary taken from the files on disk")
    elif m.boundary_from_params != m.boundary:
        m.problems.append(
            f"boundary mismatch: files on disk say {m.boundary},"
            f" parameters say {m.boundary_from_params}"
        )
        return m

    for label, idx in (("structures", original), (CONTINUED_DIR, continued)):
        if missing := gaps(idx):
            m.notes.append(
                f"{label}/ is missing {len(missing)} index(es), e.g. {missing[:5]}"
            )

    for i in continued:
        source = m.continued / f"structure_{i}.dat"
        target = m.structures / f"structure_{i + m.boundary}.dat"
        if target.exists():
            m.problems.append(f"{target.name} already exists; refusing to overwrite")
            return m
        m.moves.append((source, target))

    final = m.continued / FINAL_NAME
    if final.is_file():
        if (m.structures / FINAL_NAME).is_file():
            m.problems.append(f"{FINAL_NAME} exists in both directories")
            return m
        m.moves.append((final, m.structures / FINAL_NAME))

    if leftovers := [
        p.name
        for p in m.continued.iterdir()
        if not STRUCT_RE.fullmatch(p.name) and p.name != FINAL_NAME
    ]:
        m.notes.append(f"leaving unrecognised file(s) in place: {leftovers[:5]}")

    return m


# ── Applying ────────────────────────────────────────────────────────


def apply(m: Migration) -> int:
    """Perform one planned merge. Returns the number of files moved."""
    moved = 0
    for source, target in m.moves:
        source.rename(target)
        moved += 1

    remaining = list(m.continued.iterdir())
    if not remaining:
        m.continued.rmdir()
    return moved


def stale_input_dirs(input_roots: list[Path], continued_dirs: set[Path]) -> list[Path]:
    """Continuation input trees whose param files point at a merged-away directory."""
    stale = set()
    for root in input_roots:
        for mc_file in root.rglob("mc_params_*.json"):
            try:
                mc = json.loads(mc_file.read_text())
            except (json.JSONDecodeError, OSError):
                continue
            if Path(mc.get("checkpoint_address", "/nonexistent")).resolve() in continued_dirs:
                stale.add(mc_file.parent.parent)
    return sorted(stale)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_root", type=Path, help="Directory tree to scan for runs.")
    parser.add_argument(
        "--input-root", type=Path, action="append", default=None,
        help="Where to look for parameter files. Repeatable."
             " Defaults to <data_root>/../input.",
    )
    parser.add_argument(
        "--apply", action="store_true", help="Actually move files. Off by default."
    )
    args = parser.parse_args()

    data_root = args.data_root.resolve()
    input_roots = [p.resolve() for p in args.input_root] if args.input_root else []
    if not input_roots:
        default = (data_root.parent / "input").resolve()
        input_roots = [default] if default.is_dir() else []

    run_dirs = sorted(p.parent for p in data_root.rglob(CONTINUED_DIR) if p.is_dir())
    if not run_dirs:
        print(f"No {CONTINUED_DIR}/ directories under {data_root} — nothing to migrate.")
        return

    by_dir = mc_files_by_checkpoint(input_roots)
    if input_roots and not by_dir:
        print(f"Warning: no parameter files found under {', '.join(map(str, input_roots))}\n")

    plans = [plan(d, by_dir) for d in run_dirs]

    for m in plans:
        source = (
            "disk + parameters"
            if m.boundary_from_params is not None
            else "disk only"
        )
        print(f"{m.run_dir.relative_to(data_root)}")
        if m.ok:
            targets = [t for _, t in m.moves if STRUCT_RE.fullmatch(t.name)]
            span = (
                f"{m.boundary}..{m.boundary + len(targets) - 1}" if targets else "nothing"
            )
            print(f"  boundary {m.boundary} ({source}) -> renumber to {span}")
            print(f"  {len(m.moves)} file(s) to move")
        for note in m.notes:
            print(f"  note: {note}")
        for problem in m.problems:
            print(f"  PROBLEM: {problem}")
        print()

    good = [m for m in plans if m.ok and m.moves]
    bad = [m for m in plans if not m.ok]

    if bad:
        print(f"{len(bad)} run(s) cannot be migrated safely; none will be touched.")
        print("Resolve the problems above and re-run.")
        raise SystemExit(1)

    if not args.apply:
        total = sum(len(m.moves) for m in good)
        print(f"Dry run: {total} file(s) across {len(good)} run(s) would move.")
        print("Re-run with --apply to perform the migration.")
        return

    merged_dirs = {m.continued.resolve() for m in good}
    for m in good:
        moved = apply(m)
        print(f"{m.run_dir.relative_to(data_root)}: moved {moved} file(s)")

    print(f"\nMigrated {len(good)} run(s).")
    print(
        "Their checkpoints are now a single series in structures/, numbered as an"
        "\nuninterrupted run would have numbered them — resume from the original"
        "\nmanifest and the offset follows from the files on disk."
    )

    if stale := stale_input_dirs(input_roots, merged_dirs):
        print("\nThese continuation input trees now point at directories that no longer")
        print("exist. They are superseded; move or delete them so nothing picks them up:")
        for path in stale:
            print(f"  {path}")


if __name__ == "__main__":
    main()
