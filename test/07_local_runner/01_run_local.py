"""Vincent Ouazan-Reboul, 2026-07-30. Co-written by Claude Opus 5.

Check that the local runner stages, skips and resumes correctly.

Run 00_generate_inputs.py first.
"""

import json
from pathlib import Path

from frusa_lattice_mc.gen_param_functions import continuation, manifest
from frusa_lattice_mc.runners.local import run_manifest

ROOT = Path(__file__).parent.resolve()
MANIFEST_FILE = ROOT / "input" / "staged" / "manifest.txt"


def check(label: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(label)
    print(f"  ok: {label}")


def structures_dir(mc_file: Path) -> Path:
    return Path(json.loads(mc_file.read_text())["final_structure_address"])


def interrupt(mc_file: Path, keep: int) -> None:
    """Leave `keep` checkpoints and remove the final structure."""
    structures = structures_dir(mc_file)
    (structures / "final_structure.dat").unlink()
    for checkpoint in structures.glob("structure_*.dat"):
        if int(checkpoint.stem.removeprefix("structure_")) >= keep:
            checkpoint.unlink()


def main() -> None:
    jobs = manifest.read(MANIFEST_FILE)
    n_stages = sum(len(stages) for stages in jobs)

    print("1. clean run")
    results = run_manifest(MANIFEST_FILE, workers=2, verbose=False)
    check("no failures", all(r.returncode == 0 for r in results))
    check(f"{n_stages} stages ran", sum(r.n_ran for r in results) == n_stages)
    check(
        "every stage wrote its final structure",
        all(continuation.is_finished(s.mc_file) for stages in jobs for s in stages),
    )

    print("2. re-run is a no-op")
    results = run_manifest(MANIFEST_FILE, workers=2, verbose=False)
    check("nothing ran", sum(r.n_ran for r in results) == 0)
    check(f"{n_stages} stages skipped", sum(r.n_skipped for r in results) == n_stages)

    print("3. resume after an interruption")
    victim = jobs[0][-1]
    interrupt(victim.mc_file, keep=3)
    check("final structure is gone", not continuation.is_finished(victim.mc_file))
    check(
        "checkpoints survive",
        continuation.checkpoint_progress(
            Path(json.loads(victim.mc_file.read_text())["checkpoint_address"])
        )
        == 2,
    )

    results = run_manifest(MANIFEST_FILE, workers=2, resume=True, verbose=False)
    check("no failures", all(r.returncode == 0 for r in results))
    check("exactly one stage resumed", sum(r.n_resumed for r in results) == 1)
    check("final structure is back", continuation.is_finished(victim.mc_file))
    check(
        "continuation inputs were written",
        any((MANIFEST_FILE.parent / "_local_continued").glob("*/*/mc_params_*.json")),
    )

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
