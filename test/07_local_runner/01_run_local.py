"""Vincent Ouazan-Reboul, 2026-07-30. Co-written by Claude Opus 5.

Check that the local runner stages, skips and resumes correctly.

Run 00_generate_inputs.py first.
"""

import json
import subprocess
import sys
from pathlib import Path

from frusa_lattice_mc.gen_param_functions import continuation, manifest, slurm
from frusa_lattice_mc.runners.local import run_manifest

ROOT = Path(__file__).parent.resolve()
MANIFEST_FILE = ROOT / "input" / "staged" / "manifest.txt"


def check(label: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(label)
    print(f"  ok: {label}")


def structures_dir(mc_file: Path) -> Path:
    return Path(json.loads(mc_file.read_text())["final_structure_address"])


def checkpoint_indices(mc_file: Path) -> list[int]:
    """Every structure index on disk for this stage, in order."""
    return sorted(
        int(p.stem.removeprefix("structure_"))
        for p in structures_dir(mc_file).glob("structure_*.dat")
    )


def run_array_task(index: int) -> str:
    """One SLURM array task: `frusa-mc-task <manifest> <index>`."""
    result = subprocess.run(
        [sys.executable, "-m", "frusa_lattice_mc.runners.task",
         str(MANIFEST_FILE), str(index)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise AssertionError(f"array task {index} failed:\n{result.stdout}{result.stderr}")
    return result.stdout


def interrupt(mc_file: Path, keep: int) -> None:
    """Leave `keep` checkpoints and remove the final structure."""
    structures = structures_dir(mc_file)
    (structures / "final_structure.dat").unlink()
    for checkpoint in structures.glob("structure_*.dat"):
        if int(checkpoint.stem.removeprefix("structure_")) >= keep:
            checkpoint.unlink()


def main() -> None:
    jobs = manifest.read_manifest(MANIFEST_FILE)
    n_stages = sum(len(stages) for stages in jobs)
    n_t = json.loads(jobs[0][0].mc_file.read_text())["Nt"]

    print("1. clean run")
    results = run_manifest(MANIFEST_FILE, workers=2, verbose=False)
    check("no failures", all(r.returncode == 0 for r in results))
    check(f"{n_stages} stages ran", sum(r.n_ran for r in results) == n_stages)
    check(
        "every stage wrote its final structure",
        all(continuation.is_finished(s.mc_file) for stages in jobs for s in stages),
    )
    check(
        "every stage wrote one checkpoint per temperature",
        all(
            checkpoint_indices(s.mc_file) == list(range(n_t))
            for stages in jobs
            for s in stages
        ),
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
        any((MANIFEST_FILE.parent / "_continued").glob("*/*/*/mc_params_*.json")),
    )
    check(
        "checkpoints form one gapless series",
        checkpoint_indices(victim.mc_file) == list(range(n_t)),
    )
    check(
        "no second structures directory was created",
        not any(structures_dir(victim.mc_file).parent.glob("structures_*")),
    )

    print("4. resume again after a second interruption")
    interrupt(victim.mc_file, keep=5)
    check("only the last checkpoint is missing", checkpoint_indices(victim.mc_file) == list(range(5)))

    results = run_manifest(MANIFEST_FILE, workers=2, resume=True, verbose=False)
    check("no failures", all(r.returncode == 0 for r in results))
    check("exactly one stage resumed", sum(r.n_resumed for r in results) == 1)
    check("final structure is back", continuation.is_finished(victim.mc_file))
    check(
        "checkpoints still form one gapless series",
        checkpoint_indices(victim.mc_file) == list(range(n_t)),
    )

    print("5. both stages of one job interrupted at once")
    # jobs[1] has never been resumed, so no earlier attempt masks this one.
    for stage in jobs[1]:
        interrupt(stage.mc_file, keep=2)

    results = run_manifest(MANIFEST_FILE, workers=2, resume=True, verbose=False)
    check("no failures", all(r.returncode == 0 for r in results))
    check("both stages resumed", sum(r.n_resumed for r in results) == 2)
    check(
        "each stage kept its own gapless series",
        all(checkpoint_indices(s.mc_file) == list(range(n_t)) for s in jobs[1]),
    )
    # The two stages share a slug and differ only by run name, so their
    # continuation inputs must not land on the same path.
    attempts = sorted((MANIFEST_FILE.parent / "_continued").glob("*"))
    check(
        "their continuation inputs stayed distinct",
        len(list(attempts[-1].glob("*/*/mc_params_*.json"))) == 2,
    )

    print("6. the SLURM array task: one script starts, continues and no-ops")
    # Wipe one job entirely so the same command has to start it from scratch,
    # and interrupt the other so the same command has to continue it.
    for stage in jobs[0]:
        for leftover in structures_dir(stage.mc_file).glob("*.dat"):
            leftover.unlink()
    # Keep more than group 5's resume point: a real interruption always leaves
    # more progress than the attempt before it.
    interrupt(jobs[1][-1].mc_file, keep=4)

    fresh, continued = run_array_task(1), run_array_task(2)
    check("nothing was continued in the job that had no checkpoints", "0 continued" in fresh)
    check("the interrupted job was continued", "1 continued" in continued)
    check(
        "both jobs are complete",
        all(continuation.is_finished(s.mc_file) for stages in jobs for s in stages),
    )
    check(
        "every stage is one gapless series",
        all(
            checkpoint_indices(s.mc_file) == list(range(n_t))
            for stages in jobs
            for s in stages
        ),
    )

    again = run_array_task(1)
    check("resubmitting is a no-op", "0 run" in again and f"{len(jobs[0])} skipped" in again)
    check(
        "a no-op task leaves no empty attempt directory",
        all(any(d.iterdir()) for d in (MANIFEST_FILE.parent / "_continued").glob("*")),
    )

    print("7. the generated array script")
    script = slurm.generate_self_resuming_script(MANIFEST_FILE, len(jobs), job_name="t7")
    check("one srun per array task, not per stage", script.count("srun") == 1)
    check("array covers every manifest line", f"--array=1-{len(jobs)}" in script)
    check("passes the task id through", '"$SLURM_ARRAY_TASK_ID"' in script)

    print("8. independent stages on one line run in parallel; chains do not")
    from frusa_lattice_mc.runners.stages import independent_stages, split_independent

    # The fixture's own lines are chains — stage 2 starts from stage 1's final
    # structure — so they must survive the split untouched.
    check("a chained line is not independent", not independent_stages(jobs[0]))
    check("chained lines are left alone", len(split_independent(jobs)) == len(jobs))

    # Both first stages initialise at random, so a line holding the two of them
    # is a grouping rather than an ordering: the shape campaign 04's manifest has.
    grouped = ROOT / "input" / "staged" / "grouped.txt"
    manifest.write_stages(grouped, [[jobs[0][0], jobs[1][0]]])
    line = manifest.read_manifest(grouped)
    check("a line of random-init stages is independent", independent_stages(line[0]))
    check("it splits into one job per stage", len(split_independent(line)) == 2)

    results = run_manifest(grouped, workers=2, force=True, verbose=False)
    check("both split jobs ran", len(results) == 2 and sum(r.n_ran for r in results) == 2)
    check("no failures", all(r.returncode == 0 for r in results))
    results = run_manifest(grouped, workers=2, split=False, verbose=False)
    check("--no-split keeps them on one job", len(results) == 1)

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
