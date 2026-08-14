"""Regression tests for resuming interrupted runs.

The case that matters here is a continuation that was interrupted *before* it
wrote its first checkpoint. Its structure_index_offset then sits exactly one
ahead of the highest index on disk, which used to read as "nothing to resume
from" and silently sent the runner back to the original params — restarting the
simulation from step 0 on top of an existing dataset.
"""

import json
import shutil
from pathlib import Path

import pytest

from frusa_lattice_mc.gen_param_functions.continuation import (
    checkpoint_progress,
    continue_stage,
)
from frusa_lattice_mc.gen_param_functions.manifest import Stage
from frusa_lattice_mc.runners.stages import REFUSED_OVERWRITE_RC, run_stages


def make_run(tmp_path: Path, n_checkpoints: int, offset: int = 0) -> Stage:
    """A run with `n_checkpoints` structures on disk, as input/campaign/series/*.json.

    continue_stage takes the input root to be mc_file.parents[2], so the three
    levels below `input` are what the layout needs.
    """
    series = tmp_path / "input" / "campaign" / "series"
    series.mkdir(parents=True)
    structures = tmp_path / "data" / "series" / "run_0" / "structures"
    structures.mkdir(parents=True)

    for i in range(n_checkpoints):
        (structures / f"structure_{i}.dat").write_text("0 0\n-1 -1\n")

    model_file = series / "model_params_0.json"
    mc_file = series / "mc_params_0.json"
    model: dict[str, object] = {"lx": 4, "ly": 4}
    if offset:
        # A non-zero offset only ever occurs on a continuation, which by
        # construction initialises from the checkpoint before its first index.
        model["initialize_option"] = "from_file"
        model["state_input"] = f"{structures}/structure_{offset - 1}.dat"
    else:
        model["initialize_option"] = "random"
    model_file.write_text(json.dumps(model))
    mc = {
        "mcs_eq": 10,
        "mcs_av": 1,
        "cooling_schedule": "linear",
        "Ti": 1.0,
        "Tf": 2.0,
        "Nt": 100,
        "checkpoint_option": True,
        "checkpoint_address": f"{structures}/",
        "final_structure_address": f"{structures}/",
        "model_params_file": str(model_file),
    }
    if offset:
        mc["structure_index_offset"] = offset
    mc_file.write_text(json.dumps(mc))
    return Stage(model_file, mc_file)


def test_checkpoint_progress_counts_from_the_offset(tmp_path):
    stage = make_run(tmp_path, n_checkpoints=10)
    structures = json.loads(stage.mc_file.read_text())["checkpoint_address"]
    assert checkpoint_progress(Path(structures)) == 9
    assert checkpoint_progress(Path(structures), offset=4) == 5


def test_continue_stage_resumes_from_the_last_checkpoint(tmp_path):
    stage = make_run(tmp_path, n_checkpoints=10)
    attempt = tmp_path / "_continued" / "attempt_1"

    resumed = continue_stage(stage, attempt)

    assert resumed is not None
    mc = json.loads(resumed.mc_file.read_text())
    model = json.loads(resumed.model_file.read_text())
    # Ten structures written (0..9), so the next one to produce is index 10
    assert mc["structure_index_offset"] == 10
    assert mc["Nt"] == 90
    assert model["initialize_option"] == "from_file"
    assert model["state_input"].endswith("structure_9.dat")


def test_continuation_interrupted_before_its_first_checkpoint_is_still_resumable(
    tmp_path,
):
    """The regression: offset one ahead of the newest structure must not restart.

    A continuation resuming at index 10 whose run died before writing
    structure_10.dat leaves max_idx == 9 == offset - 1. That is not a dead end —
    structure_9.dat is still exactly the right state to restart from.
    """
    # 10 structures on disk (0..9) and a continuation that starts at 10
    stage = make_run(tmp_path, n_checkpoints=10, offset=10)
    attempt = tmp_path / "_continued" / "attempt_2"

    resumed = continue_stage(stage, attempt)

    assert resumed is not None, "must not fall back to running the original params"
    mc = json.loads(resumed.mc_file.read_text())
    assert mc["structure_index_offset"] == 10
    assert json.loads(resumed.model_file.read_text())["state_input"].endswith(
        "structure_9.dat"
    )


def test_continuation_whose_checkpoints_were_deleted_starts_over(tmp_path):
    """Wiping a run's data must not leave stale continuations chasing a gone file.

    The params survive the deletion and still name structure_{offset-1}. Handing
    them back would launch a run whose initial state does not exist.
    """
    stage = make_run(tmp_path, n_checkpoints=10, offset=10)
    structures = Path(json.loads(stage.mc_file.read_text())["checkpoint_address"])
    for leftover in structures.glob("structure_*.dat"):
        leftover.unlink()

    assert continue_stage(stage, tmp_path / "_continued" / "attempt_4") is None


def test_a_run_that_never_checkpointed_starts_over(tmp_path):
    stage = make_run(tmp_path, n_checkpoints=0)
    assert continue_stage(stage, tmp_path / "_continued" / "attempt_3") is None


def test_run_stages_refuses_to_overwrite_existing_checkpoints(tmp_path):
    """Without --resume, restarting over a populated checkpoint dir must not run."""
    stage = make_run(tmp_path, n_checkpoints=10)
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = run_stages(
        [stage],
        log_dir,
        None,  # attempt_dir None == no --resume
        tmp_path / "_continued",
        False,  # force
        Path("/nonexistent/frusa_mc"),  # never reached: the guard fires first
    )

    assert result.returncode == REFUSED_OVERWRITE_RC
    assert result.n_ran == 0
    assert result.note is not None and "refused" in result.note


def test_run_stages_still_starts_a_fresh_run(tmp_path):
    """The guard must not block a campaign that has not produced anything yet."""
    stage = make_run(tmp_path, n_checkpoints=0)
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    # `false` stands in for the simulator: it exists, ignores args, exits non-zero
    false_bin = shutil.which("false")
    assert false_bin is not None
    result = run_stages(
        [stage], log_dir, None, tmp_path / "_continued", False, Path(false_bin)
    )

    # It got past the guard and actually launched the binary
    assert result.returncode != REFUSED_OVERWRITE_RC
    assert result.note is None


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
