"""Per-file checkpoint compression must never cost a resume or a readable frame."""

import json
from pathlib import Path

import numpy as np

from frusa_lattice_mc.checkpoints import compress_checkpoints, structure_path
from frusa_lattice_mc.config import load_structure
from frusa_lattice_mc.gen_param_functions.continuation import continue_stage
from frusa_lattice_mc.gen_param_functions.manifest import Stage
from frusa_lattice_mc.runners.compress import compress_line
from frusa_lattice_mc.runners.stages import has_checkpoints

STRUCTURE = "0 1 0 1\n-1 3 -1 7\n"


def make_run(tmp_path: Path, n_checkpoints: int, finished: bool) -> Stage:
    """A run as input/campaign/series/*.json, with paths relative to tmp_path."""
    series = tmp_path / "input" / "campaign" / "series"
    series.mkdir(parents=True)
    structures = Path("data") / "series" / "run_0" / "structures"
    (tmp_path / structures).mkdir(parents=True)
    for i in range(n_checkpoints):
        (tmp_path / structures / f"structure_{i}.dat").write_text(STRUCTURE)
    if finished:
        (tmp_path / structures / "final_structure.dat").write_text(STRUCTURE)

    model_file = series / "model_params_0.json"
    mc_file = series / "mc_params_0.json"
    model_file.write_text(json.dumps({"lx": 4, "ly": 4, "initialize_option": "random"}))
    mc_file.write_text(json.dumps({
        "mcs_eq": 10, "mcs_av": 1, "cooling_schedule": "linear",
        "Ti": 1.0, "Tf": 2.0, "Nt": 100,
        "checkpoint_option": True,
        "checkpoint_address": f"{structures}/",
        "final_structure_address": f"{structures}/",
        "model_params_file": str(model_file),
    }))
    return Stage(model_file, mc_file)


def names(directory: Path) -> set[str]:
    return {p.name for p in directory.iterdir()}


def test_structure_path_prefers_plain_and_falls_back_to_gz(tmp_path):
    (tmp_path / "structure_3.dat.gz").touch()
    assert structure_path(tmp_path, 3).name == "structure_3.dat.gz"
    (tmp_path / "structure_3.dat").touch()
    assert structure_path(tmp_path, 3).name == "structure_3.dat"
    assert structure_path(tmp_path, 4).name == "structure_4.dat"
    assert structure_path(tmp_path).name == "final_structure.dat"


def test_compressed_structures_read_back_identically(tmp_path):
    (tmp_path / "structure_0.dat").write_text(STRUCTURE)
    expected = load_structure(struct_folder=tmp_path, struct_index=0)
    compress_checkpoints(tmp_path)
    assert names(tmp_path) == {"structure_0.dat.gz"}
    np.testing.assert_array_equal(
        load_structure(struct_folder=tmp_path, struct_index=0), expected
    )


def test_keep_and_final_structure_stay_plain_and_stale_tmp_goes(tmp_path):
    for i in range(3):
        (tmp_path / f"structure_{i}.dat").write_text(STRUCTURE)
    (tmp_path / "final_structure.dat").write_text(STRUCTURE)
    (tmp_path / "structure_0.dat.gz.tmp").write_text("truncated")

    compress_checkpoints(tmp_path, keep=[tmp_path / "structure_1.dat"])

    assert names(tmp_path) == {
        "structure_0.dat.gz", "structure_1.dat", "structure_2.dat.gz",
        "final_structure.dat",
    }


def test_unfinished_run_stays_resumable(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    stage = make_run(tmp_path, n_checkpoints=10, finished=False)
    structures = tmp_path / "data" / "series" / "run_0" / "structures"

    assert compress_line([stage], dry_run=False) == 9
    assert compress_line([stage], dry_run=False) == 0  # idempotent
    assert "structure_9.dat" in names(structures)
    assert has_checkpoints(stage)

    resumed = continue_stage(stage, tmp_path / "_continued" / "attempt_1")
    assert resumed is not None
    assert json.loads(resumed.mc_file.read_text())["structure_index_offset"] == 10
    state_input = json.loads(resumed.model_file.read_text())["state_input"]
    assert Path(state_input).name == "structure_9.dat"
    assert Path(state_input).is_file()


def test_finished_run_compresses_every_checkpoint(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    stage = make_run(tmp_path, n_checkpoints=10, finished=True)
    structures = tmp_path / "data" / "series" / "run_0" / "structures"

    assert compress_line([stage], dry_run=False) == 10
    assert names(structures) == {f"structure_{i}.dat.gz" for i in range(10)} | {
        "final_structure.dat"
    }
