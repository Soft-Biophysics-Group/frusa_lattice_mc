"""Round-trip tests for the manifest.

write_manifest and read_manifest are inverses: whatever shape of job goes in
comes back out. The cases that matter are the ones a campaign script can hit —
a one-stage campaign, a chained one, and reading a manifest written before the
column order was settled.
"""

import json
from pathlib import Path

import pytest
from frusa_lattice_mc.gen_param_functions.manifest import (
    Stage,
    read_manifest,
    write_manifest,
)


def make_stage(tmp_path: Path, name: str) -> Stage:
    """A stage whose two param files hold just enough to be told apart."""
    model_file = tmp_path / f"model_params_{name}.json"
    mc_file = tmp_path / f"mc_params_{name}.json"
    model_file.write_text(json.dumps({"lx": 4, "ly": 4}))
    mc_file.write_text(json.dumps({"mcs_eq": 10, "Nt": 100}))
    return Stage(model_file, mc_file)


def test_one_stage_jobs_round_trip(tmp_path):
    stages = [make_stage(tmp_path, str(i)) for i in range(3)]
    manifest = tmp_path / "model_mc_files.txt"

    write_manifest(manifest, [[stage] for stage in stages])

    assert read_manifest(manifest) == [[stage] for stage in stages]


def test_chained_jobs_keep_their_stage_order(tmp_path):
    short = [make_stage(tmp_path, f"short_{i}") for i in range(2)]
    long = [make_stage(tmp_path, f"long_{i}") for i in range(2)]
    manifest = tmp_path / "manifest.txt"

    write_manifest(manifest, [[s, l] for s, l in zip(short, long)])

    assert read_manifest(manifest) == [[short[0], long[0]], [short[1], long[1]]]


def test_columns_are_model_then_mc(tmp_path):
    stage = make_stage(tmp_path, "0")
    manifest = tmp_path / "manifest.txt"

    write_manifest(manifest, [[stage]])

    assert manifest.read_text().split() == [
        str(stage.model_file.resolve()),
        str(stage.mc_file.resolve()),
    ]


def test_reads_a_legacy_mc_first_manifest(tmp_path):
    """Campaigns written before the order was settled must stay readable."""
    stage = make_stage(tmp_path, "0")
    manifest = tmp_path / "manifest.txt"
    manifest.write_text(f"{stage.mc_file} {stage.model_file}\n")

    assert read_manifest(manifest) == [[stage]]


def test_ragged_jobs_are_refused(tmp_path):
    stages = [make_stage(tmp_path, str(i)) for i in range(3)]
    manifest = tmp_path / "manifest.txt"

    with pytest.raises(ValueError, match="stages"):
        write_manifest(manifest, [[stages[0], stages[1]], [stages[2]]])


def test_a_failed_write_leaves_the_previous_manifest_intact(tmp_path):
    """The regression: the old writer truncated on open, so a failure mid-write
    destroyed a good manifest and left an empty file behind."""
    stages = [make_stage(tmp_path, str(i)) for i in range(3)]
    manifest = tmp_path / "manifest.txt"
    write_manifest(manifest, [[stage] for stage in stages])
    good = manifest.read_bytes()

    with pytest.raises(ValueError):
        write_manifest(manifest, [[stages[0], stages[1]], [stages[2]]])

    assert manifest.read_bytes() == good
