"""
Simulation parameter dataclasses with JSON serialization.
"""

import json
from dataclasses import dataclass, asdict, replace, field
from pathlib import Path
from enum import StrEnum

from .. import contact_utils as cu
from ..json_dump import make_json_file as write_json, load_json
from .paths import (
    c_path,
    run_slug,
    input_file_paths,
    data_dir,
    structures_dir,
    energy_av_dir,
    energy_records_dir,
)


# ── Dataclasses ─────────────────────────────────────────────────────


DEFAULT_MOVE_PROBAS: dict[str, float] = {
    "swap_empty_full": 0.25,
    "rotate": 0.25,
    "rotate_and_swap_w_empty": 0.25,
    "swap_full_full": 0.25,
}


def calc_ec_ed_ratio(target_radius: float) -> float:
    """Crystal-to-defect energy ratio that stabilises a domain of `target_radius`."""
    return (2 * target_radius**2 - 2 * target_radius - 1) / (3 * target_radius**2)


class InitializeOptions(StrEnum):
    RANDOM = "random"
    FROM_FILE = "from_file"


@dataclass
class ModelParams:
    """Parameters describing the physical model (lattice, couplings, moves)."""

    # lattice
    lattice_name: str = "triangular"
    lx: int = 40
    ly: int = 40
    lz: int = 1

    # particles
    n_types: int = 1
    n_particles: list[int] = field(default_factory=lambda: [1000])

    # couplings — stored as the nested list the C++ side expects
    couplings: list | None = None

    # initialization
    initialize_option: InitializeOptions = InitializeOptions.RANDOM
    state_input: str | None = None

    # observables
    state_av_option: bool = False
    e_av_option: bool = True
    e_record_option: bool = False

    # paths filled at write time
    e_av_output: str | None = None
    e_record_output: str | None = None

    # moves
    move_probas: dict[str, float] = field(
        default_factory=lambda: dict(DEFAULT_MOVE_PROBAS)
    )

    def to_dict(self) -> dict:
        """Produce the dict that gets written to JSON for the C++ reader."""
        d = asdict(self)
        # Drop None fields — the C++ side treats absent keys as disabled
        return {k: v for k, v in d.items() if v is not None}

    def get_input_from_prev_mc(
        self, previous_mc_file: str | Path, start_step: int | None = None
    ) -> None:
        self.initialize_option = InitializeOptions.FROM_FILE

        with Path(previous_mc_file).open("r") as f:
            prev_mc = json.load(f)
        if start_step is None:
            self.state_input = (
                f"{prev_mc['final_structure_address']}final_structure.dat"
            )
        else:
            self.state_input = (
                prev_mc["checkpoint_address"] + f"structure_{start_step}.dat"
            )


class CoolingSchedules(StrEnum):
    INVERSE = "inverse"
    EXPONENTIAL = "exponential"
    LINEAR = "linear"


@dataclass
class MCParams:
    """Parameters controlling the Monte Carlo schedule."""

    mcs_eq: int = int(3.125e5)
    mcs_av: int = 100
    cooling_schedule: CoolingSchedules = CoolingSchedules.INVERSE
    Ti: float = 0.0
    Tf: float = 1.0
    Nt: int = 2000

    checkpoint_option: bool = True
    checkpoint_address: str | None = None
    final_structure_address: str | None = None

    model_params_file: str | None = None

    def to_dict(self) -> dict:
        d = asdict(self)
        return {k: v for k, v in d.items() if v is not None}


# ── Coupling helpers ────────────────────────────────────────────────


def camembert_couplings(
    e_crystal: float, crystal_to_defect_ratio: float, e_repel: float = 10.0
) -> list:
    """Build the camembert coupling map from physical parameters."""
    e_defect = e_crystal / crystal_to_defect_ratio
    return cu.get_camembert_cmap(e_crystal, e_defect, e_repel)


# ── High-level builders ─────────────────────────────────────────────


def build_params(
    root: Path,
    run_name: str,
    series_index: int,
    run_index: int,
    model: ModelParams,
    mc: MCParams,
    *,
    label: str | None = None,
    continue_from_mc_file: str | Path | None = None,
    continue_from_step: int | None = None,
) -> tuple[ModelParams, MCParams, Path, Path]:
    """
    Wire a (ModelParams, MCParams) pair to its on-disk locations for a single
    run, creating output directories as needed.

    `model` and `mc` act as templates: the function returns copies with the
    path-dependent fields filled in. `model.couplings` must be set. `label`
    tags the output slug, e.g. to record explicit parameters.

    Returns (model_params, mc_params, model_file, mc_file).
    """
    if model.couplings is None:
        raise ValueError("model.couplings must be set.")

    slug = run_slug(series_index, label)
    model_file, mc_file = input_file_paths(root, run_name, slug, run_index)
    data = data_dir(root, run_name, slug, run_index)

    # Copy the templates so per-run path wiring never leaks back to the caller.
    model = replace(model)
    mc = replace(mc)

    # -- model output dirs --
    if model.e_av_option:
        d = energy_av_dir(data)
        d.mkdir(parents=True, exist_ok=True)
        model.e_av_output = c_path(d)

    if model.e_record_option:
        d = energy_records_dir(data)
        d.mkdir(parents=True, exist_ok=True)
        model.e_record_output = c_path(d)

    if continue_from_mc_file is not None:
        model.get_input_from_prev_mc(continue_from_mc_file, continue_from_step)

    # -- mc output dirs --
    structs = structures_dir(data)
    structs.mkdir(parents=True, exist_ok=True)
    mc.checkpoint_address = c_path(structs)
    mc.final_structure_address = c_path(structs)
    mc.model_params_file = str(model_file.resolve())

    return model, mc, model_file, mc_file


def write_run(
    root: Path,
    run_name: str,
    series_index: int,
    run_index: int,
    model: ModelParams,
    mc: MCParams,
    *,
    label: str | None = None,
    continue_from_mc_file: str | Path | None = None,
    continue_from_step: int | None = None,
) -> tuple[Path, Path]:
    """Build and write parameter files for a single run. Returns (model_file, mc_file)."""
    model, mc, model_file, mc_file = build_params(
        root,
        run_name,
        series_index,
        run_index,
        model,
        mc,
        label=label,
        continue_from_mc_file=continue_from_mc_file,
        continue_from_step=continue_from_step,
    )

    # write_json creates parent directories as needed.
    write_json(model.to_dict(), model_file)
    write_json(mc.to_dict(), mc_file)

    return model_file, mc_file


def write_run_series(
    root: Path,
    run_name: str,
    series_index: int,
    model: ModelParams,
    mc: MCParams,
    n_runs: int = 20,
    *,
    label: str | None = None,
    continue_from_mc_file: str | Path | None = None,
    continue_from_step: int | None = None,
) -> list[tuple[Path, Path]]:
    """Write parameter files for a series of independent runs sharing one model/mc template."""
    return [
        write_run(
            root,
            run_name,
            series_index,
            i,
            model,
            mc,
            label=label,
            continue_from_mc_file=continue_from_mc_file,
            continue_from_step=continue_from_step,
        )
        for i in range(n_runs)
    ]


def write_manifest(path: Path, jobs: list[list[tuple[Path, Path]]]) -> None:
    """Write a manifest of the model and mc parameter files, used as input for the slurm scripts.

    Args:
        path: Path is the path to the manifest file
        jobs: list[list[tuple[Path, Path]]]. Every list element is a list associated with
        a different run. Every member of that list contains the successive stages of that run,
        in the form of (model_file, mc_file) tuples.
    """
    with path.open("w") as f:
        for job_stages in jobs:
            parts = [f"{model.resolve()} {mc.resolve()}" for model, mc in job_stages]
            f.write(" ".join(parts) + "\n")
