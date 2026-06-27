"""
Simulation parameter dataclasses with JSON serialization.
"""

import json
from dataclasses import dataclass, asdict, replace, field
from pathlib import Path
from enum import StrEnum

import contact_utils as cu
from json_dump import make_json_file as write_json, load_json
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


MOVES_WITH_SWAPS: dict[str, float] = {
    "swap_empty_full": 0.25,
    "rotate": 0.25,
    "rotate_and_swap_w_empty": 0.25,
    "swap_full_full": 0.25,
}

MOVES_WITHOUT_SWAPS: dict[str, float] = {
    "swap_empty_full": 1 / 3,
    "rotate": 1 / 3,
    "rotate_and_swap_w_empty": 1 / 3,
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
    n_particles: list[int] = field(default_factory = lambda:[1000])

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
    move_probas: dict[str, float] = field(default_factory=lambda: MOVES_WITH_SWAPS)

    def to_dict(self) -> dict:
        """Produce the dict that gets written to JSON for the C++ reader."""
        d = asdict(self)
        # Drop None fields — the C++ side treats absent keys as disabled
        return {k: v for k, v in d.items() if v is not None}

    def get_input_from_prev_mc(
        self, previous_mc_file: str | Path, start_step: int | None=None
    ) -> None:
        self.initialize_option = InitializeOptions.FROM_FILE

        with Path(previous_mc_file).open("r") as f:
            prev_mc = json.load(f)
        if start_step is None:
            self.state_input = (
                f"{prev_mc['final_structure_address']}final_structure.dat"
            )
        else:
            self.state_input = prev_mc["checkpoint_address"] + f"structure_{start_step}.dat"


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
    run_name:str,
    series_index: int,
    crystal_to_defect_ratio: float,
    n_steps_per_T: int,
    n_particles: int,
    lattice_side:int,
    run_index: int,
    *,
    allow_p_p_swaps: bool = True,
    record_option: bool = False,
    couplings: list | None = None,
    e_crystal: float = -18.7,
    e_repel: float = 10.0,
    model_overrides: dict | None = None,
    mc_overrides: dict | None = None,
    continue_from_mc_file: str | Path | None = None,
    continue_from_step : int | None = None
) -> tuple[ModelParams, MCParams, Path, Path]:
    """
    Construct a (ModelParams, MCParams) pair for a single run,
    creating output directories as needed.

    `couplings` is the flattened coupling map fed to the C++ side. When None it
    defaults to the camembert map built from `e_crystal`, `crystal_to_defect_ratio`
    and `e_repel`; pass any other map (e.g. from `contact_utils.ContactMapWrapper`)
    to run a different model.

    Returns (model_params, mc_params, model_file, mc_file).
    """
    slug = run_slug(series_index, crystal_to_defect_ratio)
    model_file, mc_file = input_file_paths(root, run_name, slug, run_index)
    data = data_dir(root, run_name, slug, run_index)

    # -- model --
    moves = MOVES_WITH_SWAPS if allow_p_p_swaps else MOVES_WITHOUT_SWAPS

    if couplings is None:
        couplings = camembert_couplings(e_crystal, crystal_to_defect_ratio, e_repel)

    model = ModelParams(
        couplings=couplings,
        move_probas=dict(moves),  # copy to avoid shared mutation
        e_av_option=True,
        e_record_option=record_option,
        n_particles=[n_particles],
        lx=lattice_side,
        ly=lattice_side,
    )

    if model.e_av_option:
        d = energy_av_dir(data)
        d.mkdir(parents=True, exist_ok=True)
        model.e_av_output = c_path(d)

    if model.e_record_option:
        d = energy_records_dir(data)
        d.mkdir(parents=True, exist_ok=True)
        model.e_record_output = c_path(d)

    if model_overrides:
        model = replace(model, **model_overrides)

    if continue_from_mc_file is not None:
        model.get_input_from_prev_mc(continue_from_mc_file, continue_from_step)

    # -- mc --
    structs = structures_dir(data)
    structs.mkdir(parents=True, exist_ok=True)

    mc = MCParams(
        mcs_eq=n_steps_per_T,
        checkpoint_address=c_path(structs),
        final_structure_address=c_path(structs),
        model_params_file=str(model_file.resolve()),
    )

    if mc_overrides:
        mc = replace(mc, **mc_overrides)

    return model, mc, model_file, mc_file


def write_run(
    root: Path,
    run_name: str,
    series_index: int,
    crystal_to_defect_ratio: float,
    n_steps_per_T: int,
    n_particles: int,
    lattice_side: int,
    run_index: int,
    *,
    allow_p_p_swaps: bool = True,
    record_option: bool = False,
    couplings: list | None = None,
    model_overrides: dict | None = None,
    mc_overrides: dict | None = None,
    continue_from_mc_file: str | Path | None = None,
    continue_from_step : int | None = None
) -> tuple[Path, Path]:
    """Build and write parameter files for a single run. Returns (model_file, mc_file)."""
    model, mc, model_file, mc_file = build_params(
        root,
        run_name,
        series_index,
        crystal_to_defect_ratio,
        n_steps_per_T,
        n_particles,
        lattice_side,
        run_index,
        allow_p_p_swaps=allow_p_p_swaps,
        record_option=record_option,
        couplings=couplings,
        model_overrides=model_overrides,
        mc_overrides=mc_overrides,
        continue_from_mc_file=continue_from_mc_file,
        continue_from_step=continue_from_step,
    )

    # write_json creates parent directories as needed.
    write_json(model.to_dict(), model_file)
    write_json(mc.to_dict(), mc_file)

    return model_file, mc_file


def write_run_series(
    root: Path,
    run_name:str,
    series_index: int,
    crystal_to_defect_ratio: float,
    n_steps_per_T: int,
    n_particles: int,
    lattice_side:int,
    n_runs: int = 20,
    *,
    allow_p_p_swaps: bool = True,
    record_option: bool = True,
    couplings: list | None = None,
    continue_from_mc_file: str | Path | None = None,
    continue_from_step : int | None = None,
    model_overrides: dict | None = None,
    mc_overrides: dict | None = None,
) -> list[tuple[Path, Path]]:
    """Write parameter files for a series of independent runs."""
    return [
        write_run(
            root,
            run_name,
            series_index,
            crystal_to_defect_ratio,
            n_steps_per_T,
            n_particles,
            lattice_side,
            i,
            allow_p_p_swaps=allow_p_p_swaps,
            record_option=record_option,
            couplings=couplings,
            continue_from_mc_file=continue_from_mc_file,
            continue_from_step=continue_from_step,
            model_overrides=model_overrides,
            mc_overrides=mc_overrides
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
