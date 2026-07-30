#! /usr/bin/env python3
# Vincent Ouazan-Reboul, 2024.05.14
# Set of utilities which the Python code will have to use frequently
# To be loaded before doing anything else

from pathlib import Path
import subprocess
import numpy as np
from numpy.typing import NDArray
from .json_dump import load_json
import os

from typing import TypedDict, cast

# Reliably make absolute paths to the right place.
# #ILovePathLib <3 <3 <3
REPO_ROOT = Path(__file__).parents[3]      # …/frusa_lattice_mc
parent_path = Path(__file__).parent.parent.parent.absolute().resolve()

# ----- OFTEN-USED PATHS -----
input_path = parent_path/"input/"
default_model_params_file: Path = input_path / "model_params.json"
default_mc_params_file: Path = input_path/"mc_params.json"

data_path = parent_path/"data/"
averages_path = data_path/"average_state/"
energy_path = data_path/"energy_moments/"
structures_path = data_path/"structures/"

exec_path = REPO_ROOT/"build/app/frusa_mc"

python_path = parent_path/"python"


# ----- STANDARDIZED FUNCTIONS TO LOAD FILES -----

# Better way to structure the json inputs / outputs in Python:
# this way your type checker knows what to expect. It does not restrict what you do at runtime
# though!
class ModelParams(TypedDict, total=False):
    lattice_name: str

    lx: int
    ly: int
    lz: int

    n_types:int
    n_particles: list[int]

    couplings: list[float]

    initialize_option:str
    state_av_option:bool
    e_av_option: bool
    e_record_option: bool

    e_record_output: str
    e_av_output: str

    move_probas: dict[str, float]

def load_model_file(model_file: str | Path = default_model_params_file) -> ModelParams:
    return cast(ModelParams, load_json(model_file))


class McParams(TypedDict):
    mcs_eq: int
    mcs_av: int

    cooling_schedule:str
    Ti: float | int
    Tf: float | int
    Nt: int

    checkpoint_option: bool
    checkpoint_address: str
    final_structure_address: str

def load_mc_file(
    mc_file: str | Path = default_mc_params_file,
) -> McParams:
    return cast(McParams, load_json(mc_file))


# ----- STATE HELPERS -----
# Superseded by LatticeState (see lattice_state.py). Kept until the remaining callers in
# plotting/ and analysis/ have migrated; delete once nothing imports them.


def load_structure(
    struct_index: int | None = None,
    struct_folder: str | Path | None = None,
    struct_file: str | Path | None = None,
) -> NDArray[np.int_]:
    """
    Load a configuration of the particles on the lattice recorded during the simulation.
    Default behavior (no parameters): loads the structure at the end of the simulation.
    If struct_index is specified: load structure with index struct_index.
    If struct_file is specified: overrides struct_index, directly fetches results in
    struct_file.
    Returned structure is a (2, Nsites) dimensional numpy array, with line 1 corresponding to particle type
    and line 2 to particle orientation.
    orientation -1 always corresponds to an empty site, irrespective of particle type.
    """
    if struct_file is not None:
        file_path = struct_file
    else:
        if struct_folder is None:
            struct_folder = structures_path
        struct_folder = Path(struct_folder)
        if struct_index is not None:
            file_path = struct_folder / f"structure_{struct_index}.dat"
        else:
            file_path = struct_folder / "final_structure.dat"

    return np.loadtxt(file_path, dtype=int)


def get_full_sites(site_orientations: NDArray[np.int_]) -> NDArray[np.int_]:
    return np.where(site_orientations[1, :] != -1)[0]


def get_full_sites_characteristics(site_orientations):
    sites = get_full_sites(site_orientations)
    types = site_orientations[0, sites]
    orientations = site_orientations[1, sites]
    return np.vstack([sites, types, orientations]).T


# ----- RUN SIMULATIONS FROM PYTHON  -----


def check_data_existence(
    struct_path: str | Path = structures_path, e_path: str | Path = energy_path
):
    struct_path = Path(struct_path)
    e_path = Path(e_path)
    files_exist = False

    if struct_path.is_dir():
        if any(struct_path.iterdir()):
            files_exist = True

    if e_path.is_dir():
        if any(e_path.iterdir()):
            files_exist = True

    return files_exist


def find_executable(explicit: str | Path | None = None) -> Path:
    """Binary location: argument, then $FRUSA_MC_BIN, then the repo's build/."""
    if explicit is not None:
        return Path(explicit)
    if env := os.environ.get("FRUSA_MC_BIN"):
        return Path(env)
    candidate = REPO_ROOT / "build/app/frusa_mc"
    if candidate.is_file():
        return candidate
    raise FileNotFoundError("No frusa_mc binary — build it or set FRUSA_MC_BIN.")


def run_simulation(
    model_file=default_model_params_file,
    mc_file=default_mc_params_file,
    overwrite=False,
):
    """
    Starts the frusa_mc program with default parameters from the root directory of the project.
    Should allow you to run the code smoothly from a python interpreter in which this file has
    been imported.

    Keyword arguments:
    - `model_file`: string or path. json model parameter file. Defaults to `input/model_params.json` from the
      root of the program.
    - `mc_file`: string or path. json Monte Carlo parameter file. Defaults to `input/mc_params.json` from the
      root of the program.
    - `overwrite`: boolean, `False` by default. If false, does not overwrite data directory
      contents if they're populated.
    """
    # execute(str(exec_path), cwd = parent_path)
    mc_params = load_mc_file(mc_file)
    struct_path = mc_params["checkpoint_address"]
    model_params = load_model_file(model_file)
    e_path = model_params["e_av_output"]
    if check_data_existence(struct_path=struct_path, e_path=e_path):
        print("At least one of the output folders is already populated!")
        if not overwrite:
            print("overwrite flag set to False: exiting.")
            return
        else:
            print("overwrite flag set to True: running anyway.")
    subprocess.run(
        [str(find_executable()), "-m", str(model_file), "-M", str(mc_file)]
    )
