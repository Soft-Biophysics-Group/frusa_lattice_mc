#! /usr/bin/env python3
# Vincent Ouazan-Reboul, 2024.05.14
# Set of utilities which the Python code will have to use frequently
# To be loaded before doing anything else

from pathlib import Path
import json
import numpy as np
import subprocess
import sys

# Reliably make absolute paths to the right place.
# #ILovePathLib <3 <3 <3
parent_path = Path(__file__).parent.parent.parent.absolute().resolve()

# ----- OFTEN-USED PATHS -----
input_path = parent_path/"input/"
default_model_params_file = input_path/"model_params.json"
default_mc_params_file = input_path/"mc_params.json"

data_path = parent_path/"data/"
averages_path = data_path/"average_state/"
energy_path = data_path/"energy_moments/"
structures_path = data_path/"structures/"

exec_path = parent_path/"build/app/frusa_mc"

python_path = parent_path/"python"

# ----- STANDARDIZED FUNCTIONS TO LOAD FILES -----
def load_model_file(model_file=default_model_params_file):
    model_file_str = str(model_file)
    with open(model_file_str, "r") as f:
        params = json.load(f)
    return params


def load_mc_file(mc_file=default_mc_params_file):
    mc_file_str = str(mc_file)
    with open(mc_file_str, "r") as f:
        params = json.load(f)
    return params


def load_structure(
    struct_index: int = -1, struct_folder: str | Path = "", struct_file: str | Path = ""
):
    """
    Load a configuration of the particles on the lattice recorded during the simulation.
    Default behavior (no parameters): loads the structure at the end of the simulation.
    If struct_index is specified: load structure with index struct_index.
    If struct_file is specified: overrides struct_index, directly fetches results in
    struct_file.
    Returned structure is a dimensional numpy array, with line 1 corresponding to particle type
    and line 2 to particle orientation.
    orientation -1 always corresponds to an empty site, irrespective of particle type.
    """
    if struct_folder == "":
        struct_folder = structures_path
    struct_folder = Path(struct_folder)
    if struct_file != "":
        file_path = struct_file
    if struct_index != -1:
        file_path = struct_folder / f"structure_{struct_index}.dat"
    else:
        file_path = struct_folder / "final_structure.dat"
    return np.loadtxt(file_path, dtype=int)

def get_full_sites(site_orientations):
    return np.where(site_orientations != -1)[0]


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
    subprocess.run([str(exec_path), "-m", str(model_file), "-M", str(mc_file)])
