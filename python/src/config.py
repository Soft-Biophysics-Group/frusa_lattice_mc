#! /usr/bin/env python3
# Vincent Ouazan-Reboul, 2024.05.14
# Set of utilities which the Python code will have to use frequently

from pathlib import Path
import json
import numpy as np
from subprocess import run, Popen

# Reliably make absolute paths to the right place.
# #ILovePathLib <3 <3 <3
parent_path = Path(__file__).parent.parent.parent.absolute().resolve()

input_path = parent_path/"input/"
default_model_params_file = input_path/"model_params.json"
default_mc_params_file = input_path/"mc_params.json"

data_path = parent_path/"data/"
averages_path = data_path/"average_state/"
energy_path = data_path/"energy_moments/"
structures_path = data_path/"structures/"

exec_path = parent_path/"build/app/frusa_mc"

##### STANDARDIZED FUNCTIONS TO LOAD FILES #####
def load_model_file(model_file = default_model_params_file):
    model_file_str = str(model_file)
    with open(model_file_str, "r") as f:
        params = json.load(f)
    return params

def load_structure(struct_index:int = -1, struct_file:str|Path = ""):
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
    if struct_file != "":
        file_path = struct_file
    if struct_index != -1:
        file_path = structures_path / f"structure_{struct_index}.dat"
    else:
        file_path = structures_path / "final_structure.dat"
    return np.loadtxt(file_path, dtype = int)

##### RUN SIMULATIONS FROM PYTHON  #####
def run_simulation():
    Popen(str(exec_path), cwd = parent_path)
