# Copyright (c) 2024 Soft Biophysics Group LPTMS
# Part of frusa_mc, released under BSD 3-Clause License.

import numpy as np
import os
from json_dump import * 

def make_dir(name):
    """
    Creates a new directory if it doesn't already exist at the 
    specified address/name.
    
    Arguments:
    
    name - (str) name and location of the new directory.
    """
    
    if not os.path.exists(name):
        os.makedirs(name)

### Define model parameters

model_params = {}

# Number of particles
model_params["N"] = 20

# Energy gap
model_params["delta"] = 1

# Initialization option
model_params["initialize_option"] = "random"

# If "initialize_option" is set to "from_file", we must specify the location of
# the input file
# model_params["state_input"] = "./..."

# Options for average collection
model_params["state_av_option"] = True
model_params["e_av_option"] = True

if model_params["state_av_option"]:
    make_dir("../data/average_state")
    model_params["state_av_output"] = "./data/average_state/"

if model_params["e_av_option"]:
    make_dir("../data/energy_moments")
    model_params["e_av_output"] = "./data/energy_moments/"

make_json_file(model_params,"../input/model_params.json")

### Define mc parameters

mc_params = {}

# Number of MC steps used for equilibration
mc_params["mcs_eq"] = 100

# Number of MC steps used for averaging
mc_params["mcs_av"] = 1

# Type of cooling schedule
mc_params["cooling_schedule"] = "exponential"

# Initial annealing temperature
mc_params["Ti"] = 0

# Final annealing temperature
mc_params["Tf"] = -5

# Number of annealing steps
mc_params["Nt"] = 10

# Option to collect state checkpoints at the end of each temperature cycle
mc_params["checkpoint_option"] = False

# If checkpoint is True, we need to provide the output address for the 
# checkpoint files
if mc_params["checkpoint_option"]:
    mc_params["checkpoint_address"] = ""

# Output location of the final state configuration (must end with "/")
mc_params["final_structure_address"] = "./"

make_json_file(mc_params,"../input/mc_params.json")
