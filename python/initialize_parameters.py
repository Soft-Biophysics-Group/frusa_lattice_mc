# Copyright (c) 2024 Soft Biophysics Group LPTMS
# Part of frusa_mc, released under BSD 3-Clause License.

import numpy as np
import json

def make_json_file(data,address):
    """
    Writes the data, provided as a python dictionary, to a .json file
    """
    with open(address,"w") as write_file:
        json.dump(data, write_file, indent=2)

#Define model parameters

model_params = {}

#Lattice dimensions
model_params["Lx"] = 100
model_params["Ly"] = 1
model_params["Lz"] = 1

#Number of particles
model_params["Np"] = 20

#Model parameters
k11 = -1
k12 = 0
k21 = 0
T_model = 0.1
model_params["couplings"] = [k11,k12,k21,T_model]

make_json_file(model_params,"../input/model_params.json")

#Define mc parameters

mc_params = {}

#Number of MC steps used for equilibration
mc_params["mcs_eq"] = 100

#Number of MC steps used for averaging
mc_params["mcs_av"] = 1

#Type of cooling schedule
mc_params["cooling_schedule"] = "exponential"

#Initial annealing temperature
mc_params["Ti"] = 0

#Final annealing temperature
mc_params["Tf"] = -5

#Number of annealing steps
mc_params["Nt"] = 10

#Option to collect state checkpoints at the end of each temperature cycle
mc_params["checkpoint_option"] = False

#If checkpoint is True, we need to provide the output address for the 
#checkpoint files
if mc_params["checkpoint_option"]:
    mc_params["checkpoint_address"] = ""

make_json_file(mc_params,"../input/mc_params.json")
