"""
Automatically generate the input files for the simulations shown in Fig. 4H.
"""
# -----------------------------------------------------------------------------------------------
# ----- PICKING PARAMETERS -----
# Energy scales
CRYSTAL_ENERGY = -15.0
DEFECT_FACTORS = [
    0.37205221,
    0.47155716,
    0.52151215,
    0.58,
]
REPEL_ENERGY = 100

# GENERATE INPUT
RUN_NAME = "fig_4h_simulations"

# Lattice dimensions
LATTICE_SIZE = 20

# Number of particles
N_PARTICLES = 1728

# ------ PICKING MONTE-CARLO PARAMETERS -----
N_STEPS_PER_T_PER_SITE = N_PARTICLES * 24
N_TEMPERATURE_STEPS = 120



# ----- IMPORTS -----
import numpy as np
from json_dump import *
from pathlib import Path
from designs.fcc import (
    CRYSTAL_INTS_CORNER,
    DEFECT_CORNER
)
from contact_utils import ContactMapWrapper
import config as cfg
from plotting import BlenderPlot
from plotting.plot_3d_utils import create_line_material

# ----- CREATING FOLDERS -----
ROOT_FOLDER = Path(__file__).parent.parent

NPARAMS = len(DEFECT_FACTORS)
NSIMS = NPARAMS

input_folder = ROOT_FOLDER / f"input/{RUN_NAME}"
input_folder.mkdir(exist_ok=True)
data_folder = ROOT_FOLDER / f"data/{RUN_NAME}"
data_folder.mkdir(exist_ok=True, parents=True)

# ----- CREATING INDIVIDUAL INPUT FILES -----
for ec_ed in DEFECT_FACTORS:
    this_run_name = f"{RUN_NAME}_eced_{ec_ed:.4g}"
    defect_energy = CRYSTAL_ENERGY / ec_ed

    model_path = input_folder / (this_run_name + "_model_params.json")
    mc_path = input_folder / (this_run_name + "_mc_params.json")

    # ----- MODEL PARAMETERS -----
    model_params = {}

    # ---------- LATTICE OPTIONS ----------
    # Options:
    # "chain", "square", "triangular", "cubic", "bcc", "fcc"
    model_params["lattice_name"] = "fcc"

    # Lattice dimensions
    model_params["lx"] = LATTICE_SIZE
    model_params["ly"] = LATTICE_SIZE  # Has to be 1 for chain
    model_params["lz"] = LATTICE_SIZE  # Has to be 1 for square & triangular

    # Number of particle types
    model_params["n_types"] = 1

    # Number of particles of each type
    model_params["n_particles"] = [N_PARTICLES]

    # ----- GENERATING ENERGY MATRIX -----
    # Separate class to handle contact energies
    cmap = ContactMapWrapper.from_lattice_name("fcc", 1, REPEL_ENERGY)
    cmap.set_contacts(CRYSTAL_INTS_CORNER, CRYSTAL_ENERGY)
    cmap.set_contacts(DEFECT_CORNER, defect_energy)

    # And process the couplings into flattened form
    model_params["couplings"] = cmap.get_formatted_couplings()
    # Putting the energy scales in the model parameters, because we can and it's convenient
    model_params["crystal_e"] = CRYSTAL_ENERGY
    model_params["defect_e"] = defect_energy
    model_params["repel_e"] = REPEL_ENERGY
    # Initialization option
    model_params["initialize_option"] = "random"

    # If "initialize_option" is set to "from_file", we must specify the location of
    # the input file
    # model_params["state_input"] = str(cfg.structures_path/"final_structure.dat")

    # Options for average collection

    model_params["state_av_option"] = False
    model_params["e_av_option"] = True
    model_params["e_record_option"] = False

    if model_params["e_record_option"]:
        model_params["e_record_output"] = str(data_folder) + "/e_record/"
        Path(model_params["e_record_output"]).resolve().mkdir(
            exist_ok=True, parents=True
        )

    if model_params["e_av_option"]:
        model_params["e_av_output"] = str(data_folder) + "/average_e/"
        Path(model_params["e_av_output"]).resolve().mkdir(exist_ok=True, parents=True)

    # Pick the probabilities of different moves. Has to sum to 1.
    # Options:
    """"
            "swap_empty_full",
            "swap_full_full",
            "rotate",
            "mutate",
            "rotate_and_swap_w_empty"
        """

    moves_dict = {}
    moves_dict["swap_empty_full"] = 1 / 3
    moves_dict["rotate"] = 1 / 3
    moves_dict["rotate_and_swap_w_empty"] = 1 / 3

    model_params["move_probas"] = moves_dict

    model_params["mc_file"] = str(mc_path)

    make_json_file(model_params, model_path)

    # ----- MONTE CARLO PARAMETERS -----
    mc_params = {}

    # Number of MC steps used for equilibration
    mc_params["mcs_eq"] = N_STEPS_PER_T_PER_SITE

    # Number of MC steps used for averaging
    mc_params["mcs_av"] = 10

    # Type of cooling schedule
    # Options: "linear", "inverse", "exponential"
    # if inverse chosen: specify 1/T as initial and final temperatures
    # if exponential chosen: specify log10(T) as initial and final temperatures

    mc_params["cooling_schedule"] = "exponential"

    # Initial annealing temperature
    mc_params["Ti"] = 2

    # Final annealing temperature
    mc_params["Tf"] = 0

    # Number of annealing steps
    mc_params["Nt"] = N_TEMPERATURE_STEPS

    # Option to collect state checkpoints at the end of each temperature cycle
    mc_params["checkpoint_option"] = True

    # If checkpoint is True, we need to provide the output address for the
    # checkpoint files

    if mc_params["checkpoint_option"]:
        mc_params["checkpoint_address"] = str(data_folder) + "/structures/"
        Path(mc_params["checkpoint_address"]).resolve().mkdir(
            exist_ok=True, parents=True
        )

    # Output location of the final state configuration (must end with "/")
    mc_params["final_structure_address"] = mc_params["checkpoint_address"]

    make_json_file(mc_params, mc_path)
