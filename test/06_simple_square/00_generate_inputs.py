"""Vincent Ouazan-Reboul, 2025/04/04

Make crystal-forming particles with a different number of particles to observe finite-size
effects.
"""

from contact_utils import ContactMapWrapper
from geometry import ParticleGeometry
from json_dump import *
from pathlib import Path
import numpy as np

# e_repel = 100.0
e_crystal = -10.0


def gen_params(n_particles: int):
    ### Key parameter: how we will name this run
    run_name = f"n_particles_{n_particles}"

    ### Define model parameters

    model_params = {}

    # ---------- LATTICE OPTIONS ----------

    # Options:
    # "chain", "square", "triangular", "cubic", "bcc", "fcc"
    model_params["lattice_name"] = "square"

    # Lattice dimensions
    model_params["lx"] = 11
    model_params["ly"] = 11  # Has to be 1 for chain
    model_params["lz"] = 11  # Has to be 1 for square & triangular

    # ---------- MODEL PARAMETERS ----------

    # Number of particle types
    model_params["n_types"] = 1

    # Number of particles of each type
    model_params["n_particles"] = [n_particles]

    # Couplings is its own beast. Should be gotten with the appropriate helper
    # function.

    cu = ContactMapWrapper.from_lattice_name(
        model_params["lattice_name"], init_energy=e_crystal
    )
    model_params["couplings"] = cu.get_formatted_couplings()

    # Initialization option
    model_params["initialize_option"] = "random"

    # If "initialize_option" is set to "from_file", we must specify the location of
    # the input file
    # model_params["state_input"] = str(cfg.structures_path/"final_structure.dat")

    # Options for average and record collection

    model_params["state_av_option"] = False
    model_params["e_av_option"] = True
    model_params["e_record_option"] = False

    if model_params["state_av_option"]:
        state_path = Path("./data/" + run_name + "/average_state")
        state_path.mkdir(parents=True, exist_ok=True)
        model_params["state_av_output"] = str(state_path.resolve()) + "/"

    if model_params["e_av_option"]:
        energy_path = Path("./data/" + run_name + "/average_energy")
        energy_path.mkdir(parents=True, exist_ok=True)
        model_params["e_av_output"] = str(energy_path.resolve()) + "/"

    if model_params["e_record_option"]:
        e_records_path = Path("./data/" + run_name + "/energy_records")
        e_records_path.mkdir(parents=True, exist_ok=True)
        model_params["e_record_output"] = str(e_records_path.resolve()) + "/"

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
    moves_dict["swap_empty_full"] = 1 / 4
    moves_dict["rotate"] = 1 / 4
    moves_dict["rotate_and_swap_w_empty"] = 1 / 4
    moves_dict["swap_full_full"] = 1 / 4

    model_params["move_probas"] = moves_dict

    make_json_file(model_params, f"input/n_particles_{n_particles}_model_params.json")

    ### Define mc parameters

    mc_params = {}

    # Number of MC steps used for equilibration
    mc_params["mcs_eq"] = 10 * model_params["n_particles"][0]

    # Number of MC steps used for averaging
    mc_params["mcs_av"] = 100

    # Type of cooling schedule
    # Options: "exponential", "linear", "inverse"
    # if exponential chosen: specify log10(T) as initial and final temperatures
    # if inverse chosen: specify 1 / T as initial and final temperatures, in decreasing order.
    mc_params["cooling_schedule"] = "exponential"

    # Initial annealing temperature
    mc_params["Ti"] = 1

    # Final annealing temperature
    mc_params["Tf"] = -2

    # Number of annealing steps
    mc_params["Nt"] = 20

    # Option to collect state checkpoints at the end of each temperature cycle
    mc_params["checkpoint_option"] = True

    # If checkpoint is True, we need to provide the output address for the
    # checkpoint files
    structures_path = Path("./data/" + run_name + "/structures")
    structures_path.mkdir(parents=True, exist_ok=True)

    if mc_params["checkpoint_option"]:
        mc_params["checkpoint_address"] = str(structures_path.resolve()) + "/"

    # Output location of the final state configuration (must end with "/")
    mc_params["final_structure_address"] = str(structures_path.resolve()) + "/"

    make_json_file(mc_params, f"input/n_particles_{n_particles}_mc_params.json")


# Create all the possible single-orientation, dimer-forming contact maps and the
# associated input.
n_particles_max = 990
step = 10
n_particles_range = range(1, n_particles_max, step)
for n_particles in n_particles_range:
    gen_params(n_particles)
