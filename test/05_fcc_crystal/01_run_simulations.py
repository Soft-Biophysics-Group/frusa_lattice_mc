"""
Vincent Ouazan-Reboul, 2025

Utility to run several MC simulations in parallel.
"""

from multiprocessing import Pool
from config import run_simulation

overwrite = True

def do_all_runs(overwrite=False, n_processes=6):
    full_run_inputs = (
        (
            f"./input/n_particles_{n_particles}_model_params.json",
            f"./input/n_particles_{n_particles}_mc_params.json",
            overwrite,
        )
        for n_particles in range(1, 1001, 10)
    )

    with Pool(processes=n_processes) as p:
        print("Got here")
        p.starmap(run_simulation, full_run_inputs)

    return

if __name__ == '__main__':
    print("Running all")
    do_all_runs(overwrite)
