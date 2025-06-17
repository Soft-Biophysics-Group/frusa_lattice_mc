from plotting import BlenderPlot
from pathlib import Path

def plot_one_run(
    n_particles,
):

    parent_folder = Path(__file__).resolve().parent
    print(parent_folder)
    run_name = f"n_particles_{n_particles}"
    struct_file = parent_folder / f"data/{run_name}/structures/final_structure.dat"
    model_file = parent_folder / f"input/{run_name}_model_params.json"

    bp = BlenderPlot.from_model_file(model_file)

    print(bp.lattice_geometry)

    # if struct_file.is_file():
    #     bp.plot_particles_from_simulation_results(struct_file=struct_file)
    #     path_to_fig = parent_folder / f"3dFigures/{run_name}.blend"
    #     bp.save(path_to_fig)

    return


for face_index in range(1, 1001, 50):
    plot_one_run(face_index)
