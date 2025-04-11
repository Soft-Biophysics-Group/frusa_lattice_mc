import plotting.plot_cubic as pc
from pathlib import Path

def plot_one_run(
    face_index,
):
    pc.clear_blender_fig()

    parent_folder = Path(__file__).resolve().parent
    print(parent_folder)
    run_name = f"face_{face_index}"
    struct_file = parent_folder / f"data/dimer_{run_name}/structures/final_structure.dat"
    model_file = parent_folder / f"input/{run_name}_model_params.json"

    if struct_file.is_file():
        pc.plot_cubes_from_simulation_results(
            struct_file=struct_file,
            model_file=model_file,
            path_to_cube=pc.path_to_numbered_cube ,
        )
        path_to_fig = parent_folder / f"3dFigures/{run_name}.blend"
        pc.save_blender_fig(path_to_fig)

    return


for run_prefix in ("hedgehog_camembert_full_crystal", "hedgehog_camembert_one_axis"):
    for face_index in range(24):
        plot_one_run(face_index)
