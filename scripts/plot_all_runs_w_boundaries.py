"""
Plots all the simulation results associated with the runs whose input files are located in the
target directory.

Usage: python plot_all_runs_w_boundaries.py {path_to_input_folder/}. e.g: python
plot_all_runs_w_boundaries input/notebook_run.
"""

from plotting import BlenderPlot
import config as cfg
from pathlib import Path
from designs.fcc import DEFECT_CORNER, CRYSTAL_CONTACTS

import sys

from plotting.plot_3d_utils import create_line_material
from typing import cast


class ModelParamsWMc(cfg.ModelParams):
    mc_file: str


def plot_all_runs(
    input_folder,
):
    parent_folder = Path(__file__).resolve().parent.parent
    print("\n")
    input_files = Path(input_folder).glob("*_model_params.json")

    for file in input_files:
        model_file = file.resolve()
        model_params = cast(ModelParamsWMc, cfg.load_model_file(model_file))

        mc_file = model_params["mc_file"]
        mc_params = cfg.load_mc_file(mc_file)

        structure_folder = Path(mc_params["final_structure_address"])

        # Hacky way to get the run name
        run_name = str(
            structure_folder.parent.relative_to(structure_folder.parent.parent)
        )
        print(run_name)

        bp = BlenderPlot.from_model_file(model_file)
        bp.load_particle(bp.particle_paths["path_to_one_axis_diagonal_colored_gray"])

        path_to_fig = parent_folder / f"3dFigures/{run_name}_w_defect.blend"
        print(path_to_fig)

        bp.plot_canonical_style(
            mc_file,
            path_to_fig,
            DEFECT_CORNER,
            crystal_contacts=CRYSTAL_CONTACTS,
            cubified=False,
            centered=True,
        )

    return


input_folder = sys.argv[1]
plot_all_runs(input_folder)
