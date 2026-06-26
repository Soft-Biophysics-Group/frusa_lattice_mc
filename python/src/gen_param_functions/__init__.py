"""Generate frusa_mc run inputs, SLURM scripts, and continuations.

Typical use:

    from gen_param_functions import write_run_series, slurm

    jobs = write_run_series(root, "00_my_run", series_index=0,
                            crystal_to_defect_ratio=0.5, n_steps_per_T=1000,
                            n_particles=500, lattice_side=40, n_runs=20)
    script = slurm.generate_script_for_prefix(root, "00_my_run", job_name="my_run")
"""

from . import slurm, continuation
from .params import (
    ModelParams,
    MCParams,
    InitializeOptions,
    CoolingSchedules,
    calc_ec_ed_ratio,
    camembert_couplings,
    build_params,
    write_run,
    write_run_series,
    write_manifest,
    load_json,
    write_json,
)
from .paths import run_slug

__all__ = [
    "slurm",
    "continuation",
    "ModelParams",
    "MCParams",
    "InitializeOptions",
    "CoolingSchedules",
    "calc_ec_ed_ratio",
    "camembert_couplings",
    "build_params",
    "write_run",
    "write_run_series",
    "write_manifest",
    "load_json",
    "write_json",
    "run_slug",
]
