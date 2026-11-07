"""Generate frusa_mc run inputs, SLURM scripts, and continuations.

Typical use:

    from gen_param_functions import ModelParams, MCParams, write_run_series, slurm

    run_name = "00_my_run"
    model = ModelParams(couplings=couplings, n_particles=[500], lx=40, ly=40)
    mc = MCParams(mcs_eq=1000, Ti=0.0, Tf=2.0, Nt=200)
    stages = write_run_series(root, run_name,
                              model=model, mc=mc, n_runs=20)
    script = slurm.generate_array_script(input_root, stages, job_name="my_run")
"""

from . import continuation, manifest, slurm
from .manifest import Stage, read_manifest, run_dir, write_manifest
from .params import (
    CoolingSchedules,
    InitializeOptions,
    MCParams,
    ModelParams,
    build_params,
    calc_ec_ed_ratio,
    camembert_couplings,
    load_json,
    write_json,
    write_multi_stage_run,
    write_multi_stage_run_series,
    write_run,
    write_run_series,
)

__all__ = [
    "CoolingSchedules",
    "InitializeOptions",
    "MCParams",
    "ModelParams",
    "Stage",
    "build_params",
    "calc_ec_ed_ratio",
    "camembert_couplings",
    "continuation",
    "load_json",
    "manifest",
    "read_manifest",
    "run_dir",
    "slurm",
    "write_json",
    "write_manifest",
    "write_multi_stage_run",
    "write_multi_stage_run_series",
    "write_run",
    "write_run_series",
]
