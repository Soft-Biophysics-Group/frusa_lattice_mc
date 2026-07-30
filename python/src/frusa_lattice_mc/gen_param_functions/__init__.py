"""Generate frusa_mc run inputs, SLURM scripts, and continuations.

Typical use:

    from gen_param_functions import ModelParams, MCParams, write_run_series, slurm

    run_name = "00_my_run"
    model = ModelParams(couplings=couplings, n_particles=[500], lx=40, ly=40)
    mc = MCParams(mcs_eq=1000, Ti=0.0, Tf=2.0, Nt=200)
    jobs = write_run_series(root, run_name, series_index=0,
                            model=model, mc=mc, n_runs=20)
    script = slurm.generate_script_for_prefix(root, run_name, job_name="my_run")
"""

from . import slurm, continuation, manifest
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
    "manifest",
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
