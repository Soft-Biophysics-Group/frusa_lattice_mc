"""
Generate a batch of frusa_mc runs on a triangular lattice and the SLURM
array script that runs them.
"""

from pathlib import Path

from frusa_lattice_mc.gen_param_functions import (
    ModelParams,
    MCParams,
    write_run_series,
    calc_ec_ed_ratio,
    camembert_couplings,
    slurm,
)

ROOT = Path(__file__).parent.resolve()
RUN_NAME = "00_triangular_demo"

# Sweep target domain radii, converted to crystal/defect energy ratios.
TARGET_RADII = [3.0, 5.0, 8.0]

mc = MCParams(mcs_eq=int(1e4), Ti=0.0, Tf=2.0, Nt=200)

jobs = []
for series_index, radius in enumerate(TARGET_RADII):
    couplings = camembert_couplings(e_crystal=-18.7, crystal_to_defect_ratio=calc_ec_ed_ratio(radius))
    model = ModelParams(couplings=couplings, n_particles=[500], lx=40, ly=40)
    jobs += write_run_series(
        ROOT,
        RUN_NAME,
        series_index=series_index,
        model=model,
        mc=mc,
        n_runs=4,
        label=f"r_{radius:g}",
    )

input_path = ROOT / "input" / RUN_NAME
print(f"Wrote {len(jobs)} runs under {input_path}")

script = slurm.generate_array_script(input_path, jobs, job_name=RUN_NAME)
script_path = ROOT / f"{RUN_NAME}.slurm"
script_path.write_text(script)
print(f"Wrote SLURM script to {script_path}")
