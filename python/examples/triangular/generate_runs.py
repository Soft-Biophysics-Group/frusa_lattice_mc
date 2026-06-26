"""
Generate a batch of frusa_mc runs on a triangular lattice and the SLURM
array script that runs them.
"""

from pathlib import Path

from gen_param_functions import write_run_series, calc_ec_ed_ratio, slurm

# Where input/ and data/ trees get written.
ROOT = Path(__file__).parent.resolve()

# Names this batch; becomes the input/<RUN_NAME>/ subtree.
RUN_NAME = "00_triangular_demo"

# Sweep over a few target domain radii, converted to crystal/defect energy ratios.
TARGET_RADII = [3.0, 5.0, 8.0]

# One (model, mc) pair per run, grouped by sweep point — the manifest/SLURM input.
jobs = []
for series_index, radius in enumerate(TARGET_RADII):
    ratio = calc_ec_ed_ratio(radius)
    jobs += write_run_series(
        ROOT,
        RUN_NAME,
        series_index=series_index,
        crystal_to_defect_ratio=ratio,
        n_steps_per_T=int(1e4),
        n_particles=500,
        lattice_side=40,
        n_runs=4,
        mc_overrides={"Ti": 0.0, "Tf": 2.0, "Nt": 200},
    )

print(f"Wrote {len(jobs)} runs under {ROOT / 'input' / RUN_NAME}")

# Emit the SLURM array script covering every run just written.
script = slurm.generate_script_for_prefix(ROOT, RUN_NAME, job_name=RUN_NAME)
script_path = ROOT / f"{RUN_NAME}.slurm"
script_path.write_text(script)
print(f"Wrote SLURM script to {script_path}")
