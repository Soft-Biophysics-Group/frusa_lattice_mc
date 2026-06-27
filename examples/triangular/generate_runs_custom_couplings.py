"""
Generate a batch of frusa_mc runs on a triangular lattice with an arbitrary
coupling map, built with `contact_utils.ContactMapWrapper`.
"""

from pathlib import Path

from contact_utils import ContactMapWrapper
from gen_param_functions import ModelParams, MCParams, write_run_series, slurm

ROOT = Path(__file__).parent.resolve()
RUN_NAME = "01_custom_couplings_demo"


def make_couplings(e_bond: float, e_repel: float) -> list:
    """Single species: bind each face to the one across the hexagon, repel otherwise."""
    cmap = ContactMapWrapper.triangular(1, init_energy=e_repel)
    for face in range(3):
        cmap[face, face + 3] = e_bond
    return cmap.get_formatted_couplings()


model = ModelParams(
    couplings=make_couplings(e_bond=-15.0, e_repel=10.0),
    n_particles=[500],
    lx=40,
    ly=40,
)
mc = MCParams(mcs_eq=int(1e4), Ti=0.0, Tf=2.0, Nt=200)

jobs = write_run_series(ROOT, RUN_NAME, series_index=0, model=model, mc=mc, n_runs=4)
print(f"Wrote {len(jobs)} runs under {ROOT / 'input' / RUN_NAME}")

script = slurm.generate_script_for_prefix(ROOT, RUN_NAME, job_name=RUN_NAME)
script_path = ROOT / f"{RUN_NAME}.slurm"
script_path.write_text(script)
print(f"Wrote SLURM script to {script_path}")
