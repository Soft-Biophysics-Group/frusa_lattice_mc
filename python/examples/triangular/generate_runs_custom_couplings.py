"""
Generate a batch of frusa_mc runs on a triangular lattice using an *arbitrary*
coupling map, rather than the built-in camembert map.

The coupling map is built with `contact_utils.ContactMapWrapper` and handed to
`write_run_series` through its `couplings=` argument. Anything that produces a
flattened coupling list the C++ side accepts works the same way.
"""

from pathlib import Path

from contact_utils import ContactMapWrapper
from gen_param_functions import write_run_series, slurm

# Where input/ and data/ trees get written.
ROOT = Path(__file__).parent.resolve()

# Names this batch; becomes the input/<RUN_NAME>/ subtree.
RUN_NAME = "01_custom_couplings_demo"


def make_couplings(e_bond: float, e_repel: float) -> list:
    """A simple single-species map: every face binds its opposite attractively,
    everything else repels. Illustrates building an arbitrary map by hand."""
    cmap = ContactMapWrapper.triangular(1, init_energy=e_repel)
    # Faces 0..5 around the hexagon; bind each face to the one across from it.
    for face in range(3):
        cmap[face, face + 3] = e_bond
    return cmap.get_formatted_couplings()


couplings = make_couplings(e_bond=-15.0, e_repel=10.0)

# `crystal_to_defect_ratio` no longer drives the physics here — it only labels the
# output directory (input/<RUN_NAME>/<series>_ec_div_ed_<ratio>/), so pass a tag.
jobs = write_run_series(
    ROOT,
    RUN_NAME,
    series_index=0,
    crystal_to_defect_ratio=0.0,
    n_steps_per_T=int(1e4),
    n_particles=500,
    lattice_side=40,
    n_runs=4,
    couplings=couplings,
    mc_overrides={"Ti": 0.0, "Tf": 2.0, "Nt": 200},
)

print(f"Wrote {len(jobs)} runs under {ROOT / 'input' / RUN_NAME}")

# Emit the SLURM array script covering every run just written.
script = slurm.generate_script_for_prefix(ROOT, RUN_NAME, job_name=RUN_NAME)
script_path = ROOT / f"{RUN_NAME}.slurm"
script_path.write_text(script)
print(f"Wrote SLURM script to {script_path}")
