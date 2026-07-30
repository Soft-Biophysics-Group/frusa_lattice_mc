"""Vincent Ouazan-Reboul, 2026-07-30. Co-written by Claude Opus 5.

Generate a handful of two-stage runs for the local runner test.
"""

import shutil
from pathlib import Path

from frusa_lattice_mc.contact_utils import ContactMapWrapper
from frusa_lattice_mc.designs.triangular import (
    CRYSTAL_CONTACTS,
    VORTEX_CAMEMBERT_CONTACTS,
)
from frusa_lattice_mc.gen_param_functions import (
    CoolingSchedules,
    MCParams,
    ModelParams,
    manifest,
    write_run,
)

ROOT = Path(__file__).parent.resolve()
MANIFEST_FILE = ROOT / "input" / "staged" / "manifest.txt"

N_RUNS = 2
CRYSTAL_ENERGY = -18.7
REPEL_ENERGY = 10.0
EC_ED_RATIO = 0.4792


def mc_params() -> MCParams:
    """Short enough to finish in a second, long enough to leave checkpoints."""
    return MCParams(
        mcs_eq=20,
        mcs_av=10,
        cooling_schedule=CoolingSchedules.LINEAR,
        Ti=1.15,
        Tf=1.150001,
        Nt=6,
    )


def main() -> None:
    for stale in (ROOT / "input", ROOT / "data"):
        shutil.rmtree(stale, ignore_errors=True)

    cmap = ContactMapWrapper.triangular(1, init_energy=REPEL_ENERGY)
    cmap.set_contacts(CRYSTAL_CONTACTS, CRYSTAL_ENERGY)
    cmap.set_contacts(VORTEX_CAMEMBERT_CONTACTS, CRYSTAL_ENERGY / EC_ED_RATIO)
    model = ModelParams(
        "triangular",
        lx=10,
        ly=10,
        n_particles=[25],
        couplings=cmap.get_formatted_couplings(),
    )

    jobs = []
    for run_index in range(N_RUNS):
        stage_1 = write_run(ROOT, "short", 0, run_index, model, mc_params(), label="lbl")
        stage_2 = write_run(
            ROOT,
            "long",
            0,
            run_index,
            model,
            mc_params(),
            label="lbl",
            continue_from_mc_file=stage_1[1],
        )
        jobs.append([stage_1, stage_2])

    manifest.write_stages(MANIFEST_FILE, jobs)
    print(f"Wrote {len(jobs)} two-stage runs; manifest at {MANIFEST_FILE}")


if __name__ == "__main__":
    main()
