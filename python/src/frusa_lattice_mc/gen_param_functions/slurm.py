"""
SLURM job script generation for frusa_mc runs.
"""

from pathlib import Path
from textwrap import dedent
import numpy as np

from .params import load_json

# Legacy import
from .manifest import write_single_stage as write_stage_manifest, Stage


def generate_array_script(
    input_root: Path,
    jobs: list[Stage],
    *,
    job_name: str = "frusa_mc",
    partition: str = "q-2sem",
    time_limit: str = "336:00:00",
    mem: str = "4gb",
    mail_user: str = "vincent.ouazan-reboul@universite-paris-saclay.fr",
    executable: str = "./src/frusa_lattice_mc/build/app/frusa_mc",
    log_dir: str = "./logs",
    file_list_path: Path | None = None,
    nodelist: str = "1-4",
) -> str:
    """Write a model_mc_files.txt manifest and return the array script."""
    input_root = Path(input_root).resolve()

    if file_list_path is None:
        file_list_path = input_root / "model_mc_files.txt"

    write_stage_manifest(file_list_path, jobs)

    n = len(jobs)

    return dedent(f"""\
        #!/bin/bash

        #SBATCH --partition={partition}
        #SBATCH --job-name={job_name}
        #SBATCH --array=1-{n}
        #SBATCH --mail-type=END,FAIL
        #SBATCH --ntasks=1
        #SBATCH --nodes=1
        #SBATCH --cpus-per-task=1
        #SBATCH --mem={mem}
        #SBATCH --time={time_limit}
        #SBATCH --mail-user={mail_user}
        #SBATCH --output={log_dir}/{job_name}_%A_%a.log
        #SBATCH --nodelist=titan-node[{nodelist}]

        FILES=$(awk -v line="$SLURM_ARRAY_TASK_ID" 'NR==line {{print $1, $2}}' {file_list_path})

        MCFILE=$(echo $FILES | awk '{{print $1}}')
        MODELFILE=$(echo $FILES | awk '{{print $2}}')

        srun {executable} \\
            -M "${{MCFILE}}" \\
            -m "${{MODELFILE}}"
    """)


def generate_self_resuming_script(
    manifest_path: Path,
    n_jobs: int,
    *,
    job_name: str = "frusa_mc",
    partition: str = "q-2sem",
    time_limit: str = "336:00:00",
    mem: str = "4gb",
    mail_user: str = "vincent.ouazan-reboul@universite-paris-saclay.fr",
    runner: str = "./.venv/bin/frusa-mc-task",
    executable: str = "./src/frusa_lattice_mc/build/app/frusa_mc",
    log_dir: str = "./logs",
    nodelist: str = "1-4",
) -> str:
    """
    An array script that starts a campaign and continues it, unchanged.

    Each array task hands one manifest line to frusa-mc-task, which skips the
    stages that finished, continues the one that was interrupted and runs the
    rest. Submit it to start; submit the same file again after a timeout to pick
    up where it stopped; once everything is done it is a no-op. Unlike
    generate_array_script_by_stages, the stage loop lives in the runner, so the
    script does not need regenerating when the number of stages changes.
    """
    return f"""\
#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --job-name={job_name}
#SBATCH --array=1-{n_jobs}
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem={mem}
#SBATCH --time={time_limit}
#SBATCH --mail-user={mail_user}
#SBATCH --output={log_dir}/{job_name}_%A_%a.log
#SBATCH --nodelist=titan-node[{nodelist}]

set -euo pipefail

# Resuming is the default: finished stages are skipped, an interrupted one
# continues from its last checkpoint, anything untouched runs from scratch.
srun {runner} \\
    {Path(manifest_path).resolve()} \\
    "$SLURM_ARRAY_TASK_ID" \\
    --executable {executable}
"""


def generate_array_script_by_stages(
    manifest_path: Path,
    n_jobs: int,
    *,
    job_name: str = "frusa_mc",
    partition: str = "q-2sem",
    time_limit: str = "336:00:00",
    mem: str = "4gb",
    mail_user: str = "vincent.ouazan-reboul@universite-paris-saclay.fr",
    executable: str = "./frusa_lattice_mc/build/app/frusa_mc",
    log_dir: str = "./logs",
    nodelist: str = "1-4",
) -> str:
    srun_lines = []
    with manifest_path.open("r") as f:
        test_line = f.readline()
    n_stages = len(test_line.split()) // 2

    for stage_idx in range(n_stages):
        col_model = stage_idx * 2 + 1
        col_mc = col_model + 1
        srun_lines.append(
            f"MCFILE=$(echo $FILES | awk '{{print ${col_mc}}}')\n"
            f"MODELFILE=$(echo $FILES | awk '{{print ${col_model}}}')\n"
            f'echo "Stage {stage_idx + 1}/{n_stages}: $MCFILE $MODELFILE"\n'
            f"srun {executable} \\\n"
            f'    -M "${{MCFILE}}" \\\n'
            f'    -m "${{MODELFILE}}"'
        )
    srun_block = "\n\n".join(srun_lines)

    return f"""\
#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --job-name={job_name}
#SBATCH --array=1-{n_jobs}
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem={mem}
#SBATCH --time={time_limit}
#SBATCH --mail-user={mail_user}
#SBATCH --output={log_dir}/{job_name}_%A_%a.log
#SBATCH --nodelist=titan-node[{nodelist}]

set -euo pipefail

FILES=$(awk -v line="$SLURM_ARRAY_TASK_ID" 'NR==line {{print}}' {manifest_path})

{srun_block}
    """
