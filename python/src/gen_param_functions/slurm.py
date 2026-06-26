"""
SLURM job script generation for frusa_mc runs.
"""

from pathlib import Path
from textwrap import dedent

from .params import load_json


def generate_array_script(
    input_root: Path,
    *,
    job_name: str = "frusa_mc",
    mc_glob: str = "**/mc_params_*.json",
    partition: str = "q-2sem",
    time_limit: str = "336:00:00",
    mem: str = "4gb",
    mail_user: str = "vincent.ouazan-reboul@universite-paris-saclay.fr",
    executable: str = "./frusa_lattice_mc/build/app/frusa_mc",
    log_dir: str = "./logs",
    file_list_path: Path | None = None,
) -> str:
    """
    Build a SLURM array job script for all MC param files under input_root.

    Writes a model_mc_files.txt manifest and returns the script as a string.
    """
    input_root = Path(input_root).resolve()

    mc_files = sorted(input_root.glob(mc_glob))
    if not mc_files:
        raise FileNotFoundError(f"No MC param files matching {mc_glob!r} in {input_root}")

    if file_list_path is None:
        file_list_path = input_root / "model_mc_files.txt"

    with file_list_path.open("w") as f:
        for mc_file in mc_files:
            mc = load_json(mc_file)
            f.write(f"{mc_file.resolve()} {mc['model_params_file']}\n")

    n = len(mc_files)

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
        #SBATCH --nodelist=titan-node[1-4]

        FILES=$(awk -v line="$SLURM_ARRAY_TASK_ID" 'NR==line {{print $1, $2}}' {file_list_path})

        MCFILE=$(echo $FILES | awk '{{print $1}}')
        MODELFILE=$(echo $FILES | awk '{{print $2}}')

        srun {executable} \\
            -M "${{MCFILE}}" \\
            -m "${{MODELFILE}}"
    """)


def generate_script_for_prefix(
    root: Path,
    script_prefix: str,
    **kwargs,
) -> str:
    """
    Convenience wrapper: find all MC files under the input directory
    for a given script prefix, then generate the SLURM script.

    Usage:
        script = generate_script_for_prefix(
            root_folder,
            "00_gen_params_constant_n_low_A",
            job_name="constant_n_low_A",
        )
        Path("run.slurm").write_text(script)
    """
    search_root = root / "input" / script_prefix
    return generate_array_script(search_root, **kwargs)

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
) -> str:
    srun_lines = []
    with manifest_path.open("r") as f:
        test_line = f.readline()
    n_stages = len(test_line.split(" ")) // 2

    for stage_idx in range(n_stages):
        col_model = stage_idx * 2 + 1
        col_mc = col_model + 1
        srun_lines.append(
            f'MCFILE=$(echo $FILES | awk \'{{print ${col_mc}}}\')\n'
            f'MODELFILE=$(echo $FILES | awk \'{{print ${col_model}}}\')\n'
            f'echo "Stage {stage_idx + 1}/{n_stages}: $MCFILE $MODELFILE"\n'
            f'srun {executable} \\\n'
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
#SBATCH --nodelist=titan-node[1-4]

set -euo pipefail

FILES=$(awk -v line="$SLURM_ARRAY_TASK_ID" 'NR==line {{print}}' {manifest_path})

{srun_block}
    """
