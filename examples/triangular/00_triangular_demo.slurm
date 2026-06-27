#!/bin/bash

#SBATCH --partition=q-2sem
#SBATCH --job-name=00_triangular_demo
#SBATCH --array=1-12
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=336:00:00
#SBATCH --mail-user=vincent.ouazan-reboul@universite-paris-saclay.fr
#SBATCH --output=./logs/00_triangular_demo_%A_%a.log
#SBATCH --nodelist=titan-node[1-4]

FILES=$(awk -v line="$SLURM_ARRAY_TASK_ID" 'NR==line {print $1, $2}' /Users/vincent/research/projects/2606_camembertEntropicControl/src/frusa_lattice_mc/python/examples/triangular/input/00_triangular_demo/model_mc_files.txt)

MCFILE=$(echo $FILES | awk '{print $1}')
MODELFILE=$(echo $FILES | awk '{print $2}')

srun ./frusa_lattice_mc/build/app/frusa_mc \
    -M "${MCFILE}" \
    -m "${MODELFILE}"
