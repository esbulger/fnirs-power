#!/bin/bash
#SBATCH --job-name=<job_name>
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<email_address>
#SBATCH --output=<log_dir>/<job_name>%A_%a.out
#SBATCH --error=<log_dir>/<job_name>%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=1G
#SBATCH --array=0-<njobs-1>%<n_parallel_jobs> # explanation: X blocks × Y subject-counts = n_jobs → indices 0-<njobs-1>; run n_parallel_jobs in parallel

# Load the required Python environment.
module load YOUR_PYTHON_OR_ANACONDA_MODULE
source /path/to/your/virtual-environment/bin/activate

# Move to the project directory so package imports and relative paths work.
cd /path/to/fnirspower/python

python ../examples/sim_cluster_based_power_parallelized.py