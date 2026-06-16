#!/bin/bash
#SBATCH --job-name=plot_power
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --output=/path/to/logs/plot_power_%A.out
#SBATCH --error=/path/to/logs/plot_power_%A.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

# Load the required Python environment.
module load YOUR_PYTHON_OR_ANACONDA_MODULE
source /path/to/your/virtual-environment/bin/activate

# Move to the project directory so package imports and relative paths work.
cd /path/to/fnirspower/python

python ../examples/sim_cluster_based_power_process.py