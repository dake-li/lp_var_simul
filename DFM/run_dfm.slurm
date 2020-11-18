#!/bin/bash
#SBATCH --job-name=lp_var_simul
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:50:00
#SBATCH --array=1-2
#SBATCH --mail-type=end
#SBATCH --mail-user=dakel@princeton.edu

export MATLAB_DGP_TYPE=G
export MATLAB_ESTIMAND_TYPE=ObsShock
export MATLAB_LAG_TYPE=4

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

cd /home/dakel/lp_var_simul/DFM
/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -r run_dfm