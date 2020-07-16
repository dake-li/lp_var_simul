#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=71:59:00
 
cd /home/dakel/Codes/DFM/SW_G
/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -r SW_DFM_G_ObsShock