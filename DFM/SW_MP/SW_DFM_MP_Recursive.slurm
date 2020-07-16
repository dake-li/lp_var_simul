#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=71:59:00
 
cd /home/dakel/Codes/DFM/SW_MP
/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -r SW_DFM_MP_Recursive