#!/bin/bash

#SBATCH --job-name=sarabande_full_params
#SBATCH --output=full_params.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jamessunseri@berkeley.edu
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH  --account=zslepian
#SBATCH  --qos=zslepian

module load conda
conda activate speed_py
python param_testing.py 4 full
echo Finished with 4PCF!

