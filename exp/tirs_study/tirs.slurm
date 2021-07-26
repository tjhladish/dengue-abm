#!/bin/bash
#SBATCH --job-name=tirs
#SBATCH --output=./auto_output/foi_impact_%A_%a.out
#SBATCH --error=./auto_output/foi_impact_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --chdir=/home/tjhladish/work/dengue/exp/tirs_study
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=8:00:00
#SBATCH --array=0-249
##SBATCH --partition=hpg2-compute

module load gcc/9.3.0 gsl

for i in `seq 1 150`;
do
    ./abc_sql abc_tirs_foi_exp.json --simulate
done
