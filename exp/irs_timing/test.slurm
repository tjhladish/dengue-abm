#!/bin/bash
#SBATCH --job-name=test
#SBATCH --error=./auto_output/test_%A_%a-vv.err
#SBATCH --output=./auto_output/test_%A_%a-vv.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --workdir=/home/tjhladish/work/dengue/exp/irs_timing
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=0:10:00

module load gcc/5.2.0 gsl
touch asdf
