#!/bin/bash
#SBATCH --job-name=mpeak-dengue
#SBATCH --output=./auto_output/mpeak_%A_%a.out
#SBATCH --error=./auto_output/mpeak_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi-b
#SBATCH --workdir=/home/tjhladish/work/dengue/exp/mpeak_incidence_response
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=24:00:00
#SBATCH --array=0-2999
#SBATCH --partition=hpg2-compute

module load gcc/7.3.0 gsl

for i in `seq 1 20`;
do
    ./abc_sql abc_mpeak_incidence_response.json --simulate
done
