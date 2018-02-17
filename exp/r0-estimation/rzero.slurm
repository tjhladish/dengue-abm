#!/bin/bash
#SBATCH --job-name=irs-dengue
#SBATCH --output=./auto_output/rzero-sensitivity_%A_%a.out
#SBATCH --error=./auto_output/rzero-sensitivity_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi-b
#SBATCH --workdir=/home/tjhladish/work/dengue/r0-estimation
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=24:00:00
#SBATCH --array=0-1999
#SBATCH --partition=hpg2-compute

module load gcc/5.2.0
for i in `seq 1 10`;
do
    ./rzero abc_rzero.json --simulate
done
