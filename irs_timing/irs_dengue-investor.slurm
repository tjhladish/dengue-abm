#!/bin/bash
#SBATCH --job-name=irs-dengue
#SBATCH --output=./auto_output/refit-stopping_%A_%a.out
##SBATCH --output=/ufrc/longini/tjhladish/irs_daily_output-irs_refit/refit0_%A_%a.out
#SBATCH --error=./auto_output/refit-stopping_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --workdir=/home/tjhladish/work/dengue/irs_timing
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=24:00:00
#SBATCH --array=0-499

module load gcc/5.2.0 gsl
for i in `seq 1 50`;
do
    ./abc_sql abc_stopping-effect.json --simulate
done
