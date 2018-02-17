#!/bin/bash
#SBATCH --job-name=irs-dengue
#SBATCH --output=./auto_output/refit2-results_%A_%a.out
##SBATCH --output=/ufrc/longini/tjhladish/irs_daily_output-irs_refit2/refit0_%A_%a.out
#SBATCH --error=./auto_output/refit2-results_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi-b
#SBATCH --workdir=/home/tjhladish/work/dengue/irs_timing_refit2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=24:00:00
#SBATCH --array=0-105
#SBATCH --partition=hpg2-compute

module load gcc/5.2.0 gsl
for i in `seq 1 1`;
do
    ./abc_sql abc_irs-summer_winter-foi.json --simulate
done
