#!/bin/bash
#SBATCH --job-name=irs-dengue
#SBATCH --output=./auto_output/refit-stopping-repeat_%A_%a.out
##SBATCH --output=/ufrc/longini/tjhladish/irs_weekly_output2/reg_%A_%a-vv.out
#SBATCH --output=./auto_output/cont_%A_%a.out
#SBATCH --error=./auto_output/cont_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi-b
#SBATCH --workdir=/home/tjhladish/work/dengue/exp/abc-irs_refit2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=24:00:00
#SBATCH --array=0-999
#SBATCH --partition=hpg2-compute

module load gcc/7.3.0 gsl

for i in `seq 1 1`;
do
    ./abc_sql run_posterior.json --simulate
done
