#!/bin/bash
#SBATCH --job-name=irs-dengue
#--output=/ufrc/longini/tjhladish/irs_weekly_output2/reg_%A_%a-vv.out
#SBATCH --output=./auto_output/post_%A_%a.out
#SBATCH --error=./auto_output/post_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --workdir=/home/tjhladish/work/dengue/exp/abc-irs_refit
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=24:00:00

#SBATCH --qos=epi
#SBATCH --array=0-56
#SBATCH --partition=hpg2-compute

module load gcc/5.2.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/home/tjhladish/work/AbcSmc/gsl_local/lib:

for i in `seq 1 1`;
do
    ./abc_sql run_posterior.json --simulate
    #./abc_sql abc-irs_refit.json --simulate
done
