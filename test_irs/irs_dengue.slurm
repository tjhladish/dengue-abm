#!/bin/bash
#SBATCH --job-name=who-dengue
#SBATCH --output=./auto_output/reg_%A_%a.out
#SBATCH --error=./auto_output/reg_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --workdir=/home/tjhladish/work/dengue/test_irs
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=5:00:00
#SBATCH --array=0-499

module load gcc/5.2.0 gsl
for i in `seq 1 3`;
do
    ./abc_sql abc_config.json --simulate
done
