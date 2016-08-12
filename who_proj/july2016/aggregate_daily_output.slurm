#!/bin/bash
#SBATCH --job-name=who-agg
#SBATCH --output=./agg_output/reg_%A_%a.out
#SBATCH --error=./agg_output/reg_%A_%a.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=tjhladish@gmail.com

#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --workdir=/home/tjhladish/work/dengue/who_proj/july2016
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=2:00:00
#SBATCH --array=0-69

module load gcc/5.2.0 gsl
scenario=`sed $PBS_ARRAYID'q;d' scenarios`
./process_who who-80-rev-forecast-attempt2.sqlite /ufrc/longini/tjhladish/who-jul-2016 /ufrc/longini/tjhladish/who-jul-2016-aggregated $scenario
