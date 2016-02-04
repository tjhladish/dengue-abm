#!/bin/bash
logdir="february2016/auto_output"
echo "seed,year,vax_age,seropositive,pop" > seroprevalence.csv 
grep "seed, year, vaccine cohort age, seroprevalence" $logdir/*.err-* | cut -d' ' -f7 >> seroprevalence.csv
grep "SCENARIO" $logdir/*.err-* | cut -d' ' -f2-16 > seroscenarios.csv
grep "COHORT" $logdir/*.err-* | cut -d, -f2- > cohort.csv
