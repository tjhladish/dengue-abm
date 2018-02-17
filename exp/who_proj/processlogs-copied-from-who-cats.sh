#!/bin/bash
logs="july2016/good_logs"
echo "seed,year,vax_age,seroprevalence" > seroprevalence.csv 
grep "seed, year, vaccine cohort age, seroprevalence" $logs | cut -d' ' -f7 >> seroprevalence.csv
grep "SCENARIO" $logs | cut -d' ' -f2-16 > seroscenarios.csv
grep "COHORT" $logs | cut -d, -f2- > cohort.csv
