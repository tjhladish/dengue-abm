#!/bin/bash
rm processed.csv
rm scenarios.csv
touch processed.csv
touch scenarios.csv
ls *.err-* | xargs cat | grep "seed, year, vaccine cohort age, seroprevalence" >> processed.csv
ls *.err-* | xargs cat | grep "SCENARIO" >> scenarios.csv
cut -d' ' -f2-17 scenarios.csv > seroscenarios.csv
cut -d: -f2 processed.csv > reprocessed.csv
echo "seed,year,vax_age,seroprevalence" > seroprevalence.csv 
sed 's/^ *//g' < reprocessed.csv >> seroprevalence.csv
