#!/bin/bash
rm processed.csv
touch processed.csv
ls *.err-* | xargs cat | grep "seed, year, vaccine cohort age, seroprevalence" >> processed.csv
cut -d: -f2 processed.csv > reprocessed.csv
echo "seed,year,vax_age,seroprevalence" > seroprevalence.csv 
sed 's/^ *//g' < reprocessed.csv >> seroprevalence.csv
