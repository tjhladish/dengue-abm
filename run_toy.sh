#valgrind --tool=callgrind ./model -randomseed $SEED \
SEED=5500
./model -randomseed $SEED \
        -locfile ./pop-toy/locations-toy.txt \
        -netfile ./pop-toy/network-toy.txt \
        -popfile ./pop-toy/population-toy.txt \
        -immfile ./pop-toy/immunity-toy.txt \
        -dailyeipfile ./pop-toy/eip-toy.txt \
        -mosquitomove 0.15 \
        -mosquitomovemodel weighted \
        -mosquitoteleport 0.0 \
        -betapm 0.5 \
        -betamp 0.5 \
        -mosquitomultipliers 52 7 0.05 7 0.04 7 0.05 7 0.04 7 0.03 7 0.04 7 0.05 7 0.03 7 0.02 7 0.02 7 0.03 7 0.03 7 0.05 7 0.04 7 0.05 7 0.04 7 0.06 7 0.07 7 0.08 7 0.09 7 0.11 7 0.15 7 0.18 7 0.19 7 0.24 7 0.28 7 0.38 7 0.46 7 0.45 7 0.61 7 0.75 7 0.97 7 0.91 7 1.00 7 0.94 7 0.85 7 0.79 7 0.71 7 0.65 7 0.65 7 0.42 7 0.30 7 0.26 7 0.27 7 0.11 7 0.10 7 0.11 7 0.12 7 0.09 7 0.08 7 0.04 8 0.07 \
        -primarypathogenicity 1.0 0.25 1.0 0.25 \
        -secondaryscaling 1.0 1.0 1.0 1.0 \
        -mosquitocapacity 65 \
        -daysimmune 730 \
        -runlength 365 \
        -dailyexposed 1 1 1 1 \
        -ves 0.7 \
        -yearlypeopleoutputfile ./toy-output/people-output-toy-multiseason-randomseed$SEED-y \
        -dailyoutputfile ./toy-output/daily-output-toy-multiseason-randomseed$SEED.csv \
        > ./toy-output/output-toy-multiseason-randomseed$SEED.csv
