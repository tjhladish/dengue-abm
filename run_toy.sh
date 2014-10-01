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
