#valgrind --tool=callgrind ./model -randomseed 5489 \
./model -randomseed 5489 \
        -locfile locations-bangphae.txt \
        -netfile network-bangphae.txt \
        -popfile population-bangphae.txt \
        -immfile immunity-bangphae.txt \
        -probfile swap_probabilities-bangphae.txt \
        -mosquitomove 0.15 \
        -mosquitomovemodel weighted \
        -mosquitoteleport 0.0 \
        -betapm 0.1 \
        -betamp 0.25 \
        -mosquitomultipliers 12 31 0.11 28 0.09 31 0.20 30 0.30 31 0.83 30 1.00 31 0.50 31 0.37 30 0.38 31 0.26 30 0.30 31 0.19 \
        -primarypathogenicity 1.0 0.25 1.0 0.25 \
        -secondaryscaling 1.0 1.0 1.0 1.0 \
        -mosquitocapacity 42 \
        -daysimmune 120 \
        -runlength 730 \
        -dailyexposed 2 2 2 2 \
        -ves 0.7 \
        -yearlypeoplefile ./test/people-output-bangphae-multiseason-randomseed5489-y \
        -dailyfile ./test/daily-output-bangphae-multiseason-randomseed5489.csv \
        > ./test/output-bangphae-multiseason-randomseed5489.csv
