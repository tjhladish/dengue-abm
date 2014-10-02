#valgrind --tool=callgrind ./model -randomseed $SEED \
SEED=5501
./model -randomseed $SEED \
        -locfile ./pop-yucatan/locations-yucatan.txt \
        -netfile ./pop-yucatan/network-yucatan.txt \
        -popfile ./pop-yucatan/population-yucatan.txt \
        -immfile ./pop-yucatan/immunity-yucatan.txt \
        -probfile ./pop-yucatan/swap_probabilities-yucatan.txt \
        -annualintroscoef 0.1 \
        -mosquitomove 0.15 \
        -mosquitomovemodel weighted \
        -mosquitoteleport 0.0 \
        -betapm 0.1 \
        -betamp 0.25 \
        -primarypathogenicity 1.0 0.25 1.0 0.25 \
        -secondaryscaling 1.0 1.0 1.0 1.0 \
        -mosquitocapacity 65 \
        -daysimmune 730 \
        -runlength 2920 \
        -simulateannualserotypes \
        -dailyeipfile ./pop-yucatan/merida_eip.out \
        -yearlypeopleoutputfile ./junk/people-output-yucatan-multiseason-randomseed$SEED-y \
        -dailyoutputfile ./junk/daily-output-yucatan-multiseason-randomseed$SEED.csv \
        > ./junk/output-yucatan-multiseason-randomseed$SEED.csv
