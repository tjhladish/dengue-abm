#valgrind --tool=callgrind ./model -randomseed $SEED \
for ((SEED = 5500; SEED<5510; SEED++))
do
./model -randomseed $SEED \
        -locfile ./pop-yucatan/locations-yucatan.txt \
        -netfile ./pop-yucatan/network-yucatan.txt \
        -popfile ./pop-yucatan/population-yucatan.txt \
        -immfile ./pop-yucatan/immunity-yucatan.txt \
        -probfile ./pop-yucatan/swap_probabilities-yucatan.txt \
        -annualintrosfile ./pop-yucatan/paho_cases \
        -annualintroscoef 1e-6 \
        -mosquitomove 0.15 \
        -mosquitomovemodel weighted \
        -mosquitoteleport 0.0 \
        -betapm 0.1 \
        -betamp 0.25 \
        -mosquitomultipliers 52 7 0.05 7 0.04 7 0.05 7 0.04 7 0.03 7 0.04 7 0.05 7 0.03 7 0.02 7 0.02 7 0.03 7 0.03 7 0.05 7 0.04 7 0.05 7 0.04 7 0.06 7 0.07 7 0.08 7 0.09 7 0.11 7 0.15 7 0.18 7 0.19 7 0.24 7 0.28 7 0.38 7 0.46 7 0.45 7 0.61 7 0.75 7 0.97 7 0.91 7 1.00 7 0.94 7 0.85 7 0.79 7 0.71 7 0.65 7 0.65 7 0.42 7 0.30 7 0.26 7 0.27 7 0.11 7 0.10 7 0.11 7 0.12 7 0.09 7 0.08 7 0.04 8 0.07 \
        -primarypathogenicity 1.0 0.25 1.0 0.25 \
        -secondaryscaling 1.0 1.0 1.0 1.0 \
        -mosquitocapacity 65 \
        -daysimmune 730 \
        -runlength 2920 \
        -dailyexposed 1 1 1 0 \
        -ves 0.7 \
        -yearlypeoplefile ./junk/people-output-yucatan-multiseason-randomseed$SEED-y \
        -dailyfile ./junk/daily-output-yucatan-multiseason-randomseed$SEED.csv \
        > ./junk/output-yucatan-multiseason-randomseed$SEED.csv
done
