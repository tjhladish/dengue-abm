#valgrind --tool=callgrind ./model -randomseed $SEED \
SEED=5503
./model -randomseed $SEED \
        -locfile ./pop-yucatan/locations-yucatan.txt \
        -netfile ./pop-yucatan/network-yucatan.txt \
        -popfile ./pop-yucatan/population-yucatan.txt \
        -probfile ./pop-yucatan/swap_probabilities-yucatan.txt \
        -annualintroscoef 0.033 \
        -normalizeserotypeintros \
        -mosquitomove 0.15 \
        -mosquitomovemodel weighted \
        -mosquitoteleport 0.0 \
        -mosquitodistribution constant \
        -betapm 0.15 \
        -betamp 0.25 \
        -primarypathogenicity 1.0 0.25 1.0 0.25 \
        -secondaryscaling 1.0 1.0 1.0 1.0 \
        -mosquitocapacity 100 \
        -daysimmune 730 \
        -runlength 12785 \
        -simulateannualserotypes \
        -dailyeipfile ./pop-yucatan/merida_eip.out \
        -monthlyoutput \
        -yearlyoutput \
        >  ./testrun/80mos-randomseed$SEED.out #\
        #2> ./testrun/80mos-randomseed$SEED.err
        #-yearlypeopleoutputfile ./junk/people-output-yucatan-multiseason-randomseed$SEED-y \
        #-dailyoutputfile ./junk/daily-output-yucatan-multiseason-randomseed$SEED.csv \
