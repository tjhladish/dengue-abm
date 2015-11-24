serotype=$1

if [[ -n "$serotype" ]]; then

    touch sero${serotype}.log

    for i in `seq 1 10`;
    do
        ./simulate_serotypes abc_config_denv${serotype}.json --process --simulate -n 1000000 --serotype ${serotype} 2>> sero${serotype}.log
    done

    ./simulate_serotypes abc_config_denv${serotype}.json --process 2>> sero${serotype}.log
else
    echo "argument error"
fi
