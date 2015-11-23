mpiexec -n 10 simulate_serotypes abc_config_denv1.json 2> sero.log 
echo "1 done\n"
mpiexec -n 10 simulate_serotypes abc_config_denv2.json 2>> sero.log 
echo "2 done\n"
mpiexec -n 10 simulate_serotypes abc_config_denv3.json 2>> sero.log 
echo "3 done\n"
mpiexec -n 10 simulate_serotypes abc_config_denv4.json 2>> sero.log 
echo "4 done\n"
