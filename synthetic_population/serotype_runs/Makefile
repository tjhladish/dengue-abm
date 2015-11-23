CC=mpicxx
CFLAGS = -g -std=c++11
#CFLAGS = -O2 -std=c++11
ABCDIR = $(HOME)/work/AbcSmc

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp   
GSL_LIB = -lm -L$$TACC_GSL_LIB/ -L$$HPC_GSL_LIB/ -lgsl -lgslcblas

libabc_mpi:  
	$(MAKE) -C $(ABCDIR) -f Makefile all_ubuntu_mpi

gen_sero: 
	mpicxx -O2 -std=c++11 generate_serotype_runs.cpp -o gen_sero -lm -lgsl -lgslcblas

sero_mpi: libabc_mpi
	mpicxx -O2 -std=c++11 -D USING_MPI fit_serotype_data.cpp -o simulate_serotypes -I/home/tjhladish/work/AbcSmc -L/home/tjhladish/work/AbcSmc -labc -ljsoncpp -lm -lgsl -lgslcblas

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	rm -f abc abc_mpi
