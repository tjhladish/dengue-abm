CC=g++-4.9
#CFLAGS = -g -std=c++11
CFLAGS = -O2 -std=c++11
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$(SQLDIR)

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -lgsl -lgslcblas -lpthread -ldl

sero: libabc
	$(CC) $(CFLAGS) $(INCLUDE) -Wno-deprecated-declarations fit_serotype_data.cpp -o simulate_serotypes -I/home/tjhladish/work/AbcSmc $(ABC_LIB) $(GSL_LIB)

generator: yucatan_serotype_generator.h
	$(CC) $(CFLAGS) yucatan_serotype_generator.cpp -o sero_generator $(GSL_LIB)

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile all_no_mpi

gen_sero:
	$(CC) -O2 -std=c++11 generate_serotype_runs.cpp -o gen_sero -lm -lgsl -lgslcblas

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	rm -f abc abc_mpi
