#CC=mpicxx
CC=g++
#CFLAGS = -g -std=c++11
CFLAGS = -O2 -std=c++11 -pedantic -Wall
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
INCLUDES = -I$(ABCDIR) -I$(DENDIR) -I$$TACC_GSL_INC $$HPC_GSL_INC .
SQLDIR = $(ABCDIR)/sqdb

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o 
GSL_LIB = -lm -L$$TACC_GSL_LIB/ -L$$HPC_GSL_LIB/ -lgsl -lgslcblas -lpthread -ldl

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile all_no_mpi

rzero: dengue main.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -I$(SQLDIR) main.cpp -o rzero -L$(DENDIR) $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	rm -f rzero
