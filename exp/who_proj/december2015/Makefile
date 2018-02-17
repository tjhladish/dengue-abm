CC=g++
#CC=mpicxx
#CFLAGS = -g -std=c++11
CFLAGS = -O2 -std=c++11
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
IMMDIR = $(HOME)/work/dengue/synthetic_population/initial_immunity
IMMOBJ = $(IMMDIR)/ImmunityGenerator.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$(IMMDIR)
ifdef TACC_GSL_INC
INCLUDE += -I$$TACC_GSL_INC
endif
ifdef HPC_GSL_INC
INCLUDE += -I$$HPC_GSL_INC
endif

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -L$$TACC_GSL_LIB/ -L$$HPC_GSL_LIB/ -lgsl -lgslcblas -lpthread -ldl

default: all_no_mpi

all_no_mpi: CFLAGS += -std=c++0x -Wall
all_no_mpi: libabc_no_mpi abc_sql

libabc_no_mpi:  
	$(MAKE) -C $(ABCDIR) -f Makefile all_no_mpi

dengue_no_mpi:
	$(MAKE) -C $(DENDIR) -f Makefile all_no_mpi

abc_sql: libabc_no_mpi dengue_no_mpi main.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main.cpp -o abc_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

process_who: libabc_no_mpi process_daily_output.cpp
	$(CC) $(CFLAGS) -fopenmp $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc abc_mpi
