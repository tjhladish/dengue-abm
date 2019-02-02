-include local.mk

SHELL = /bin/bash
G++VER = $(shell command -v g++-4.9)

ifndef G++VER
CPP = g++
else
CPP = g++-4.9
endif

#CFLAGS = -g -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
CFLAGS = -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue
GSL_PATH = $(HOME)/work/AbcSmc/gsl_local
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
IMMDIR = $(HOME)/work/dengue/synthetic_population/initial_immunity
IMMOBJ = $(IMMDIR)/ImmunityGenerator.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$(IMMDIR) -I$(GSL_PATH)/include/
#ifdef TACC_GSL_INC
#INCLUDE += -I$$TACC_GSL_INC
#endif
#ifdef HPC_GSL_INC
#INCLUDE += -I$$HPC_GSL_INC
#endif

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -lpthread -ldl

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile

abc_sql: libabc dengue main.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main.cpp -o abc_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_sql-foi_gen: libabc dengue main.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main-generate_alt_foi_immunity_files.cpp -o abc_sql-foi_gen $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_sql-ivn_lag: libabc dengue main.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main-ivn_lag.cpp -o abc_sql-ivn_lag $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

process_who: libabc process_daily_output.cpp
	#$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)
	$(CPP) $(CFLAGS) -fopenmp $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc abc_sql
