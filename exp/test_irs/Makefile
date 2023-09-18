SHELL=/bin/bash
G++VER := $(shell command -v g++-4.9)

ifndef G++VER
CPP:=g++
else
CPP:=g++-4.9
endif

#CFLAGS = -g -std=c++11 -Wall
CFLAGS = -O2 -std=c++11 -Wall -Wno-deprecated-declarations --pedantic
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue-abm
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

default: libabc abc_sql 

libabc:  
	$(MAKE) -C $(ABCDIR) -f Makefile

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile

abc_sql: libabc dengue main.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main.cpp -o abc_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

process_who: libabc process_daily_output.cpp
	#$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)
	$(CPP) $(CFLAGS) -fopenmp $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc abc_sql
