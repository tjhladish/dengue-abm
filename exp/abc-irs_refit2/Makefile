SHELL=/bin/bash
G++VER := $(shell command -v g++-4.9)

ifndef G++VER
CXX:=g++
else
CXX:=g++-4.9
endif

#CFLAGS = -g -std=c++11 -pedantic -Wall -Wextra -Wno-deprecated-declarations
CFLAGS = -O2 -std=c++11 -pedantic -Wall -Wextra -Wno-deprecated-declarations
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue-abm
GSL_PATH = $(HOME)/work/AbcSmc/gsl_local
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
IMMDIR = $(DENDIR)/synthetic_population/initial_immunity
SERDIR = $(DENDIR)/synthetic_population/serotype_runs
IMMOBJ = $(IMMDIR)/ImmunityGenerator.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$(IMMDIR) -I$(GSL_PATH)/include/ -I$(ABCDIR)/jsoncpp/include
#ifdef TACC_GSL_INC
#INCLUDE += -I$$TACC_GSL_INC
#endif
#ifdef HPC_GSL_INC
#INCLUDE += -I$$HPC_GSL_INC
#endif

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -lpthread -ldl

default: libabc abc_sql

libabc:  
	$(MAKE) -C $(ABCDIR) -f Makefile

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile model

abc_sql: libabc dengue main.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) -I$(SERDIR) main.cpp -o abc_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_thai: libabc dengue main_sql-thai.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_sql-thai.cpp -o abc_thai $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_deb: libabc dengue main_deb.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_deb.cpp -o main_deb $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc abc_sql
