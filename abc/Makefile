SHELL=/bin/bash
G++VER := $(shell command -v g++-4.9)

ifndef G++VER
CPP:=g++
else
CPP:=g++-4.9
endif

#CFLAGS = -g -std=c++11 -Wall
CFLAGS = -O2 -std=c++11 -pedantic -Wall
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
IMMDIR = $(HOME)/work/dengue/synthetic_population/initial_immunity
SERDIR = $(DENDIR)/synthetic_population/serotype_runs
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
	$(MAKE) -C $(DENDIR) -f Makefile model

abc_sql: libabc dengue main_sql.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) -I$(SERDIR) main_sql.cpp -o abc_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_thai: libabc dengue main_sql-thai.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_sql-thai.cpp -o abc_thai $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_deb: libabc dengue main_deb.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_deb.cpp -o main_deb $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc abc_sql
