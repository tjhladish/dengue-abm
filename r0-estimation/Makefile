SHELL=/bin/bash
G++VER := $(shell command -v g++-4.9)

ifndef G++VER
CPP:=g++
else
CPP:=g++-4.9
endif

#CFLAGS = -g -std=c++11
CFLAGS = -O2 -std=c++11 -pedantic -Wall#-Wextra
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$$TACC_GSL_INC $$HPC_GSL_INC
SQLDIR = $(ABCDIR)/sqdb

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o 
GSL_LIB = -lm -L$$TACC_GSL_LIB/ -L$$HPC_GSL_LIB/ -lgsl -lgslcblas -lpthread -ldl

rzero: dengue main.cpp 
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main.cpp -o rzero -L$(DENDIR) $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	rm -f rzero
