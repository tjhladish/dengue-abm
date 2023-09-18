SHELL=/bin/bash
G++VER := $(shell command -v g++-4.9)

ifndef G++VER
CPP:=g++
else
CPP:=g++-4.9
endif

#CFLAGS = -g -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
CFLAGS = -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue-abm
GSL_PATH = $(HOME)/work/AbcSmc/gsl_local
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$(GSL_PATH)/include/
ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -lpthread -ldl

default: libabc rzero

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

rzero: dengue main.cpp 
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main.cpp -o rzero -L$(DENDIR) $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	rm -f rzero
