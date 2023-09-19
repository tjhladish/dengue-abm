-include local.mk

#CFLAGS = -g -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
CFLAGS = -O2 -std=c++17 -Wall -Wextra -Wno-deprecated-declarations --pedantic
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue-abm
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
IMMDIR = $(DENDIR)/synthetic_population/initial_immunity
IMMOBJ = $(IMMDIR)/ImmunityGenerator.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$(IMMDIR) -I$(ABCDIR)/jsoncpp/include
#ifdef TACC_GSL_INC
#INCLUDE += -I$$TACC_GSL_INC
#endif
#ifdef HPC_GSL_INC
#INCLUDE += -I$$HPC_GSL_INC
#endif

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -lgsl -lgslcblas -lpthread -ldl

default: libabc rzero

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

rzero: dengue main.cpp 
	$(CXX) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main.cpp -o rzero -L$(DENDIR) $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f rzero
