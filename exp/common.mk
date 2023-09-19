CPP = g++
#CFLAGS = -g -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
CFLAGS = -O2 -std=c++17 -Wall -Wextra -Wno-deprecated-declarations --pedantic
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue-abm
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
IMMDIR = $(DENDIR)/synthetic_population/initial_immunity
IMMOBJ = $(IMMDIR)/ImmunityGenerator.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$(IMMDIR) -I$(ABCDIR)/jsoncpp/include

INCLUDE += -I$(ABCDIR)/include -I$(ABCDIR)/lib -I$(ABCDIR)/lib/sqdb/include
INCLUDE += -I$(ABCDIR)/lib/CCRC32/include -I$(ABCDIR)/lib/PLS/include
INCLUDE += -I$(ABCDIR)/lib/PLS/lib/eigen -I$(ABCDIR)/lib/jsoncpp/include

#ifdef TACC_GSL_INC
#INCLUDE += -I$$TACC_GSL_INC
#endif
#ifdef HPC_GSL_INC
#INCLUDE += -I$$HPC_GSL_INC
#endif

ABC_LIB := -L$(ABCDIR)/build -L$(ABCDIR)/build/PLS -labc -lpls -ljsoncpp -lm
GSL_LIB = -lm -lgsl -lgslcblas -lpthread -ldl
