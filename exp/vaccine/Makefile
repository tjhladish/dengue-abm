CC=g++
CFLAGS = -O2 -std=c++11 -pedantic -Wall
DENDIR = $(HOME)/work/dengue-abm
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
IMMDIR = $(DENDIR)/synthetic_population/initial_immunity
IMMOBJ = $(IMMDIR)/ImmunityGenerator.o
INCLUDES = -I$(DENDIR) -I$(IMMDIR)

GSL_LIB = -lm -lgsl -lgslcblas

default: vaccinate

dengue_no_mpi:
	$(MAKE) -C $(DENDIR) -f Makefile all_no_mpi

immunity_no_mpi:
	$(MAKE) -C $(IMMDIR) -f Makefile all_no_mpi

vaccinate: main.cpp dengue_no_mpi immunity_no_mpi
	$(CC) $(CFLAGS) $(INCLUDES) main.cpp -o vaccinate -L$(DENDIR) $(DENOBJ) $(IMMOBJ) $(GSL_LIB)

clean:
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f vaccinate
