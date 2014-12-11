CC=mpicxx
#CFLAGS = -g -std=c++11
CFLAGS = -O2 -std=c++11
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
IMMDIR = $(HOME)/work/dengue/synthetic_population/initial_immunity
IMMOBJ = $(IMMDIR)/ImmunityGenerator.o
INCLUDES = -I$(ABCDIR) -I$(DENDIR) -I$(IMMDIR) -I$$TACC_GSL_INC $$HPC_GSL_INC

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp   
GSL_LIB = -lm -L$$TACC_GSL_LIB/ -L$$HPC_GSL_LIB/ -lgsl -lgslcblas

default: all_no_mpi

all_no_mpi: CFLAGS += -std=c++0x -Wall
all_no_mpi: libabc_no_mpi abc 

all_mpi: CC = mpicxx
all_mpi: CFLAGS += -cxx=icc -w0 -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
all_mpi: libabc_mpi abc_mpi 

libabc_no_mpi:  
	$(MAKE) -C $(ABCDIR) -f Makefile all_no_mpi

libabc_mpi:  
	$(MAKE) -C $(ABCDIR) -f Makefile all_mpi

dengue_no_mpi:
	$(MAKE) -C $(DENDIR) -f Makefile all_no_mpi

dengue_mpi:
	$(MAKE) -C $(DENDIR) -f Makefile all_mpi

immunity_mpi:
	$(MAKE) -C $(IMMDIR) -f Makefile all_mpi

abc: main.cpp
	$(CC) $(CFLAGS) main.cpp -o abc $(INCLUDES) $(ABC_LIB) $(GSL_LIB)

abc_mpi: libabc_mpi dengue_mpi immunity_mpi main_mpi.cpp
	$(CC) $(CFLAGS) -D USING_MPI main_mpi.cpp -o abc_mpi $(INCLUDES) $(DENOBJ) $(IMMOBJ) $(ABC_LIB) $(GSL_LIB)

vaccinate_mpi: libabc_mpi dengue_mpi immunity_mpi main_mpi.cpp 
	$(CC) $(CFLAGS) -D USING_MPI main_mpi.cpp -o vaccinate_abc_mpi $(INCLUDES) -L$(DENDIR) $(DENOBJ) $(IMMOBJ) $(ABC_LIB) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc abc_mpi
