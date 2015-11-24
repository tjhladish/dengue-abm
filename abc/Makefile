CC=g++-4.9
#CC=mpicxx
#CFLAGS = -g -std=c++11
CFLAGS = -O2 -std=c++11
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

default: all_no_mpi

all_no_mpi: CFLAGS += -std=c++11 -Wall
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
	$(CC) $(CFLAGS) $(INCLUDE) main.cpp -o abc $(ABC_LIB) $(GSL_LIB)

abc_mpi: libabc_mpi dengue_mpi immunity_mpi main_mpi.cpp
	$(CC) $(CFLAGS) -D USING_MPI $(INCLUDE) main_mpi.cpp -o abc_mpi $(DENOBJ) $(IMMOBJ) $(ABC_LIB) $(GSL_LIB)

abc_sql: libabc_no_mpi dengue_no_mpi main_sql.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) -I$(SERDIR) main_sql.cpp -o abc_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_thai: libabc_no_mpi dengue_no_mpi main_sql-thai.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_sql-thai.cpp -o abc_thai $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_deb: libabc_no_mpi dengue_no_mpi main_deb.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_deb.cpp -o main_deb $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc abc_mpi
