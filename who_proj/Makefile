CC=g++
#CC=mpicxx
#CFLAGS = -g -std=c++11
CFLAGS = -O2 -std=c++11 -Wall
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
IMMDIR = $(HOME)/work/dengue/synthetic_population/initial_immunity
IMMOBJ = $(IMMDIR)/ImmunityGenerator.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$(IMMDIR) -I$(SQLDIR)
ifdef TACC_GSL_INC
INCLUDE += -I$$TACC_GSL_INC
endif
ifdef HPC_GSL_INC
INCLUDE += -I$$HPC_GSL_INC
endif

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o -lpthread -ldl
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
	$(CC) $(CFLAGS) main.cpp -o abc $(INCLUDE) $(ABC_LIB) $(GSL_LIB)

abc_brazil: libabc_no_mpi dengue_no_mpi main_brazil.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_brazil.cpp -o abc_brazil $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_thai: libabc_no_mpi dengue_no_mpi main_thai.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_thai.cpp -o abc_thai $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

test_sql: libabc_no_mpi dengue_no_mpi test_sql.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) test_sql.cpp -o test_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

process_who: libabc_no_mpi process_daily_output.cpp
	$(CC) $(CFLAGS) -fopenmp $(INCLUDE) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc abc_mpi process_who
