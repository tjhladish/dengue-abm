CC=g++
CFLAGS = -O2 -std=c++11 -Wall -fopenmp
ABCDIR = $(HOME)/work/AbcSmc
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(SQLDIR) 

ABC_LIB = -lsqdb $(ABCDIR)/sqlite3.o -lpthread -ldl

libabc_no_mpi:  
	$(MAKE) -C $(ABCDIR) -f Makefile all_no_mpi

process_who: libabc_no_mpi process_daily_output.cpp
	$(CC) $(CFLAGS) $(INCLUDE) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	rm -f process_who
