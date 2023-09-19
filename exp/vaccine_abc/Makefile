-include ../common.mk

SERDIR = $(DENDIR)/synthetic_population/serotype_runs

default: libabc abc_sql

libabc:
	$(MAKE) -C $(ABCDIR)/build -f Makefile abc

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile

abc_sql: libabc dengue main_sql.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) -I$(SERDIR) main_sql.cpp -o abc_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

test_sql: libabc_no_mpi dengue_no_mpi test_sql.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) test_sql.cpp -o test_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc test_sql
