-include ../common.mk
-include local.mk

default: libabc abc_sql

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile

abc_sql: libabc dengue main.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main.cpp -o abc_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

process_who: libabc process_daily_output.cpp
	#$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)
	$(CPP) $(CFLAGS) -fopenmp $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc abc_sql
