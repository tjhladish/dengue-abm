-include ../../exp/common.mk

sero: libabc
	$(CXX) $(CFLAGS) $(INCLUDE) -Wno-deprecated-declarations fit_serotype_data.cpp -o simulate_serotypes -I/home/tjhladish/work/AbcSmc $(ABC_LIB) $(GSL_LIB)

generator: yucatan_serotype_generator.h
	$(CXX) $(CFLAGS) yucatan_serotype_generator.cpp -o sero_generator $(GSL_LIB)

libabc:
	$(MAKE) -C $(ABCDIR)/build -f Makefile abc

gen_sero:
	$(CXX) -O2 -std=c++11 generate_serotype_runs.cpp -o gen_sero -lm -lgsl -lgslcblas

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	rm -f simulate_serotypes sero_generator gen_sero
