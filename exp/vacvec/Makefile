SHELL=/bin/bash
G++VER := $(shell command -v g++-4.9)

ifndef G++VER
CPP:=g++
else
CPP:=g++-4.9
endif

#CFLAGS = -g -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
CFLAGS = -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
ABCDIR = $(HOME)/work/AbcSmc
DENDIR = $(HOME)/work/dengue
GSL_PATH = $(HOME)/work/AbcSmc/gsl_local
DENOBJ = $(DENDIR)/Person.o $(DENDIR)/Location.o $(DENDIR)/Mosquito.o $(DENDIR)/Community.o $(DENDIR)/Parameters.o $(DENDIR)/Utility.o
IMMDIR = $(HOME)/work/dengue/synthetic_population/initial_immunity
IMMOBJ = $(IMMDIR)/ImmunityGenerator.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$(IMMDIR) -I$(GSL_PATH)/include/
#ifdef TACC_GSL_INC
#INCLUDE += -I$$TACC_GSL_INC
#endif
#ifdef HPC_GSL_INC
#INCLUDE += -I$$HPC_GSL_INC
#endif

ABC_LIB = -L$(ABCDIR) -L$(DENDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -lpthread -ldl

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile

abc_sql: libabc dengue main.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main.cpp -o abc_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_sql-ivn_lag: libabc dengue main.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main-ivn_lag.cpp -o abc_sql-ivn_lag $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

process_who: libabc process_daily_output.cpp
	#$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)
	$(CPP) $(CFLAGS) -fopenmp $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(DENDIR) -f Makefile clean
	$(MAKE) -C $(IMMDIR) -f Makefile clean
	rm -f abc abc_sql

-include local.mk

baseline.rds: parse.R $(REFDB)
	Rscript $^ $@

intervention.rds: parse.R $(REFDB)
	Rscript $^ $@

lagged_intervention.rds: lagparse.R $(LAGDB)
	Rscript $^ $@

effectiveness.rds: eff.R baseline.rds intervention.rds
	Rscript $^ $@

lag_effectiveness.rds: eff.R baseline.rds lagged_intervention.rds
	Rscript $^ $@

comboeff.rds: comboeff.R effectiveness.rds
	Rscript $^ $@

effstats.rds: effstats.R comboeff.rds
	Rscript $^ $@

digests: effectiveness.rds lag_effectiveness.rds comboeff.rds effstats.rds

fig_2.png: fig_2.R effstats.rds
	Rscript $^ $@

# old version
#fig_1.png: fig_1.R baseline.rds intervention.rds effstats.rds
#	Rscript $^ $@

lag_plot.png: lag_plot.R baseline.rds lagged_intervention.rds lag_effectiveness.rds
	Rscript $^ $@

fig_1.png: fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ $@

# show only the incidence panels of fig 1
fig_1_inc.png: fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ only=inc $@

# show only the effectiveness panels of fig 1
fig_1_eff.png: fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ only=eff $@

# show only the vaccine panels of fig 1
fig_1_vac.png: fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ only=vac $@

# show only the vector control panels of fig 1
fig_1_vec.png: fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ only=vec $@

combo_fig_mech_facet.png: combo_fig_mech_facet.R effstats.rds
	Rscript $^ $@

combo_fig_coverage_facet.png: combo_fig_coverage_facet.R effstats.rds
	Rscript $^ $@

SI_IQ_fig.png: IQ_fig.R effectiveness.rds
	Rscript $^ $@

SI_ceff_fig.png: ceff_fig.R effectiveness.rds comboeff.rds
	Rscript $^ $@

.PHONY: figures

MANUSCRIPTFIGS := alt_fig_1.png fig_1.png fig_2.png combo_fig_mech_facet.png combo_fig_coverage_facet.png SI_IQ_fig.png

figures: $(MANUSCRIPTFIGS)

clean_figs: $(MANUSCRIPTFIGS)
	rm $^
