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

baseline.rds intervention.rds: parse.R $(REFDB) projref.R
	Rscript $^ $@

foi_baseline.rds foi_intervention.rds: parse.R $(FOIDB) projref.R
	Rscript $^ $@

lag_intervention.rds: parse.R $(LAGDB) projref.R
	Rscript $^ $@

rootdigests: baseline.rds intervention.rds

foidigests: foi_baseline.rds foi_intervention.rds

effectiveness.rds: eff.R baseline.rds intervention.rds
	Rscript $^ $@

foi_effectiveness.rds: eff.R foi_baseline.rds foi_intervention.rds
	Rscript $^ $@

lag_effectiveness.rds: eff.R baseline.rds lag_intervention.rds
	Rscript $^ $@

effectiveness: effectiveness.rds foi_effectiveness.rds lag_effectiveness.rds

comboeff.rds: comboeff.R effectiveness.rds
	Rscript $^ $@

foi_comboeff.rds: comboeff.R foi_effectiveness.rds
	Rscript $^ $@

lag_comboeff.rds: lag_comboeff.R lag_effectiveness.rds effectiveness.rds
	Rscript $^ $@

comboeff: comboeff.rds foi_comboeff.rds lag_comboeff

effstats.rds: effstats.R comboeff.rds
	Rscript $^ $@

foi_effstats.rds: effstats.R foi_comboeff.rds
	Rscript $^ $@

lag_effstats.rds: effstats.R lag_comboeff.rds
	Rscript $^ $@

effstat: effstats.rds foi_effstats.rds lag_effstats.rds

# not the real fig 2
fig_2.png: fig_2.R effstats.rds
	Rscript $^ $@

# old version
old_fig_1.png: old_fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ $@

lag_plot.png: lag_plot.R lag_effstats.rds
	Rscript $^ $@

foi_plot.png: foi_plot.R foi_baseline.rds foi_intervention.rds foi_effectiveness.rds
	Rscript $^ $@

foi_fig_2.png: foi_fig_2.R foi_effstats.rds effstats.rds
	Rscript $^ $@

fig_1.png fig_1.tiff: fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ $@

# show only the * panels of fig 1
#  if % = vac / vec, then show just the vaccine-only (vector control only, resp) results
#  if % = inc / eff, then show just the incidence (effectiveness, resp) results
fig_1_vac.png fig_1_vec.png fig_1_inc.png fig_1_eff.png: fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ only=$(patsubst fig_1_%.png,%,$@) $@

clean_fig_1_slices: | fig_1_*.png
	rm $|

combo_fig_mech_facet.png: combo_fig_mech_facet.R effstats.rds
	Rscript $^ $@

# real fig 2
combo_fig_coverage_facet.png: combo_fig_coverage_facet.R effstats.rds
	Rscript $^ $@

SI_IQ_fig.png: IQ_fig.R effectiveness.rds
	Rscript $^ $@

SI_ceff_fig.png: ceff_fig.R effectiveness.rds comboeff.rds
	Rscript $^ $@

.PHONY: figures

MANUSCRIPTFIGS := $(addsuffix .png,fig_1 combo_fig_coverage_facet foi_fig_2 SI_IQ_fig SI_ceff_fig)

figures: $(MANUSCRIPTFIGS)

clean_figs: $(MANUSCRIPTFIGS)
	rm $^

# NEW DB IS *ivn_lag2.sqlite
