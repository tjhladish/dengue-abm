-include ../common.mk

default: abc_sql

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

dengue:
	$(MAKE) -C $(DENDIR) -f Makefile

abc_sql: libabc dengue main.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main.cpp -o abc_sql $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_sql-foi_gen: libabc dengue main.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main-generate_alt_foi_immunity_files.cpp -o abc_sql-foi_gen $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

abc_sql-ivn_lag: libabc dengue main.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main-ivn_lag.cpp -o abc_sql-ivn_lag $(DENOBJ) $(ABC_LIB) $(GSL_LIB)

process_who: libabc process_daily_output.cpp
	#$(CXX) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)
	$(CXX) $(CFLAGS) -fopenmp $(INCLUDE) -I$(SQLDIR) process_daily_output.cpp -o process_who $(ABC_LIB) -L$(ABCDIR) $(GSL_LIB)

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

effectiveness.rds: eff.R baseline.rds intervention.rds
	Rscript $^ $@

comboeff.rds: comboeff.R effectiveness.rds
	Rscript $^ $@

effstats.rds: effstats.R comboeff.rds
	Rscript $^ $@

digests: baseline.rds intervention.rds effectiveness.rds comboeff.rds effstats.rds

fig_2.png: fig_2.R effstats.rds
	Rscript $^ $@

fig_1.png: fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ $@

alt_fig_1.png: alt_fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ $@

combo_fig_mech_facet.png: combo_fig_mech_facet.R effstats.rds
	Rscript $^ $@

combo_fig_coverage_facet.png: combo_fig_coverage_facet.R effstats.rds
	Rscript $^ $@

SI_IQ_fig.png: IQ_fig.R effectiveness.rds
	Rscript $^ $@

.PHONY: figures

MANUSCRIPTFIGS := alt_fig_1.png fig_1.png fig_2.png combo_fig_mech_facet.png combo_fig_coverage_facet.png SI_IQ_fig.png

figures: $(MANUSCRIPTFIGS)

clean_figs: $(MANUSCRIPTFIGS)
	rm $^
