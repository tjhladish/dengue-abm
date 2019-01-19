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

# INJESTING DB

# specific output database dependencies for each digest target
baseline.rds intervention.rds baseline.csv intervention.csv: $(REFDB)
foi_baseline.rds foi_intervention.rds foi_baseline.csv foi_intervention.csv: $(FOIDB)
lag_intervention.rds lag_intervention.csv: $(LAGDB)
ALLDIG := $(addsuffix .rds,baseline intervention foi_baseline foi_intervention lag_intervention)

baseline.sql intervention.sql foi_baseline.sql foi_intervention.sql lag_intervention.sql: sql.R projref.R
	Rscript $^ $@

%.csv: %.sql
	sqlite3 -csv $(word 2,$^) < $< > $@

# create data.table translations of baseline / intervention results
# includes the specific digest<-database dependencies (above)
# via make dependency combination rules
#
# baseline schema:
# key(particle, replicate, year), seronegative + imm1-imm4 (all %), s, c.s (ints)
#
# intervention schema, same, but with:
# key(vc, vac, vc_coverage, vaccine, catchup) added (on the left) to cols
#
# foi_% schema: same as baseline + intervention, w/ key(foi) added on the left
# lag_% schema: same as intervention, w/ key(vac_first) added on the left
baseline.rds intervention.rds foi_baseline.rds foi_intervention.rds lag_intervention.rds: parse.R projref.R
	Rscript $^ $@

# NB: make -j safe
digests: $(ALLDIG)
clean_digests:
	rm -f $(ALLDIG)

# EFFECTIVENESS CALCULATIONS

# specific input digest dependencies for each effectiveness target
effectiveness.rds: baseline.rds intervention.rds
foi_effectiveness.rds: foi_baseline.rds foi_intervention.rds
lag_effectiveness.rds: baseline.rds lag_intervention.rds
ALLEFF := $(addsuffix .rds,effectiveness foi_effectiveness lag_effectiveness)

# TODO document c = cumulative
# TODO rename ind = naive?
# TODO breaking comment bars in R?

# create data.table effectiveness computations
#
# basic (effectiveness.rds) schema:
# key(vc, vac, vc_coverage, vaccine, catchup, particle, replicate, year)
# bcases (baseline cases), icases (int. cases), c.bcases, c.icases (cum. cases of base and int.)
# eff ((bcases - icases)/bcases), c.eff ((c.bcases - c.icases)/c.bcases),
# logm (log(icases/bcases)), c.logm(log(c.icases/c.bcases))
#
# foi_% schema: same as basic, w/ key(foi) added on the left
# lag_% schema: same as basic, w/ key(vac_first) added on the left
effectiveness.rds foi_effectiveness.rds lag_effectiveness.rds: eff.R
	Rscript $^ $@

# make -j safe
effectiveness: $(ALLEFF)
clean_effectiveness:
	rm -f $(ALLEFF)

# COMBINATION EFFECTIVENESS CALCULATIONS
# combination effectiveness calcs compare naive combinations of interventions
# (naive == each intervention is acting independently)
# with observed combined outcomes
# schema:
# key'd according to underlying effectiveness, then
# combo.eff, c.combo.eff
# vec.eff, c.vec.eff, vac.eff, c.vac.eff
# ind.eff, c.ind.eff
comboeff.rds: effectiveness.rds
foi_comboeff.rds: foi_effectiveness.rds

comboeff.rds foi_comboeff.rds: comboeff.R
	Rscript $^ $@

# lag effectiveness has to be treated specially
lag_comboeff.rds: lag_comboeff.R lag_effectiveness.rds effectiveness.rds
	Rscript $^ $@

# NB: make -j safe
comboeff: comboeff.rds foi_comboeff.rds lag_comboeff.rds

## EFFECTIVENESS STATS CALCS
# collapses the particle*replicate keys for effectiveness results
# schema is key of underlying effectiveness / comboeff.rds +
# variable (the combo.eff non-key cols) +
# lo.lo (2.5%), lo (25%), med (50%), hi (75%), hi.hi (97.5%) quantiles

effstats.rds: comboeff.rds
foi_effstats.rds: foi_comboeff.rds
lag_effstats.rds: lag_comboeff.rds
ALLEFS := $(addsuffix .rds,effstats foi_effstats lag_effstats)

effstats.rds foi_effstats.rds lag_effstats.rds: effstats.R
	Rscript $^ $@

# NB: make -j safe
effstat: $(ALLEFS)
clean_effstat:
	rm -f $(ALLEFS)

## FIGURES

# TODO: Figure 1 generation - should be somewhat automated? notes:
#  four panel diagram, as column (possibly row? if row, first three panels over fourth)
#   * p1: night time lights + translation into network;
#     show night lights + admin boundaries, but have a vertical (horiz if row)
#     slice removed / transparent to "reveal" network
#   * p2: zoom of network from "revealed" segment, differentiating location types
#     by color + shape (basic shapes or maybe silhouettes? solid color)
#     show mosquito movement probability weights (network line width? alpha? bidirectional?)
#     also, show human movement lines? e.g., work / school connections - color same as those
#     locations?  Maybe self loops for stay-at-home?
#   * p3: expand vertices (eliminating or leaving largely diminished edges) to show people
#     this part doesn't need to be automatic; perhaps just the background elements (+mosquitos?) +
#     a table of data to inform making the people puppets to insert manually
#     in addition to showing demographic heterogeneity, should also show some transmission states
#     e.g., empty mosquito silhouette (S, just one, show location mosquito count),
#     yellow (E, multiple mos), red (I, multiple) mosquitos +
#     humans with empty outline (full S), increasingly filled outline (exposure ordinality; solid
#     black = completely immune), yellow (E), red (I), blue (transient immunity?)
#     maybe show S moquito "perched" on vertix (as if feeding), with E/I mosquitos inside
#   * p4: last panel showing state machine, colored to match human / mosquito colors in panel 3

NIGHTLIGHT := some.csv # some csv? shapefile?
ADMINBOUNDS := admin.shp # seems like shapefile
LOCATIONS := vertices.csv # some csv? shapefile?
CONNECTIONS := edges.csv # some csv? shapefile?
TARRECT := tar.json # something specifying lat / long extent of zoom, other annotations (e.g., vertical reveal?)

fig1p1.png: diagram_panel_wideview.R $(NIGHTLIGHT) $(ADMINBOUNDS) $(LOCATIONS) $(CONNECTIONS) $(TARRECT)
	Rscript $^ $@

fig1p2.png: diagram_panel_zoom.R $(LOCATIONS) $(CONNECTIONS) $(TARRECT) # $(ADMINBOUNDS)? in case the zoom rect crosses admin lines?
	Rscript $^ $@

fig1p3back.png: diagram_panel_model_background.R $(LOCATIONS) $(CONNECTIONS) $(TARRECT) # $(ADMINBOUNDS)? in case the zoom rect crosses admin lines?
	Rscript $^ $@

fig1.png: fig1p1.png fig1p2.png fig1p3.png fig1p4.png
	@echo $^

# figure 2: show prereq results (single interventions), raw main results (combo effectiveness) 
raw_effectiveness_plot.png: raw_effectiveness_plot.R effstats.rds
	Rscript $^ $@



# not the real fig 2
fig_2.png: fig_2.R effstats.rds
	Rscript $^ $@

# old version
old_fig_1.png: old_fig_1.R baseline.rds intervention.rds effstats.rds
	Rscript $^ $@

lag_plot.png: lag_plot.R lag_effstats.rds effstats.rds
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

MANUSCRIPTFIGS := $(addsuffix .png,fig_1 combo_fig_coverage_facet foi_fig_2 lag_plot SI_IQ_fig SI_ceff_fig)

figures: $(MANUSCRIPTFIGS)

clean_figs: $(MANUSCRIPTFIGS)
	rm $^

# christmas figures
# fig 1 - diagram, probably not managed via make

# fig 2: single intervention incidence, annual effectiveness, cumulative effectiveness
# row 1: incidence
# row 2: ann. effectiveness
# row 3: cum. effectiveness
# one col; show only "best" versions of intervention (w/ catchup, max vc)
# distinguish vaccine by line type = solid, none; dashes, EDV; dotted, Dengvaxia-like
# complementary SI fig:
# multiple cols: one for each vaccine type (none, EDV, Dengvaxia-like)
# IQs

new_SI_fig_2.png new_SI_fig_2.tiff new_fig_2.png new_fig_2.tiff: new_fig_2.R baseline.rds intervention.rds effectiveness.rds
	Rscript $^ $@
	# 1 col, 3 rows - incidence, annual effectiveness, cumulative effectiveness
	# both vaccines w/ catchup, 75% vc - IQs?

new_fig_3.png new_fig_3.tiff: new_fig_3.R effectiveness.rds comboeff.rds
	Rscript $^ $@
	# 2 rows, 3 cols - annual effectiveness of combination interventions, with background ref lines
	# for single interventions; rows are vaccine type, cols are vc level

new_fig_4.png new_fig_4.tiff: new_fig_4.R effectiveness.rds
	Rscript $^ $@
	# repeat fig 3, but with shaded effect regions?
	# for single interventions; rows are vaccine type, cols are vc level

new_fig_5.png new_fig_5.tiff: new_fig_5.R lag_effectiveness.rds
	Rscript $^ $@
	# lagged interventions