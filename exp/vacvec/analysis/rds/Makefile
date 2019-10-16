vpath %.rds .:..
vpath %.R .:..

-include local.mk

POSTJSON := ../../../abc-irs_refit2/abc-irs_refit2.json

wrap = $(addsuffix $(3),$(addprefix $(1),$(2)))

default: effstat foi_v_mpop.rds posterior.rds

projref.rda: projref.R
	Rscript $^ $@

# INJESTING DB

# create data.table translations of baseline / intervention results
# includes the specific digest<-database dependencies (first few lines)
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


# specific output database dependencies for each digest target
#baseline.rds intervention.rds raw_baseline.csv raw_intervention.csv: $(REFDB)
#test_intervention.rds: $(TSTDB)
#foi_intervention.rds raw_foi_intervention.csv foi_baseline.rds raw_foi_baseline.csv: $(FOIDB)
#foi_test_intervention.rds: $(TSTFOIDB)
#lag_intervention.rds raw_lag_intervention.csv: $(LAGDB)
#testsens_intervention.rds: $(TESTDB)

#baseline.sql intervention.sql foi_baseline.sql foi_intervention.sql lag_intervention.sql testsens_intervention.sql: sql.R projref.rda
#	Rscript $^ $@

#raw_%.csv: %.sql
#	sqlite3 -csv $(word 2,$^) < $< > $@

#digest_%.csv: tocsv.R %.rds
#	Rscript $^ $@

#baseline.rds intervention.rds foi_intervention.rds foi_test_intervention.rds lag_intervention.rds testsens_intervention.rds test_intervention.rds: parse.R utils.R projref.rda
#	Rscript $^ $@

#foi_baseline.rds: parse.R utils.R projref.rda | baseline.rds
#	Rscript $^ tmp_$@
#	Rscript -e "require(data.table); saveRDS({ ref <- readRDS('tmp_$@'); setkeyv(rbind(ref,readRDS('$|')[, foi := 1.0 ]), key(ref)) }, '$@')"
#	rm tmp_$@
	

#ALLDIG := $(addsuffix .rds,baseline intervention foi_baseline foi_intervention lag_intervention testsens_intervention)
# NB: make -j safe
#digests: $(ALLDIG)
#clean_digests:
#	rm -f $(ALLDIG)

#mrg_intervention.rds: substitute.R intervention.rds test_intervention.rds
#	Rscript $^ $@

#sub_intervention.rds: mrg_intervention.rds

#mrg_foi_intervention.rds: substitute.R foi_intervention.rds foi_test_intervention.rds
#	Rscript $^ $@

#sub_foi_intervention.rds: mrg_foi_intervention.rds

#WHICHINT := mrg_intervention
#WHICHFOI := mrg_foi_intervention

# EFFECTIVENESS CALCULATIONS

# specific input digest dependencies for each effectiveness target
#effectiveness.rds: $(addsuffix .rds,baseline $(WHICHINT))
#sub_effectiveness.rds: $(addsuffix .rds,baseline sub_intervention)
#foi_effectiveness.rds: $(addsuffix .rds,foi_baseline $(WHICHFOI))
#sub_foi_effectiveness.rds: $(addsuffix .rds,foi_baseline sub_foi_intervention)

#lag_effectiveness.rds: $(addsuffix .rds,baseline lag_intervention)
#testsens_effectiveness.rds: $(addsuffix .rds,foi_baseline testsens_intervention)
#ALLEFF := $(addsuffix .rds,effectiveness foi_effectiveness lag_effectiveness testsens_effectiveness)

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
#effectiveness.rds foi_effectiveness.rds sub_foi_effectiveness.rds lag_effectiveness.rds testsens_effectiveness.rds: eff.R
#	Rscript $^ $@

# make -j safe
#effectiveness: $(ALLEFF)
#clean_effectiveness:
#	rm -f $(ALLEFF)

# COMBINATION EFFECTIVENESS CALCULATIONS
# combination effectiveness calcs compare naive combinations of interventions
# (naive == each intervention is acting independently)
# with observed combined outcomes
# schema:
# key'd according to underlying effectiveness, then
# combo.eff, c.combo.eff
# vec.eff, c.vec.eff, vac.eff, c.vac.eff
# ind.eff, c.ind.eff
#comboeff.rds: effectiveness.rds
#sub_comboeff.rds: sub_effectiveness.rds effectiveness.rds
#sub_foi_comboeff.rds: sub_foi_effectiveness.rds
#foi_comboeff.rds: foi_effectiveness.rds
#testsens_comboeff.rds: testsens_effectiveness.rds effectiveness.rds
#ALLCMB := $(addsuffix .rds,comboeff foi_comboeff lag_comboeff testsens_comboeff)

#ALLCMB := $(addsuffix .rds,comboeff sub_comboeff foi_comboeff sub_foi_comboeff lag_comboeff)

#comboeff.rds sub_comboeff.rds foi_comboeff.rds sub_foi_comboeff.rds: comboeff.R
#	Rscript $^ $@

#testsens_comboeff.rds: ts_comboeff.R
#	Rscript $^ $@

# lag effectiveness has to be treated specially
#lag_comboeff.rds: lag_comboeff.R lag_effectiveness.rds effectiveness.rds
#	Rscript $^ $@

# NB: make -j safe
#comboeff: $(ALLCMB)
#clean_comboeff:
#	rm -f $(ALLCMB)

## EFFECTIVENESS STATS CALCS
# collapses the particle*replicate keys for effectiveness results
# schema is key of underlying effectiveness / comboeff.rds +
# variable (the combo.eff non-key cols) +
# lo.lo (2.5%), lo (25%), med (50%), hi (75%), hi.hi (97.5%) quantiles

#effstats.rds: comboeff.rds
#foi_effstats.rds: foi_comboeff.rds
#lag_effstats.rds: lag_comboeff.rds
# sub_effstats.rds: sub_comboeff.rds
# sub_foi_effstats.rds: sub_foi_comboeff.rds
# ALLEFS := $(addsuffix .rds,effstats foi_effstats sub_effstats sub_foi_effstats lag_effstats)
# 
# effstats.rds sub_effstats.rds foi_effstats.rds sub_foi_effstats.rds lag_effstats.rds: effstats.R utils.R projref.rda
#testsens_effstats.rds: testsens_comboeff.rds
#ALLEFS := $(addsuffix .rds,effstats foi_effstats lag_effstats testsens_effstats)

#effstats.rds foi_effstats.rds lag_effstats.rds testsens_effstats.rds: effstats.R utils.R projref.rda
#	Rscript $^ $@

# NB: make -j safe
#effstat: $(ALLEFS)
#clean_effstat:
#	rm -f $(ALLEFS)

#clean: clean_digests clean_effectiveness clean_comboeff clean_effstat
#	rm projref.rda

base_single_intervention.rds: base_single_intervention.R projref.rda effstats.rds
	Rscript $^ $@

foi_v_mpop.rds: special.R utils.R $(FOISIDEDB) $(FOISIDETGZ)
	Rscript $^ $@

posterior.rds: posterio.R utils.R $(ABCDB) $(POSTJSON)
	Rscript $^ $@

alt_eff_interventions.rds: alt_parse.R utils.R projref.rda $(ALTEFFDB) all_interventions.rds
	Rscript $^ $@

alt_effectiveness.rds: eff.R all_baseline.rds alt_eff_interventions.rds
	Rscript $^ $@

alt_comboeff.rds: comboeff.R alt_effectiveness.rds
	Rscript $^ $@

all_baseline.rds all_interventions.rds: join_them_all.R $(REFDB) $(FOIDB) $(LAGDB) $(TSTDB) $(TSTFOIDB) $(TESTDB)
	Rscript $^ $@

all_effectiveness.rds: eff.R all_baseline.rds all_interventions.rds
	Rscript $^ $@

nolag_comboeff.rds: comboeff.R all_effectiveness.rds
	Rscript $^ $@

lag_comboeff.rds: lag_comboeff.R all_effectiveness.rds
	Rscript $^ $@

%_effstats.rds: effstats.R utils.R projref.rda %_comboeff.rds
	Rscript $^ $@

effstat: nolag_effstats.rds lag_effstats.rds alt_effstats.rds