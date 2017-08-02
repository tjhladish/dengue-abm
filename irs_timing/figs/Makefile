DBROOT := ~/Dropbox
DBMAIN := $(DBROOT)/who/fig1_data

# usage: $(call first,middles...,last)
wrap = $(addsuffix $(3),$(addprefix $(1),$(2)))

default: $(DBMAIN)/fig1.png

# rules for making fig-ready data

# this rule should come first; other rules are adding specific prerequisites, and need to come later
$(DBMAIN)/%.rds: %.R
	Rscript $^ $@
	# all prerequisites (starts w/ the same-name R file from this rule, others added below) + target

$(DBMAIN)/obs.rds: $(DBMAIN)/weekly_confirmed_cases_1995-2015_wide.csv

$(DBMAIN)/sim.rds: ./daily_output-refit

$(DBMAIN)/eip.rds: $(DBMAIN)/seasonal_EIP_24hr.out

$(DBMAIN)/R0.rds: $(DBMAIN)/rzero-refit.sqlite

$(DBMAIN)/mos.rds: $(DBMAIN)/NOAA_yucatan_daily.csv

$(DBMAIN)/baseline.rds: $(DBROOT)/who/irs_figs/irs_timing-refit0.sqlite

$(DBMAIN)/interventions.rds: $(DBROOT)/who/irs_figs/irs_timing-refit0.sqlite $(DBROOT)/who/irs_insecticide-durability-effect.sqlite

$(DBMAIN)/stat.eff10.rds: $(DBMAIN)/baseline.rds $(DBMAIN)/interventions.rds

$(DBMAIN)/coverage.rds $(DBMAIN)/duration.rds $(DBMAIN)/durability.rds: sensitivity.R $(DBMAIN)/stat.eff10.rds
	Rscript $^ $(DBMAIN)/coverage.rds $(DBMAIN)/duration.rds $(DBMAIN)/durability.rds
	# override generic rule; these have shared features, so build them all at once

TABSRCS := $(call wrap,$(DBMAIN)/,coverage duration durability,.rds)
SRCS := $(call wrap,$(DBMAIN)/,obs sim eip R0 mos,.rds) $(TABSRCS)

$(DBMAIN)/fig1.png: main_fig.R $(SRCS)
	Rscript $^ $@

$(DBMAIN)/tab1.csv: main_tab.R $(TABSRCS)
	Rscript $^ $@

$(DBMAIN)/offset_peak.csv: offset.R $(DBMAIN)/tab1.csv
	Rscript $^ $@
