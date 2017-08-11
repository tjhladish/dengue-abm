DBROOT := ~/Dropbox
DBMAIN := $(DBROOT)/who/fig1_data

# usage: $(call first,middles...,last)
wrap = $(addsuffix $(3),$(addprefix $(1),$(2)))

default: supfigs mainfigs

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

$(DBMAIN)/baseline.cases.rds: $(DBMAIN)/stopping-baseline.rds

$(DBMAIN)/interventions.cases.rds: $(DBMAIN)/stopping-interventions.rds

$(DBMAIN)/averted.cases.rds: $(DBMAIN)/baseline.cases.rds $(DBMAIN)/interventions.cases.rds

$(DBMAIN)/stat.eff10.rds: $(DBMAIN)/baseline.rds $(DBMAIN)/interventions.rds

$(DBMAIN)/coverage.rds $(DBMAIN)/duration.rds $(DBMAIN)/durability.rds: sensitivity.R $(DBMAIN)/stat.eff10.rds
	Rscript $^ $(DBMAIN)/coverage.rds $(DBMAIN)/duration.rds $(DBMAIN)/durability.rds
	# override generic rule; these have shared features, so build them all at once

TABSRCS := $(call wrap,$(DBMAIN)/,coverage duration durability,.rds)
SRCS := $(call wrap,$(DBMAIN)/,obs sim eip R0 mos,.rds) $(TABSRCS)

$(DBMAIN)/fig1.png: main_fig.R $(SRCS)
	Rscript $^ $@

$(DBMAIN)/stopping-baseline.rds $(DBMAIN)/stopping-interventions.rds: $(DBROOT)/who/irs_stopping-effect_rerun.sqlite

$(DBMAIN)/stopping-eff.rds $(DBMAIN)/stopping-sero.rds $(DBMAIN)/stopping-sero-SI.rds: $(DBMAIN)/stopping-baseline.rds $(DBMAIN)/stopping-interventions.rds

$(DBMAIN)/fig2.png: main_fig2.R $(DBMAIN)/stopping-eff.rds $(DBMAIN)/stopping-sero.rds
	Rscript $^ $@

$(DBMAIN)/fig2_alt.png: main_fig2_alt.R $(DBMAIN)/stopping-eff.rds $(DBMAIN)/stopping-sero.rds
	Rscript $^ $@

mainfigs: $(call wrap,$(DBMAIN)/,fig1 fig2_alt,.png)

supfigs: $(call wrap,$(DBMAIN)/,irs_timing-realigned-SI vc_4-panel-impact_10yr vc_4-panel-impact_50yr stacked_immunity pro_vs_re_active-prevalence.local foi_effect,.png)

$(DBMAIN)/SI-fig-sero.png: SI-fig-sero.R $(DBMAIN)/stopping-sero-SI.rds
	Rscript $^ $@

$(DBMAIN)/irs_timing-realigned-SI.png: timing_plot_realigned-supplement.R $(DBMAIN)/stat.eff10.rds $(DBMAIN)/R0.rds $(DBMAIN)/mos.rds
	Rscript $^ $@

$(DBMAIN)/vc_4-panel-impact_10yr.png $(DBMAIN)/vc_4-panel-impact_50yr.png: annual_effectiveness_plot.R $(DBMAIN)/baseline.cases.rds $(DBMAIN)/interventions.cases.rds $(DBMAIN)/averted.cases.rds $(DBMAIN)/stopping-eff.rds
	Rscript $^ $@

$(DBMAIN)/foi_data.rds: $(DBMAIN)/irs_timing-summer_winter-foi.sqlite

$(DBMAIN)/foi_effect.png: foi_plot.R $(DBMAIN)/foi_data.rds
	Rscript $^ $@

$(DBMAIN)/stacked_immunity.png: annual_immunity_plot-SI.R $(DBMAIN)/stopping-baseline.rds $(DBMAIN)/stopping-interventions.rds
	Rscript $^ $@

CMPRES = incidence.intro incidence.local prevalence.intro prevalence.local

# override production of these .rds, since they all use same script
$(call wrap,$(DBMAIN)/campaign-timing-,$(CMPRES),.rds): campaign-timing.R $(DBMAIN)/irs_timing-summer_winter-foi.sqlite $(DBROOT)/daily-irs_refit-simple-nofilenames-uniq.out
	Rscript $^ $@

$(DBMAIN)/pro_vs_re_active-%.png: pro_vs_re_active.R $(DBMAIN)/campaign-timing-%.rds
	Rscript $^ $@

$(DBMAIN)/tab1.csv: main_tab.R $(TABSRCS)
	Rscript $^ $@

$(DBMAIN)/offset_peak.csv: offset.R $(DBMAIN)/tab1.csv
	Rscript $^ $@
