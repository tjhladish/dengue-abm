# defines
REFS ?= data_locations.mk
include $(REFS)

# usage: $(call first,middles...,last)
wrap = $(addsuffix $(3),$(addprefix $(1),$(2)))

#default: $(DBMAIN)/fig1.png
default: supfigs mainfigs $(DBMAIN)/check_runs.png

clean:
	rm $(DBMAIN)/*.rds

.SECONDARY:

# rules for making fig-ready data

$(DBMAIN):
	mkdir $(DBMAIN)

# this rule should come first; other rules are adding specific prerequisites, and need to come later
$(DBMAIN)/%.rds: %.R | $(DBMAIN)
	Rscript $^ $@

########################
# main text figure data

$(DBMAIN)/obs.rds: $(DBMAIN)/weekly_confirmed_cases_1995-2015_wide.csv

$(DBMAIN)/sim.rds: ./daily_output-refit

$(DBMAIN)/eip.rds: $(DBMAIN)/seasonal_EIP_24hr.out

$(DBMAIN)/R0.rds: $(DBMAIN)/rzero-refit.sqlite

$(DBMAIN)/mos.rds: $(DBMAIN)/NOAA_yucatan_daily.csv

$(DBMAIN)/baseline.rds: $(SQMAIN)

$(DBMAIN)/interventions.rds: $(SQMAIN)

$(DBMAIN)/stat.eff10.rds: $(DBMAIN)/baseline.rds $(DBMAIN)/interventions.rds

$(DBMAIN)/coverage.rds $(DBMAIN)/duration.rds $(DBMAIN)/durability.rds: sensitivity.R $(DBMAIN)/stat.eff10.rds
	Rscript $^ $(DBMAIN)/coverage.rds $(DBMAIN)/duration.rds $(DBMAIN)/durability.rds
	# override generic rule; these have shared features, so build them all at once

$(DBMAIN)/confidence-intervals.rds: $(SQMAIN)

#########################
# conserved rules

$(DBMAIN)/%-eff.rds: alt-eff.R $(DBMAIN)/%-baseline.rds $(DBMAIN)/%-interventions.rds
	Rscript $^ $@

$(DBMAIN)/stopping-sero.rds $(DBMAIN)/stopping-sero-SI.rds: $(DBMAIN)/stopping-baseline.rds $(DBMAIN)/stopping-interventions.rds

$(DBMAIN)/foi-sero.rds $(DBMAIN)/foi-sero-SI.rds: $(DBMAIN)/foi-baseline.rds $(DBMAIN)/foi-interventions.rds

$(DBMAIN)/%.cases.rds: cases.R $(DBMAIN)/%.rds
	Rscript $^ $@

$(DBMAIN)/%-averted.rds: averted.cases.R $(DBMAIN)/%-baseline.cases.rds $(DBMAIN)/%-interventions.cases.rds
	Rscript $^ $@


#########################
# SI 4-panel figure data

$(DBMAIN)/stopping-baseline.rds $(DBMAIN)/stopping-interventions.rds: $(SQMAIN)

######################
# foi analysis data

$(DBMAIN)/foi-baseline.rds $(DBMAIN)/foi-interventions.rds: foi.R $(FOI)
	Rscript $^ $@

$(DBMAIN)/4panel_foi_20yr.png: foi_annual_effectiveness_plot.R $(DBMAIN)/foi-baseline.cases.rds $(DBMAIN)/foi-interventions.cases.rds $(DBMAIN)/foi-averted.rds $(DBMAIN)/foi-eff.rds
	Rscript $^ $@


TABSRCS := $(call wrap,$(DBMAIN)/,coverage duration durability,.rds)
SRCS := $(call wrap,$(DBMAIN)/,obs sim eip R0 mos,.rds) $(TABSRCS)

$(DBMAIN)/fig1.png: main_fig.R $(SRCS)
	Rscript $^ $@

$(DBMAIN)/fig2.png: main_fig2.R $(DBMAIN)/stopping-eff.rds $(DBMAIN)/stopping-sero.rds
	Rscript $^ $@

$(DBMAIN)/fig2_alt.png: main_fig2_alt.R $(DBMAIN)/stopping-eff.rds $(DBMAIN)/stopping-sero.rds
	Rscript $^ $@

mainfigs: $(call wrap,$(DBMAIN)/,fig1 fig2_alt,.png)
fig1: $(call wrap,$(DBMAIN)/,fig1,.png)
fig2: $(call wrap,$(DBMAIN)/,fig2_alt,.png)

supfigs: $(call wrap,$(DBMAIN)/,irs_timing-realigned-SI irs_timing-realigned-SI-alt vc_4-panel-impact_10yr vc_4-panel-impact_50yr stacked_immunity pro_vs_re_active-prevalence.local foi_effect,.png)

$(DBMAIN)/SI-fig-sero.png: SI-fig-sero.R $(DBMAIN)/stopping-sero-SI.rds
	Rscript $^ $@

$(DBMAIN)/irs_timing-realigned-SI.png: timing_plot_realigned-supplement.R $(DBMAIN)/stat.eff10.rds $(DBMAIN)/R0.rds $(DBMAIN)/mos.rds
	Rscript $^ $@

$(DBMAIN)/irs_timing-realigned-SI-alt.png: timing_plot_realigned-supplement-alt.R $(DBMAIN)/stat.eff10.rds $(DBMAIN)/R0.rds $(DBMAIN)/mos.rds
	Rscript $^ $@

$(DBMAIN)/vc_4-panel-impact_10yr.png $(DBMAIN)/vc_4-panel-impact_50yr.png: annual_effectiveness_plot.R $(DBMAIN)/stopping-baseline.cases.rds $(DBMAIN)/stopping-interventions.cases.rds $(DBMAIN)/stopping-averted.rds $(DBMAIN)/stopping-eff.rds
	Rscript $^ $@

$(DBMAIN)/foi_effect.png: foi_plot.R $(DBMAIN)/foi-eff.rds
	Rscript $^ $@

$(DBMAIN)/stacked_immunity.png: annual_immunity_plot-SI.R $(DBMAIN)/stopping-baseline.rds $(DBMAIN)/stopping-interventions.rds
	Rscript $^ $@

CMPRES = incidence.intro incidence.local prevalence.intro prevalence.local

# override production of these .rds, since they all use same script
$(call wrap,$(DBMAIN)/campaign-timing-,$(CMPRES),.rds): campaign-timing.R $(FOI) $(DAILY)
	Rscript $^ $@

$(DBMAIN)/pro_vs_re_active-%.png: pro_vs_re_active.R $(DBMAIN)/campaign-timing-%.rds
	Rscript $^ $@

$(DBMAIN)/tab1.csv: main_tab.R $(TABSRCS)
	Rscript $^ $@

$(DBMAIN)/offset_peak.csv: offset.R $(DBMAIN)/tab1.csv
	Rscript $^ $@

$(DBMAIN)/confidence-intervals.rds: $(SQMAIN) $(DBROOT)/irs_insecticide-durability-effect_intro-fix.sqlite

$(DBMAIN)/check_runs.png: check_runs.R $(SQMAIN)
	Rscript $^ $@
