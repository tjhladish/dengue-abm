# defines
REFS ?= data_locations.mk
include $(REFS)

# usage: $(call first,middles...,last)
wrap = $(addsuffix $(3),$(addprefix $(1),$(2)))

#default: $(DBMAIN)/fig1.png
default: supfigs mainfigs $(DBMAIN)/check_runs.png

clean:
	rm $(DBMAIN)/*.rds

check_requires:
  # some command to grep all the requires
  # feed that into R

.SECONDARY:

# rules for making fig-ready data

$(DBMAIN):
	mkdir $(DBMAIN)

# this rule should come first; other rules are adding specific prerequisites, and need to come later
$(DBMAIN)/%.rds: %.R | $(DBMAIN)
	Rscript $^ $@

########################
numbers:
	Rscript discontinue_VC.R $(DBMAIN)/stopping-eff.rds


########################
# main text figure data

$(DBMAIN)/obs.rds: $(DBMAIN)/weekly_confirmed_cases_1995-2015_wide.csv

#$(DBMAIN)/sim.rds: ./daily_output-refit
$(DBMAIN)/sim.rds: $(DBMAIN)/daily_cases_1995-2015_refit2.out

$(DBMAIN)/eip.rds: $(DBMAIN)/seasonal_EIP_24hr.out

$(DBMAIN)/R0.rds: $(DBMAIN)/rzero-refit2.sqlite

$(DBMAIN)/mos.rds: $(DBMAIN)/NOAA_yucatan_daily.csv

$(DBMAIN)/baseline.rds: $(SQMAIN)

$(DBMAIN)/interventions.rds: $(SQMAIN) $(SQREV)

$(DBMAIN)/extra-interventions.rds: $(SQEXT)

$(DBMAIN)/stat.eff10.rds: $(DBMAIN)/baseline.rds $(DBMAIN)/interventions.rds
	
$(DBMAIN)/extra-stat.eff10.rds: stat.eff10.R $(DBMAIN)/baseline.rds $(DBMAIN)/extra-interventions.rds
	Rscript $^ $@

SENSITIVITES := coverage duration durability efficacy

$(call wrap,$(DBMAIN)/,$(SENSITIVITES),.rds): sensitivity.R $(DBMAIN)/stat.eff10.rds
	Rscript $^ $(call wrap,$(DBMAIN)/,$(SENSITIVITES),.rds)
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


TABSRCS := $(call wrap,$(DBMAIN)/,coverage duration durability efficacy,.rds)
SRCS := $(call wrap,$(DBMAIN)/,obs sim eip R0 mos,.rds) $(TABSRCS)

$(DBMAIN)/fig1.png $(DBMAIN)/fig1.tiff: main_fig.R $(SRCS)
	Rscript $^ $@

$(DBMAIN)/fig1-SI.png: $(DBMAIN)/fig1.png

$(DBMAIN)/fig1-SI.tiff: $(DBMAIN)/fig1-SI.png

$(DBMAIN)/fig2.png: main_fig2.R $(DBMAIN)/stopping-eff.rds $(DBMAIN)/stopping-sero.rds
	Rscript $^ $@

$(DBMAIN)/fig2_alt.png $(DBMAIN)/fig2_alt.tiff: main_fig2_alt.R $(DBMAIN)/stopping-eff.rds $(DBMAIN)/stopping-sero.rds
	Rscript $^ $@

mainfigs: $(call wrap,$(DBMAIN)/,fig1 fig2_alt,.png)
fig1: $(call wrap,$(DBMAIN)/,fig1,.png)
fig2: $(call wrap,$(DBMAIN)/,fig2_alt,.png)

supfigs: $(call wrap,$(DBMAIN)/,irs_timing-realigned-SI irs_timing-realigned-SI-alt vc_4-panel-impact_10yr vc_4-panel-impact_50yr stacked_immunity pro_vs_re_active-prevalence.local foi_effect coverage_illustration,.png)

$(DBMAIN)/SI-fig-sero.png: SI-fig-sero.R $(DBMAIN)/stopping-sero-SI.rds
	Rscript $^ $@

$(DBMAIN)/irs_timing-realigned-SI.png: timing_plot_realigned-supplement.R $(DBMAIN)/stat.eff10.rds $(DBMAIN)/extra-stat.eff10.rds $(DBMAIN)/R0.rds $(DBMAIN)/mos.rds
	Rscript $^ $@

$(DBMAIN)/irs_timing-realigned-SI-alt.png: timing_plot_realigned-supplement-alt.R $(DBMAIN)/stat.eff10.rds $(DBMAIN)/R0.rds $(DBMAIN)/mos.rds
	Rscript $^ $@

$(DBMAIN)/vc_4-panel-impact_10yr.png $(DBMAIN)/vc_4-panel-impact_50yr.png: annual_effectiveness_plot.R $(DBMAIN)/stopping-baseline.cases.rds $(DBMAIN)/stopping-interventions.cases.rds $(DBMAIN)/stopping-averted.rds $(DBMAIN)/stopping-eff.rds
	Rscript $^ $@

$(DBMAIN)/foi_effect.png: foi_plot.R $(DBMAIN)/foi-eff.rds
	Rscript $^ $@

$(DBMAIN)/coverage_illustration.png: coverage_illustration.R $(call wrap,$(DBMAIN)/,eip R0 mos,.rds)
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

$(DBMAIN)/tab1.tex: main_tab_tex.R $(TABSRCS)
	Rscript $^ > $@

$(DBMAIN)/tab2.tex: main_tab2_tex.R $(DBMAIN)/foi-baseline.cases.rds $(DBMAIN)/foi-interventions.cases.rds $(DBMAIN)/foi-averted.rds
	Rscript $^ > $@

$(DBMAIN)/offset_peak.csv: offset.R $(DBMAIN)/tab1.csv
	Rscript $^ $@

$(DBMAIN)/tabX.tex: main_tabX_tex.R $(DBMAIN)/stopping-eff.rds $(DBMAIN)/stopping-sero.rds
	Rscript $^ > $@

$(DBMAIN)/confidence-intervals.rds: $(SQMAIN) $(DBROOT)/irs_insecticide-durability-effect_intro-fix.sqlite

$(DBMAIN)/check_runs.png: check_runs.R $(SQMAIN)
	Rscript $^ $@
