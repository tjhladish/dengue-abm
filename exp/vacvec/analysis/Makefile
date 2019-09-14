FIGDIR := fig
TEXDIR := tex

vpath %.rds rds
vpath %.rda rds
vpath %.png $(FIGDIR)
vpath %.tiff $(FIGDIR)
vpath %.gif $(FIGDIR)
vpath %.tex $(TEXDIR)

-include local.mk

wrap = $(addsuffix $(3),$(addprefix $(1),$(2)))

default: figures

.PHONY: parse showfigs

parse:
	+$(MAKE) -C rds

cleanparse:
	+$(MAKE) -C rds clean

%.rds %.rda:
	+$(MAKE) -C rds $@

RFIG = Rscript $^ $|/$@

$(FIGDIR):
	mkdir $@

figref.rda: figref.R utils.R projref.rda
	Rscript $^ $@

IMGTYPES := png tiff eps pdf

$(addprefix fig_2.,$(IMGTYPES)): fig_raw_effectiveness.R figref.rda nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

$(addprefix fig_3.,$(IMGTYPES)): fig_combo_eff_illustration.R figref.rda nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

$(addprefix fig_4.,$(IMGTYPES)): fig_full_combo_interaction.R figref.rda nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

$(addprefix fig_4_alt.,$(IMGTYPES)): fig_alt_combo_interaction.R figref.rda alt_effstats.rds | $(FIGDIR)
	$(RFIG)

$(addprefix fig_5.,$(IMGTYPES)): fig_lagged_interventions.R figref.rda lag_effstats.rds nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

$(addprefix fig_6.,$(IMGTYPES)): fig_foi_alternatives.R figref.rda nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

# single intervention IQs
$(addprefix SIfig_2.,$(IMGTYPES)): SIfig_raw_effectiveness.R figref.rda nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

# combination dels from naive effectiveness
$(addprefix SIfig_3.,$(IMGTYPES)): SIfig_combo_gaps.R figref.rda nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

# combination IQs
$(addprefix SIfig_4.,$(IMGTYPES)): SIfig_all_combos.R figref.rda nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

# lagged IQs
$(addprefix SIfig_5.,$(IMGTYPES)): SIfig_lagged.R figref.rda lag_effstats.rds nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

# foi study IQs
$(addprefix SIfig_6.,$(IMGTYPES)): SIfig_foi.R figref.rda nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

# foi study IQs
$(addprefix SIfig_7.,$(IMGTYPES)): SIfig_foi_v_mpop.R figref.rda foi_v_mpop.rds | $(FIGDIR)
	$(RFIG)

$(addprefix SIfig_8.,$(IMGTYPES)): SIfig_foi_both_cydtdv.R figref.rda nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

$(addprefix SIfig_9a.,$(IMGTYPES)): SIfig_test_sens.R figref.rda nolag_effstats.rds | $(FIGDIR)
	$(RFIG)

SIfig_9b.png SIfig_9c.png SIfig_9d.png: SIfig_9a.png

SIfig_9b.eps SIfig_9c.eps SIfig_9d.eps: SIfig_9a.eps

SIfig_9b.pdf SIfig_9c.pdf SIfig_9d.pdf: SIfig_9a.pdf

# lagged showing single interventions
non_fig_1.png non_fig_1.tiff: fig_lagged_interventions.R figref.rda lag_effstats.rds effstats.rds | $(FIGDIR)
	$(RFIG)


MTFIGS := $(addprefix fig_,$(addsuffix .pdf,2 3 4 5 6))
SIFIGS := $(addprefix SIfig_,$(addsuffix .pdf,2 3 4 5 6 7 8 9a 9b 9c 9d))

maintextfigs: $(MTFIGS)

sifigs: $(SIFIGS)

figures: maintextfigs sifigs

ALLFIGS := $(MTFIGS) $(SIFIGS)

showfigs:
	open $(addprefix fig/,$(ALLFIGS))

clean_figures:
	rm -rf $(FIGDIR)

$(TEXDIR):
	mkdir $@

numerical_results.tex: texgen.R defs.template effstats.rds | $(TEXDIR)
	$(RFIG)

tab_pars.tex: tex_posterior_params.R utils.R posterior.rds parkey.json | $(TEXDIR)
	$(RFIG)

tab_mets.tex: tex_posterior_metrics.R utils.R posterior.rds metkey.json | $(TEXDIR)
	$(RFIG)

textabs: tab_pars.tex tab_mets.tex

pullstats: pullstats.R nolag_effstats.rds lag_effstats.rds foi_v_mpop.rds
	@Rscript $^
