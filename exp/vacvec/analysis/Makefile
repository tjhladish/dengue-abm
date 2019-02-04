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
	make -C rds

cleanparse:
	make -C rds clean

%.rds %.rda:
	make -C rds $@

RFIG = Rscript $^ $|/$@

$(FIGDIR):
	mkdir $@

figref.rda: figref.R utils.R projref.rda
	Rscript $^ $@

fig_2.png fig_2.tiff: fig_raw_effectiveness.R figref.rda effstats.rds | $(FIGDIR)
	$(RFIG)

fig_3.png fig_3.tiff: fig_combo_eff_illustration.R figref.rda effstats.rds | $(FIGDIR)
	$(RFIG)

fig_4.png fig_4.tiff: fig_full_combo_interaction.R figref.rda effstats.rds | $(FIGDIR)
	$(RFIG)

fig_5.png fig_5.tiff: fig_lagged_interventions.R figref.rda lag_effstats.rds effstats.rds | $(FIGDIR)
	$(RFIG)

fig_6.png fig_6.tiff: fig_foi_alternatives.R figref.rda foi_effstats.rds effstats.rds
	$(RFIG)

# single intervention IQs
SIfig_2.png SIfig_2.tiff: SIfig_raw_effectiveness.R figref.rda effstats.rds | $(FIGDIR)
	$(RFIG)

# combination dels from naive effectiveness
SIfig_3.png SIfig_3.tiff: SIfig_combo_gaps.R figref.rda effstats.rds | $(FIGDIR)
	$(RFIG)

# combination IQs
SIfig_4.png SIfig_4.tiff: SIfig_all_combos.R figref.rda effstats.rds | $(FIGDIR)
	$(RFIG)

# lagged IQs
SIfig_5.png SIfig_5.tiff: SIfig_lagged.R figref.rda lag_effstats.rds effstats.rds | $(FIGDIR)
	$(RFIG)

# foi study IQs
SIfig_6.png SIfig_6.tiff: SIfig_foi.R figref.rda foi_effstats.rds effstats.rds
	$(RFIG)

# foi study IQs
SIfig_7.png SIfig_7.tiff: SIfig_foi_v_mpop.R figref.rda foi_v_mpop.rds
	$(RFIG)


MTFIGS := fig_2.png fig_3.png fig_4.png fig_5.png fig_6.png
SIFIGS := SIfig_2.png SIfig_3.png SIfig_4.png SIfig_5.png SIfig_6.png SIfig_7.png

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

posterior_par.tex: tex_posterior_params.R posterior.rds | $(TEXDIR)
	$(RFIG)

posterior_met.tex: tex_posterior_metrics.R posterior.rds | $(TEXDIR)
	$(RFIG)