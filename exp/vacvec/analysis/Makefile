vpath %.rds rds
vpath %.rda rds
vpath %.png fig
vpath %.tiff fig
vpath %.gif fig

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

FIGDIR := fig
RFIG = Rscript $^ $(FIGDIR)/$@

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

MTFIGS := fig_2.png fig_3.png fig_4.png fig_5.png fig_6.png
SIFIGS := SIfig_2.png SIfig_4.png SIfig_5.png SIfig_6.png

maintextfigs: $(MTFIGS)

sifigs: $(SIFIGS)

figures: maintextfigs sifigs

ALLFIGS := $(MTFIGS) $(SIFIGS)

showfigs:
	open $(addprefix fig/,$(ALLFIGS))

clean_figures:
	rm -rf $(FIGDIR)

## FIGURES

# TODO: automate figure 1 generation:
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
#
# NIGHTLIGHT := some.csv # some csv? shapefile?
# ADMINBOUNDS := admin.shp # seems like shapefile
# LOCATIONS := vertices.csv # some csv? shapefile?
# CONNECTIONS := edges.csv # some csv? shapefile?
# TARRECT := tar.json # something specifying lat / long extent of zoom, other annotations (e.g., vertical reveal?)
# 
# fig1p1.png: diagram_panel_wideview.R $(NIGHTLIGHT) $(ADMINBOUNDS) $(LOCATIONS) $(CONNECTIONS) $(TARRECT)
# 	Rscript $^ $@
# 
# fig1p2.png: diagram_panel_zoom.R $(LOCATIONS) $(CONNECTIONS) $(TARRECT) # $(ADMINBOUNDS)? in case the zoom rect crosses admin lines?
# 	Rscript $^ $@
# 
# fig1p3back.png: diagram_panel_model_background.R $(LOCATIONS) $(CONNECTIONS) $(TARRECT) # $(ADMINBOUNDS)? in case the zoom rect crosses admin lines?
# 	Rscript $^ $@
# 
# fig1.png: fig1p1.png fig1p2.png fig1p3.png fig1p4.png
# 	@echo $^

numerical_results.tex: texgen.R defs.template effstats.rds
	Rscript $^ $@