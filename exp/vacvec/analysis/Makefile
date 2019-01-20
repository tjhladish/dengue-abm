vpath %.rds rds
vpath %.rda rds
vpath %.png fig
vpath %.tiff fig
vpath %.gif fig

-include local.mk

wrap = $(addsuffix $(3),$(addprefix $(1),$(2)))

default: parse

.PHONY: parse

parse:
	make -C rds

cleanparse:
	make -C rds clean

%.rds:
	make -C rds $@

fig:
	mkdir $@

figref.rda: figref.R utils.R projref.rda
	Rscript $^ $@

alt_fig_2.png fig_2.png fig_2.tiff: fig_raw_effectiveness.R figref.rda effstats.rds | fig
	Rscript $^ fig/$@

# TODO add incidence as a top row? move legend?
SIfig_2.png SIfig_2.tiff: SIfig_raw_effectiveness.R figref.rda effstats.rds | fig
	Rscript $^ fig/$@

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

# lag_plot.png: lag_plot.R lag_effstats.rds effstats.rds
# 	Rscript $^ $@
#
# foi_plot.png: foi_plot.R foi_baseline.rds foi_intervention.rds foi_effectiveness.rds
# 	Rscript $^ $@
#
# foi_fig_2.png: foi_fig_2.R foi_effstats.rds effstats.rds
# 	Rscript $^ $@
#
# fig_1.png fig_1.tiff: fig_1.R baseline.rds intervention.rds effstats.rds
# 	Rscript $^ $@
#
# # show only the * panels of fig 1
# #  if % = vac / vec, then show just the vaccine-only (vector control only, resp) results
# #  if % = inc / eff, then show just the incidence (effectiveness, resp) results
# fig_1_vac.png fig_1_vec.png fig_1_inc.png fig_1_eff.png: fig_1.R baseline.rds intervention.rds effstats.rds
# 	Rscript $^ only=$(patsubst fig_1_%.png,%,$@) $@
#
# clean_fig_1_slices: | fig_1_*.png
# 	rm $|
#
# combo_fig_mech_facet.png: combo_fig_mech_facet.R effstats.rds
# 	Rscript $^ $@
#
# # real fig 2
# combo_fig_coverage_facet.png: combo_fig_coverage_facet.R effstats.rds
# 	Rscript $^ $@
#
# SI_IQ_fig.png: IQ_fig.R effectiveness.rds
# 	Rscript $^ $@
#
# SI_ceff_fig.png: ceff_fig.R effectiveness.rds comboeff.rds
# 	Rscript $^ $@
#
# .PHONY: figures
#
# MANUSCRIPTFIGS := $(addsuffix .png,fig_1 combo_fig_coverage_facet foi_fig_2 lag_plot SI_IQ_fig SI_ceff_fig)
#
# figures: $(MANUSCRIPTFIGS)
#
# clean_figs: $(MANUSCRIPTFIGS)
# 	rm $^
#
# # christmas figures
# # fig 1 - diagram, probably not managed via make
#
# # fig 2: single intervention incidence, annual effectiveness, cumulative effectiveness
# # row 1: incidence
# # row 2: ann. effectiveness
# # row 3: cum. effectiveness
# # one col; show only "best" versions of intervention (w/ catchup, max vc)
# # distinguish vaccine by line type = solid, none; dashes, EDV; dotted, Dengvaxia-like
# # complementary SI fig:
# # multiple cols: one for each vaccine type (none, EDV, Dengvaxia-like)
# # IQs
#
# new_SI_fig_2.png new_SI_fig_2.tiff new_fig_2.png new_fig_2.tiff: new_fig_2.R baseline.rds intervention.rds effectiveness.rds
# 	Rscript $^ $@
# 	# 1 col, 3 rows - incidence, annual effectiveness, cumulative effectiveness
# 	# both vaccines w/ catchup, 75% vc - IQs?
#
# new_fig_3.png new_fig_3.tiff: new_fig_3.R effectiveness.rds comboeff.rds
# 	Rscript $^ $@
# 	# 2 rows, 3 cols - annual effectiveness of combination interventions, with background ref lines
# 	# for single interventions; rows are vaccine type, cols are vc level
#
# new_fig_4.png new_fig_4.tiff: new_fig_4.R effectiveness.rds
# 	Rscript $^ $@
# 	# repeat fig 3, but with shaded effect regions?
# 	# for single interventions; rows are vaccine type, cols are vc level
#
# new_fig_5.png new_fig_5.tiff: new_fig_5.R lag_effectiveness.rds
# 	Rscript $^ $@
# 	# lagged interventions
