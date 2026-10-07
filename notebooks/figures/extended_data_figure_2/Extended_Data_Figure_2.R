# Extended Data Figure 2

# This script requires IGVFFI3804AVJR.csv.gz which is available from
# https://data.igvf.org/tabular-files/IGVFFI3804AVJR/ and expected at
# data/input/biobank/IGVFFI3804AVJR.csv.gz. Pending a re-export of that
# file, it currently instead reads an AoU-only stand-in (see
# load_biobank_or_table.R) -- swap the read below back to a plain
# read_csv() once IGVFFI3804AVJR.csv.gz is available again.

# Load libraries
library(tidyverse)
library(patchwork)
library(ggh4x)
library(extrafont)

source('../../../src/lib/figures/load_biobank_or_table.R')

# Font setup
loadfonts(device = 'all')

# Figure theme
point_in_mm = 0.3527778
nature_theme <- theme_linedraw() +
  theme(
    text = element_text(family = 'Arial', size = 7),
    line = element_line(linewidth = 1*point_in_mm),
    geom = element_geom(
      linewidth = 1*point_in_mm,
      borderwidth = 1*point_in_mm,
      pointsize = 3*point_in_mm
    ),
    #spacing = unit(3, 'points'),
    axis.ticks = element_line(linewidth = 0.5*point_in_mm),
    #axis.title = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7, face = 'bold', color = 'black', margin = margin(3,3,3,3)),
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    ggh4x.facet.nestline = element_line(linewidth = 1*point_in_mm, color="black"),
    plot.tag = element_text(face = 'bold')
  )

# Load main table -- keep_consequence = TRUE so panel a (condensed_assay_plot)
# can overlay a Missense-only series on the aggregate ("All") one; every
# other use below restricts back to Consequence == 'All' to stay unaffected.
or_df <- load_biobank_or_table(
  "../../../data/input/biobank/AoU-OR-estimates_2026-10-06_merged.tsv.gz",
  keep_consequence = TRUE
)

# Broad gene-phenotype classes

gene_groups_df <- tribble(
  ~Gene, ~Disease, ~D2, ~D3,
  "BAP1", "Cancer", "Cancer", "Cancer",
  "BARD1", "Cancer", "Cancer", "Cancer",
  "BRCA1", "Cancer", "Cancer", "Cancer",
  "BRCA2", "Cancer", "Cancer", "Cancer",
  "CALM3", "Cardiovascular", "Cardio- vascular", "Cardiovascular",
  "CHEK2", "Cancer", "Cancer", "Cancer",
  "G6PD", "Metabolic", "Meta- bolic", "Meta- bolic",
  "GCK", "Metabolic", "Meta- bolic", "Meta- bolic",
  "KCNE1", "Cardiovascular", "Cardio- vascular", "Cardiovascular",
  "KCNH2", "Cardiovascular", "Cardio- vascular", "Cardiovascular",
  "KCNQ4", "Hearing loss", "Hearing loss", "Hearing loss", #Rare disease",
  "LDLR", "Metabolic", "Meta- bolic", "Meta- bolic",
  "MSH2", "Cancer", "Cancer", "Cancer",
  "OTC", "Metabolic", "Meta- bolic", "Meta- bolic",
  "PALB2", "Cancer", "Cancer", "Cancer",
  "PTEN", "Cancer", "Cancer", "Cancer",
  "RAD51C", "Cancer", "Cancer", "Cancer",
  "RAD51D", "Cancer", "Cancer", "Cancer",
  "SCN5A", "Cardiovascular", "Cardio- vascular", "Cardiovascular",
  "TARDBP", "Rare disease", "Rare disease", "Rare disease",
  "TP53", "Cancer", "Cancer", "Cancer",
  "TSC2", "Cancer", "Cancer", "Cancer"
)

# Predictor plot component

# Order levels
vep_levels = c(
  "≤ -4", "≤ -3", "≤ -2", "≤ -1",
  #"0",
  "≥ +1", "≥ +2", "≥ +3", "≥ +4"
)

# Select predictor results
predictor_plot_df <- or_df %>%
  filter(
    Consequence == 'All',
    `Cases with variants` > 0,
    Classifier %in% c(
      'Gene-specific',
      'Genome-wide aggregation'
    ),
    Classification %in% vep_levels
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    significance = if_else(
      (LogOR_LI > 0) | (LogOR_UI < 0),
      'Significant at 95%',
      'Not significant at 95%'
    ),
    Classification = factor(
      Classification,
      levels = vep_levels
    ),
    Classifier = Classifier %>% fct(
      levels = c(
        'Gene-specific',
        #'Domain aggregation',
        'Genome-wide aggregation'
      )
    )
    ,
    Dataset = fct(
      Dataset,
      levels = c('REVEL', 'AlphaMissense', 'MutPred2')
    )
  )

# Colors
color_mapping <- c(
  `Gene-specific` = "#CC6015", 
  `Domain aggregation` = "#D2A624",
  `Genome-wide aggregation` = "#544439"
)


# Plotting function
make_predictor_plot <- function(
    in_df,
    scales = 'free_x',
    legend_position = 'bottom'
)
  ggplot(
    in_df,
    aes(
      x=`Odds Ratio`,
      xmin=OR_LI,
      xmax=OR_UI,
      y=Classification,
      color = Classifier,
      shape = significance
    )
  ) +
  facet_grid(
    cols = vars(Gene),
    rows = vars(Dataset),
    scales = scales
  ) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  geom_pointrange(position = position_dodge(width=0.5), fill = 'white') +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks"
  ) +
  scale_shape_manual(
    values = c(
      'Not significant at 95%' = 21,
      'Significant at 95%' = 16),
    breaks = c('Significant at 95%'),
    guide = guide_legend(position = legend_position)
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_color_manual(
    values = color_mapping,
    guide = guide_legend(
      position = legend_position,
      override.aes = aes(shape = 21, fill = 'white')
    )
  )

# Make plot frame
figure_predictor_plot_df <- predictor_plot_df %>%
  filter(
    Gene %in% c('BRCA1', 'BRCA2', 'MSH2', 'TP53', 'TSC2'), # Have gene-specific calibration
  )

# Make plot for final figure
figure_predictor_plot <- figure_predictor_plot_df %>%
  make_predictor_plot(legend_position = 'bottom') +
  geom_blank(aes(x=deframe(
    predictor_plot_df %>%
      filter(Gene == 'MSH2') %>%
      summarize(max(OR_UI))
  )))


# Assay plot component

# Select assays for plot
assay_classification_levels = c(
  "NORMAL", "ABNORMAL",
  "≤ -8", "≤ -7", "≤ -6", "≤ -5", "≤ -4", "≤ -3", "≤ -2", "≤ -1",
  "0",
  "≥ +1", "≥ +2", "≥ +3", "≥ +4", "≥ +5", "≥ +6", "≥ +7", "≥ +8"
)

assay_plot_df <-
  or_df %>%
  filter(
    # Missense is kept alongside All so condensed_assay_plot_df below can
    # overlay a missense-only series on the aggregate one -- every other
    # per-consequence breakdown (Stop gained/Deletion/Splicing/...) isn't
    # needed here and is dropped.
    Consequence %in% c('All', 'Missense'),
    Classifier %in% c('ExCALIBR', 'StandardizedClass'),
    `Cases with variants` > 0,
    Classification %in% assay_classification_levels
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    significance = if_else(
      (LogOR_LI > 0) | (LogOR_UI < 0),
      'Significant at 95%',
      'Not significant at 95%'
    ),
    Classification = factor(
      Classification,
      levels = assay_classification_levels
    )
  )

condensed_assay_datasets <- c(
  'BAP1_Waters_2024',
  'BARD1_IGVF',
  #'BRCA1_Adamovich_2022_Cisplatin_Resistance',
  #'BRCA1_Adamovich_2022_HDR',
  'BRCA1_Findlay_2018',
  'BRCA2_Hu_2024',
  #'BRCA2_Sahu_2025_SGE',
  'GCK_Gersing_2023_complementation',
  'KCNH2_Jiang_2022',
  'KCNQ4_Zheng_2022_current_homozygous',
  'MSH2_Jia_2021',
  'PALB2_IGVF',
  'RAD51C_Olvera-León_2024',
  'RAD51D_IGVF',
  'SCN5A_Ma_2024',
  'TP53_Fayer_2021_meta',
  #'TP53_Boettcher_2019'
  #'TP53_Fortuno_2021',
  #'TP53_Giacomelli_2018_combined_score'
  'TSC2_IGVF'
)

# Sample sizes shown above each gene panel: missense vs. other (derived as
# All minus Missense -- ExCALIBR has no direct "non-missense" consequence
# row the way "Functional class" does) carrier counts among cases/controls.
# The ExCALIBR score bins are cumulative, so there's no single row that
# already holds "every carrier of this consequence" -- "<= -1"/">= +1" are
# the least-restrictive (broadest) bin in each direction, so summing just
# those two covers everyone with a non-zero score. Carrier counts are
# parsed from the source data with privacy-censored "<= 20" cells treated
# as 20 (see parse_censored_count() in load_biobank_or_table.R) -- fine for
# this descriptive annotation, not precise enough for anything computed.
# Short-scale (17.8K-style) number formatting -- these panels are narrow
# (5 across), and full comma-formatted digit strings don't fit.
count_fmt <- scales::label_number(scale_cut = scales::cut_short_scale())

assay_dataset_counts <- or_df %>%
  filter(
    Dataset %in% condensed_assay_datasets,
    Classifier == 'ExCALIBR',
    Consequence %in% c('All', 'Missense'),
    Classification %in% c('≤ -1', '≥ +1')
  ) %>%
  group_by(Dataset, Gene, Consequence) %>%
  summarise(
    `Carrier cases` = sum(`Carrier cases`),
    `Carrier controls` = sum(`Carrier controls`),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = Consequence,
    values_from = c(`Carrier cases`, `Carrier controls`)
  )

condensed_assay_plot_df <-
  assay_plot_df %>%
  filter(
    Classifier == 'ExCALIBR',
    Dataset %in% condensed_assay_datasets,
    Classification != "0"
  ) %>%
  mutate(
    Consequence = factor(Consequence, levels = c('All', 'Missense'))
  ) %>%
  left_join(gene_groups_df)

# Per-gene left edge for the annotation's left-aligned header below: `x =
# -Inf` would be the natural way to pin it to the panel's left border, but
# with scale_x_log10() that silently drops the layer (log10(-Inf) is NaN,
# not -Inf -- confirmed with a minimal repro), unlike `x = Inf` for the
# right-aligned counts, which works fine (log10(Inf) = Inf). Anchoring at
# each panel's own leftmost plotted value instead sidesteps that and lands
# in the same place anyway.
assay_panel_x_min <- condensed_assay_plot_df %>%
  group_by(Gene) %>%
  summarise(panel_x_min = min(OR_LI, na.rm = TRUE), .groups = 'drop')

# One row per gene per annotation line: a left-aligned "Cases / Controls"
# header sharing its row with the right-aligned Other (All minus Missense)
# count, and a second row below it for the right-aligned Missense count in
# blue (matching its series in the plot). This lands below the native
# per-gene strip (gene name, bold) and above the panel's own plotted data --
# not above the strip/Disease header -- because that strip is given extra
# bottom margin below (reserving blank room for this content within its own
# cell) rather than this layer having to clear the whole strip height, as a
# strip and the panel directly below it are flush with no gap by default.
# `vjust` (more negative = further above the panel) only needs to be large
# enough to clear into that reserved margin via coord_cartesian(clip =
# 'off') below -- actual values tuned against this plot's real panel/strip
# size, not transferable to a plot with different dimensions.
assay_count_annotations_df <- bind_rows(
  assay_dataset_counts %>%
    transmute(
      Gene,
      part = 'header',
      hjust = 0,
      vjust = -2.0,
      # Smaller than the counts (4pt vs 5pt) -- "Cases/Controls" sharing a
      # row with a wide count (KCNQ4's "28.318K / 210.2K") collided at
      # equal size; this is a caption next to data, not data itself, so
      # shrinking it a bit is a reasonable way to buy the needed room back.
      size = 4,
      color = 'black',
      text = 'Cases/Controls'
    ),
  assay_dataset_counts %>%
    transmute(
      Gene,
      part = 'other',
      hjust = 1,
      vjust = -2.0,
      size = 5,
      color = 'black',
      text = sprintf(
        '%s / %s',
        count_fmt(`Carrier cases_All` - `Carrier cases_Missense`),
        count_fmt(`Carrier controls_All` - `Carrier controls_Missense`)
      )
    ),
  assay_dataset_counts %>%
    transmute(
      Gene,
      part = 'missense',
      hjust = 1,
      vjust = -0.6,
      size = 5,
      color = '#1D7AAB',
      text = sprintf(
        '%s / %s',
        count_fmt(`Carrier cases_Missense`),
        count_fmt(`Carrier controls_Missense`)
      )
    )
) %>%
  left_join(assay_panel_x_min, by = 'Gene') %>%
  # facet_nested_wrap facets on both Disease and Gene -- without Disease
  # here too, ggplot can't route each row to its one matching panel and
  # broadcasts it into every panel instead (confirmed: this is exactly what
  # happened before this join was added).
  left_join(gene_groups_df, by = 'Gene') %>%
  mutate(x = if_else(part == 'header', panel_x_min, Inf))

# Limits for most panels
assay_plot_common_limits <- condensed_assay_plot_df %>%
  filter(!(Gene %in% c("BRCA1", "MSH2", "KCNH2", "TSC2", "GCK"))) %>%
  summarise(
    OR_LI = min(OR_LI),
    OR_UI = max(OR_UI)
  ) %>% deframe()

# Build assay plot
condensed_assay_plot <- ggplot(
  condensed_assay_plot_df,
  aes(
    x=`Odds Ratio`,
    xmin=OR_LI,
    xmax=OR_UI,
    y=Classification,
    shape = significance,
    color = Consequence
  )
) +
  facet_nested_wrap(
    facets = vars(Disease, Gene),
    nrow = 3,
    scales = 'free_x',
    # Default solo_line = FALSE only draws the nest line under a Disease
    # header spanning 2+ genes (e.g. "Cancer"), not a single-gene one (e.g.
    # "Cardiovascular" above just KCNH2) -- that asymmetry misaligned the
    # single-gene headers (no line reserving/centering their text the same
    # way). Drawing it for every Disease header, solo or not, equalizes them.
    solo_line = TRUE,
    # Disease keeps the default strip style; the per-gene strip gets extra
    # bottom margin, reserving blank room below the gene name (within its
    # own cell) for the Cases/Controls annotation below, so that content
    # lands below the gene name and above the panel -- not above the
    # Disease header, and not inside the panel either.
    strip = strip_nested(
      text_x = list(
        element_text(size = 7, face = 'bold'),
        element_text(size = 7, face = 'bold', margin = margin(t = 3, r = 3, b = 18, l = 3))
      ),
      # strip_nested()'s by_layer_x defaults to FALSE, which does NOT treat
      # the text_x list above as one element per nesting depth (Disease,
      # then Gene) -- without this, the margin meant for just the Gene
      # level was applied inconsistently, stretching the Disease level's
      # own cell and leaving dead space below "Cancer" etc. above its own
      # separator line.
      by_layer_x = TRUE
    )
  ) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  geom_pointrange(position = position_dodge(width=0.5), fill='white') +
  # Sample-size annotation (Cases/Controls header + counts), positioned
  # below each panel's native gene-name strip -- in the room that strip's
  # own enlarged bottom margin reserves for it above -- rather than inside
  # the panel, via coord_cartesian(clip = 'off') below, which lets this
  # layer's content overflow past the panel's own border into that margin.
  geom_text(
    data = assay_count_annotations_df,
    aes(x = x, y = Inf, label = text, hjust = hjust, vjust = vjust, size = I(size / .pt), color = I(color)),
    inherit.aes = FALSE,
    # Without this, ggplot merges this layer into the Consequence color
    # legend below (it maps color too, even via I()) and draws geom_text's
    # placeholder key glyph -- literally the letter "a" -- on top of the
    # All/Missense keys.
    show.legend = FALSE
  ) +
  coord_cartesian(clip = 'off') +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks",
  ) +
  scale_shape_manual(
    values = c(
      'Not significant at 95%' = 21,
      'Significant at 95%' = 16),
    breaks = c('Significant at 95%'),
    guide = guide_legend(position = 'bottom')
  ) +
  scale_color_manual(
    values = c('All' = 'black', 'Missense' = '#1D7AAB'),
    guide = guide_legend(
      position = 'bottom',
      # Matches panel b's own color legend (make_predictor_plot above):
      # without this, the key glyph falls back to a solid dot instead of
      # the open (fill = 'white') circle every point in this panel uses.
      override.aes = aes(shape = 21, fill = 'white')
    )
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed')

# Compose figure together
fig_exd2 <- (condensed_assay_plot + nature_theme +
    # Extra room above each row of panels for the sample-size annotation
    # (geom_text + coord_cartesian(clip = 'off') above) to occupy without
    # overlapping the panel row above it.
    theme(panel.spacing.y = unit(4, 'mm'))) +
  (figure_predictor_plot + nature_theme + theme(plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))) +
  plot_annotation(tag_levels='a') +
  plot_layout(
    design = c(
      area(1,1,6,6),
      area(7,1,9,6)
    )
  )

# Show figure
print(fig_exd2)

# Save figure
ggsave(
  '../../../data/output/figures/assets/extended_data_figure_2.pdf',
  fig_exd2,
  width = 160, # Max 183
  height = 247,
  units = 'mm',
  device = pdf, # JS 20260714
  # device = cairo_pdf)
  create.dir = TRUE)
