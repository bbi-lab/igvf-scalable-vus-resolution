# Figure 6b

# Reads the biobank carrier odds-ratio table via load_biobank_or_table.R,
# currently the checked-in All of Us stand-in
# (data/input/biobank/AoU-OR-estimates_2026-10-07.tsv.gz). See
# docs/biobank_analysis.md for the full workflow and why this will be
# swapped back to a plain read_csv() against IGVFFI3804AVJR.csv.gz once a
# fresh re-export of that file is available.

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
    legend.background = element_blank(),
    ggh4x.facet.nestline = element_line(linewidth = 1*point_in_mm, color="black"),
    plot.tag = element_text(face = 'bold')
  )

# Load main table
or_df <- load_biobank_or_table("../../../data/input/biobank/AoU-OR-estimates_2026-10-07.tsv.gz")

# Magic number: y-placement of the horizontal line separating the "pathogenic" and "benign" zones in this plot
PATH_BENIGN_BOUNDARY <- 11.5

# Broad gene-phenotype classes

gene_groups_df <- tribble(
  ~Gene, ~Disease, ~D2, ~D3,
  "BAP1", "Cancer", "Cancer", "Cancer",
  "BARD1", "Cancer", "Cancer", "Cancer",
  "BRCA1", "Cancer", "Cancer", "Cancer",
  "BRCA2", "Cancer", "Cancer", "Cancer",
  "CALM3", "Cardiovascular", "Cardio- vascular", "Cardio- vascular",
  "CHEK2", "Cancer", "Cancer", "Cancer",
  # D3's Metabolic rows are left as the plain "Metabolic" (not pre-split
  # like "Meta- bolic" below) -- str_wrap()'s default whitespace_only = TRUE
  # (confirmed against the installed stringr) only ever breaks at an
  # existing space, so a single word with no internal space, like this one,
  # is never wrapped regardless of panel width. Unlike Cardiovascular (see
  # "Cardio- vascular" below), this never needed splitting: Metabolic's
  # header is wide enough for "Metabolic" unwrapped regardless of how many
  # gene panels it spans.
  "G6PD", "Metabolic", "Meta- bolic", "Metabolic",
  "KCNE1", "Cardiovascular", "Cardio- vascular", "Cardio- vascular",
  "KCNH2", "Cardiovascular", "Cardio- vascular", "Cardio- vascular",
  "KCNQ4", "Hearing loss", "Hearing loss", "Hearing loss", #Rare disease",
  "LDLR", "Metabolic", "Meta- bolic", "Metabolic",
  "MSH2", "Cancer", "Cancer", "Cancer",
  "OTC", "Metabolic", "Meta- bolic", "Metabolic",
  "PALB2", "Cancer", "Cancer", "Cancer",
  "RAD51C", "Cancer", "Cancer", "Cancer",
  "SCN5A", "Cardiovascular", "Cardio- vascular", "Cardio- vascular",
  "TARDBP", "Rare disease", "Rare disease", "Rare disease",
  "TP53", "Cancer", "Cancer", "Cancer"
)


points_levels = c(
  "≤ -16", "≤ -15", "≤ -14", "≤ -13", "≤ -12", "≤ -11", "≤ -10", "≤ -9",
  "≤ -8", "≤ -7", "≤ -6", "≤ -5", "≤ -4", "≤ -3", "≤ -2", "≤ -1",
  #"0",
  "≥ +1", "≥ +2", "≥ +3", "≥ +4", "≥ +5", "≥ +6", "≥ +7", "≥ +8",
  "≥ +9", "≥ +10", "≥ +11", "≥ +12", "≥ +13", "≥ +14", "≥ +15", "≥ +16"
)

combined_points_plot_df <-
  or_df %>%
  filter(
    Dataset == 'Total points',
    `Cases with variants` > 0,
    Classification %in% points_levels,
    !Gene %in% c("TSC2", "GCK", "PTEN", "RAD51D")
  ) %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    significance = case_when(
      `p-value (Holm adjusted)` < 0.01 ~ 'Significant at 99% (Holm adjusted)',
      `p-value (Holm adjusted)` < 0.05 ~ 'Significant at 95% (Holm adjusted)',
      TRUE ~ 'Not significant at 95%'
    ),
    Classification = factor(
      Classification,
      levels = points_levels
    )
  )

combined_points_condensed_plot_df <-
  combined_points_plot_df %>%
  inner_join(
    combined_points_plot_df %>%
      group_by(Gene) %>%
      summarize(keep = any(significance != 'Not significant at 95%')) %>%
      filter(keep) %>%
      select(Gene)
  ) %>%
  left_join(gene_groups_df)

combined_points_condensed_plot <- ggplot(
  combined_points_condensed_plot_df,
  aes(
    x = `Odds Ratio`,
    xmin = OR_LI,
    xmax = OR_UI,
    y = Classification,
    fill = significance
  )
) +
  facet_nested(
    cols = vars(D3, Gene),
    scales = 'free',
    labeller = as_labeller(
      function(s) s %>% str_replace_all('_', ' ') %>% str_wrap(width=7)
    ),
    # Italicize only the inner (Gene) strip layer -- gene symbols are
    # conventionally italicized, D3's broad disease-category labels aren't.
    # by_layer_x = TRUE maps text_x's two elements to cols = vars(D3, Gene)'s
    # two nesting layers in order; NULL for the D3 layer means "use the
    # theme" (cfig_theme's strip.text below), same as every other property
    # element_text(face=...) here doesn't set.
    strip = strip_nested(
      by_layer_x = TRUE,
      text_x = list(NULL, element_text(face = 'bold.italic'))
    )
  ) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  geom_pointrange(position = position_dodge(width=0.5), shape = 21) +
  scale_x_log10(
    labels = scales::label_number(drop0trailing=TRUE),
    minor_breaks = NULL,
    guide = "axis_logticks",
  ) +
  scale_fill_manual(
    values = c(
      'Not significant at 95%' = 'white',
      'Significant at 95% (Holm adjusted)' = 'grey',
      'Significant at 99% (Holm adjusted)' = 'black'),
    breaks = c('Significant at 95% (Holm adjusted)', 'Significant at 99% (Holm adjusted)')
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_hline(yintercept = PATH_BENIGN_BOUNDARY, linetype = 'dotted') +
  geom_blank(aes(x = 6))

# Theme adjustments
cfig_theme <- nature_theme +
  theme(
    strip.text = element_text(
      size = 6,
      face = 'bold',
      color = 'black',
      margin = margin(0,0,0,0),
      family = "Arial"
    ),
    strip.clip = 'off',
    panel.spacing = unit(2, 'pt'),
    axis.title.x = element_text(hjust = 0.2),
    legend.position = 'bottom',
    legend.justification = 0.65,
    legend.box.spacing = unit(-18, 'pt')
  )

fig6b_plot = combined_points_condensed_plot +
  cfig_theme +
  theme(
    axis.text.x = element_text(angle = 60),
    # BARD1 (left-most panel, the only one that prints y-axis tick labels)
    # has up to 24 categories -- at this figure's 1.55in panel height
    # (see the ggsave height comment below) nature_theme's inherited 6pt
    # leaves too little room per row, and the <=/>= glyphs render
    # overlapping between adjacent rows at sizes above 4pt. Set to 5pt
    # anyway per explicit request, accepting that overlap.
    axis.text.y = element_text(size = 5)
  )

# Show plot (interactive only -- under Rscript this would otherwise open
# the default graphics device and leave behind a stray Rplots.pdf)
if (interactive()) print(fig6b_plot)

# ggh4x::force_panelsizes(rows = unit(1.55, 'in')) would be the direct way to
# pin the panel height, but it renders the y-axis's <=/>= glyphs
# overlapping/stacked across rows at this panel height (confirmed with a
# minimal repro -- a rendering bug in the two-pass layout it forces, not
# something specific to this plot). Hitting the same panel height by instead
# tuning the device height avoids that pass entirely: every other element
# (strips, axis text/title, legend) has a fixed absolute size, so growing or
# shrinking the overall height changes only the panel -- 66.6mm total here
# measured out to a 1.55in-tall panel against the fixed ~27.2mm the rest of
# the figure takes up.
ggsave(
  '../../../data/output/figures/assets/figure_6/figure_6b.pdf',
  fig6b_plot,
  width = 183, # Max 183
  height = 66.6,
  units = 'mm',
  device = cairo_pdf,
  create.dir = TRUE)
ggsave(
  '../../../data/output/figures/assets/figure_6/figure_6b.svg',
  fig6b_plot,
  width = 183, # Max 183
  height = 66.6,
  units = 'mm',
  create.dir = TRUE)
