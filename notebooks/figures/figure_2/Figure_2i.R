# Figure 2i

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

# Add a separate "Functionally Abnormal (<consequence>)" row alongside the
# unbroken-down "Functionally Abnormal (All)" row for each consequence
# listed here, so consequence-specific estimates (often the ones driving
# classification decisions) are visible next to the overall one. Set to
# c() for the simpler 2-category (Normal/Abnormal) version instead.
ABNORMAL_BREAKDOWN_CONSEQUENCES <- c('Missense', 'Truncating')

# Load main table
or_df <- load_biobank_or_table(
  "../../../data/input/biobank/AoU-OR-estimates_2026-10-07.tsv.gz",
  keep_consequence = length(ABNORMAL_BREAKDOWN_CONSEQUENCES) > 0
)

# Filter to IGVF functional assays
func_df <- or_df %>%
  filter(
    str_ends(Dataset, '_IGVF'),
    Classifier == 'StandardizedClass',
    Classification %in% c('NORMAL', 'ABNORMAL'),
    `Cases with variants` > 0
  )

if (length(ABNORMAL_BREAKDOWN_CONSEQUENCES) > 0) {
  func_df <- func_df %>%
    filter(
      Consequence == 'All' | (Classification == 'ABNORMAL' & Consequence %in% ABNORMAL_BREAKDOWN_CONSEQUENCES)
    ) %>%
    mutate(
      Classification = factor(
        case_when(
          Classification == 'NORMAL' ~ 'Functionally Normal',
          Consequence == 'All' ~ 'Functionally Abnormal (All)',
          TRUE ~ str_c('Functionally Abnormal (', Consequence, ')')
        ),
        levels = c(
          'Functionally Normal',
          'Functionally Abnormal (All)',
          str_c('Functionally Abnormal (', ABNORMAL_BREAKDOWN_CONSEQUENCES, ')')
        )
      )
    )
} else {
  # load_biobank_or_table() already restricted to Consequence == 'All'
  # (and dropped that column) when keep_consequence = FALSE above.
  func_df <- func_df %>%
    mutate(
      Classification = factor(
        str_c('Functionally ', str_to_title(Classification)),
        levels = c('Functionally Normal', 'Functionally Abnormal')
      )
    )
}

func_df <- func_df %>%
  filter(Gene != 'TSC2') %>% # TSC2 excluded -- remove this filter to re-include
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    Gene = factor(
      Gene,
      levels = c('BARD1', 'PALB2', 'RAD51D', 'XRCC2', 'CTCF', 'SFPQ')
    )
  )

# Calc limits for small plots
limits_df <- func_df %>%
  summarise(
    OR_LI = min(OR_LI),
    OR_UI = max(OR_UI)
  )

# Make plot
make_func_plot <- function(plot_df, x_limits)
  ggplot(
    plot_df,
    aes(
      x=`Odds Ratio`,
      xmin=OR_LI,
      xmax=OR_UI,
      y=Classification
    )
  ) +
  facet_grid(cols = vars(Gene), scales = 'free_x') +
  geom_pointrange() +
  geom_errorbar(width = 0.2) +
  scale_x_log10(
    #labels = scales::label_number(accuracy = 0.1),
    minor_breaks = NULL,
    #minor_breaks = scales::minor_breaks_log(),
    guide = "axis_logticks"
  ) +
  geom_blank(aes(x = x_limits)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10)) +
  geom_vline(xintercept = 1, linetype = 'dashed')

fig2i_plot <- make_func_plot(func_df, deframe(limits_df))


# Show plot
print(fig2i_plot + nature_theme)


# Save plot -- each added breakdown category needs ~15mm more vertical
# room per panel than the base 2-category version (30mm) to keep the
# wrapped y-axis labels from colliding.
ggsave(
  '../../../data/output/figures/assets/figure_2/figure_2i.pdf',
  fig2i_plot + nature_theme,
  width = 100,
  height = 30 + 15 * length(ABNORMAL_BREAKDOWN_CONSEQUENCES),
  units = 'mm',
  family = 'Arial',
  device = cairo_pdf,
  create.dir = TRUE
)
