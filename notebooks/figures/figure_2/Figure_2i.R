# Figure 2i

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
    # Gene is italicized at the facet level instead (make_func_plot below,
    # via facet_grid2()'s `strip` argument) -- not here. nature_theme is a
    # *complete* theme (built from theme_linedraw()), and `plot +
    # <complete theme>` replaces the plot's whole existing theme rather than
    # merging into it; a theme(strip.text.x = ...) tweak set on the plot
    # before nature_theme is applied would simply vanish (confirmed with a
    # minimal repro). A facet's own `strip` lives outside the theme object,
    # so it survives.
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
#
# Missense is hidden by default -- every gene's Functionally Abnormal
# (Missense) row is currently empty (no non-NA odds ratio) -- flip
# INCLUDE_MISSENSE_BREAKDOWN to TRUE once that changes.
INCLUDE_MISSENSE_BREAKDOWN <- FALSE
ABNORMAL_BREAKDOWN_CONSEQUENCES <- c(
  if (INCLUDE_MISSENSE_BREAKDOWN) 'Missense',
  'Truncating'
)

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
  filter(Gene != 'RAD51D') %>% # RAD51D excluded -- remove this filter to re-include
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
  # facet_nested() (ggh4x) with a single facet var behaves like plain
  # facet_grid() -- used here instead so the italics (`strip` below) live in
  # the same strip_nested()-based mechanism as Figure_6b.R and panel a of
  # Extended_Data_Figure_4.R, both confirmed to render correctly under this
  # figure's own cairo_pdf ggsave device below. facet_grid2()'s
  # strip_themed() (ggh4x's other themed-strip mechanism) was tried first
  # and worked under the base pdf() device (Extended_Data_Figure_4.R panel
  # b) but not under cairo_pdf here -- switched to the mechanism with a
  # confirmed cairo_pdf track record rather than chase that further.
  facet_nested(
    cols = vars(Gene),
    scales = 'free_x',
    strip = strip_nested(text_x = list(element_text(face = 'bold.italic')))
  ) +
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


# Show plot (interactive only -- under Rscript this would otherwise open
# the default graphics device and leave behind a stray Rplots.pdf)
if (interactive()) print(fig2i_plot + nature_theme)


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
