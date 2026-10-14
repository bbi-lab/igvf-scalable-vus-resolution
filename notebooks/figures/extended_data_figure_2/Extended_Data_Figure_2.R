# Extended Data Figure 2

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

# Some magic numbers:
# y-placements of the horizontal lines separating the "pathogenic" and "benign" zones,
VEP_PATH_BENIGN_BOUNDARY <- 4.5
ASSAY_PATH_BENIGN_BOUNDARY <- 8.5

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
    # Half the default vertical gaps around a bottom-positioned legend:
    # legend.box.spacing (panel-to-legend, default 11pt) and legend.margin's
    # top/bottom (default 5.5pt each, left/right left unchanged).
    legend.box.spacing = unit(5.5, 'pt'),
    legend.margin = margin(2.75, 5.5, 2.75, 5.5),
    ggh4x.facet.nestline = element_line(linewidth = 1*point_in_mm, color="black"),
    plot.tag = element_text(face = 'bold')
  )

# Load main table -- keep_consequence = TRUE so panel a (condensed_assay_plot)
# can overlay a Missense-only series on the aggregate ("All") one; every
# other use below restricts back to Consequence == 'All' to stay unaffected.
or_source_df <- load_biobank_or_table(
  "../../../data/input/biobank/AoU-OR-estimates_2026-10-07.tsv.gz",
  keep_consequence = TRUE
)

# TSC2_IGVF has a non-NA odds ratio only in the older 2026-10-06_merged
# export; under the current 2026-10-07 export it is NA, like every other
# gene's <= 20-carrier ExCALIBR bins (see SHOW_LOW_CARRIER_STYLE below), so
# flipping this on currently adds nothing but blank data. Off by default;
# pairs with SHOW_LOW_CARRIER_STYLE = TRUE. Affects gene_groups_df,
# condensed_assay_datasets, and figure_predictor_plot_df's gene filter below.
INCLUDE_TSC2 <- FALSE

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
  # "GCK", "Metabolic", "Meta- bolic", "Meta- bolic", # excluded -- no longer has any non-NA odds ratio (see condensed_assay_datasets)
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
) %>%
  filter(INCLUDE_TSC2 | Gene != "TSC2")

# Transformations that are common to both the predictor plot (panel b) and the assay plot (panel a)
or_df <- or_source_df %>%
  mutate(
    `Odds Ratio` = exp(LogOR),
    OR_LI = exp(LogOR_LI),
    OR_UI = exp(LogOR_UI),
    significance = case_when(
      `p-value (Holm adjusted)` < 0.01 ~ 'Significant at 99% (Holm adjusted)',
      `p-value (Holm adjusted)` < 0.05 ~ 'Significant at 95% (Holm adjusted)',
      TRUE ~ 'Not significant at 95%'
    )
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
    Consequence == 'Missense',
    Classifier %in% c(
      'Gene-specific',
      'Genome-wide aggregation'
    ),
    Classification %in% vep_levels
  ) %>%
  mutate(
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
  # facet_grid2() (ggh4x), not plain facet_grid() -- its `strip` argument
  # below italicizes the Gene (column) strips without touching Dataset
  # (row), and -- unlike a theme(strip.text.x = ...) call on this plot --
  # survives `figure_predictor_plot + nature_theme` further down. nature_theme
  # starts from theme_linedraw(), a *complete* theme; ggplot2 has a plot
  # (not theme) add a complete theme by replacing its whole existing theme
  # rather than merging into it, which would silently drop a
  # theme(strip.text.x = ...) set here (confirmed with a minimal repro --
  # the strip vanished). A facet's own `strip` lives outside the theme
  # object entirely, so it isn't subject to that replacement.
  facet_grid2(
    cols = vars(Gene),
    rows = vars(Dataset),
    scales = scales,
    strip = strip_themed(text_x = element_text(face = 'bold.italic'))
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
      'Significant at 95% (Holm adjusted)' = 16,
      'Significant at 99% (Holm adjusted)' = 18),
    breaks = c('Significant at 95% (Holm adjusted)', 'Significant at 99% (Holm adjusted)'),
    guide = guide_legend(position = legend_position)
  ) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  geom_hline(yintercept = VEP_PATH_BENIGN_BOUNDARY, linetype = 'dotted') +
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
    # Have gene-specific calibration; TSC2 gated by INCLUDE_TSC2 above.
    Gene %in% c('BRCA1', 'BRCA2', 'MSH2', 'TP53', if (INCLUDE_TSC2) 'TSC2'),
  )

# Make plot for final figure
figure_predictor_plot <- figure_predictor_plot_df %>%
  make_predictor_plot(legend_position = 'bottom') +
  geom_blank(aes(x=deframe(
    predictor_plot_df %>%
      filter(Gene == 'MSH2') %>%
      summarize(max(OR_UI))
  ))) +
  # Pin every per-gene/per-dataset panel's x-axis (scales = 'free_x' above)
  # to span at least 0.9-3.0, same broadcast-constant trick as the MSH2
  # geom_blank just above: a scalar aes(x=...) with no `data` argument
  # inherits figure_predictor_plot_df, so it's recycled across every row and
  # lands one invisible point per existing Gene/Dataset facet. A panel whose
  # real OR_LI/OR_UI already exceeds this range is unaffected; one narrower
  # than it gets extended out to exactly 0.9 and/or 3.0.
  geom_blank(aes(x = 0.9)) +
  geom_blank(aes(x = 3.0))


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
  'CHEK2_Gebbia_2024',
  #'GCK_Gersing_2023_complementation', # excluded -- no longer has any non-NA odds ratio
  'KCNE1_Muhammad_2024_potassium_flux',
  'KCNH2_Jiang_2022',
  # 'KCNH2_O_Neill_2024_surface_expression',
  'KCNQ4_Zheng_2022_current_homozygous',
  # LDLR_Tabet_2025_presence_VLDL excluded -- removed from the dataset
  # entirely (see docs/variant_classification.md). LDLR_Tabet_2025_uptake
  # excluded here too -- this version shows abundance only.
  'LDLR_Tabet_2025_abundance',
  'MSH2_Jia_2021',
  'PALB2_IGVF',
  'PTEN_Mighell_2018',
  'RAD51C_Olvera-León_2024',
  'RAD51D_IGVF',
  'SCN5A_Ma_2024',
  'TP53_Fayer_2021_meta',
  #'TP53_Boettcher_2019'
  #'TP53_Fortuno_2021',
  #'TP53_Giacomelli_2018_combined_score'
  if (INCLUDE_TSC2) 'TSC2_IGVF'
)

# Almost every dataset above reports *both* an ExCALIBR score-interval
# classification and a separate, coarser StandardizedClass (Functional
# class: NORMAL/ABNORMAL) one -- e.g. BAP1_Waters_2024 has both. Panel a
# should keep using each dataset's ExCALIBR points as before; StandardizedClass
# is only the right (and only available) choice for a dataset like
# TP53_Fayer_2021_meta that has no ExCALIBR calibration at all. Naively
# filtering Classifier %in% c('ExCALIBR', 'StandardizedClass') pulled in
# *both* for every other gene too -- double-counting carriers in the
# sample-size annotations and adding spurious extra points alongside the
# real ExCALIBR ones. This resolves the choice per dataset instead.
condensed_assay_datasets_without_excalibr <- setdiff(
  condensed_assay_datasets,
  or_df %>%
    filter(Classifier == 'ExCALIBR', Dataset %in% condensed_assay_datasets) %>%
    distinct(Dataset) %>%
    pull(Dataset)
)

# Sample sizes shown above each gene panel: whether that's each dataset's
# total cohort ("Total cases"/"Total controls", constant across every row of
# a Dataset -- not broken down by assay result at all) or just the subset of
# that cohort who carry a classified variant from this gene's assay
# ("Carrier cases"/"Carrier controls") is controlled by SHOW_TOTAL_COUNTS
# below. Total is simpler and immune to the privacy-censoring issues carrier
# counts have (see the carrier-mode comment further down) -- it's the
# default -- but carrier counts remain available by flipping the flag.
SHOW_TOTAL_COUNTS <- TRUE

# Whether to draw the "N cases"/"N controls" sample-size strip above each
# gene panel at all (see assay_count_annotations_df and the geom_text()
# layer further down). Off by default; set to TRUE to bring it back.
SHOW_SAMPLE_SIZE_ANNOTATIONS <- FALSE

# Whether a point's "Carrier cases"/"Carrier controls" is privacy-censored
# (AoU reports only "<= 20" rather than an exact count whenever the true
# value is 1-20) is a direct signal of an imprecise OR -- these are the bins
# that otherwise show as insignificant/unstable and (per ExCALIBR's own
# figure, which this mirrors) can be visually set apart rather than dropped.
# Off by default; set to TRUE to draw those points with a thinner, dashed
# CI (point/marker size is unchanged -- see low_carrier below and
# geom_errorbar()/geom_pointrange() further down). The current data source
# (2026-10-07, see or_df above) blanks LogOR/LogOR_LI/LogOR_UI for every
# low-carrier row, so turning this on currently has nothing to draw --
# only the older 2026-10-06_merged export still carried real values there.
SHOW_LOW_CARRIER_STYLE <- FALSE

# Plain comma-grouped integer up to 9,999 (e.g. "1,410"); above that, one
# decimal place of thousands with a "K" suffix (e.g. "12.9K") up to
# 99,999, and no decimal place at all from 100,000 up (e.g. "120K") -- a
# tenths digit stops being meaningful noise at that point. Unlike
# scales::label_number(scale_cut = ...), which switches to "K" at 1,000 and
# varies its decimal count by magnitude (e.g. "1.410K", "284.49K").
count_fmt <- function(x) {
  dplyr::case_when(
    x >= 100000 ~ paste0(scales::label_number(accuracy = 1)(x / 1000), 'K'),
    x > 9999 ~ paste0(scales::label_number(accuracy = 0.1)(x / 1000), 'K'),
    TRUE ~ scales::label_comma(accuracy = 1)(x)
  )
}

assay_dataset_counts <- if (SHOW_TOTAL_COUNTS) {
  # "Total cases"/"Total controls" (see load_biobank_or_table.R) are each
  # dataset's whole cohort, the same for every row of that Dataset -- so
  # this is just one row per dataset, no binning/summing/censoring to deal
  # with. Columns are still named `Carrier ... upper/known/n_censored` to
  # match the carrier-mode branch below, so format_censored_count() and
  # assay_count_annotations_df don't need to know which mode is active:
  # n_censored = 0 always takes format_censored_count()'s bare-number path.
  or_df %>%
    filter(Dataset %in% condensed_assay_datasets) %>%
    distinct(Dataset, Gene, `Total cases`, `Total controls`) %>%
    transmute(
      Dataset, Gene,
      `Carrier cases upper` = `Total cases`,
      `Carrier controls upper` = `Total controls`,
      `Carrier cases known` = `Total cases`,
      `Carrier controls known` = `Total controls`,
      `Carrier cases n_censored` = 0L,
      `Carrier controls n_censored` = 0L
    )
} else {
  # Total carriers among cases/controls with any classified (non-"0")
  # variant from that gene's assay. Not broken down by consequence (e.g.
  # missense vs. other) -- "All" and "Missense" aren't disjoint (a carrier
  # can have variants of both kinds in the same assay), so an
  # All-minus-Missense subtraction risks double-counting/undercounting
  # people rather than giving an exact non-missense count. The ExCALIBR
  # score bins are cumulative, so there's no single row that already holds
  # "every carrier of this consequence" -- "≤ -1"/"≥ +1" are the
  # least-restrictive (broadest) bin in each direction, so summing just
  # those two covers everyone with a non-zero score.
  #
  # Carrier counts are parsed from the source data with privacy-censored
  # "≤ 20" cells treated as 20 (see parse_censored_count() in
  # load_biobank_or_table.R) -- an upper bound on that cell, not its real
  # value. A sum that includes one or more censored cells is therefore only
  # known to fall somewhere between (sum of the exact cells + 1 per censored
  # cell, since a censored cell is never really 0) and (sum of the exact
  # cells + 20 per censored cell) -- shown as that range (e.g. "971-990") in
  # assay_count_annotations_df, unless the range degenerates to a single
  # censored cell against an otherwise-empty sum (no exact contribution at
  # all), in which case the true total is simply ≤ 20, same as the source
  # cell itself.
  or_df %>%
    filter(
      Dataset %in% condensed_assay_datasets,
      (Classifier == 'ExCALIBR' & !(Dataset %in% condensed_assay_datasets_without_excalibr)) |
        (Classifier == 'StandardizedClass' & Dataset %in% condensed_assay_datasets_without_excalibr),
      Consequence == 'All',
      # '≤ -1'/'≥ +1' are ExCALIBR's broadest bins (see comment above); for
      # StandardizedClass (TP53_Fayer_2021_meta's OddsPath calibration) every
      # carrier is NORMAL or ABNORMAL, so that pair is the equivalent "everyone
      # with a classification" set.
      Classification %in% c('≤ -1', '≥ +1', 'NORMAL', 'ABNORMAL')
    ) %>%
    group_by(Dataset, Gene) %>%
    summarise(
      # Upper bound: every censored cell counted at its substituted value (20,
      # see above) -- equal to the exact total when nothing was censored.
      `Carrier cases upper` = sum(`Carrier cases`),
      `Carrier controls upper` = sum(`Carrier controls`),
      # The portion of the sum known exactly, excluding any censored cell.
      `Carrier cases known` = sum(if_else(`Carrier cases censored`, 0L, `Carrier cases`)),
      `Carrier controls known` = sum(if_else(`Carrier controls censored`, 0L, `Carrier controls`)),
      `Carrier cases n_censored` = sum(`Carrier cases censored`),
      `Carrier controls n_censored` = sum(`Carrier controls censored`),
      .groups = 'drop'
    )
}

# A StandardizedClass (OddsPath) dataset only ever reports a binary
# NORMAL/ABNORMAL call, not a score -- but every variant in one of those
# classes gets the exact same point value, fixed by that dataset's own
# OddsPath calibration strength (e.g. TP53_Fayer_2021_meta's Evidence Code
# Abnormal is PS3_moderate, i.e. +2 for every ABNORMAL carrier). Thresholds
# mirror annotate_OP_points() in Variant_Classification_analysis.ipynb.
oddspath_points <- read_csv(
  '../../../data/output/mave_calibration/OddsPath_calibrations.csv.gz',
  show_col_types = FALSE
) %>%
  transmute(
    Dataset,
    OddsNormal = suppressWarnings(as.numeric(OddsNormal)),
    OddsAbnormal = suppressWarnings(as.numeric(OddsAbnormal))
  ) %>%
  mutate(
    normal_magnitude = case_when(
      is.na(OddsNormal) ~ NA_integer_,
      OddsNormal < 0.053 ~ 4L,
      OddsNormal < 0.23 ~ 2L,
      OddsNormal < 0.48 ~ 1L,
      TRUE ~ NA_integer_
    ),
    abnormal_magnitude = case_when(
      is.na(OddsAbnormal) ~ NA_integer_,
      OddsAbnormal > 350 ~ 8L,
      OddsAbnormal > 18.7 ~ 4L,
      OddsAbnormal > 4.3 ~ 2L,
      OddsAbnormal > 2.1 ~ 1L,
      TRUE ~ NA_integer_
    )
  )

condensed_assay_plot_df <-
  assay_plot_df %>%
  filter(
    # Most datasets here report both Classifier values -- ExCALIBR is kept
    # for those; StandardizedClass only takes over for a dataset like
    # TP53_Fayer_2021_meta that has no ExCALIBR calibration at all (see
    # condensed_assay_datasets_without_excalibr above). Without this
    # per-dataset split, every other gene would get both its ExCALIBR points
    # and its StandardizedClass NORMAL/ABNORMAL call plotted together.
    (Classifier == 'ExCALIBR' & !(Dataset %in% condensed_assay_datasets_without_excalibr)) |
      (Classifier == 'StandardizedClass' & Dataset %in% condensed_assay_datasets_without_excalibr),
    Dataset %in% condensed_assay_datasets,
    Classification != "0"
  ) %>%
  left_join(oddspath_points, by = 'Dataset') %>%
  mutate(
    magnitude = case_when(
      Classifier == 'StandardizedClass' & Classification == 'NORMAL' ~ normal_magnitude,
      Classifier == 'StandardizedClass' & Classification == 'ABNORMAL' ~ abnormal_magnitude,
      TRUE ~ NA_integer_
    )
  ) %>%
  # ExCALIBR's own point bins are cumulative -- '>= +1' includes every
  # carrier in '>= +2' (and anyone scoring exactly +1). A StandardizedClass
  # dataset's ABNORMAL call is a single fixed point value (e.g. +2 above),
  # which means nobody from that assay ever scores exactly +1 -- so its
  # '>= +1' row covers exactly the same carriers as its '>= +2' row, and the
  # two should be identical rather than '>= +1' being left blank. Expand
  # each StandardizedClass row into one copy per cumulative bin from 1 up to
  # its own magnitude (ExCALIBR rows get a single no-op copy via the
  # coalesce to 1).
  rowwise() %>%
  mutate(bin = list(seq_len(coalesce(magnitude, 1L)))) %>%
  ungroup() %>%
  unnest(bin) %>%
  mutate(
    Classification = case_when(
      Classifier == 'StandardizedClass' & Classification == 'NORMAL' ~ paste0('≤ -', bin),
      Classifier == 'StandardizedClass' & Classification == 'ABNORMAL' ~ paste0('≥ +', bin),
      TRUE ~ as.character(Classification)
    ),
    Classification = factor(Classification, levels = assay_classification_levels)
  ) %>%
  # A StandardizedClass row whose calibration didn't clear even a Supporting
  # threshold (magnitude NA) has no bin to land on -- drop it rather than
  # plot an unlabeled point.
  filter(!is.na(Classification)) %>%
  select(-normal_magnitude, -abnormal_magnitude, -magnitude, -bin) %>%
  mutate(
    Consequence = factor(Consequence, levels = c('All', 'Missense')),
    # See SHOW_LOW_CARRIER_STYLE above -- `Carrier cases`/`Carrier controls`
    # are never an exact nonzero count under 21 (AoU always reports those as
    # "<= 20"), so the censored flag already is the "<= 20" test.
    low_carrier = `Carrier cases censored` | `Carrier controls censored`,
    # Dataset names (e.g. "RAD51C_Olvera-León_2024") are too long to fit one
    # panel-width strip at this font size and ran into neighboring panels
    # unwrapped. str_wrap() only breaks at whitespace, and these names have
    # none -- the underscore swap gives it break points (also reads better
    # than a bare underscored identifier as a heading); done here, after
    # every upstream Dataset %in%/join use above, so those keep matching the
    # real identifier. Must match the same transform applied to Dataset in
    # assay_count_annotations_df below, since facet routing for that layer
    # is by exact string match.
    # Width 20, not 18 -- at 18, KCNE1_Muhammad_2024_potassium_flux's name
    # wrapped to 3 lines (every other dataset wraps to at most 2), and
    # facet_nested_wrap shares one strip-row height across every panel row,
    # not just this one's -- that single 3-line strip was eating into the
    # panel-height budget for every row, squeezing the Classification axis's
    # 18 tick labels into overlapping text even in rows with short headers.
    Dataset = str_wrap(gsub('_', ' ', Dataset), width = 20)
  ) %>%
  left_join(gene_groups_df)

# One row per gene per annotation line: a left-aligned "N cases" and a
# right-aligned "N controls", sharing a single row above each panel (no
# longer split by consequence -- see assay_dataset_counts above). This
# lands below the native per-gene strip (gene name, bold) and above the
# panel's own plotted data -- not above the strip/Disease header -- because
# that strip is given extra bottom margin below (reserving blank room for
# this content within its own cell) rather than this layer having to clear
# the whole strip height, as a strip and the panel directly below it are
# flush with no gap by default. `vjust` (more negative = further above the
# panel) only needs to be large enough to clear into that reserved margin
# via coord_cartesian(clip = 'off') below -- actual value tuned against
# this plot's real panel/strip size, not transferable to a plot with
# different dimensions.
# Formats a (possibly partially censored) sum: a bare number when nothing
# was censored; "≤ <upper>" when the only uncertainty is a single censored
# cell with no other (exact) contribution, so the true total is, like the
# source cell itself, simply at most 20; otherwise a "<lower>-<upper>"
# range, since a second nonzero contribution (exact or itself censored)
# means "≤ 20" would no longer be a tight -- or even correct -- bound (see
# assay_dataset_counts above). A range whose ends round to the same
# count_fmt() display (e.g. 24,792-24,811, both "24.8K") collapses to that
# one value instead of the redundant-looking "24.8K-24.8K".
format_censored_count <- function(known, n_censored, upper) {
  lower <- known + n_censored
  lower_fmt <- count_fmt(lower)
  upper_fmt <- count_fmt(upper)
  dplyr::case_when(
    n_censored == 0 ~ upper_fmt,
    n_censored == 1 & known == 0 ~ paste0('≤ ', upper_fmt),
    lower_fmt == upper_fmt ~ lower_fmt,
    TRUE ~ paste0(lower_fmt, '-', upper_fmt)
  )
}

assay_count_annotations_df <- bind_rows(
  assay_dataset_counts %>%
    transmute(
      Gene, Dataset,
      hjust = 0,
      text = sprintf(
        '%s cases',
        format_censored_count(`Carrier cases known`, `Carrier cases n_censored`, `Carrier cases upper`)
      )
    ),
  assay_dataset_counts %>%
    transmute(
      Gene, Dataset,
      hjust = 1,
      text = sprintf(
        '%s controls',
        format_censored_count(`Carrier controls known`, `Carrier controls n_censored`, `Carrier controls upper`)
      )
    )
) %>%
  # facet_nested_wrap facets on both Disease and Dataset -- without Disease
  # here too, ggplot can't route each row to its one matching panel and
  # broadcasts it into every panel instead (confirmed: this is exactly what
  # happened before this join was added). Dataset itself already comes from
  # assay_dataset_counts above, so only Disease (keyed by Gene) needs joining.
  left_join(gene_groups_df, by = 'Gene') %>%
  mutate(
    # condensed_assay_plot below plots pre-log10'd Odds Ratio/CI columns on
    # a plain continuous scale (see that scale's own comment for why) --
    # which, unlike scale_x_log10(), handles -Inf/Inf natively, so both
    # ends of this annotation can sit exactly at their panel's true border.
    x = if_else(hjust == 0, -Inf, Inf),
    # Must match the Dataset wrap applied in condensed_assay_plot_df above --
    # this layer is routed to its panel by exact Dataset string match.
    Dataset = str_wrap(gsub('_', ' ', Dataset), width = 20),
    # BAP1/MSH2/RAD51D's cases count is a range (e.g. "1,399-1,418"), and in
    # these narrower (5-across) panels that text is too wide for the
    # default 5pt without colliding with the controls text on the same
    # row -- shrink just these three enough to clear it.
    size = if_else(Gene %in% c('BAP1', 'MSH2', 'RAD51D'), 4.4, 5)
  )

# Limits for most panels
assay_plot_common_limits <- condensed_assay_plot_df %>%
  filter(!(Gene %in% c("BRCA1", "MSH2", "KCNH2"))) %>%
  summarise(
    OR_LI = min(OR_LI),
    OR_UI = max(OR_UI)
  ) %>% deframe()

# Build assay plot
# x/xmin/xmax are pre-log10'd here (and scale_x_continuous() below, not
# scale_x_log10(), does the display/guide work) rather than leaving the
# transform to the scale: scale_x_log10() maps -Inf to NaN rather than
# -Inf, which silently drops any layer placed there -- confirmed with a
# minimal repro -- which ruled out using x = -Inf to pin the left-aligned
# "cases" annotation (further down) exactly to each panel's true left
# border the same way x = Inf already does for "controls" on the right.
# Pre-transforming sidesteps that: a plain continuous scale passes -Inf
# through unchanged, same as Inf, so both ends of that annotation can now
# sit exactly at their panel's border with no approximation needed.
condensed_assay_plot <- ggplot(
  condensed_assay_plot_df,
  aes(
    x = log10(`Odds Ratio`),
    xmin = log10(OR_LI),
    xmax = log10(OR_UI),
    y=Classification,
    shape = significance,
    color = Consequence
  )
) +
  facet_nested_wrap(
    # Dataset, not Gene -- one panel per assay, so a gene with multiple
    # assays (e.g. LDLR's uptake/abundance) gets one panel each instead of
    # being collapsed/overlaid into a single gene-level panel.
    facets = vars(Disease, Dataset),
    nrow = 3,
    scales = 'free_x',
    # Dataset names (e.g. "RAD51C Olvera-León 2024") are much longer than
    # the gene symbols the strip used to show. A labeller-based wrap
    # (labeller(Dataset = label_wrap_gen(...))) left the text unwrapped and
    # overrunning into neighboring panels -- confirmed this is specific to
    # facet_nested_wrap/strip_nested, not a width problem, since the same
    # labeller wraps fine on a plain facet_wrap(). Wrapping is done directly
    # on the Dataset values instead (see condensed_assay_plot_df above),
    # which renders correctly here.
    # Default solo_line = FALSE only draws the nest line under a Disease
    # header spanning 2+ genes (e.g. "Cancer"), not a single-gene one (e.g.
    # "Cardiovascular" above just KCNH2) -- that asymmetry misaligned the
    # single-gene headers (no line reserving/centering their text the same
    # way). Drawing it for every Disease header, solo or not, equalizes them.
    solo_line = TRUE,
    # Disease keeps the default strip style; the per-assay strip gets extra
    # bottom margin, reserving blank room below the assay name (within its
    # own cell) for the Cases/Controls annotation below, so that content
    # lands below the assay name and above the panel -- not above the
    # Disease header, and not inside the panel either.
    strip = strip_nested(
      text_x = list(
        element_text(size = 7, face = 'bold'),
        # Plain bold, not bold.italic -- the Dataset/assay level isn't a
        # gene symbol, so it doesn't get the italic convention Gene used to.
        element_text(
          size = 7, face = 'bold',
          margin = margin(t = 3, r = 3, b = if (SHOW_SAMPLE_SIZE_ANNOTATIONS) 9 else 3, l = 3)
        )
      ),
      # strip_nested()'s by_layer_x defaults to FALSE, which does NOT treat
      # the text_x list above as one element per nesting depth (Disease,
      # then Dataset) -- without this, the margin meant for just the Dataset
      # level was applied inconsistently, stretching the Disease level's
      # own cell and leaving dead space below "Cancer" etc. above its own
      # separator line.
      by_layer_x = TRUE
    )
  ) +
  # low_carrier & SHOW_LOW_CARRIER_STYLE (rather than low_carrier alone)
  # collapses to a single FALSE value for every row when the flag is off, so
  # the scale_*_manual() calls below are no-ops and every point renders
  # exactly as before. Point/marker size is left unmapped -- every circle is
  # the same size regardless of low_carrier, only the CI is styled.
  # geom_pointrange() draws its own xmin-xmax connecting line on top of
  # geom_errorbar()'s (the point/range geoms overlap by design, for the
  # caps) -- linetype has to be mapped on both, or geom_pointrange's always-
  # solid line shows through underneath geom_errorbar()'s dashed one.
  geom_errorbar(
    aes(linetype = low_carrier & SHOW_LOW_CARRIER_STYLE),
    width = 0.5, position = position_dodge(width=0.5)
  ) +
  geom_pointrange(
    aes(linewidth = low_carrier & SHOW_LOW_CARRIER_STYLE, linetype = low_carrier & SHOW_LOW_CARRIER_STYLE),
    position = position_dodge(width=0.5), fill='white'
  ) +
  # Pin every per-Disease/Dataset panel's x-axis (facet_nested_wrap scales =
  # 'free_x' above) to span at least Odds Ratio 0.9-3.0, same
  # broadcast-constant geom_blank trick as figure_predictor_plot uses for its
  # MSH2 panel above: a scalar aes(x=...) with no `data` argument inherits
  # condensed_assay_plot_df, so it's recycled across every row and lands one
  # invisible point per existing Disease/Dataset facet. x is log10'd here (see
  # this plot's own x aes above) since scale_x_continuous(), not
  # scale_x_log10(), does the display work. A panel whose real OR_LI/OR_UI
  # already exceeds this range is unaffected; one narrower than it gets
  # extended out to exactly 0.9 and/or 3.0.
  geom_blank(aes(x = log10(0.9))) +
  geom_blank(aes(x = log10(3.0))) +
  # Sample-size annotation ("N cases" / "N controls"), positioned below
  # each panel's native gene-name strip -- in the room that strip's own
  # enlarged bottom margin reserves for it above -- rather than inside the
  # panel, via coord_cartesian(clip = 'off') below, which lets this layer's
  # content overflow past the panel's own border into that margin.
  (if (SHOW_SAMPLE_SIZE_ANNOTATIONS) geom_text(
    data = assay_count_annotations_df,
    # size is per-gene (see assay_count_annotations_df) -- I() bypasses the
    # default area-based size scale, treating these as literal point sizes.
    aes(x = x, y = Inf, label = text, hjust = hjust, size = I(size / .pt)),
    vjust = -0.9,
    color = 'black',
    # geom_text() doesn't inherit nature_theme's text family -- unlike
    # theme elements (axis/strip/legend text), a geom's font is set per
    # layer, not from the plot theme.
    family = 'Arial',
    inherit.aes = FALSE,
    # Without this, ggplot merges this layer into the Consequence color
    # legend below and draws geom_text's placeholder key glyph -- literally
    # the letter "a" -- on top of the All/Missense keys.
    show.legend = FALSE
  ) else NULL) +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(
    # x/xmin/xmax are already log10'd (see aes() above) -- breaks are
    # computed the same way scale_x_log10()'s own default would (nice
    # breaks in the untransformed space, via scales::breaks_log(), then
    # re-log10'd to place them on this now-linear axis), and
    # guide_axis_logticks(prescale.base = 10) draws the usual log-style
    # major/minor ticks for data that's pre-transformed rather than
    # transformed by the scale itself.
    # name: without this, the axis title falls back to deparsing the aes()
    # expression itself (literally "log10('Odds Ratio')") instead of the
    # plain column name a non-computed aes() mapping would have shown.
    name = 'Odds Ratio',
    breaks = function(lims) log10(scales::breaks_log()(10^lims)),
    labels = function(x) scales::label_number(drop0trailing = TRUE)(10^x),
    minor_breaks = NULL,
    guide = guide_axis_logticks(prescale.base = 10)
  ) +
  scale_shape_manual(
    values = c(
      'Not significant at 95%' = 21,
      'Significant at 95% (Holm adjusted)' = 16,
      'Significant at 99% (Holm adjusted)' = 18),
    breaks = c('Significant at 95% (Holm adjusted)', 'Significant at 99% (Holm adjusted)'),
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
  # Low-carrier styling (see SHOW_LOW_CARRIER_STYLE above): dashed, thinner
  # CI only -- point/marker size is intentionally untouched (and unmapped)
  # so every circle stays exactly the size it was before this option
  # existed. 'FALSE' = 1*point_in_mm matches nature_theme's own default
  # linewidth (confirmed via ggplot_build() against the unmodified script --
  # nature_theme's element_geom(pointsize = ...) does not end up governing
  # GeomPointrange's point size here, so don't reuse it as a reference
  # value). guide = 'none' on both -- this is a data-quality cue, not a new
  # legend.
  scale_linetype_manual(values = c('FALSE' = 'solid', 'TRUE' = '22'), guide = 'none') +
  scale_linewidth_manual(values = c('FALSE' = 1*point_in_mm, 'TRUE' = 0.6*point_in_mm), guide = 'none') +
  geom_vline(xintercept = log10(1), linetype = 'dashed') +
  geom_hline(yintercept = ASSAY_PATH_BENIGN_BOUNDARY, linetype = 'dotted')

# Compose figure together
fig_exd4 <- (condensed_assay_plot + nature_theme +
    # facet_nested_wrap() (ggh4x, used here) renders panel.border at roughly
    # double the linewidth facet_grid() (panel b, below) does for the same
    # theme value -- confirmed by rendering each panel alone and measuring
    # border pixel width at a fixed DPI. Halving it here (nature_theme's own
    # panel.border is unchanged, since panel b renders it correctly) matches
    # the two panels' border weight.
    theme(panel.border = element_rect(fill = NA, colour = 'black', linewidth = 0.25))) +
  (figure_predictor_plot + nature_theme + theme(plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))) +
  plot_annotation(tag_levels='a') +
  plot_layout(
    design = c(
      area(1,1,6,6),
      area(7,1,9,6)
    )
  )

# Show figure
print(fig_exd4)

# Save figure
ggsave(
  '../../../data/output/figures/assets/extended_data_figure_2/extended_data_figure_2.pdf',
  fig_exd4,
  # Bumped from 160 towards the 183 max -- panel a now has enough assay
  # panels (nrow = 3) to need up to 6 columns per row instead of the
  # original 5, and the extra width gives each narrower panel's 2-line
  # strip text room so neighboring panels' headings don't run into each
  # other.
  width = 183,
  # Bumped from the original 247 in two steps. First to 270: panel a's
  # assay-name strips now wrap onto 2 lines, which otherwise ate into the
  # fixed design = area(1,1,6,6) panel region and squeezed the
  # Classification axis's 18 levels into overlapping tick labels. Then to
  # 340, switching to cairo_pdf below: the base pdf() device was silently
  # substituting a narrow Hershey-vector approximation for the
  # Classification axis's "≥"/"≤" labels (confirmed by comparing
  # the two devices' output directly) rather than embedding the real
  # Arial glyph, so tick labels fit at 270mm; cairo_pdf's correctly-sized
  # (taller) glyphs need the extra room, confirmed with an isolated
  # facet_wrap() repro that reproduced the overlap at 270mm-equivalent
  # spacing and resolved it once each row had enough height.
  height = 340,
  units = 'mm',
  family = 'Arial',
  # cairo_pdf, not pdf -- the base pdf() device doesn't actually render
  # the Classification axis's "≥"/"≤" tick labels as those
  # Unicode glyphs at all; see the height comment above for how this was
  # confirmed.
  device = cairo_pdf,
  create.dir = TRUE)
