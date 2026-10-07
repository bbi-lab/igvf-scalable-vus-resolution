# Shared loader for the biobank-cohort odds-ratio table used by Figure 2i,
# Figure 6b, and Extended Data Figure 2 (func_df/combined_points_plot_df/
# predictor_plot_df/assay_plot_df). Pending a real IGVFFI3804AVJR.csv.gz
# re-export from Zenodo, these three scripts instead read an AoU-only export
# (data/input/biobank/AoU-OR-estimates_*.tsv.gz, optionally concatenated
# with a supplemental per-dataset file like excalibr-TP53-OR-estimates_*.tsv.gz
# -- see docs/data.md) that uses a different column layout and vocabulary.
# This function translates it back to the original schema those three
# scripts were written against (Dataset, Gene, Classifier, Classification,
# LogOR, LogOR_LI, LogOR_UI, `Cases with variants`), so none of them need to
# change beyond the read_csv() -> load_biobank_or_table() call itself. If a
# real IGVFFI3804AVJR.csv.gz export becomes available again, point the
# scripts back at read_csv() directly and retire this file.
#
# Expects the caller to have already loaded: dplyr, readr.

# Renames `x` wherever it matches a name in `map`, leaving everything else
# untouched -- unlike dplyr::recode(), this doesn't require every input
# value to be enumerated up front.
rename_vocab <- function(x, map) {
  matched <- match(x, names(map))
  dplyr::if_else(is.na(matched), x, unname(map[matched]))
}

# Parses a `Carrier cases`/`Carrier controls` string into a count, treating
# a privacy-censored "<= 20" as exactly 20 (the censoring threshold, and the
# most defensible single-number estimate without the real value). Only for
# descriptive display (e.g. sample-size annotations) -- never for filtering
# or statistics, where that approximation would be inappropriate.
parse_censored_count <- function(x) {
  as.integer(stringr::str_remove(x, '^.*\\s'))
}

# By default this restricts to the unbroken-down "All" row per Dataset/
# Gene/Classifier/Classification, reproducing the one row per bin the
# original file had. Pass keep_consequence = TRUE (Figure_2i.R's
# truncating/non-truncating split) to instead keep every per-consequence
# row (Missense/Stop gained/Splicing/Deletion/their combination/All) with
# its own `Consequence` column, and do that filtering yourself.
load_biobank_or_table <- function(path, keep_consequence = FALSE) {
  grouping_to_classifier <- c(
    'Functional class' = 'StandardizedClass'
    # Every other `Grouping/Calibration` value (ExCALIBR, Gene-specific,
    # Genome-wide aggregation, Domain aggregation, ExCALIBR + REVEL
    # single-gene) already matches the original `Classifier` vocabulary.
  )
  range_to_classification <- c(
    'Functionally abnormal' = 'ABNORMAL',
    'Functionally normal' = 'NORMAL'
    # Every other `Group/Range` value (the ExCALIBR/points/predictor-score
    # bins, e.g. "≤ -1", "≥ +1") already matches the original
    # `Classification` vocabulary.
  )
  dataset_rename <- c(
    'Combined points' = 'Total points'
  )

  df <- readr::read_tsv(path, show_col_types = FALSE)
  if (!keep_consequence) {
    df <- df %>% dplyr::filter(Consequence == 'All')
  }
  df %>%
    dplyr::transmute(
      Dataset = rename_vocab(`Dataset/Predictor`, dataset_rename),
      Gene = Gene,
      Classifier = rename_vocab(`Grouping/Calibration`, grouping_to_classifier),
      Classification = rename_vocab(`Group/Range`, range_to_classification),
      Consequence = if (keep_consequence) Consequence else NULL,
      LogOR = LogOR,
      LogOR_LI = LogOR_LI,
      LogOR_UI = LogOR_UI,
      # `Carrier cases` is privacy-censored to the string "<= 20" whenever
      # the true count is 1-20 (real zeros are left as "0", unambiguous and
      # not privacy-sensitive). Every script that reads this table only
      # ever uses `Cases with variants` for a `> 0` presence filter, never
      # the actual count -- so collapsing it to a 0/1 indicator on whether
      # the raw string is literally "0" preserves that filter exactly
      # without fabricating a fake number for the censored rows.
      `Cases with variants` = dplyr::if_else(`Carrier cases` == '0', 0L, 1L),
      # Parsed approximate counts (censored -> 20) for descriptive sample-
      # size annotations only -- see parse_censored_count() above.
      `Carrier cases` = parse_censored_count(`Carrier cases`),
      `Carrier controls` = parse_censored_count(`Carrier controls`)
    )
}
