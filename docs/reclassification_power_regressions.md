# Reclassification power regressions since the original submission

Gene-by-gene comparison of reclassification power (rank-based plot,
Extended Data Figure 5b, `% VUS remaining`) between the original submission
and the current pipeline, plus root-cause findings for the genes that moved
significantly. One real code bug found and fixed (RAD51C). Two genes
(KCNE1, PTEN) traced to a loss of exCALIBR calibration strength that isn't
a code issue and isn't currently understood. One gene (TPK1) traced to a
different upstream MaveDB score set. Two improvements (DDX3X, CHEK2)
already explained by prior fixes. Several genes still open.

## Genes with a significant shift

`% VUS remaining` is the rank-based plot's metric: of a gene's assayed
missense VUS, the share whose combined evidence (MAVE functional points +
REVEL/AM/MP2 predictor points) still lands in the ACMG-style Uncertain
band (0-5 total points) instead of resolving to a Benign or Pathogenic
bucket. Lower is better -- more of the gene's VUS got reclassified.
`Overall` is the mean of REVEL/AM/MP2 for that gene. "Original" and
"current" are both computed the same way, directly from
`Supplementary_Data_5.preprint.xlsx` (as-submitted) and the current
`Supplementary_Data_5.xlsx`, so the columns are directly comparable.

Ranks are out of the same 39 genes in both columns (1 = best reclassification
power, i.e. lowest `% VUS remaining`); `LDLR` is excluded from ranking since
it isn't in the original submission.

| Gene | Orig. rank | Cur. rank | Δ rank | Orig. Overall | Cur. Overall | Δ Overall | Orig. REVEL | Cur. REVEL | Δ REVEL | Direction | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| KCNE1 | 7th | 17th | +10 | 10.6% | 28.0% | +17.4 pp | 11.1% | 29.5% | +18.4 pp | Worse | Root cause found -- loss of exCALIBR calibration strength, not reverting (see below) |
| SCN5A | 36th | 39th | +3 | 83.3% | 100.0% | +16.7 pp | 83.3% | 100.0% | +16.7 pp | Worse | Not yet investigated |
| PTEN | 17th | 24th | +7 | 30.5% | 44.8% | +14.3 pp | 30.3% | 46.2% | +15.9 pp | Worse | Root cause found -- lost non-2018 ExCALIBR calibration plus weaker 2018-vintage calibration (see below) |
| NDUFAF6 | 21st | 25th | +4 | 33.5% | 46.1% | +12.6 pp | 33.9% | 47.3% | +13.4 pp | Worse | Not yet investigated |
| HMBS | 18th | 23rd | +5 | 31.2% | 40.7% | +9.5 pp | 36.0% | 49.6% | +13.6 pp | Worse | Not explained by controls; not yet investigated further |
| TPK1 | 27th | 32nd | +5 | 59.4% | 65.8% | +6.4 pp | 62.4% | 74.7% | +12.3 pp | Worse | Root cause found -- different upstream MaveDB score set (see below) |
| BRCA2 | 12th | 4th | -8 | 15.7% | 6.1% | -9.6 pp | 18.7% | 5.9% | -12.8 pp | Better | Noted, not chased -- smaller magnitude |
| CTCF | 32nd | 28th | -4 | 68.5% | 58.9% | -9.6 pp | 48.2% | 40.0% | -8.2 pp | Better | Not yet investigated |
| TSC2 | 24th | 19th | -5 | 43.2% | 29.2% | -14.0 pp | 47.0% | 32.9% | -14.1 pp | Better | Not yet investigated |
| CHEK2 | 19th | 9th | -10 | 32.0% | 14.1% | -17.9 pp | 30.9% | 13.0% | -17.9 pp | Better | Already explained -- CHEK2 QC-flag transcript-prefix fix, no action needed |
| DDX3X | 31st | 18th | -13 | 66.9% | 28.5% | -38.4 pp | 45.3% | 11.1% | -34.2 pp | Better | Already explained -- multi-range ExCALIBR fix (`f1f8514`), no action needed |


| RAD51C | 9th | 11th | +2 | 12.3% | 16.3% | +4.0 pp | 13.6% | 19.5% | +5.9 pp | Worse (fixed) | Fixed this session (Unicode bug); numbers above are post-fix -- pre-fix current was 43.2% REVEL |

Rank and percentage mostly agree on magnitude, but not always -- KCNE1 moves
10 rank positions on a 17-point swing while TPK1 moves 5 positions on a
smaller 6-point swing, because the gene list is much more tightly packed in
TPK1's part of the distribution (60-70% VUS remaining) than in KCNE1's
(10-30%). Rank is the more eyeball-friendly number; percentage is the one
that's actually additive/comparable across genes.

Everything not listed above moved less than ~9 percentage points either
direction, in line with noise.

## Per-gene status

**RAD51C -- fixed.** `load_excalibr_calibrations.py` derived each dataset's
name from its calibration JSON filename's stem. macOS decomposes
non-ASCII filenames to NFD on extraction, so `RAD51C_Olvera-León`'s
accented o came back as `o` + a combining acute accent, which doesn't
match the NFC-typed name in `score_sets.tsv`. The merge against `pp`
failed silently -- no error, just `NaN` -- so `Fxn_points` fell through to
the zero default for all 808 of RAD51C's VUS variants. Fixed at the source
(normalize to NFC when reading the filename) and defensively in the
notebook's merge cell. Restores real functional evidence to 780 of 808
variants; VUS-remaining drops from 43.2/62.1/62.5% to 19.4/16.2/13.2%
(REVEL/AM/MP2). Changes staged, not committed yet.

**KCNE1 -- root cause found, not reverting.** All three KCNE1 datasets
calibrated worse on the current exCALIBR run than at original submission:
`trafficking` lost all 16 calibration points outright (full calibration →
empty), `potassium_flux` lost its two weakest positive-direction points and
switched control-pooling method (`benign_method` avg → benign), and
`trafficking_WT_background_DN` lost its entire negative-direction range and
dropped from 3 control categories to 2. Restoring the original calibration
values (from the original Supplementary Data 4) drops VUS-remaining from
28.4/25.9/22.2% to 11.1/11.1/8.6% on the 81 shared variants, plus surfaces
14 more variants that should be excluded outright for conflicting evidence
across assays but currently aren't, because `trafficking`'s always-zero
value can't trigger the opposite-sign conflict check. Dedup strategy isn't
a factor -- only 1 of 95 duplicated variants picks a different dataset
between `v1` and `abs_max`, and even that one nets identical output.
ClinVar controls aren't the cause either -- KCNE1 gained controls (35 →
53), didn't lose any. This is a real loss of calibration strength on the
exCALIBR side -- not a code regression, and the cause isn't known
(exCALIBR is external to this pipeline). Decision is to keep the current
calibration rather than restore the original. No code change.

**TPK1 -- root cause found.** Old and current draw `TPK1_Weile_2017` from
two different MaveDB score sets. Old uses `urn:mavedb:00000001-d-2`
("TPK1 imputed and refined") -- 4,860 variants, amino-acid resolution,
explicitly built by "Random-Forest-based machine learning ... to impute
missing values and refine low-confidence measurements." Current
`score_sets.tsv` points only at `urn:mavedb:00001251-a-1`
("TPK1 DMS-TileSeq") -- 9,542 variants, nucleotide resolution, the raw
Weile et al. 2017 screen with no imputation. That single switch explains
both the resolution difference (old's 4,584 condensed rows are all
`nucleotide_or_aa == 'aa'`; current's 7,639 are all `'nt'`) and a clean
147-nt / 49-residue void in current data (c.355-c.501, aa 119-167) that old
appears to cover: old's values there are ML-imputed, not measured, and the
raw screen genuinely has no data for that stretch. Confirmed this isn't a
position-mapping artifact -- checked whether old's rows in that range are
the same physical variants as some shifted position in current data, and
genomic coordinates rule it out (43 kb apart, different exons). ClinVar
controls dropped 58 → 43 as a side effect of the score-set switch. Open
question: why ExCALIBR's calibration for `TPK1_Weile_2017` comes back
completely empty -- 7,639 real variants and 43 controls is enough to
attempt a calibration, so the empty result isn't obviously explained by the
score-set switch alone.

The project's own curation sheet (the v2 sheet built by the curation team)
independently treats `TPK1_Weile_2017` and `TPK1_Weile_2017_imputed` as two
separate datasets, and excludes the imputed one from the pillar project.
That matches the score-set finding above and confirms the raw/imputed
split is a known, intentional distinction on the curation side -- not
something this investigation is surfacing for the first time.

**PTEN -- root cause found.** Two compounding calibration losses, found by
comparing the original `Supplementary_Data_4.v1.xlsx` to the current
calibration sheet. First, PTEN's non-2018-vintage ExCALIBR calibration is
gone entirely: the original submission had four calibration rows
(`PTEN_Matreyek_2018`, `PTEN_Matreyek_2018_clinvar_2018`,
`PTEN_Mighell_2018`, `PTEN_Mighell_2018_clinvar_2018`); only the two
`_clinvar_2018` rows exist now. PTEN is one of three genes
(`priority_genes_2018 = ['BRCA1', 'PTEN', 'MSH2']`) that prefers the
2018-vintage calibration and falls back to the 2025 one when 2018 is
missing -- that fallback no longer exists. Second, the remaining
`_clinvar_2018` calibrations themselves weakened, with the same signature
as KCNE1: `n_c` went from a bare `3c` to `3c_benign`
(`PTEN_Matreyek_2018_clinvar_2018`) and `3c_avg`
(`PTEN_Mighell_2018_clinvar_2018`). Two unrelated genes losing calibration
strength through the identical `n_c` pattern points to a broad exCALIBR-side
event rather than something specific to either gene. ClinVar controls
aren't a factor -- PTEN gained controls (378 → 489), didn't lose any.
Dedup strategy isn't a factor either -- 48 of 628 variants pick a different
dataset between `v1` and `abs_max`, but none of those 48 change class.
Restoring the original calibration on the 597 shared variants drops
VUS-remaining from 44.7/41.2/44.2% to 31.3/30.5/33.7% (REVEL/AM/MP2); 80
REVEL variants flip class, 72 `Uncertain → Likely Benign` and 8
`Uncertain → Likely Pathogenic`. No code change -- like KCNE1, this is on
the exCALIBR side and the cause isn't known.

**HMBS -- open.** Real, significant drop in reclassification power vs. the
original submission. ClinVar controls aren't the cause -- control set is
identical (119 → 119, zero lost, zero gained). No working hypothesis yet.

**DDX3X, CHEK2 -- already explained, no new work needed.** Both
improvements trace to fixes already in the commit history: DDX3X to
`f1f8514` (ExCALIBR calibrations with disjoint multi-range point values --
DDX3X was the motivating case), CHEK2 to the transcript-prefix fix for
QC-flag matching.

**BRCA2 -- noted, not chased.** Smaller shift than the others (Overall
15.7% → 6.1%). No working hypothesis.

**SCN5A, NDUFAF6, CTCF, TSC2 -- not yet looked at.** SCN5A's one dataset
(`SCN5A_Glazer_2020`) is already known to have an empty ExCALIBR
calibration (see below), which plausibly explains its jump to 100% VUS
remaining -- not yet confirmed the way TPK1 was. NDUFAF6, CTCF, and TSC2
have no finding attached at all.

## Systemic finding: datasets with zero functional evidence

Audited all 91 MAVE datasets for 100% zero `Fxn_points`. Twelve came back
positive, in three distinct categories:

1. **Unicode join bug (fixed):** `RAD51C_Olvera-León_2024` only. No other
   dataset shares the failure signature -- the other affected gene names
   (`Grønbæk-Thygesen`, `Bolognesi_Faure`) use characters that don't have
   the same NFC/NFD ambiguity.
2. **Deliberate exclusion, not a bug:** `SFPQ_IGVF`, documented in
   `docs/variant_classification.md` as dropped for insufficient controls.
3. **ExCALIBR itself returned empty `point_ranges`** (join is fine, the
   calibration content is just empty, and `uncalibratable_reason` wasn't
   populated to say why): `ASPA_Grønbæk-Thygesen_2024_abundance`,
   `CARD11_Meitlis_2020_SGE_Ibrutinib_GoF`,
   `CARD11_Meitlis_2020_SGE_LoF`, `CBS_Sun_2020_high_B6`,
   `CRX_Shepherdson_2024`, `KCNE1_Muhammad_2024_trafficking`,
   `KCNQ4_Zheng_2022_v12_homozygous`, `SCN5A_Glazer_2020`,
   `TARDBP_Bolognesi_Faure_2019`, `TPK1_Weile_2017`.

Genes left with **no working ExCALIBR dataset at all** (single affected
dataset, no fallback): CARD11, CRX, TARDBP, TPK1 (TPK1's coverage gap is
explained above; the empty calibration itself is still open, and the other
three haven't been looked at). SFPQ is the same situation but expected.
RAD51C was in this group before the fix.

Genes with a **working fallback dataset** (so less severe): ASPA, CBS,
KCNE1, KCNQ4, SCN5A.

## Open items

- CARD11, CRX, TARDBP: same "empty ExCALIBR calibration" symptom as TPK1.
  Haven't checked whether they share TPK1's cause (wrong/imputed score set)
  or something else.
- TPK1: coverage gap explained (score-set switch); why ExCALIBR's
  calibration is empty despite adequate data is still open.
- PTEN: root cause found (lost calibration, see above); like KCNE1, cause
  on the exCALIBR side is unknown and no decision yet on whether to
  restore the original calibration.
- HMBS: significant power loss, ClinVar controls ruled out, no working
  hypothesis yet.
- SCN5A: likely explained by its empty ExCALIBR calibration (same symptom
  as CARD11/CRX/TARDBP/TPK1) -- not yet confirmed the way TPK1 was.
- NDUFAF6, CTCF, TSC2: real double-digit swings, no investigation started.
- RAD51C fix: staged, needs a commit decision.
- KCNE1 and PTEN both show the exact same `n_c` calibration-weakening
  signature -- worth checking whether other genes share it too, since it
  may be a single broad exCALIBR-side event rather than two coincidences.
