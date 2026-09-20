# Rose oil comparability analysis

Analysis code for a study of whether volatile-composition measurements of *Rosa*
essential oils and flowers, published by different laboratories, can be compared
with one another.

The corpus is a set of open-access primary studies whose composition tables were
extracted into a single long-format table of sample-by-compound measurements.
This repository contains the code that analyses that table. The literature
search, full-text retrieval, table extraction, compound-identifier resolution and
manual curation that produced it are a separate stage and are not included here;
the tables they produced are the input.

## What the analysis does

**Within-study contrasts.** Every study is searched for a factor whose levels were
measured side by side in the same laboratory on the same instrument. Eight factors
are enumerated — cultivar, species, developmental stage, organ, geographic origin,
post-harvest processing, experimental treatment, and chromatographic column class —
and levels are read both from the declared metadata fields and from sample codes in
the column headers. A contrast counts as genuine only when at least two levels are
present, every level shares at least three compounds that all levels measured, and
the platform and column polarity are constant across the levels; contrasts failing
the last condition are dropped rather than flagged. A separate field records when a
factor covaries with another factor, and another when its levels coincide with the
source table.

**Composition tests.** For each genuine contrast the compositions are closed to unit
sum, zero-replaced multiplicatively with δ set to 0.65 of the smallest positive value
in the column, and centred-log-ratio transformed; contrasts whose replaced mass would
exceed 30% of the total are refused. The factor is then tested by PERMANOVA on
Aitchison distances with Bray–Curtis as a comparator. Labels are permuted freely, and
the smallest p-value the design can produce under exact enumeration is reported
alongside the p-value itself, so that a floor is never mistaken for a null result.
Per-compound effects are reported for the twelve strongest drivers of each contrast,
in CLR units for two-level contrasts and as η² or a level range where there are more
levels, together with the difference on the original closed percentages.

**Variance partition and separability.** The corpus is cut into complete-case blocks —
chosen by three different criteria, so the blocks overlap and do not cover every
sample — and for each block the code asks whether a factor is separable from
laboratory identity at all. Five verdicts are distinguished: separable within study,
separable within study only, study-level only (nested within study, which the code
explicitly does not count as separability), structurally inseparable, and no
contrast. The verdict is recorded next to the tests rather than used to gate them:
marginal tests are run for every factor, including factors that are perfectly
confounded with study, and Cramér's V against study identity is reported so that
such a result can be recognised.

**Direction agreement.** Where two laboratories measured the same pair of levels on
the same axis, the signed per-compound effects are compared. Agreement is counted over
shared compounds, tested against a binomial null of one half, corrected within four
declared families, and recomputed at four effect-size thresholds so that agreement
driven by compounds that barely moved can be seen for what it is.

**Prediction baselines.** A leave-study-out evaluation asks how much of a held-out
study's composition can be predicted from the recorded metadata. Studies with fewer
than three samples are never held out. Every score is referenced to the mean of the
training folds rather than to the mean of the held-out study, because referencing to
the held-out mean would hide exactly the effect being measured. Group-mean baselines
are reported with the number of held-out samples whose level was actually present in
training, and the per-fold differences are tested with an exact two-sided sign test
with ties dropped. A linear probe recovers study identity, column polarity and
platform from the same feature space, its balanced accuracy reported against a
baseline of five label permutations.

**Within-study centring.** The probe is repeated on features from which each study's
own centroid has been subtracted. Centring requires no labels and can therefore be
applied to a held-out study without leakage, so the difference between the two probe
scores measures how much of the recoverable identity survives an operator that a
downstream user could actually apply.

**Reporting census.** Each study is scored against 32 reporting items covering
identification basis, retention indices, internal standards and response factors,
instrument and column configuration, extraction conditions, quantification mode,
dispersion, and the biological metadata needed to reuse the numbers. Items are graded
essential, recommended or desirable. The census is computed from the article text
where the item is a methods statement and from the curated tables where it is a
property of the data, and the re-parse is checked field by field against the
curated analytical table. Two costs are reported separately: how far the studies
that report an item differ from those that do not, and how many cross-study
comparisons are lost to silence at each rung of the ladder.

**Commensurability ladder.** Every pair of studies that shares a biological level is
walked up a cumulative sequence of six requirements, from "the two measured the same
kind of quantity" to "both named the same retention-index convention". A pair is
poolable when all requirements are met, incompatible when any is met with a stated
difference, and undecidable when any is simply not reported. Because the verdict
depends on how attributes are compared, the whole ladder is recomputed under all 32
combinations of five coding switches and over nine rung sets, and the full surface is
written out; the console headline quotes one named coding, and the spread across the
others is read from the table. An optimistic counterfactual then imputes the
unreported attributes from the empirical distribution among the studies that did
report them, giving an upper bound on what better reporting alone would buy.

**Published claims retested.** Five families of composition claim of the kind routinely made in this
literature — conformity to the ISO reference intervals, ranking two species by a
single marker, a characteristic profile for a growing origin, an effect of processing
route, and a correlation between two compounds — are tested on the pooled corpus and
then re-tested with laboratory held constant, to see which survive. Where no
within-laboratory contrast exists, the claim is recorded as untestable rather than
confirmed.

## Requirements

Python 3.9 or later, with:

```
numpy
scipy
scikit-learn
```

```
pip install -r requirements.txt
```

`09_reporting_census.py` additionally needs `pdfplumber`, but only when the article
text cache does not already cover a study's PDF; it is imported lazily and is not in
`requirements.txt`. Install it with `pip install pdfplumber` if you are rebuilding
the cache from PDFs.

## Data

### Obtaining the input data

The input tables are not distributed with this code. They are derived from the
composition tables of published studies, every one of which is open access.

`sources.tsv` lists all 60 studies with their DOI, PMC accession, year and licence.
`fetch_sources.py` reads that list and downloads the JATS full text and the
supplementary archives from Europe PMC into `$ROSE_DATA/fulltext/` and
`$ROSE_DATA/supplementary/`, which is the layout the scripts expect:

```
python3 fetch_sources.py
```

It needs no key and no account. Studies without a PMC record are printed with their
DOI so they can be retrieved from the publisher. The endpoints it uses are

```
https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:<doi>&format=json
https://www.ebi.ac.uk/europepmc/webservices/rest/<pmcid>/fullTextXML
https://www.ebi.ac.uk/europepmc/webservices/rest/<pmcid>/supplementaryFiles
```

and compound identifiers are resolved against PubChem at
`https://pubchem.ncbi.nlm.nih.gov/rest/pug/`.

Turning the downloaded articles into the tables below is the extraction and curation
stage, which is not part of this release.

The code reads its inputs from one directory and writes every result to another. Both
are set by environment variable and default to `data/` and `results/` inside the
repository:

```
export ROSE_DATA=/path/to/the/input/tables
export ROSE_RESULTS=/path/to/write/results
```

The scripts derive these paths from their own location, so they can be run from any
working directory. Importing the package creates the results directory.

The input directory is expected to contain:

| file | what it holds |
|---|---|
| `measurements_long.tsv` | one row per (sample, compound) cell, joined to its sample and study metadata; this is the table most of the analyses read |
| `samples.tsv` | one row per measurement column: species, cultivar, organ, stage, origin, material state, platform, column phase and polarity, and whether the column is in composition scope |
| `studies.tsv` | one row per primary study: identifiers, year, tier, analytical platform, country, and in-scope sample and observation counts |
| `observations_relative_pct.tsv` | one row per measurement reported as a relative percentage |
| `observations_absolute_conc.tsv` | the same for absolute concentrations |
| `observations_other.tsv` | the same for measurements in any other unit |
| `sample_material_types.tsv` | the material each measurement column describes: essential oil, fresh flower headspace, extract, and so on |
| `study_analytical.tsv` | the analytical method fields parsed from each article, against which the re-parse is checked |
| `compound_cache.json` | resolved compound identifiers, keyed by reported name |
| `smiles_cache.json` | SMILES strings for the resolved compounds |
| `article_text_cache.json` | the extracted text of each open-access PDF, keyed by file name |
| `fulltext/` | the JATS XML of each study available from Europe PMC, named by PMC accession |
| `oa_pdf/`, `oa_pdf_recovered/` | the open-access PDFs |

The same set of tables also contains `compounds.tsv` and `exclusions.tsv`, which
document the corpus but are not read here.

A missing input file raises `FileNotFoundError` naming the path; a missing column
raises `KeyError` naming the column. One case fails more quietly and is worth knowing
about: if the PDF directories are absent, studies without cached text are scored as
not reporting the text-derived items rather than as unscorable.
`10_commensurability_ladder.py` stops with a message naming the script to run first.

## Running

The scripts are numbered in the order they were developed and are listed in that
order below.

```
python3 analysis/01_within_study_contrasts.py
python3 analysis/02_contrast_results.py
python3 analysis/03_variance_partition.py
python3 analysis/04_species_cross_study.py
python3 analysis/05_direction_agreement.py
python3 analysis/06_prediction_baselines.py
python3 analysis/07_baseline_folds.py
python3 analysis/08_invariance_probe.py
python3 analysis/09_reporting_census.py
python3 analysis/10_commensurability_ladder.py
python3 analysis/11_published_claims.py
```

Most of them are independent of one another and read only the input tables.
One dependency is strict: `10` reads the article text cache written by `09`. Note
that `02` does not read `01`: it rebuilds the contrast inventory from the measurement
table, and `01` exists to write that inventory out as a table.

Every script writes tab-separated tables and, where the result is nested, JSON. Some
tables carry a hand-written caveat column recording what is wrong with a particular
study's design; those sentences are constants in the code, not generated text.

### Reproducibility

Every permutation, bootstrap and simulation draws from a generator seeded with a
fixed constant written in the source, so a rerun on the same input reproduces the
same numbers exactly. The constants are `20260912` for the composition tests and the
variance partition, `20260913` for the species comparison, the direction analysis,
the prediction baselines and the claim retests, and `20260914` for the ladder
counterfactual. Two tables carry their seed and permutation count as columns
(`species_cross_study_tests.tsv`, `ladder_counterfactual.tsv`); for the rest the seed
is the constant in the script that wrote them.

`09_reporting_census.py` draws every random quantity from one generator in file
order, so a row's value depends on how many tests preceded it; the script reproduces
exactly as a whole, but a single test cannot be rerun in isolation.

`04_species_cross_study.py` reads an optional `ROSE_N_PERMUTATIONS` (default 9999).
Changing it changes the p-values in `species_cross_study_tests.tsv`, which is why
that table records the value it used.

### Outputs

| script | writes |
|---|---|
| 01 | `within_study_contrasts.tsv` |
| 02 | `contrast_results.tsv`, `contrast_compound_effects.tsv`, `contrast_direction_agreement.tsv` |
| 03 | `sample_compound_matrix.tsv`, `sample_compound_matrix_absolute.tsv`, `variance_partition_tests.tsv`, `variance_separability.tsv` |
| 04 | `species_design.tsv`, `species_marker_profiles.tsv`, `species_cross_study_tests.tsv` |
| 05 | `direction_effects.tsv`, `direction_cross_study_pairs.tsv`, `direction_model.tsv`, `direction_summary.tsv`, `direction_effect_thresholds.tsv`, `direction_probes.json` |
| 06 | `baseline_results.tsv`, `probe_results.json` |
| 07 | `baseline_per_fold_r2.tsv`, `baseline_fold_tests.tsv` |
| 08 | `invariant_results.tsv`, `invariant_probes.json` |
| 09 | `reporting_census.tsv`, `reporting_study_items.tsv`, `reporting_vs_comparability.tsv`, `reporting_damascena_studies.tsv`, `reporting_item_cost.tsv`, `reporting_yield.tsv`, `reporting_summary.json`, and `article_text_cache.json` for any PDF it had to parse |
| 10 | `ladder_sensitivity.tsv`, `ladder_pair_detail.tsv`, `ladder_exclusion_reasons.tsv`, `ladder_counterfactual.tsv`, `ladder_by_axis.tsv` |
| 11 | `false_conclusions.tsv`, `main_figure_data.tsv` |

## Layout

```
comparability/            the library
  config.py               input and output directories
  samples.py              measurement columns and their metadata
  tables.py               JATS article parsing
  analytical.py           analytical method fields from article text
  methods.py              country, region and cultivar patterns for methods prose
  contrasts.py            enumeration and scoring of within-study contrasts
  contrast_analysis.py    composition tests on those contrasts
  stats.py                closure, CLR, distances, PERMANOVA, effect sizes
  separability.py         whether a factor is separable from laboratory identity
  species_analysis.py     the cross-study species comparison
  reporting.py            the reporting census
  meta.py                 reference intervals and marker profiles
  model/
    data.py               the modelling matrix
    baselines.py          metadata baselines
    evaluate.py           leave-study-out folds and scoring
    invariant.py          within-study centring
    direction.py          cross-study direction agreement
analysis/                 one script per analysis
```

## Three decisions in the code that change the answers

**Compounds are joined on the InChIKey skeleton in the matrix analyses and on the
full InChIKey in the marker analyses.** The full key distinguishes stereoisomers; the
skeleton does not. Geraniol and nerol are the case that matters here: they share a
skeleton but carry separate reference intervals. The contrast, variance and
prediction analyses join on the skeleton, because most studies do not report
stereochemistry and the full key would split one reported compound into several; the
species, reference-interval and marker analyses join on the full key, because there
the distinction is the point. Compounds with no resolved skeleton fall back to a
join on the cleaned compound name. The choice is fixed per analysis in the source,
not exposed as a parameter.

**"Not detected" and "not measured" are different cells.** A compound a study looked
for and did not find is a zero; a compound a study never reported is missing.
Conflating them turns a short compound list into a claim about absence. The long
format keeps four states — quantified, not detected, detected but not quantified,
not measured — and the matrix writers report the four counts per sample and encode
the last two distinctly. The modelling matrix collapses them into one value array
and one boolean array, in which a detected-but-not-quantified cell is marked present
with no number; compound filters count that array, so a compound seen only as a
trace can pass a filter it contributes no values to.

**No single coding of the ladder is presented as the answer.** All 32 combinations of
the five coding switches are computed over nine rung sets and written to
`ladder_sensitivity.tsv`. One coding is named in the source as the default and is
what the console prints; every other coding is in the table, and a result that
survives only one coding can be recognised as such.

## What this code does not do

It analyses the curated tables; it does not rebuild them. The reference intervals
used in the claim retests are hard-coded constants transcribed from open secondary
sources rather than from the standard itself, and the code says so in the caveat
column it writes. Figure and table generation is not included.

Several corpus-specific facts are hard-coded in the library: per-study caveats and
axis corrections keyed by accession, and a list of studies whose tables are not
volatile compositions. Run against a different corpus, those branches are silently
inert.

