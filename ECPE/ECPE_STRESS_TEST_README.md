# ECPE stress testing for the current scoring system

This standalone R program tests the **current Somers Dxy version** of your scoring system against `edmdata::items_ecpe`: 2,922 examinees and 28 binary grammar items. It also reads the 28-by-3 expert Q matrix. This dataset covers ECPE grammar items; the resulting weighted scores are experimental, not official ECPE scores.

The implementation matches `AndersH3/Scoring/item_analysis.r`, Git blob `a336d1b3d9964bd5f2e0f90f214e15d540be827d`, retrieved on 1 September 2026. The older rank-cut agreement algorithm is not used.

## Install and run

Requires R 4.1 or newer. Only one additional package is required:

```bash
Rscript -e 'install.packages("edmdata", repos="https://cloud.r-project.org")'
Rscript ecpe_scoring_stress_test.R --out ecpe_results
```

Run these commands in the directory containing the downloaded program. On your Fedora setup, use the toolbox where you already run R. The program does not install software automatically, modify your existing scoring program, or upload results.

The output directory must be new or empty. Each run records its settings, package versions, input checksums, and random seed.

For a fast functional check:

```bash
Rscript ecpe_scoring_stress_test.R --quick --out ecpe_quick
```

Quick mode uses 30 bootstraps, 10 replicates per sensitivity condition and 30 permutations. These are too few for stable tail estimates or precise p-values.

For more Monte Carlo precision:

```bash
Rscript ecpe_scoring_stress_test.R \
  --bootstrap 2000 \
  --repeats 500 \
  --permutations 5000 \
  --seed 20260901 \
  --out ecpe_extended
```

Default counts are 500 bootstraps, 100 replicates per sensitivity condition, and 500 permutations. All random sampling runs sequentially for reproducibility. The full default run has already been executed for the supplied results.

Other options:

```bash
Rscript ecpe_scoring_stress_test.R --self-test
Rscript ecpe_scoring_stress_test.R --help
Rscript ecpe_scoring_stress_test.R --score-type total --out ecpe_total
Rscript ecpe_scoring_stress_test.R --sizes 50,100,250,500,1000 --out ecpe_sizes
```

`--score-type rest` is the default and matches the current source. Rest-versus-total diagnostics and both permutation nulls are produced in either mode. Sample sizes above the available calibration size are omitted; the entire calibration sample is always included as a zero-resampling-error endpoint.

## Exact formula

For item j, n calibration examinees and c correct responses:

1. Compute the corresponding score as the sum of the other items (default), or the full raw total if explicitly requested.
2. Calculate C as the fraction of correct-versus-incorrect pairs ordered in the expected direction, giving a half point to tied corresponding scores.
3. Set `Dxy = 2*C - 1` and `Cplus = max(0, Dxy)`.
4. Set `pad1 = (c+1)/(n+1)` and `pad2 = (n-n*Cplus+1)/(n+1)`.
5. Set `raw_weight = (1-log(pad1))*(1-log(pad2))` using natural logarithms.
6. Divide raw weights by their minimum to obtain operational weights with minimum 1.
7. The person score is `sum(response*weight)/sum(weight)`, on a 0–1 scale.

The fractional value `n*Cplus` is never rounded. Empirical C and Dxy remain undefined for constant items. For weighting only, these items receive Cplus=0. Negative Dxy likewise gives Cplus=0: the rarity factor remains and there is no negative-discrimination penalty. These choices reproduce your source program.

The stress-test kernel counts pairs on the small integer score support instead of repeatedly loading Hmisc. This is mathematically equivalent to `Hmisc::somers2(corresponding_score, item_response)` with half-credit ties. If Hmisc is already installed, the built-in self-tests check against it too.

## Experimental design

A seeded random split assigns **2,045 examinees to calibration** and **877 to audit** by default. All repeated calibration experiments estimate weights solely from calibration responses, then score the original audit responses. A separate five-fold cross-fit provides a score for every examinee using weights estimated without that examinee's fold.

| Test | What it examines |
|---|---|
| Full-data item and person audit | Difficulty, C, Dxy, tie frequency, factors, raw/operational/normalized weights, weighted scores, and rank reversals |
| Calibration-person bootstrap | Sampling uncertainty in item weights and resulting audit scores; all items and their normalization are recalculated together |
| Calibration size | Stability at n=50, 100, 250, 500, 1,000, 2,000 and 2,045; sampling without replacement within each replicate |
| Population selection | Calibration using lower-score, upper-score or unselected examinees; common calibration size and untouched audit sample |
| Response corruption | 1%, 5% and 10% cell flips, or replaced random/all-zero/all-one respondent rows in calibration |
| Independent-column permutation | Null discrimination while preserving every item's observed correct count; comparison of rest-score and total-score self-inclusion effects |
| Item deletion | Score and ranking changes when each item is removed; separates deletion from recalibration of remaining weights |
| Added artificial items | All-zero, all-one, rare 1%-correct random, 50%-correct random, reversed, and duplicated items |
| Replicated calibration rows | Repeating exactly the same observations 2, 5, 10 or 20 times to isolate the sample-size effect of +1 padding |
| Split-half consistency | Raw and weighted consistency across 100 random 14/14 item splits, with weights calibrated independently of audit persons |
| Q-matrix coverage | Whether weighting changes coverage of morphosyntactic, cohesive and lexical rules; skills can overlap |

The artificial-item scenarios are one reproducible draw each, not a repeated sampling study. They append items to the original form. Reversed and duplicated items use the item with the largest calibration weight. Random 50% responses are an explicit stress condition, not an estimate of ECPE guessing probability.

## Read these outputs first

| File in the output directory | Contents |
|---|---|
| `REPORT.md` | Run summary, formula, key results and interpretation limits |
| `diagnostic_plots.pdf` | Five pages: scores/weights, bootstrap intervals, calibration sizes, null/consistency diagnostics and contamination |
| `item_weights_full.csv` | Complete item-level scoring audit on all 2,922 examinees |
| `person_scores.csv` | Raw, full-fit and cross-fitted scores, average-rank percentiles, fold and calibration/audit membership |
| `bootstrap_item_intervals.csv` | Pointwise 95% percentile intervals for normalized, raw and operational weights and Dxy |
| `bootstrap_audit_score_intervals.csv` | Calibration-induced score uncertainty for the 877 audit examinees |
| `sample_size_summary.csv` | Monte Carlo summaries of score/rank sensitivity by calibration size |
| `population_shift_summary.csv` | Results under raw-score-selected calibration populations |
| `contamination_summary.csv` | Sensitivity to calibration-response corruption |
| `permutation_item_tests.csv` | One-sided permutation tests for positive association, with BH-adjusted p-values separately for rest and total |
| `augmented_item_stress.csv` | Score/rank effects of constant, rare, reversed and duplicate items |
| `item_deletion.csv` | Effects of omitting each original item |
| `sample_replication.csv` | Changes caused solely by the effective n in the padding formula |
| `split_half_consistency.csv` | Raw and weighted split-half correlations and Spearman–Brown adjustments |
| `run_state.rds` | Bootstrap draws, audit scores, column/person orders, split IDs, fold IDs and initial/final RNG states |
| `settings.R`, `sessionInfo.txt`, `input_md5.csv`, `output_manifest.csv` | Reproducibility record; the manifest excludes itself |

Additional CSVs contain individual replicate results, rest-versus-total comparisons, original response/Q matrices and within-fold diagnostics. `program_used.R` is a copy of the executed program when run through Rscript or R's `--file` option.

## Interpreting the numbers

- **Score RMSE** and absolute score changes use the 0–1 scale. Multiply by 100 for percentage-score points. RMSE 0.01 is one percentage-score point.
- **Rank shifts** are already expressed in percentile percentage points, using average ranks for ties. They compare positions within the same audit sample.
- **Top-decile overlap** is between 0 and 1. Examinees tied at the boundary receive fractional membership, avoiding arbitrary row-order decisions.
- **Weight share** is raw weight divided by the sum of raw weights. It is the useful quantity for comparisons between calibrations: the minimum-to-one operational scaling cancels from person scores.
- **Effective item count** is `1/sum(share^2)`. It measures weight concentration only; it is not the number of statistically independent items or a reliability measure.
- **share_L1** is the sum of absolute differences between normalized weights, between 0 and 2.
- **Raw-score reversals** count pairs where the person with fewer correct answers receives a strictly larger weighted score. Raw-score ties are excluded; weighted differences of at most 1e-12 are treated as equal for this count.
- The `p025` and `p975` columns in sensitivity summaries are the empirical spread across generated scenarios. They are not confidence intervals for the reported mean.

The bootstrap reference is the original calibration fit, not a known true proficiency score. The bootstrap intervals quantify calibration uncertainty conditional on fixed observed item responses; they do not measure all test-score error or latent ability uncertainty. Intervals are pointwise, and truncation at Dxy=0 can make their coverage irregular.

High agreement with raw scores and high consistency do not establish that the weighting system is better. ECPE has multiple grammar skills, so even Spearman–Brown-adjusted random-half correlations should be treated as descriptive. This program does not fit an IRT/CDM model, estimate demographic fairness, apply official pass thresholds, or claim external criterion validity.

The selected-score population experiments can themselves induce or alter associations through selection. The augmented-item and deletion scenarios alter the test form and denominator. Duplicate rows add no information. The reports explicitly distinguish these effects from ordinary calibration uncertainty.

## Validation performed for this delivery

- Ran the program under R 4.3.3 with edmdata 1.3.0 on the actual 2,922-by-28 matrix.
- Passed 120 exhaustive small-sample pair/rank comparisons plus checks for perfect/reversed discrimination, tied scores, constants, input rejection, row/item ordering, score endpoints and weight normalization.
- Compared 174 item/dataset/mode combinations with the current source's pair-count and weighting functions and the retrieved upstream Hmisc `somers2` function. Maximum C difference was 1.67e-15; maximum raw-weight difference was 4.45e-15.
- Checked 100 reversal-count cases against explicit pair matrices.
- Verified output checksums and expected row/replicate counts; reopened the saved run state in a fresh R process; visually inspected all five plot pages.
- Completed the full default run: 500 calibration bootstraps, 100 replicates per sensitivity condition and 500 permutations, with all 28 items retained.

The Hmisc package itself was not installed in the validation environment; the standalone upstream function was used for the separate source-equivalence check. Runtime self-test output distinguishes this from an installed-package check.

## Sources and authorship

- [ECPE response data](https://tmsalab.github.io/edmdata/reference/items_ecpe.html)
- [ECPE expert Q matrix](https://tmsalab.github.io/edmdata/reference/qmatrix_ecpe.html)
- [Your current scoring source](https://github.com/AndersH3/Scoring/blob/main/item_analysis.r)
- [Hmisc concordance implementation](https://github.com/harrelfe/Hmisc/blob/master/R/somers2.s)

The program exports `citation("edmdata")` into each result directory. The underlying ECPE references include Templin and Hoffman (2013), DOI 10.1111/emip.12010, and Templin and Bradshaw (2014), DOI 10.1007/s11336-013-9362-0.

The scoring-system specification is Anders Hellström's. OpenAI ChatGPT generated this stress-test design, R implementation and guide, and executed the recorded validation. No new empirical superiority claim is implied.
