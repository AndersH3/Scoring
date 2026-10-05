# Scoring — Experimental Item-Weighted Test Scoring

Research code and empirical stress tests for an experimental **item-weighted scoring system for dichotomous test items**. The current implementation combines item rarity/difficulty with positive item discrimination measured by **Somers' Dxy**, then uses the resulting item weights to compute normalized person scores.

The repository contains:

- the reference R implementation, `item_analysis.r`;
- ECPE empirical stress tests and a full technical report;
- PISA 2012 booklet-pattern and repeated-item stability analyses;
- generated result tables and figures;
- a GitHub Pages project browser at **https://scoring.hellstrom.pw/**.

> **Status:** research software. The scoring method is experimental. The included empirical analyses investigate implementation behavior, sensitivity, and internal stability; they do not establish external validity, fairness, causal interpretation, or superiority over established psychometric scoring models.

## Current scoring method

For item `j`, let:

- `n` = number of test takers;
- `c_j` = number answering item `j` correctly;
- `x_ij` = 0/1 response of person `i` to item `j`;
- `T_i,-j` = the default **item-rest score**, i.e. the number of other items answered correctly.

The program computes Harrell's Somers' Dxy using `Hmisc::somers2()`:

    Dxy_j = 2 * C_j - 1

where `C_j` is the concordance probability.

Only **positive discrimination** receives an additional weight reward:

    Cplus_j = max(0, Dxy_j)

Thus:

- `Dxy <= 0` gives `Cplus = 0`;
- `Dxy = 1` gives `Cplus = 1`;
- negative discrimination is not assigned a negative weight or penalty beyond receiving no positive-discrimination reward.

### Two padded factors

The difficulty/rarity pair is:

    rarity_pair_j = (c_j, n)

with

    pad1_j = (c_j + 1) / (n + 1)

The discrimination pair is:

    concordance_pair_j = (n * Cplus_j, n)

with

    pad2_j = (n - n * Cplus_j + 1) / (n + 1)

The fractional value `n * Cplus_j` is deliberately **not rounded**.

The raw item weight is:

    raw_weight_j =
        [1 - log(pad1_j)] *
        [1 - log(pad2_j)]

using the natural logarithm.

Operational weights are rescaled so that the smallest item weight is 1:

    weight_j = raw_weight_j / min_k(raw_weight_k)

This final rescaling is purely multiplicative. It preserves all item-weight ratios and does not change normalized person scores.

### Person score

The normalized weighted score for person `i` is:

    S_i = sum_j(x_ij * weight_j) / sum_j(weight_j)

so the ordinary output score is on a 0–1 scale, with a corresponding percentage also reported.

## Why item-rest scoring is the default

The default corresponding score for item discrimination is:

    T_i,-j = sum_{k != j} x_ik

rather than the total score including the target item itself.

This avoids the direct part-whole self-inclusion created when an item's own response is included in its discrimination criterion. Full total-score calibration remains available with:

```bash
--score-type total
```

The default is:

```bash
--score-type rest
```

## Constant and non-positive-discrimination items

For an all-zero or all-one item, empirical Somers' Dxy and concordance are undefined. The reference implementation retains the item and, **for weighting only**, uses a neutral pre-remap concordance value `C = 0.5`, which maps to:

    Cplus = 0

The reported empirical Somers statistics remain undefined rather than being presented as estimated values.

This design choice is important. In particular, an extremely rare but uninformative item can still receive a substantial **rarity factor** even when its discrimination factor receives no reward. The ECPE stress tests explicitly examine this behavior.

## Reference implementation

The main program is:

- [`item_analysis.r`](item_analysis.r)

### Input format

The program accepts a CSV-like response matrix containing complete binary 0/1 item responses.

Typical format:

```csv
id,item1,item2,item3,item4
p001,1,0,1,1
p002,0,1,1,0
p003,1,1,1,1
```

The first column can be used as a person identifier. With the default `--id-column auto`, the program attempts to recognize an ID column; `--no-id-column` treats every column as an item.

The core scoring program currently requires a **complete binary response matrix**: missing or non-binary item values are rejected.

### Required R packages

```r
install.packages(c("data.table", "optparse", "Hmisc"))
```

Optional research diagnostics use:

```r
install.packages(c("boot", "coin"))
```

- `boot` is required only when bootstrap replicates are requested.
- `coin` is required only when permutation replicates are requested.

### Basic run

```bash
Rscript item_analysis.r responses.csv
```

Default outputs are:

```text
item_weights.csv
person_scores.csv
```

You can choose the files explicitly:

```bash
Rscript item_analysis.r responses.csv \
  --items-output my_item_weights.csv \
  --persons-output my_person_scores.csv
```

### Useful options

```bash
# Use item-rest discrimination (default)
Rscript item_analysis.r responses.csv --score-type rest

# Use full total score instead
Rscript item_analysis.r responses.csv --score-type total

# Add diagnostic columns and internal validation
Rscript item_analysis.r responses.csv --diagnostics --validation

# Add research pair-count diagnostics
Rscript item_analysis.r responses.csv --research-diagnostics

# Bootstrap item-statistic/weight uncertainty
Rscript item_analysis.r responses.csv --bootstrap-reps 1000

# Monte Carlo permutation inference
Rscript item_analysis.r responses.csv --permutation-reps 5000

# Run built-in mathematical/integration tests
Rscript item_analysis.r --self-test
```

Run:

```bash
Rscript item_analysis.r --help
```

for the complete command-line interface.

## Main outputs

### `item_weights.csv`

The item-level audit includes, among other fields:

- counts correct/incorrect;
- Somers' Dxy and concordance probability;
- remapped positive concordance;
- rarity and concordance padded pairs;
- difficulty and discrimination weight factors;
- raw and operational item weights;
- warnings/status fields;
- optional pair-count diagnostics;
- optional bootstrap uncertainty;
- optional permutation p-values.

### `person_scores.csv`

The person-level output includes:

- number correct;
- unweighted proportion correct;
- weighted raw score;
- normalized weighted score;
- normalized weighted percentage;
- item-response columns;
- item contribution columns;
- normalization/status diagnostics.

Because operational item weights differ from raw formula weights only by a common multiplicative constant, the program checks that normalized person scores are unchanged by that rescaling.

## Validation features

The reference implementation contains several fail-fast consistency checks, including:

- the identity `Dxy = 2*C - 1`;
- concordant/discordant/tied pair-count identities;
- Mann-Whitney-U equivalence;
- optional rank-formula cross-checks;
- input validation;
- positive finite weight checks;
- weight-normalization invariance;
- optional frozen-output regression checks;
- built-in self-tests for padding, remapping, perfect concordance, reversal, ties, constants, and integration behavior.

Use:

```bash
Rscript item_analysis.r --self-test
```

and, for a real dataset:

```bash
Rscript item_analysis.r responses.csv --validation
```

## ECPE study

The [`ECPE/`](ECPE/) subdirectory contains the largest dedicated stress test of the current scoring system.

It uses the historical `edmdata::items_ecpe` grammar-response matrix:

- **2,922 examinees**;
- **28 binary grammar items**.

These data do **not** represent the complete current ECPE examination, and the weighted scores are not official ECPE scores.

Key files:

| File | Purpose |
| --- | --- |
| [`ECPE/README.md`](ECPE/README.md) | Concise subproject guide |
| [`ECPE/ECPE_STRESS_TEST_README.md`](ECPE/ECPE_STRESS_TEST_README.md) | Detailed formula, execution, outputs, validation and interpretation guide |
| [`ECPE/ECPE_scoring_report.pdf`](ECPE/ECPE_scoring_report.pdf) | Full technical report |
| [`ECPE/ECPE_scoring_report.tex`](ECPE/ECPE_scoring_report.tex) | XeLaTeX report source |
| [`ECPE/ecpe_scoring_stress_test.r`](ECPE/ecpe_scoring_stress_test.r) | Standalone stress-test program |
| [`ECPE/stress_test_results.md`](ECPE/stress_test_results.md) | Selected completed-run results |
| [`ECPE/diagnostic_plots.pdf`](ECPE/diagnostic_plots.pdf) | Diagnostic figures |
| [`ECPE/item_weights.csv`](ECPE/item_weights.csv) | Completed-run item weights |
| [`ECPE/person_scores.csv`](ECPE/person_scores.csv) | Completed-run person scores |

The ECPE analysis examines:

- bootstrap calibration uncertainty;
- calibration sample size;
- population-selection shifts;
- response corruption;
- null permutations;
- item deletion;
- added constant, rare, random, reversed and duplicated items;
- replicated calibration rows;
- split-half consistency;
- Q-matrix skill coverage;
- cross-fitted person scoring.

### ECPE evidence boundary

The completed stress test found that the implementation is reproducible and that weighted and raw score rankings are generally close on the historical ECPE grammar data. The observed split-half consistency improvement is very small, while artificial-item experiments expose sensitivity to rare uninformative items and duplicated content.

Those results **do not establish** that the proposed weighting method measures proficiency better than raw scoring, IRT, Rasch-family models, cognitive-diagnostic models, or other established psychometric methods.

## PISA analyses

The [`PISA/`](PISA/) subdirectory applies the scoring method to PISA 2012 U.S. mathematics response data and examines stability under the PISA booklet/design structure.

The booklet-stability analysis reconstructs the **13 largest structural response patterns** in the `edmdata::items_pisa12_us_math` matrix and recalibrates repeated items across the major booklet-pattern groups.

The analysis distinguishes:

- ordinary form-specific item-rest calibration;
- a common-anchor rest score based on items shared across patterns containing the target item;
- raw formula weight stability;
- operational minimum-to-one weight stability;
- item difficulty and Somers' Dxy changes across booklet patterns.

Key files include:

| File | Purpose |
| --- | --- |
| [`PISA/report.pdf`](PISA/report.pdf) | PISA scoring report |
| [`PISA/pisa_booklet_stability_analysis.r`](PISA/pisa_booklet_stability_analysis.r) | Booklet/repeated-item stability analysis |
| [`PISA/pisa_repeated_item_stability.csv`](PISA/pisa_repeated_item_stability.csv) | Repeated-item stability summaries |
| [`PISA/pisa_repeated_item_pairwise.csv`](PISA/pisa_repeated_item_pairwise.csv) | Pairwise cross-pattern comparisons |
| [`PISA/pisa_all_pattern_item_calibrations.csv`](PISA/pisa_all_pattern_item_calibrations.csv) | Pattern-specific item calibrations |
| [`PISA/pisa_major_pattern_summary.csv`](PISA/pisa_major_pattern_summary.csv) | Major structural booklet-pattern groups |
| [`PISA/pisa_repeated_item_weight_stability.png`](PISA/pisa_repeated_item_weight_stability.png) | Weight-stability visualization |
| [`PISA/item_weights.csv`](PISA/item_weights.csv) | Example/calibration item weights |
| [`PISA/person_scores.csv`](PISA/person_scores.csv) | Corresponding person scores |

The PISA analyses are intended as **stability and transportability diagnostics**, not as an official PISA scoring procedure.

## Repository structure

```text
Scoring/
├── item_analysis.r          # reference scoring implementation
├── ECPE/                    # ECPE stress test, report, tables and figures
├── PISA/                    # PISA analyses, report, tables and figures
├── index.html               # GitHub Pages front end
├── repository-index.js      # generated/site repository index
├── analytics.js             # site analytics helper
├── CNAME                    # scoring.hellstrom.pw
├── exclude.txt              # website/index exclusions
└── LICENSE                  # AGPL-3.0
```

## Research questions this repository can address

The code is useful for questions such as:

- How strongly do the proposed weights depart from unit/raw scoring?
- Are rare items given disproportionate influence?
- How much does positive discrimination change the rarity-only component?
- How stable are weights across calibration samples?
- How stable are repeated-item weights across different test forms/booklets?
- How sensitive are rankings and scores to item deletion, duplication, contamination or calibration population shifts?
- Does item-rest calibration behave differently from total-score calibration?
- How much uncertainty is introduced by estimating the weights from finite samples?

It should not, by itself, be used to answer questions such as:

- Does the weighted score measure a latent trait more validly?
- Is it fair across demographic groups?
- Is it invariant across populations?
- Does it outperform IRT or Rasch scoring on external criteria?
- Does a larger weight imply that an item is intrinsically more important?
- Can the weighted score be interpreted as an official ECPE or PISA proficiency score?

Those require additional model-based and external validation evidence.

## Important methodological limitations

1. **Endogenous calibration criterion.** Item discrimination is measured against a score built from the same test. Item-rest scoring reduces direct self-inclusion but does not create an external criterion.
2. **Rarity is rewarded by design.** A rare item can receive a substantial difficulty factor even if it has no positive discrimination. This is a substantive scoring assumption, not a statistical theorem.
3. **Truncation at zero.** All non-positive Dxy values receive the same discrimination factor of 1 after remapping. The method therefore distinguishes positive discrimination but does not penalize negative discrimination.
4. **Sample dependence.** Difficulty and discrimination are empirical properties of the calibration sample, so weights can change across populations and forms.
5. **Minimum-to-one normalization is cosmetic.** It aids interpretation of item weights but does not alter normalized scores.
6. **No latent-variable model is implied.** The method is not an IRT/Rasch/CDM model and its item weights are not item parameters in those senses.
7. **Reliability is not validity.** High rank agreement, split-half consistency, or bootstrap stability cannot by themselves demonstrate better construct measurement.
8. **External validation remains necessary.** Criterion validity, fairness, population invariance, and practical decision performance require independent evidence.

## Reproducibility

The repository retains both code and many completed-run artifacts so that formulas, diagnostics and downstream results can be inspected without rerunning every analysis.

The ECPE stress-test workflow additionally records:

- settings;
- package versions;
- input checksums;
- random seed/state;
- output manifests;
- selected completed-run tables and plots.

The PISA directory contains both executable R source and generated calibration/stability tables.

## Project website

A browsable view of repository subprojects and files is published at:

**https://scoring.hellstrom.pw/**

The site is served from this repository through GitHub Pages using the included `CNAME`, `index.html`, and JavaScript index files.

## License

This repository is licensed under the **GNU Affero General Public License v3.0 (AGPL-3.0)**. See [LICENSE](LICENSE).

## Citation and reuse

If you use the scoring method or code in research, cite this repository and also cite the original sources for any empirical datasets used in the corresponding subproject.

The ECPE subproject contains detailed source references and an explicit AI-use disclosure. The PISA directory includes links to OECD/PISA source material and related references.

## Author / project direction

The scoring-system specification and project direction are by **Anders Hellström**. See the subproject documentation for dataset-specific acknowledgements, source citations, validation records, and AI-use disclosures.
