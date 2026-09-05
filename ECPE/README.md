# ECPE scoring study

This folder contains the ECPE application and stress test of Anders Hellström's
experimental item-weighted scoring system. The empirical data are the historical
`edmdata` ECPE grammar response matrix: 2,922 examinees by 28 binary items. They do
not represent the complete current ECPE examination, and the resulting scores are
not official ECPE scores.

## Main files

- [`ECPE_scoring_report.pdf`](ECPE_scoring_report.pdf) — the complete 54-page
  technical report.
- [`ECPE_scoring_report.tex`](ECPE_scoring_report.tex) and
  [`references.bib`](references.bib) — editable XeLaTeX source and bibliography.
- [`ECPE_report_source.zip`](ECPE_report_source.zip) — compact report source,
  figures, tables, code, selected data and reproduction records.
- [`ecpe_scoring_stress_test.r`](ecpe_scoring_stress_test.r) — standalone R stress
  test for the current Somers' Dxy/item-rest implementation.
- [`ecpe_scoring_stress_test_bundle.zip`](ecpe_scoring_stress_test_bundle.zip) —
  compact program, guide and completed default-run output bundle.
- [`ECPE_STRESS_TEST_README.md`](ECPE_STRESS_TEST_README.md) — detailed execution,
  formula, output and interpretation guide.
- [`diagnostic_plots.pdf`](diagnostic_plots.pdf) and the CSV files — selected
  completed-run diagnostics and reproducible result tables.

## Run the stress test

```sh
Rscript -e 'install.packages("edmdata", repos="https://cloud.r-project.org")'
Rscript ecpe_scoring_stress_test.r --out ecpe_results
```

The default run uses 500 calibration bootstraps, 100 replicates per sensitivity
condition and 500 permutations. Use a new or empty output directory.

The compact ZIP archives omit the 4.2 MiB `run_state.rds` simulation-state file to
fit the repository connector's per-file transfer limit. The readable report,
programs, figures and result tables are included; rerunning the stress test recreates
fresh simulation state. Each archive documents its exact omissions.

## Evidence boundaries

The study verifies the implementation and finds stable calibration on the original
28 items. Weighted and raw rankings are close, but the observed split-half
consistency improvement is very small. Artificial-item tests expose sensitivity to
rare uninformative items and duplicated content. These findings do not establish
external validity, fairness, improved proficiency measurement or an official ECPE
interpretation.

## AI-use disclosure

The scoring-system specification and project direction are Anders Hellström's.
OpenAI ChatGPT substantially contributed to the stress-test design, R implementation,
validation, supplementary analyses, report writing, typesetting and quality checks.
The report contains a prominent short disclosure at the beginning and a detailed
disclosure in its final appendix.
