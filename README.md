# AML PRFT transcriptomic-state analysis

This repository provides research-use code, derived result tables, locked model metadata, figure-supporting summaries, supplementary-table support files, and repository manifests supporting the manuscript entitled:

**Integrated transcriptomic analysis identifies a six-gene proteostasis-associated ferroptosis-tolerance state in acute myeloid leukemia**

## Repository status

This repository supports a retrospective public transcriptomic and pharmacogenomic reanalysis of acute myeloid leukemia (AML). It provides analysis scripts, locked Formula A metadata, derived result tables, figure-supporting summaries, supplementary-table support files, data dictionaries, code manifests, and repository manifests for transparency, reviewer access, and reproducibility checking.

This repository is not a clinical decision-support tool. It does not provide treatment recommendations, patient-level response prediction, or evidence of clinical deployment readiness. The analyses should be interpreted as retrospective transcriptomic and ex vivo pharmacogenomic associations.

Project repository: https://github.com/chenmenghao525-sys/AML-PRFT-transcriptomic-state

Zenodo archive: https://doi.org/10.5281/zenodo.21214731

## Contents

- `code/`: analysis and packaging scripts retained for traceability. Scripts are provided for review and reuse, but the repository does not promise one-command full reproduction from raw public repositories.
- `docs/`: project notes and submission-support documentation retained for transparency and reviewer access.
- `figures/final/`: compatibility copy of final Figure 1–Figure 7 files copied from the submission-ready figure set.
- `figures/main/`: final main figures (`Figure 1.pdf` through `Figure 7.pdf`, with PNG backups where available).
- `figures/supplementary_figures/`: final supplementary figures (`Supplementary Figure S1.pdf` through `Supplementary Figure S8.pdf`, with PNG backups where available).
- `metadata/`: repository manifests, data dictionary, code manifest, and locked Formula A metadata.
- `metadata/formula_A/`: locked six-gene Formula A coefficients and formula notes.
- `results/`: derived result tables needed for reproducibility checks and manuscript traceability.
- `supplementary_tables/`: final Supplementary Tables S1–S15. Supplementary Table S10 and Supplementary Table S15 are valid tables and are retained.

## Locked model

The locked model is Formula A only:

```text
CLCN5, ARHGEF5, TRIM32, ITGB2, SAT1, ACOX2
```

The locked Formula A coefficients and formula notes are provided in `metadata/formula_A/`. The Formula A gene order and coefficients should not be replaced by machine-learning ranking outputs, PPI-derived feature strategies, pharmacogenomic analyses, or processed single-cell table outputs.
