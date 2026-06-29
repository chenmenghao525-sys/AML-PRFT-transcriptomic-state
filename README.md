# AML PRFT Human Genomics

This repository provides the research-use code, derived result tables, locked model metadata, figure-supporting summaries, and repository manifests supporting the AML PRFT Human Genomics manuscript.

## Repository status

This repository is the public research repository supporting the manuscript entitled "A proteostasis-associated ferroptosis-tolerance transcriptomic state in acute myeloid leukemia is associated with adverse prognosis, myeloid immune programs and ex vivo pharmacogenomic patterns". It provides analysis scripts, locked model metadata, derived result tables, figure-supporting summaries, supplementary-table support files, and repository manifests for transparency and reviewer access. It is not a clinical decision-support tool and does not provide treatment recommendations.

Project repository: https://github.com/chenmenghao525-sys/AML-PRFT-Human-Genomics

No Zenodo DOI has been minted at this stage. No GitHub release is required for this repository state.

## Contents

- `code/`: analysis and packaging scripts retained for traceability. Scripts are provided for review and reuse, but the repository does not promise one-command full reproduction from raw public repositories.
- `docs/`: project notes and submission-support documentation retained for traceability.
- `figures/final/`: compatibility copy of final Figure 1-Figure 7 files copied from the final upload package.
- `figures/main/`: final main figures copied from the final upload package (`Figure 1.pdf` through `Figure 7.pdf`, with PNG backups where available).
- `figures/supplementary_figures/`: final supplementary figures copied from the final upload package (`Supplementary Figure S1.pdf` through `Supplementary Figure S8.pdf`, with PNG backups where available).
- `metadata/`: repository manifests, data dictionary, code manifest, and locked Formula A metadata.
- `metadata/formula_A/`: locked six-gene Formula A coefficients and formula notes.
- `results/`: derived result tables needed for reproducibility checks and manuscript traceability.
- `supplementary_tables/`: final Supplementary Tables S1-S15 copied from the final upload package. Supplementary Table S10 and Supplementary Table S15 are valid tables and are retained.

## Scientific boundaries

- The locked model is Formula A only: CLCN5, ARHGEF5, TRIM32, ITGB2, SAT1, and ACOX2.
- BeatAML analyses are interpreted as ex vivo pharmacogenomic or drug AUC associations only.
- Single-cell evidence is limited to processed-table or score-level support; no de novo single-cell analysis is claimed here.
- PPI results are used for network-level prioritization and do not establish mechanism proof.
- The repository does not claim clinical treatment response, clinical efficacy, or treatment recommendation.

## Data and code availability

Public transcriptomic and pharmacogenomic datasets remain available from their original repositories, including TCGA, GEO, BeatAML, and other sources cited in the manuscript. Raw public datasets are not redistributed here.

Processed result tables, locked model metadata, figure-supporting summaries, and code/metadata manifests generated for this study are available in this repository: https://github.com/chenmenghao525-sys/AML-PRFT-Human-Genomics

## License status

A conservative license note is included for review. The repository is intended for research transparency and reproducibility support. The license note should be reviewed by the authors before any archival DOI workflow.


Public-release boundaries:

- No patient-level identifiable information is included.
- No raw controlled-access clinical data are redistributed.
- The repository is intended for research transparency and reproducibility support.
- No Zenodo DOI has been minted at this stage.
