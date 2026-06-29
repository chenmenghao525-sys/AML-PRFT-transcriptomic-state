# AML PRFT Human Genomics

This repository provides the research-use code, derived result tables, locked model metadata, figure-supporting summaries, and repository manifests supporting the AML PRFT Human Genomics manuscript.

## Repository status

This repository is prepared for private reviewer/editor access or public release after author confirmation. It is not a clinical decision-support tool and does not provide treatment recommendations.

Project repository: https://github.com/chenmenghao525-sys/AML-PRFT-Human-Genomics

No Zenodo DOI or GitHub release has been created in this working copy.

## Contents

- `code/`: analysis and packaging scripts retained for traceability. Scripts are provided for review and reuse, but the repository does not promise one-command full reproduction from raw public repositories.
- `docs/`: project notes and submission-support documentation retained for traceability.
- `figures/final/`: final key figure files retained for repository traceability.
- `figures/main/`: final main Figure 3-Figure 6 files retained from the public candidate package.
- `figures/supplementary_figures/`: supplementary figure support files retained for S1-S8; retired supplementary figure candidates are excluded from the public working tree.
- `metadata/`: repository manifests, data dictionary, code manifest, and locked Formula A metadata.
- `metadata/formula_A/`: locked six-gene Formula A coefficients and formula notes.
- `results/`: derived result tables needed for reproducibility checks and manuscript traceability.
- `supplementary_tables/`: Supplementary Tables S1-S15. Supplementary Table S10 is a valid table and is retained.

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

A conservative license note is included for review. A formal public license should be confirmed by the authors before changing repository visibility or archiving the repository.
