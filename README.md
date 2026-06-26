# PRFT transcriptomic state in AML: processed data and analysis scripts

## Project overview

This repository staging package supports the manuscript describing a proteostasis-associated ferroptosis-tolerance (PRFT) transcriptomic state in acute myeloid leukemia (AML). It contains staged processed tables, selected analysis scripts, locked Formula A metadata, and supplementary tables intended for repository deposition after author review.

## Study purpose

The study uses public AML transcriptomic cohorts, BeatAML ex vivo pharmacogenomic data, processed single-cell module-score tables, and exploratory SigCom LINCS / CMap-L1000 perturbational outputs to evaluate a locked PRFT transcriptomic state. The repository is intended to make processed data tables and scripts inspectable for reviewers and readers.

## Data sources

The analyses use public or reused sources including TCGA-LAML, GEO validation cohorts, BeatAML, a public single-cell AML dataset, SigCom LINCS, and CMap/L1000 resources. Raw external datasets are not redistributed in this staging package unless authors verify that redistribution is allowed. The final Data Availability statement should list verified accessions, source URLs, versions or access dates.

## Folder structure

- `data/`: processed PRFT score tables and derived result evidence tables.
- `scripts/`: selected scripts supporting reproduction of processed analyses and major tables.
- `metadata/`: locked Formula A files and formula provenance notes.
- `supplementary/`: supplementary tables, figure-design tables, and source-data maps.
  
## Reproducibility note

This repository provides processed tables, locked Formula A metadata, supplementary result tables, and selected analysis scripts to support inspection and partial reproduction of the major derived analyses reported in the manuscript. Large public raw datasets and third-party-derived matrices are not redistributed here. Some scripts may require users to download the corresponding public datasets and adjust local input paths according to the folder structure described above. This repository is therefore intended as a transparent processed-data and script archive rather than a fully containerized one-command reproduction environment.

## Processed data description

The staged tables include PRFT risk scores, external validation summaries, BeatAML AUC association tables, processed single-cell module-score summaries, SigCom/LINCS input and output tables, and BeatAML-SigCom consistency tables. Large public raw data and selected large third-party-derived matrices require author review before redistribution.

## Scripts description

The staged scripts support selected preprocessing, PRFT scoring, external validation, BeatAML analysis, single-cell module localization, and SigCom/CMap preparation or curation. They are provided to support transparency of processed analyses and major tables. They are not yet guaranteed to provide a one-command complete reproduction from raw external downloads on a fresh computer.

## How to reproduce key tables

Use the scripts in `scripts/` together with the public source datasets and the staged metadata files. Before public deposition, authors should replace local absolute paths with relative paths, document required package versions, and verify that each script can run from the repository root or from a clearly documented working directory.

## Locked Formula A note

Locked Formula A was not modified in repository staging. The staged metadata files preserve the locked six-gene Formula A source, gene list, coefficients, and formula note. Supplementary Table S6 should remain the locked Formula A coefficient table unless the authors explicitly approve final supplementary renumbering.

## SigCom/LINCS note

SigCom LINCS and CMap/L1000 files are included as exploratory perturbational signature reversal resources. These outputs should be interpreted as computational prioritization and hypothesis generation, not as evidence that any perturbagen has therapeutic activity in AML.

## BeatAML AUC interpretation note

Higher BeatAML AUC values indicate lower ex vivo sensitivity / greater relative resistance. BeatAML outputs in this repository staging package should not be described as patient-level treatment guidance.

## Citation note

The final repository record should be cited using its verified repository DOI or persistent identifier after deposition. Public third-party datasets and software should also be cited in the manuscript/reference manager after human verification.  

## Contact

For questions about this repository after deposition, contact: 469309679@qq.com
