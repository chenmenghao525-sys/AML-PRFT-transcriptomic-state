# PRFT transcriptomic state in AML: processed data and selected-script archive

## Project title

Proteostasis-associated ferroptosis-tolerance transcriptomic state in acute myeloid leukemia.

## Repository purpose

This repository provides processed tables, locked Formula A metadata, supplementary result tables, machine-learning benchmarking evidence, and selected analysis scripts to support transparency, result inspection, and partial reproduction of the major derived analyses reported in the manuscript.

This repository should be interpreted as a processed-data and selected-script archive rather than a fully containerized one-command reproduction environment.

## Folder structure

- `data/`: processed PRFT score tables, external-validation summaries, BeatAML pharmacogenomic association tables, processed single-cell summaries, SigCom/LINCS input or output tables, and BeatAML-SigCom consistency tables.
- `metadata/`: locked Formula A files, formula provenance notes, and model-selection comparison metadata.
- `metadata/model_selection/`: comparison table for the locked Formula A baseline and top-ranked survival-learning configurations.
- `scripts/`: selected scripts supporting preprocessing, PRFT scoring, external validation, BeatAML analysis, single-cell module localization, and SigCom/CMap preparation or curation.
- `supplementary/`: supplementary tables, source-data maps, figure-design tables, and manuscript-supporting evidence tables.
- `supplementary/ml_benchmarking/`: Phase3A-fix survival-learning benchmarking tables.

## Reproducibility note

The files are staged to make processed results inspectable and to support partial reproduction of major derived analyses when paired with the public source datasets and the required R or Python packages. Raw third-party datasets are not redistributed here unless redistribution is permitted and confirmed by the authors.

## Processed data description

The processed tables include PRFT risk scores, external-validation performance summaries, BeatAML ex vivo pharmacogenomic association summaries, processed single-cell module-score summaries, SigCom/LINCS perturbational signature tables, and BeatAML-SigCom consistency evidence. Column-level interpretation should be checked against the final data dictionary and supplementary tables before public release.

## Selected scripts description

The selected scripts document the analysis logic used to generate or curate major processed outputs. They are provided for transparency and partial reproduction, not as a guaranteed end-to-end raw-data pipeline on a fresh computer. Before public release, package versions and any repository-specific run instructions should be verified by the authors.

## Machine-learning benchmarking evidence note

The machine-learning benchmarking files document 150 planned model combinations and 125 successfully fitted configurations across multiple survival-learning frameworks. These files support the manuscript statement of more than 100 machine-learning model configurations and should not be interpreted as 100 distinct machine-learning algorithms.

## Formula A lock note

Locked Formula A was not modified during repository staging. The repository preserves the locked six-gene Formula A gene list and coefficients. The final manuscript should remain aligned with the locked gene list, coefficients, cutoff rule, and model-selection comparison files.

## Data source note

The study reuses public AML transcriptomic cohorts, BeatAML ex vivo pharmacogenomic data, processed single-cell AML resources, SigCom LINCS, and CMap/L1000 resources. Final accessions, source URLs, versions, access dates, and third-party reuse permissions must be verified by the authors before public release.

## Citation / DOI placeholder

DOI will be inserted after Zenodo archiving. The final repository record should be cited using its verified DOI or persistent identifier after deposition.

## Contact

For questions about this repository after deposition, contact: [CONTACT EMAIL TO CONFIRM]
