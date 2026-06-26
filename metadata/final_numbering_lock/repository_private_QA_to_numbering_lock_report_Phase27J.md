# Phase27J repository private QA to numbering-lock report

## Scope

Phase27J performed a static repository-side figure, supplementary table, supplementary figure, and wording-consistency check for the private-upload candidate.

No analysis rerun was performed. No DEG, WGCNA, ML benchmarking, BeatAML, SigCom/LINCS, CMap/L1000, single-cell, cutoff-selection, model-training, or model-selection analysis was rerun. Locked Formula A was not modified. Formula B was not restored. No data files were modified.

The GitHub repository must remain private until public-release steps are explicitly approved. No public release, GitHub release, Zenodo DOI, or manuscript submission was performed in Phase27J.

## Files inspected

- `README.md`
- `metadata/Repository_upload_manifest_Phase27I_final.csv`
- `metadata/Data_dictionary_Phase27I_final.csv`
- `metadata/Code_manifest_Phase27I_final.csv`
- `metadata/README_formula_lock.md`
- `metadata/model_selection/phase3A_fix_original_6gene_vs_top_models.csv`
- `supplementary/Main_Figure_Original_Source_Map.csv`
- `supplementary/Phase17C_SigCom_BeatAML_evidence_tier_table.csv`
- `supplementary/Supplementary_Table_S6_LOCKED_formula_A_coefficients.csv`
- `supplementary/Supplementary_Table_Sx_SigCom_BeatAML_final_Phase17C.csv`
- `supplementary/Supplementary_Figure_Sx_SigCom_BeatAML_final_design_Phase17C.csv`
- `supplementary/ml_benchmarking/Supplementary_Table_S8__phase3A_fix_model_ranking_composite_score.csv`
- `supplementary/ml_benchmarking/phase3A_fix_model_plan_150.csv`
- `supplementary/ml_benchmarking/phase3A_fix_model_performance_success_only.csv`
- `supplementary/ml_benchmarking/phase3A_fix_model_performance_all.csv`
- `supplementary/ml_benchmarking/phase3A_fix_model_ranking_composite_score.csv`

## Numbering findings

- Main Figures 1-7 are supported by `supplementary/Main_Figure_Original_Source_Map.csv`, with all mapped panels marked for main-figure use. These figure labels require manuscript cross-check because much of the source artwork is author-held outside the repository.
- Supplementary Table S6 is preserved as `supplementary/Supplementary_Table_S6_LOCKED_formula_A_coefficients.csv`.
- Supplementary Table S8 is preserved as `supplementary/ml_benchmarking/Supplementary_Table_S8__phase3A_fix_model_ranking_composite_score.csv`.
- The SigCom/BeatAML supplementary table currently named `Supplementary_Table_Sx_SigCom_BeatAML_final_Phase17C.csv` is proposed as Supplementary Table S9, pending author confirmation and manuscript cross-check.
- The SigCom/BeatAML supplementary figure currently named `Supplementary_Figure_Sx_SigCom_BeatAML_final_design_Phase17C.csv` is proposed as Supplementary Figure S1, pending author confirmation and manuscript cross-check.
- `supplementary/Phase17C_SigCom_BeatAML_evidence_tier_table.csv` is retained as a repository support and wording-risk table, not assigned a new manuscript supplementary table number in this lock.

## Row-count checks

- `supplementary/Supplementary_Table_S6_LOCKED_formula_A_coefficients.csv`: 6 rows, matching the locked six-gene Formula A coefficient table.
- `supplementary/ml_benchmarking/phase3A_fix_model_plan_150.csv`: 150 planned combinations.
- `supplementary/ml_benchmarking/phase3A_fix_model_performance_success_only.csv`: 125 successfully fitted configurations.
- `supplementary/ml_benchmarking/Supplementary_Table_S8__phase3A_fix_model_ranking_composite_score.csv`: 126 ranking rows.
- `metadata/model_selection/phase3A_fix_original_6gene_vs_top_models.csv`: 17 Formula A versus top-model rows.
- `supplementary/Supplementary_Table_Sx_SigCom_BeatAML_final_Phase17C.csv`: 37 rows.
- `supplementary/Supplementary_Figure_Sx_SigCom_BeatAML_final_design_Phase17C.csv`: 4 figure-design rows.

The 126 S8 ranking rows should not be described as 126 successfully fitted configurations. The safe interpretation remains 150 planned combinations, 125 successfully fitted configurations, and 126 ranking rows including the locked Formula A reference model.

## Wording and risk scan

Phase27J targeted scan found no repository occurrences of:

- `Formula B`
- `old DEG`
- `126 successfully`
- `successfully fitted 126`
- `126 successfully fitted`
- `fully reproducible`
- `raw-data reproduction`
- `raw data reproduction`
- `causal treatment recommendation`

Occurrences of `therapeutic efficacy`, `patient-level treatment guidance`, and related efficacy terms were reviewed as cautionary or prohibited-overclaim wording, not as claims of efficacy or clinical guidance. The README explicitly states that BeatAML outputs should not be described as patient-level treatment guidance. SigCom/LINCS and CMap/L1000 files remain framed as exploratory perturbational signature-reversal and hypothesis-generating evidence.

## Consistency conclusion

No high-risk repository-manuscript numbering issue was detected in the static Phase27J audit. The only items requiring author confirmation are final manuscript-side numbering choices for the proposed SigCom/BeatAML Supplementary Table S9 and Supplementary Figure S1, plus final cross-check that main Figures 1-7 match the author-held manuscript figure set.
