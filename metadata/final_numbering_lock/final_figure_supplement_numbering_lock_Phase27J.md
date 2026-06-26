# Phase27J final figure and supplement numbering lock

Static repository-side numbering lock for the Human Genomics pre-submission repository candidate.

No analysis rerun was performed. Locked Formula A genes, coefficients, cutoff, and grouping rules were not modified. Formula B was not restored. This file records repository-manuscript numbering consistency only.

## Main figures

| Proposed final label | Proposed final title | Evidence source file | Status | Caution note |
|---|---|---|---|---|
| Figure 1 | Overall workflow, study design, PRFT transcriptomic state framework, and fixed six-gene overview | `supplementary/Main_Figure_Original_Source_Map.csv` | needs manuscript cross-check | Supported by Figure 1 panels A-D in the source map; source artwork is mostly author-held outside this public data archive. |
| Figure 2 | PRFT-high differential expression and WGCNA module context | `supplementary/Main_Figure_Original_Source_Map.csv` | needs manuscript cross-check | Supported by Figure 2 panels A-D; do not treat this lock as permission to rerun DEG or WGCNA. |
| Figure 3 | Candidate filtering, Cox/LASSO derivation context, and six-gene coefficient display | `supplementary/Main_Figure_Original_Source_Map.csv` | needs manuscript cross-check | Supported by Figure 3 panels A-D; Locked Formula A must remain unchanged. |
| Figure 4 | PRFT validation, survival/ROC summaries, and clinical model context | `supplementary/Main_Figure_Original_Source_Map.csv`; `metadata/model_selection/phase3A_fix_original_6gene_vs_top_models.csv` | needs manuscript cross-check | Data dictionary links Formula A/top-model comparison to Figure 4 and Supplementary Table S8. |
| Figure 5 | Machine-learning benchmarking, PRFT-high state recognition, SHAP, and PPI ranking | `supplementary/Main_Figure_Original_Source_Map.csv`; `supplementary/ml_benchmarking/Supplementary_Table_S8__phase3A_fix_model_ranking_composite_score.csv` | needs manuscript cross-check | Describe ML as 150 planned combinations and 125 successfully fitted configurations; 126 ranking rows include the locked Formula A reference model. |
| Figure 6 | Bulk pathway and processed single-cell PRFT localization summaries | `supplementary/Main_Figure_Original_Source_Map.csv` | needs manuscript cross-check | Supported by Figure 6 panels A-D; do not introduce new single-cell methods. |
| Figure 7 | BeatAML pharmacogenomic association and interpretation summary | `supplementary/Main_Figure_Original_Source_Map.csv`; `data/beataml_*`; `data/BeatAML_SigCom_consistency_check_Phase17B_v2.csv` | needs manuscript cross-check | Higher BeatAML AUC means lower ex vivo sensitivity / greater relative resistance; not patient-level treatment guidance. |

## Supplementary tables

| Proposed final label | Proposed final title | Current file name | Evidence source file | Status | Caution note |
|---|---|---|---|---|---|
| Supplementary Table S6 | Locked Formula A gene coefficients | `Supplementary_Table_S6_LOCKED_formula_A_coefficients.csv` | `metadata/LOCKED_PRFT_six_gene_formula_A_coefficients.csv`; `metadata/README_formula_lock.md` | locked | Must remain the locked six-gene Formula A coefficient table. |
| Supplementary Table S8 | Phase3A-fix machine-learning benchmarking ranking and composite-score table | `Supplementary_Table_S8__phase3A_fix_model_ranking_composite_score.csv` | `supplementary/ml_benchmarking/phase3A_fix_model_plan_150.csv`; `supplementary/ml_benchmarking/phase3A_fix_model_performance_success_only.csv`; `metadata/model_selection/phase3A_fix_original_6gene_vs_top_models.csv` | locked | S8 has 126 ranking rows; manuscript-facing text should state 125 successfully fitted configurations plus the locked Formula A reference model where needed. |
| Supplementary Table S9 | Exploratory SigCom LINCS / CMap-L1000 perturbational reversal and BeatAML consistency summary | `Supplementary_Table_Sx_SigCom_BeatAML_final_Phase17C.csv` | `supplementary/Phase17C_SigCom_BeatAML_evidence_tier_table.csv`; `data/SigCom_LINCS_reversal_results_curated_Phase17B_v2.csv`; `data/BeatAML_SigCom_consistency_check_Phase17B_v2.csv` | pending author confirmation | Proposed as S9 because S6 and S8 are already fixed; manuscript must confirm no existing Supplementary Table S9 conflict. Interpret as exploratory and hypothesis-generating only. |

## Supplementary figures

| Proposed final label | Proposed final title | Current file name | Evidence source file | Status | Caution note |
|---|---|---|---|---|---|
| Supplementary Figure S1 | Exploratory SigCom LINCS / CMap-L1000 perturbational reversal and BeatAML consistency design | `Supplementary_Figure_Sx_SigCom_BeatAML_final_design_Phase17C.csv` | `supplementary/Phase17C_SigCom_BeatAML_evidence_tier_table.csv`; `data/CMap_LINCS_input_signature_LOCKED_A_Phase17A_lite.csv`; `data/SigCom_LINCS_reversal_results_curated_Phase17B_v2.csv`; `data/BeatAML_SigCom_consistency_check_Phase17B_v2.csv` | pending author confirmation | Proposed as S1 only within this repository-side lock. Final manuscript must confirm no existing supplementary-figure numbering conflict. Do not describe as therapeutic efficacy or clinical treatment evidence. |

## Repository support files not assigned manuscript numbering

| File | Repository role | Status | Caution note |
|---|---|---|---|
| `supplementary/Main_Figure_Original_Source_Map.csv` | Main-figure source map for Figures 1-7 | locked as repository support file | Not itself a manuscript figure or supplementary table. |
| `supplementary/Phase17C_SigCom_BeatAML_evidence_tier_table.csv` | Evidence-tier and wording-risk table for SigCom/BeatAML interpretation | locked as repository support file | Use to police exploratory wording; not assigned a new supplementary table number here. |
| `supplementary/ml_benchmarking/phase3A_fix_model_plan_150.csv` | ML planned-combination support table | locked as repository support file | Supports 150 planned model combinations. |
| `supplementary/ml_benchmarking/phase3A_fix_model_performance_success_only.csv` | ML success-only support table | locked as repository support file | Supports 125 successfully fitted configurations. |
| `supplementary/ml_benchmarking/phase3A_fix_model_performance_all.csv` | ML all-attempted support table | locked as repository support file | Documents planned and attempted model records. |
| `supplementary/ml_benchmarking/phase3A_fix_model_ranking_composite_score.csv` | Unprefixed copy/source of the S8 ranking table | locked as repository support file | Keep consistent with S8 copy; no manuscript-facing renumbering needed unless authors choose to rename later. |

## Locked interpretation guardrails

- Supplementary Table S6 remains the locked Formula A coefficient table.
- Supplementary Table S8 remains the ML benchmarking ranking/composite table.
- SigCom LINCS, CMap/L1000, and BeatAML convergence outputs remain exploratory, perturbational, and hypothesis-generating only.
- No file reviewed in Phase27J should be used to imply therapeutic efficacy, patient-level treatment guidance, or causal treatment recommendation.
- Repository wording remains a processed-data and selected-script archive, not a fully reproducible one-command raw-data reproduction package.
- Formula B and old DEG wording were not detected by the Phase27J targeted scan.
- Manuscript-facing ML wording must not say 126 successfully fitted configurations.
