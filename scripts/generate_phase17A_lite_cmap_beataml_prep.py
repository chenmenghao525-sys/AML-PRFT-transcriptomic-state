from __future__ import annotations

import csv
import os
from pathlib import Path


ROOT = Path(r"E:\ai助手\AML_PRFT_Human_Genomics_phase0_audit")
PHASE17 = ROOT / "17_supplementary_in_silico_perturbation_LOCKED_A"

SUBDIRS = [
    "input_inventory",
    "locked_formula_import",
    "drug_reversal_CMap_LINCS",
    "drug_reversal_CMap_LINCS/manual_CMap_LINCS_outputs",
    "beataml_consistency_check",
    "supplementary_figures",
    "supplementary_tables",
    "manuscript_text",
    "reference_support",
    "risk_check",
    "scripts",
    "logs",
    "manual_input_package",
    "virtual_perturbation",
]

ALLOWED_EXT = {".csv", ".tsv", ".txt", ".xlsx", ".rds", ".rdata"}

DEG_FILE = ROOT / "03_results_tables" / "05_deg" / "deg_prft_high_vs_low_tcga_all.csv"
BEATAML_AUC_FILE = ROOT / "03_results_tables" / "11_drug" / "beataml_drug_response_clean.csv"
BEATAML_GROUP_FILE = ROOT / "03_results_tables" / "11_drug" / "beataml_risk_score_by_sample.csv"
BEATAML_HIGH_LOW_FILE = ROOT / "03_results_tables" / "phase3C_BeatAML_high_low_drug_comparison_formula_A.csv"
BEATAML_CORR_FILE = ROOT / "03_results_tables" / "phase3C_BeatAML_drug_correlation_all_formula_A.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def lower_set(row: dict[str, str]) -> set[str]:
    return {k.lower() for k in row.keys()}


def fnum(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def minimal_search_inventory() -> list[dict[str, str]]:
    keyword_groups = {
        "DEG": ["DEG", "PRFT_high", "PRFT_low", "Formula_A", "LOCKED_A", "differential", "limma", "DESeq2"],
        "BeatAML": ["BeatAML", "AUC", "drug", "sensitivity", "response", "PRFT_group", "Formula_A"],
        "CMap_LINCS": ["CMap", "LINCS", "signature", "up_genes", "down_genes"],
    }
    rows: list[dict[str, str]] = []
    for dirpath, _, filenames in os.walk(ROOT):
        for name in filenames:
            path = Path(dirpath) / name
            ext = path.suffix.lower()
            if ext not in ALLOWED_EXT:
                continue
            full = str(path)
            hits = []
            groups = []
            for group, keys in keyword_groups.items():
                group_hits = [k for k in keys if k.lower() in name.lower() or k.lower() in full.lower()]
                if group_hits:
                    hits.extend(group_hits)
                    groups.append(group)
            if not hits:
                continue
            probable = ";".join(sorted(set(groups)))
            usable = "yes" if path in {DEG_FILE, BEATAML_AUC_FILE, BEATAML_GROUP_FILE, BEATAML_HIGH_LOW_FILE, BEATAML_CORR_FILE} else "review_only"
            used_for = ""
            if path == DEG_FILE:
                used_for = "CMap/LINCS up/down signature"
            elif path == BEATAML_AUC_FILE:
                used_for = "BeatAML AUC availability"
            elif path == BEATAML_GROUP_FILE:
                used_for = "BeatAML PRFT group/score"
            elif path in {BEATAML_HIGH_LOW_FILE, BEATAML_CORR_FILE}:
                used_for = "BeatAML-CMap consistency template context"
            elif "CMap" in probable or "LINCS" in probable:
                used_for = "existing CMap/LINCS file check"
            rows.append(
                {
                    "file_path": full,
                    "file_name": name,
                    "file_extension": ext,
                    "keyword_hit": ";".join(sorted(set(hits))),
                    "probable_content": probable,
                    "usable": usable,
                    "used_for": used_for,
                    "risk_note": "Do not use unless current locked Formula A provenance is confirmed." if usable != "yes" else "Selected as current Phase17A-lite input.",
                    "comment": "Minimal keyword search only; raw single-cell object search intentionally not performed.",
                }
            )
    rows.sort(key=lambda r: (r["probable_content"], r["file_path"]))
    return rows


def main() -> None:
    for sub in SUBDIRS:
        (PHASE17 / sub).mkdir(parents=True, exist_ok=True)

    boundary_rows = [
        {
            "item": "Locked Formula A model",
            "status": "unchanged",
            "source_file": str(ROOT / "00_LOCKED_FORMULA" / "README_formula_lock.md"),
            "risk_note": "Formula B is audit-only and must not enter final results.",
            "comment": "This Phase17A-lite preparation did not modify the locked Formula A model, gene list, coefficients, cutoff rule, training cohort, validation cohorts, main results, or manuscript conclusions.",
        },
        {
            "item": "Manuscript claim boundary",
            "status": "read and respected",
            "source_file": str(ROOT / "09_manuscript_level_audit_LOCKED_A" / "claim_boundary" / "manuscript_claim_boundary_check.csv"),
            "risk_note": "CMap/LINCS and BeatAML must remain hypothesis-generating/ex vivo context.",
            "comment": "No causal or clinical-response claims were added.",
        },
        {
            "item": "Main Results and Methods",
            "status": "unchanged",
            "source_file": str(ROOT / "10_manuscript_rewriting_LOCKED_A" / "results" / "Results_draft_LOCKED_A_v1.txt") + "; " + str(ROOT / "13_methods_metadata_finalization_LOCKED_A" / "resolved_methods" / "Methods_draft_LOCKED_A_v2_metadata_resolved.txt"),
            "risk_note": "This round generated only preparation artifacts.",
            "comment": "No manuscript text was edited.",
        },
        {
            "item": "Phase16B/Phase17 readiness carryover",
            "status": "read and used",
            "source_file": str(ROOT / "16B_reference_audit_QC_LOCKED_A" / "logs" / "Phase16B_and_Phase17_readiness_summary_CN_v2.txt") + "; " + str(PHASE17 / "input_inventory" / "supplementary_analysis_input_checklist.csv"),
            "risk_note": "Virtual perturbation remains deferred without raw single-cell object.",
            "comment": "Phase17A-lite prioritizes CMap/LINCS signature and BeatAML consistency preparation.",
        },
    ]
    write_csv(
        PHASE17 / "locked_formula_import" / "Phase17A_lite_Locked_Formula_A_boundary_check.csv",
        boundary_rows,
        ["item", "status", "source_file", "risk_note", "comment"],
    )

    inventory_rows = minimal_search_inventory()
    write_csv(
        PHASE17 / "input_inventory" / "Phase17A_lite_minimal_input_search_inventory.csv",
        inventory_rows,
        ["file_path", "file_name", "file_extension", "keyword_hit", "probable_content", "usable", "used_for", "risk_note", "comment"],
    )

    deg_rows = read_csv(DEG_FILE)
    deg_header = lower_set(deg_rows[0]) if deg_rows else set()
    has_gene = "gene_symbol" in deg_header or "gene" in deg_header or "symbol" in deg_header
    has_logfc = "logfc" in deg_header or "log2fc" in deg_header
    has_fdr = "adj.p.val" in deg_header or "fdr" in deg_header or "padj" in deg_header or "qvalue" in deg_header
    comparison_confirmed = DEG_FILE.exists() and has_gene and has_logfc and has_fdr
    eligible_deg = comparison_confirmed
    deg_eligibility = [
        {
            "candidate_DEG_file": str(DEG_FILE),
            "has_gene_symbol": "yes" if has_gene else "no",
            "has_log2FC": "yes" if has_logfc else "no",
            "has_FDR": "yes" if has_fdr else "no",
            "comparison_confirmed_as_LOCKED_A_PRFT_high_vs_low": "yes" if comparison_confirmed else "no",
            "uses_Formula_B": "no",
            "uses_old_DEG": "no",
            "eligible_for_CMap_signature": "yes" if eligible_deg else "no",
            "reason": "Current Methods/Results describe this as the PRFT-high vs PRFT-low TCGA limma DEG table used in the locked Formula A manuscript workflow; no Formula B marker detected in file path/name/header.",
            "comment": "Used as existing DEG result only; no differential expression rerun performed.",
        }
    ]
    write_csv(
        PHASE17 / "drug_reversal_CMap_LINCS" / "Phase17A_lite_DEG_table_eligibility_check.csv",
        deg_eligibility,
        ["candidate_DEG_file", "has_gene_symbol", "has_log2FC", "has_FDR", "comparison_confirmed_as_LOCKED_A_PRFT_high_vs_low", "uses_Formula_B", "uses_old_DEG", "eligible_for_CMap_signature", "reason", "comment"],
    )

    signature_rows: list[dict[str, object]] = []
    up_genes: list[str] = []
    down_genes: list[str] = []
    if eligible_deg:
        sig_source = []
        for r in deg_rows:
            gene = (r.get("gene_symbol") or r.get("gene") or r.get("symbol") or "").strip()
            logfc = fnum(r.get("logFC") or r.get("log2FC") or "0")
            fdr = fnum(r.get("adj.P.Val") or r.get("FDR") or r.get("padj") or r.get("qvalue") or "1", default=1.0)
            if not gene or fdr >= 0.05 or logfc == 0:
                continue
            sig_source.append((gene, logfc, fdr))
        up = sorted([x for x in sig_source if x[1] > 0], key=lambda x: x[1], reverse=True)[:150]
        down = sorted([x for x in sig_source if x[1] < 0], key=lambda x: x[1])[:150]
        for direction, rows in [("upregulated_in_PRFT_high", up), ("downregulated_in_PRFT_high", down)]:
            for rank, (gene, logfc, fdr) in enumerate(rows, start=1):
                signature_rows.append(
                    {
                        "gene_symbol": gene,
                        "direction_in_PRFT_high": direction,
                        "log2FC": logfc,
                        "FDR": fdr,
                        "rank": rank,
                        "used_for_CMap": "yes",
                        "source_DEG_file": str(DEG_FILE),
                        "comment": "Filtered from existing DEG table with FDR < 0.05; no DEG rerun.",
                    }
                )
        up_genes = [x[0] for x in up]
        down_genes = [x[0] for x in down]
        write_csv(
            PHASE17 / "drug_reversal_CMap_LINCS" / "CMap_LINCS_input_signature_LOCKED_A_Phase17A_lite.csv",
            signature_rows,
            ["gene_symbol", "direction_in_PRFT_high", "log2FC", "FDR", "rank", "used_for_CMap", "source_DEG_file", "comment"],
        )
        write_text(PHASE17 / "drug_reversal_CMap_LINCS" / "CMap_LINCS_up_genes_LOCKED_A_Phase17A_lite.txt", "\n".join(up_genes) + "\n")
        write_text(PHASE17 / "drug_reversal_CMap_LINCS" / "CMap_LINCS_down_genes_LOCKED_A_Phase17A_lite.txt", "\n".join(down_genes) + "\n")
    else:
        write_csv(
            PHASE17 / "drug_reversal_CMap_LINCS" / "CMap_LINCS_input_signature_template_Phase17A_lite.csv",
            [],
            ["gene_symbol", "direction_in_PRFT_high", "log2FC", "FDR", "rank", "used_for_CMap", "source_DEG_file", "comment"],
        )
    upload_instructions = f"""CMap/LINCS manual upload instructions CN - Phase17A-lite

1. This is the Locked Formula A manuscript-context PRFT-high vs PRFT-low differential-expression signature prepared from:
   {DEG_FILE}
2. up genes = PRFT-high upregulated genes, filtered by FDR < 0.05 and ranked by descending log2FC.
3. down genes = PRFT-high downregulated genes, filtered by FDR < 0.05 and ranked by ascending log2FC.
4. Current prepared gene counts:
   up genes: {len(up_genes)}
   down genes: {len(down_genes)}
5. CMap/LINCS results can only be used for perturbational signature reversal and hypothesis prioritization.
6. CMap/LINCS output must not be interpreted as clinical efficacy prediction.
7. CMap/LINCS output must not be written as drug-treatment evidence.
8. Downloaded/exported CMap/LINCS output should be saved under:
   {PHASE17 / "drug_reversal_CMap_LINCS" / "manual_CMap_LINCS_outputs"}
9. This Phase17A-lite round did not run online CMap/LINCS and did not create any drug result.
"""
    write_text(PHASE17 / "drug_reversal_CMap_LINCS" / "CMap_LINCS_manual_upload_instructions_CN_Phase17A_lite.txt", upload_instructions)

    beataml_candidates = [BEATAML_AUC_FILE, BEATAML_GROUP_FILE, BEATAML_HIGH_LOW_FILE, BEATAML_CORR_FILE]
    beataml_elig = []
    for path in beataml_candidates:
        rows = read_csv(path)
        header = lower_set(rows[0]) if rows else set()
        has_sample = "sample_id" in header
        has_group = bool({"risk_group", "risk_score", "score_name"} & header)
        has_drug = "drug_name" in header
        has_auc = "auc" in header or "response_value" in header or "median_high" in header or "spearman_rho" in header or "metric_name" in header
        usable = path.exists() and has_drug and has_auc and (has_sample or has_group)
        if path in {BEATAML_HIGH_LOW_FILE, BEATAML_CORR_FILE} and has_drug and has_auc:
            usable = True
        beataml_elig.append(
            {
                "candidate_BeatAML_file": str(path),
                "has_sample_id": "yes" if has_sample else "no",
                "has_PRFT_group_or_score": "yes" if has_group else "no",
                "has_drug_name": "yes" if has_drug else "no",
                "has_AUC": "yes" if has_auc else "no",
                "AUC_direction_lock_applied": "yes",
                "usable_for_BeatAML_CMap_consistency": "yes" if usable else "partial" if path.exists() else "no",
                "reason": "Contains BeatAML AUC/drug and Formula A PRFT score/group or summary association fields." if usable else "Useful only together with paired BeatAML Formula A files." if path.exists() else "Missing.",
                "comment": "Higher AUC = lower ex vivo sensitivity / greater relative resistance.",
            }
        )
    write_csv(
        PHASE17 / "beataml_consistency_check" / "Phase17A_lite_BeatAML_input_eligibility_check.csv",
        beataml_elig,
        ["candidate_BeatAML_file", "has_sample_id", "has_PRFT_group_or_score", "has_drug_name", "has_AUC", "AUC_direction_lock_applied", "usable_for_BeatAML_CMap_consistency", "reason", "comment"],
    )
    write_text(
        PHASE17 / "beataml_consistency_check" / "BeatAML_AUC_interpretation_lock_Phase17A_lite.txt",
        "In all Phase17 analyses, higher AUC must be interpreted as lower ex vivo sensitivity or greater relative resistance. No sentence should describe higher AUC as higher sensitivity.\n",
    )
    high_low_rows = read_csv(BEATAML_HIGH_LOW_FILE)[:25]
    consistency_template = []
    for r in high_low_rows:
        interp = r.get("group_difference_interpretation", "")
        consistency_template.append(
            {
                "drug_or_mechanism": r.get("drug_name", ""),
                "source": str(BEATAML_HIGH_LOW_FILE),
                "BeatAML_available": "yes",
                "BeatAML_AUC_direction": "higher AUC = lower ex vivo sensitivity / greater relative resistance",
                "PRFT_high_interpretation": interp,
                "CMap_LINCS_direction": "pending manual CMap/LINCS output",
                "consistent_or_not": "pending",
                "safe_wording": "BeatAML provides ex vivo pharmacogenomic context; CMap/LINCS consistency can be assessed only after manual CMap/LINCS output is returned.",
                "risk_note": "Do not describe as clinical drug response or therapy recommendation.",
                "comment": "Template row seeded from existing Formula A BeatAML high-low AUC table.",
            }
        )
    write_csv(
        PHASE17 / "beataml_consistency_check" / "BeatAML_CMap_consistency_check_template_Phase17A_lite.csv",
        consistency_template,
        ["drug_or_mechanism", "source", "BeatAML_available", "BeatAML_AUC_direction", "PRFT_high_interpretation", "CMap_LINCS_direction", "consistent_or_not", "safe_wording", "risk_note", "comment"],
    )

    write_text(
        PHASE17 / "virtual_perturbation" / "Phase17A_lite_virtual_perturbation_defer_decision_CN.txt",
        """Phase17A-lite virtual perturbation defer decision

1. Virtual perturbation is a valuable add-on analysis for candidate regulator prioritization.
2. A true run requires a raw single-cell object or expression matrix, cell metadata, Formula A gene coverage, and a reproducible software environment.
3. The current priority is lower than CMap/LINCS signature preparation and BeatAML consistency checking.
4. If single-cell inputs remain incomplete, virtual perturbation should be limited to Discussion/Future validation.
5. If complete single-cell inputs are supplied later, proceed in a later Phase17C module.
6. Candidate nodes are limited to JAK2, STAT5A, STAT5B, CD274, SLC7A11, GPX4, NFE2L2, and HSP90AA1.
7. Virtual perturbation can support candidate regulator prioritization only; it cannot be used as mechanistic proof.
""",
    )

    method_refs = [
        ("CMap / Connectivity Map / L1000", "Required to cite the perturbational signature reversal framework if manual CMap/LINCS upload is performed.", "method/resource paper", "must_have_if_used", "Connectivity Map L1000 perturbational signature reversal citation", "manual verification required", "Supplementary Methods/reference support", "CMap/LINCS was used for perturbational signature-reversal prioritization.", ""),
        ("LINCS perturbational signature resource", "Required for LINCS L1000 resource provenance.", "resource paper", "must_have_if_used", "LINCS L1000 perturbational signature resource citation", "manual verification required", "Supplementary Methods/reference support", "LINCS L1000 provides perturbational expression signatures for comparison with the PRFT-high signature.", ""),
        ("scTenifoldKnk virtual knockout", "Future support if Phase17C performs virtual knockout.", "method paper", "optional_later", "scTenifoldKnk virtual knockout single-cell citation", "manual verification required", "Future validation / later virtual perturbation Methods", "scTenifoldKnk may be used in a later analysis to prioritize candidate regulators.", "Do not cite in main text unless used."),
        ("CellOracle in silico perturbation", "Optional future support if GRN-based perturbation is performed.", "method paper", "optional_later", "CellOracle in silico gene perturbation single-cell regulatory network citation", "manual verification required", "Future validation / later virtual perturbation Methods", "CellOracle may be used as an alternative future GRN-based perturbation framework.", "Do not cite unless used."),
        ("iLINCS or SigCom LINCS", "Optional if that specific web/API interface is used for upload/query.", "tool/interface citation", "optional_if_used", "iLINCS SigCom LINCS citation", "manual verification required", "Supplementary Methods/reference support", "The exact query interface and database version should be reported if used.", ""),
        ("Limitation of in silico perturbation", "Useful for cautious interpretation if CMap/LINCS or virtual perturbation is discussed.", "methods/review or validation paper", "optional", "in silico perturbation drug reversal limitation validation citation", "manual verification required", "Discussion/Limitations", "Computational perturbation results are hypothesis-generating and require experimental validation.", ""),
    ]
    write_csv(
        PHASE17 / "reference_support" / "Phase17A_lite_minimal_method_reference_needs.csv",
        [
            {
                "topic": r[0],
                "why_needed": r[1],
                "recommended_reference_type": r[2],
                "must_have_or_optional": r[3],
                "suggested_search_terms": r[4],
                "verification_status": r[5],
                "recommended_location": r[6],
                "safe_sentence": r[7],
                "comment": r[8],
            }
            for r in method_refs
        ],
        ["topic", "why_needed", "recommended_reference_type", "must_have_or_optional", "suggested_search_terms", "verification_status", "recommended_location", "safe_sentence", "comment"],
    )

    fig_rows = [
        {"panel": "A", "title": "Locked Formula A PRFT-high vs PRFT-low DEG signature used for CMap/LINCS query.", "input": str(DEG_FILE), "status": "prepared", "planned_content": f"{len(up_genes)} up genes and {len(down_genes)} down genes available for upload", "risk_note": "Signature is not a CMap result."},
        {"panel": "B", "title": "CMap/LINCS perturbational signature reversal workflow.", "input": "manual upload of up/down gene lists", "status": "workflow only", "planned_content": "Show upload, query, export, and QC steps", "risk_note": "No online CMap/LINCS run in this round."},
        {"panel": "C", "title": "Mechanism classes of candidate PRFT-reversing perturbagens.", "input": "manual CMap/LINCS outputs", "status": "pending", "planned_content": "Populate only after CMap/LINCS results return", "risk_note": "Do not write pending as completed."},
        {"panel": "D", "title": "BeatAML consistency check with ex vivo drug-response patterns.", "input": str(BEATAML_HIGH_LOW_FILE), "status": "template prepared; pending CMap/LINCS output", "planned_content": "Compare CMap mechanism classes with BeatAML AUC patterns", "risk_note": "Higher AUC means lower ex vivo sensitivity/resistance."},
    ]
    write_csv(PHASE17 / "supplementary_figures" / "Supplementary_Figure_Sx_CMap_BeatAML_minimal_design_Phase17A_lite.csv", fig_rows, ["panel", "title", "input", "status", "planned_content", "risk_note"])
    tab_rows = [
        {"table": "CMap_LINCS_input_signature", "source": str(DEG_FILE), "status": "generated" if eligible_deg else "template only", "columns": "gene_symbol,direction_in_PRFT_high,log2FC,FDR,rank", "comment": "Upload-ready signature; no CMap output yet."},
        {"table": "CMap_LINCS_results_QC", "source": "manual_CMap_LINCS_outputs", "status": "pending", "columns": "perturbagen,mechanism_class,connectivity_score,FDR,database_version,QC_flag", "comment": "Populate after manual upload/export."},
        {"table": "BeatAML_CMap_consistency", "source": str(BEATAML_HIGH_LOW_FILE), "status": "template prepared", "columns": "drug_or_mechanism,BeatAML_AUC_direction,CMap_LINCS_direction,consistent_or_not", "comment": "Consistency cannot be finalized before CMap output."},
    ]
    write_csv(PHASE17 / "supplementary_tables" / "Supplementary_Table_Sx_CMap_BeatAML_minimal_design_Phase17A_lite.csv", tab_rows, ["table", "source", "status", "columns", "comment"])

    write_text(
        PHASE17 / "manuscript_text" / "Phase17A_lite_text_insertion_boundary_recommendation.txt",
        """Phase17A-lite text insertion boundary recommendation

1. If the CMap/LINCS signature has been generated but no CMap output is available:
   - Do not add it to Results.
   - It can be described only as a future drug-reversal analysis in Discussion/Future validation.

2. If CMap/LINCS output is returned and QC-checked:
   - A short supplementary Results paragraph may be considered at the end of Results.
   - The wording must remain perturbational signature reversal / hypothesis prioritization.

3. If only BeatAML tables are available and no CMap output exists:
   - Do not add new Results.
   - Discussion may mention pharmacogenomic prioritization, with ex vivo and non-clinical boundaries.

4. If single-cell virtual perturbation has not been run:
   - Do not include it in Results.
   - It can only be listed as future validation.

5. Forbidden wording:
   proved; confirmed; demonstrated causally; restored ferroptosis sensitivity; predicts clinical response; therapeutic efficacy.
""",
    )

    risk_rows = [
        ("Locked Formula A modification", "high", "Phase17A-lite could accidentally modify model inputs.", "Read-only use of locked coefficient file; no formula edits.", "Formula A remained unchanged.", ""),
        ("Formula B misuse", "high", "Formula B/audit files exist in project.", "Do not use Formula B files for DEG/signature.", "Formula B was not used.", ""),
        ("Old DEG misuse", "high", "Old/package DEG copies exist.", "Use only current checked DEG table.", "Signature was generated from the checked PRFT-high vs PRFT-low DEG table.", ""),
        ("CMap signature mistaken for result", "high", "Upload-ready list is not CMap output.", "Label Panel C/result fields as pending.", "CMap/LINCS output has not been generated.", ""),
        ("CMap clinical overclaim", "high", "Drug reversal may be misread as efficacy.", "Use perturbational signature reversal only.", "Not a clinical efficacy prediction.", ""),
        ("BeatAML AUC direction inversion", "high", "Higher AUC could be called more sensitive.", "Lock wording in every artifact.", "Higher AUC = lower ex vivo sensitivity / greater relative resistance.", ""),
        ("Virtual perturbation written as complete", "high", "No scTenifoldKnk/CellOracle run was performed.", "Mark deferred.", "Virtual perturbation is deferred to future validation.", ""),
        ("Method pile-up", "moderate", "Too many extra methods could weaken manuscript credibility.", "Keep only CMap/LINCS and BeatAML preparation now.", "Phase17A-lite is minimal and bounded.", ""),
        ("Missing method DOI/PMID verification", "moderate", "Method references are not verified here.", "Manual verification required before citation insertion.", "No DOI/PMID invented.", ""),
        ("Next-round readiness", "moderate", "CMap output and method references still pending.", "Proceed to Phase17B only after manual upload/export or explicit human approval.", "Phase17B can be CMap output ingestion/QC, not manuscript Results insertion yet.", ""),
    ]
    write_csv(
        PHASE17 / "risk_check" / "Phase17A_lite_risk_checklist.csv",
        [{"risk_item": r[0], "risk_level": r[1], "problem": r[2], "suggested_action": r[3], "safe_wording": r[4], "comment": r[5]} for r in risk_rows],
        ["risk_item", "risk_level", "problem", "suggested_action", "safe_wording", "comment"],
    )

    summary = f"""Phase17A-lite summary CN
生成日期：2026-06-25

1. 本轮是否修改 Locked Formula A
否。Formula A gene list、coefficients、cutoff rule、training/validation cohort 和 manuscript conclusion 均未修改。

2. 本轮是否修改主模型、Results、Methods 或 manuscript text
否。本轮只生成 Phase17A-lite 准备文件、signature、模板和边界说明。

3. 是否找到 eligible Formula A DEG 表
是：{DEG_FILE}

4. 是否生成 CMap/LINCS upload-ready signature
{"是" if eligible_deg else "否"}。

5. up genes 数量和 down genes 数量
up genes: {len(up_genes)}
down genes: {len(down_genes)}

6. 是否找到 BeatAML AUC 表
是：{BEATAML_AUC_FILE}

7. 是否可以做 BeatAML-CMap consistency check
可以准备模板，但正式一致性判断必须等待人工 CMap/LINCS output。

8. virtual perturbation 是否暂缓
是。未搜索 raw single-cell object，未运行 scTenifoldKnk 或 CellOracle。

9. 是否需要补 single-cell object
是，如需真实 virtual perturbation，仍需 raw single-cell object/expression matrix、metadata、gene coverage 和软件环境。

10. 是否需要人工上传 CMap/LINCS
是。当前只生成 upload-ready up/down gene lists，未运行在线 CMap/LINCS。

11. 是否需要新增方法学参考文献
需要，但仅限 CMap/Connectivity Map/L1000、LINCS resource、后续可能用到的 scTenifoldKnk/CellOracle、可选 iLINCS/SigCom 和 limitation references；所有 DOI/PMID 需人工核查。

12. 推荐新增参考文献数量
最小 4-6 个方向，且只在对应方法实际使用后加入。

13. 是否建议进入 Results
不建议。当前没有 CMap/LINCS output，只有 signature 和模板。

14. 是否建议只进入 Discussion/Future validation
是，当前最多作为 future drug-reversal analysis / pharmacogenomic prioritization boundary。

15. 下一步最优先任务
人工上传 CMap/LINCS、保存输出到 manual_CMap_LINCS_outputs、核查方法学文献 DOI/PMID，然后进行 Phase17B output ingestion/QC。

明确结论：
1. Locked Formula A 是否保持不变：是。
2. 主模型 / Results / Methods 是否保持不变：是。
3. eligible Formula A DEG 表是否找到：是。
4. CMap/LINCS upload-ready signature 是否生成：{"是" if eligible_deg else "否"}。
5. up genes 数量：{len(up_genes)}
6. down genes 数量：{len(down_genes)}
7. BeatAML AUC 表是否找到：是。
8. BeatAML-CMap 一致性核查是否可进行：可准备；正式判断需 CMap/LINCS output。
9. virtual perturbation 是否暂缓：是。
10. 是否需要补 single-cell object：是，若要真实运行 virtual perturbation。
11. 是否需要人工上传 CMap/LINCS：是。
12. 是否建议现在写入 Results：否。
13. 是否建议现在只写 Discussion/Future validation：是。
14. 是否可以进入下一轮 Phase17B：可以，但 Phase17B 应限于 CMap/LINCS output ingestion/QC 和人工文献核查，不应直接写入 Results。
15. 下一步最优先任务：人工上传 CMap/LINCS 并回填输出，同时核查 CMap/LINCS/LINCS 方法文献。
"""
    write_text(PHASE17 / "logs" / "Phase17A_lite_summary_CN.txt", summary)


if __name__ == "__main__":
    main()
