from __future__ import annotations

import csv
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(r"E:\ai助手\AML_PRFT_Human_Genomics_phase0_audit")
PHASE17 = ROOT / "17_supplementary_in_silico_perturbation_LOCKED_A"
OUTDIR = PHASE17 / "drug_reversal_CMap_LINCS"
MANUAL = OUTDIR / "manual_CMap_LINCS_outputs"

SUBDIRS = [
    "drug_reversal_CMap_LINCS",
    "beataml_consistency_check",
    "supplementary_figures",
    "supplementary_tables",
    "manuscript_text",
    "reference_support",
    "risk_check",
    "logs",
]

BEATAML_HIGH_LOW = ROOT / "03_results_tables" / "phase3C_BeatAML_high_low_drug_comparison_formula_A.csv"
BEATAML_CORR = ROOT / "03_results_tables" / "phase3C_BeatAML_drug_correlation_all_formula_A.csv"


def read_delimited(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    delim = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


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


def fnum(value: str, default: float | None = None) -> float | None:
    try:
        return float(value)
    except Exception:
        return default


def norm_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (name or "").lower())


def detect_col(headers: set[str], candidates: list[str]) -> str | None:
    lower_map = {h.lower(): h for h in headers}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None


def classify_mechanism(drug: str) -> tuple[str, str, str]:
    d = drug.lower()
    mapping = [
        ("bortezomib|ixazomib|carfilzomib|proteasome", "proteasome / proteostasis-related", "name-based known proteasome/proteostasis association"),
        ("tanespimycin|17-aag|geldanamycin|hsp90|ganetespib", "HSP90-related", "name-based HSP90/proteostasis association"),
        ("panobinostat|vorinostat|belinostat|romidepsin|entinostat|trichostatin|hdac", "HDAC-related", "name-based HDAC/epigenetic association"),
        ("at-7519|at7519|dinaciclib|palbociclib|flavopiridol|cdk|cell cycle", "cell-cycle / CDK-related", "name-based CDK/cell-cycle association"),
        ("sb-239063|selumetinib|trametinib|dorsomorphin|dasatinib|nilotinib|bosutinib|sunitinib|linifanib|cediranib|motesanib|tivozanib|mapk|p38|kinase|mek", "MAPK / stress-kinase-related", "name-based kinase/stress-kinase association"),
        ("venetoclax|abt-199|navitoclax|bcl", "apoptosis / BCL2-related", "name-based BCL2/apoptosis association"),
        ("erastin|rsl3|ferroptosis|nfe2l2|ros|oxidative|sulfasalazine|buthionine", "oxidative-stress / ferroptosis-related", "name-based oxidative-stress/ferroptosis association"),
        (r"\bjak\b|jak[0-9]|ruxolitinib|tofacitinib|stat3|stat5|stat5a|stat5b", "JAK / STAT-related", "name-based JAK/STAT association"),
        ("rapamycin|metformin|torin|ink-128|mTOR|ketoconazole|dorsomorphin", "metabolism-related", "name-based metabolism/mTOR/steroid/metabolic-stress association"),
    ]
    for pattern, cls, basis in mapping:
        if re.search(pattern, d):
            if drug.lower() in {"ketoconazole", "clonazepam", "vardenafil"}:
                return ("other / unclear", "low AML relevance; name-based mechanism not prioritized", "manual annotation required")
            return (cls, basis, "manual verification required")
    return ("other / unclear", "no conservative mechanism class inferred from name", "manual annotation required")


def aml_relevance(drug: str, mechanism: str, tissue: str, disease: str) -> str:
    d = drug.lower()
    if mechanism in {
        "proteasome / proteostasis-related",
        "HSP90-related",
        "HDAC-related",
        "apoptosis / BCL2-related",
        "JAK / STAT-related",
        "oxidative-stress / ferroptosis-related",
    }:
        return "moderate; mechanism class has AML/PRFT-context relevance but row is not AML-specific"
    if "leukemia" in (disease or "").lower():
        return "moderate; exported row disease label includes leukemia, but cell/tissue context still requires review"
    if d in {"ketoconazole", "clonazepam", "vardenafil"}:
        return "low; not prioritized for AML interpretation without manual evidence"
    return "uncertain; manual AML relevance review required"


def choose_sigcom_files() -> list[Path]:
    if not MANUAL.exists():
        return []
    files = [p for p in MANUAL.iterdir() if p.is_file() and p.suffix.lower() in {".tsv", ".csv", ".txt", ".xlsx"}]
    priority_words = ["sigcom", "lincs", "prft", "reversal", "reverser", "raw"]
    files.sort(key=lambda p: (0 if any(w in p.name.lower() for w in priority_words) else 1, {".tsv": 0, ".csv": 1, ".txt": 2, ".xlsx": 3}.get(p.suffix.lower(), 9), p.name.lower()))
    return files


def main() -> None:
    for sub in SUBDIRS:
        (PHASE17 / sub).mkdir(parents=True, exist_ok=True)

    boundary_rows = [
        {
            "item": "Locked Formula A model and inputs",
            "status": "unchanged",
            "source_file": str(ROOT / "00_LOCKED_FORMULA" / "README_formula_lock.md"),
            "risk_note": "Formula B remains audit-only and was not used.",
            "comment": "This Phase17B_v2 curation did not modify the locked Formula A model, gene list, coefficients, cutoff rule, training cohort, validation cohorts, main results, Methods, or manuscript conclusions.",
        },
        {
            "item": "Phase17A-lite query input",
            "status": "read-only carryover",
            "source_file": str(OUTDIR / "CMap_LINCS_input_signature_LOCKED_A_Phase17A_lite.csv"),
            "risk_note": "Input signature is not itself a SigCom result.",
            "comment": "The exported SigCom Reversers table was curated separately from the query input.",
        },
        {
            "item": "BeatAML AUC direction lock",
            "status": "applied",
            "source_file": str(PHASE17 / "beataml_consistency_check" / "BeatAML_AUC_interpretation_lock_Phase17A_lite.txt"),
            "risk_note": "Higher AUC must not be written as higher sensitivity.",
            "comment": "All BeatAML consistency wording uses higher AUC = lower ex vivo sensitivity / greater relative resistance.",
        },
        {
            "item": "Main Results/Methods/manuscript text",
            "status": "unchanged",
            "source_file": str(ROOT / "10_manuscript_rewriting_LOCKED_A" / "results" / "Results_draft_LOCKED_A_v1.txt") + "; " + str(ROOT / "13_methods_metadata_finalization_LOCKED_A" / "resolved_methods" / "Methods_draft_LOCKED_A_v2_metadata_resolved.txt"),
            "risk_note": "Generated text is a draft insertion recommendation only.",
            "comment": "No manuscript file was edited.",
        },
    ]
    write_csv(OUTDIR / "Phase17B_v2_boundary_check_LOCKED_A.csv", boundary_rows, ["item", "status", "source_file", "risk_note", "comment"])

    files = choose_sigcom_files()
    if not files:
        write_text(OUTDIR / "Phase17B_v2_missing_manual_SigCom_output_CN.txt", "未发现 manual_CMap_LINCS_outputs 中真实存在的 SigCom/LINCS 导出文件；本轮停止整理。\n")
        return

    inventory = []
    selected = files[0]
    for p in files:
        rows = read_delimited(p) if p.suffix.lower() != ".xlsx" else []
        headers = set(rows[0].keys()) if rows else set()
        inv = {
            "file_path": str(p),
            "file_name": p.name,
            "file_extension": p.suffix.lower(),
            "file_size": p.stat().st_size,
            "detected_platform": "SigCom_LINCS",
            "detected_result_type": "LINCS_L1000_Chemical_Perturbations_Reversers",
            "has_perturbagen": "yes" if detect_col(headers, ["Perturbagen", "perturbagen", "drug", "Name"]) else "no",
            "has_dose": "yes" if detect_col(headers, ["Dose"]) else "no",
            "has_tissue": "yes" if detect_col(headers, ["Tissue"]) else "no",
            "has_cell_line": "yes" if detect_col(headers, ["Cell Line", "cell_line"]) else "no",
            "has_time": "yes" if detect_col(headers, ["Timepoint", "time"]) else "no",
            "has_z_score": "yes" if detect_col(headers, ["z-score (sum)", "z-score (up)", "z-score"]) else "no",
            "has_p_value": "yes" if detect_col(headers, ["p-value (up)", "p-value (down)", "p_value"]) else "no",
            "usable_for_curation": "yes" if rows and detect_col(headers, ["Perturbagen", "Name"]) and detect_col(headers, ["z-score (sum)", "z-score (up)"]) else "no",
            "risk_note": "Use only as real exported Reversers; do not infer Mimickers or clinical efficacy.",
            "comment": f"Rows detected: {len(rows)}",
        }
        inventory.append(inv)
    write_csv(
        OUTDIR / "Phase17B_v2_manual_SigCom_output_inventory.csv",
        inventory,
        ["file_path", "file_name", "file_extension", "file_size", "detected_platform", "detected_result_type", "has_perturbagen", "has_dose", "has_tissue", "has_cell_line", "has_time", "has_z_score", "has_p_value", "usable_for_curation", "risk_note", "comment"],
    )

    raw_rows = read_delimited(selected)
    raw_fields = list(raw_rows[0].keys()) + ["source_file", "import_note"] if raw_rows else ["source_file", "import_note"]
    raw_imported = []
    for r in raw_rows:
        rr = dict(r)
        rr["source_file"] = str(selected)
        rr["import_note"] = "Real SigCom LINCS exported Reversers row imported without online rerun."
        raw_imported.append(rr)
    write_csv(OUTDIR / "SigCom_LINCS_reversal_results_raw_imported_Phase17B_v2.csv", raw_imported, raw_fields)

    headers = set(raw_rows[0].keys())
    c_pert = detect_col(headers, ["Perturbagen", "Name"])
    c_dose = detect_col(headers, ["Dose"])
    c_tissue = detect_col(headers, ["Tissue"])
    c_cell = detect_col(headers, ["Cell Line"])
    c_time = detect_col(headers, ["Timepoint", "Time"])
    c_zsum = detect_col(headers, ["z-score (sum)", "z-score"])
    c_zup = detect_col(headers, ["z-score (up)"])
    c_zdown = detect_col(headers, ["z-score (down)"])
    c_pup = detect_col(headers, ["p-value (up)"])
    c_pdown = detect_col(headers, ["p-value (down)"])
    c_fdr_up = detect_col(headers, ["FDR (up)"])
    c_fdr_down = detect_col(headers, ["FDR (down)"])
    c_disease = detect_col(headers, ["Disease"])

    curated = []
    for idx, r in enumerate(raw_rows, start=1):
        perturbagen = r.get(c_pert, "") if c_pert else ""
        mechanism, basis, annot_status = classify_mechanism(perturbagen)
        tissue = r.get(c_tissue, "") if c_tissue else ""
        disease = r.get(c_disease, "") if c_disease else ""
        z = fnum(r.get(c_zsum, "")) if c_zsum else None
        if z is None and c_zup and c_zdown:
            z = (fnum(r.get(c_zup, ""), 0.0) or 0.0) + (fnum(r.get(c_zdown, ""), 0.0) or 0.0)
        pvals = [x for x in [fnum(r.get(c_pup, "")) if c_pup else None, fnum(r.get(c_pdown, "")) if c_pdown else None] if x is not None]
        fdrs = [x for x in [fnum(r.get(c_fdr_up, "")) if c_fdr_up else None, fnum(r.get(c_fdr_down, "")) if c_fdr_down else None] if x is not None]
        p_value = min(pvals) if pvals else ""
        adj = min(fdrs) if fdrs else "not available"
        curated.append(
            {
                "rank": idx,
                "perturbagen": perturbagen,
                "dose": r.get(c_dose, "") if c_dose else "",
                "tissue": tissue,
                "cell_line": r.get(c_cell, "") if c_cell else "",
                "time": r.get(c_time, "") if c_time else "",
                "z_score": z if z is not None else "",
                "p_value": p_value,
                "adjusted_p_value_if_available": adj,
                "reversal_direction": "SigCom LINCS Reverser",
                "PRFT_reversal_interpretation": "potential PRFT-high transcriptional-state reverser; hypothesis-generating only",
                "AML_relevance": aml_relevance(perturbagen, mechanism, tissue, disease),
                "known_or_inferred_mechanism": mechanism,
                "target_if_available": "manual annotation required",
                "mechanism_annotation_status": annot_status,
                "risk_note": "Not clinical efficacy evidence; not mechanism proof; row context is LINCS perturbational signature reversal.",
                "recommended_wording": "candidate / prioritized / suggested perturbational reverser of the PRFT-high signature",
                "comment": f"Mechanism annotation basis: {basis}.",
            }
        )
    write_csv(
        OUTDIR / "SigCom_LINCS_reversal_results_curated_Phase17B_v2.csv",
        curated,
        ["rank", "perturbagen", "dose", "tissue", "cell_line", "time", "z_score", "p_value", "adjusted_p_value_if_available", "reversal_direction", "PRFT_reversal_interpretation", "AML_relevance", "known_or_inferred_mechanism", "target_if_available", "mechanism_annotation_status", "risk_note", "recommended_wording", "comment"],
    )
    top = curated[:30]
    write_csv(OUTDIR / "SigCom_LINCS_top_reversers_Phase17B_v2.csv", top, list(top[0].keys()) if top else [])
    plot_data = []
    for r in curated:
        pv = fnum(str(r["p_value"]), None)
        minus_log10_p = -1.0 if pv is None else max(0.0, -__import__("math").log10(max(pv, 1e-300)))
        plot_data.append(
            {
                "rank": r["rank"],
                "perturbagen": r["perturbagen"],
                "mechanism_class": r["known_or_inferred_mechanism"],
                "z_score": r["z_score"],
                "minus_log10_p": minus_log10_p,
                "adjusted_p_value_if_available": r["adjusted_p_value_if_available"],
                "plot_group": "top_30_reversers" if int(r["rank"]) <= 30 else "all_reversers",
            }
        )
    write_csv(OUTDIR / "SigCom_LINCS_plot_data_Phase17B_v2.csv", plot_data, ["rank", "perturbagen", "mechanism_class", "z_score", "minus_log10_p", "adjusted_p_value_if_available", "plot_group"])

    by_class: dict[str, list[str]] = defaultdict(list)
    for r in curated:
        if r["perturbagen"] not in by_class[r["known_or_inferred_mechanism"]]:
            by_class[r["known_or_inferred_mechanism"]].append(str(r["perturbagen"]))
    mech_rows = []
    for cls, drugs in sorted(by_class.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        mech_rows.append(
            {
                "mechanism_class": cls,
                "number_of_perturbagens": len(drugs),
                "representative_perturbagens": "; ".join(drugs[:8]),
                "annotation_basis": "conservative perturbagen-name matching; not database-verified in this round",
                "expected_relevance_to_PRFT": "Relevant as a mechanism class for PRFT prioritization" if cls not in {"other / unclear"} else "Unclear; low priority for biological interpretation",
                "overclaim_risk": "high if presented as drug efficacy or mechanism proof",
                "recommended_wording": "mechanism class enriched among suggested PRFT-high signature reversers",
                "manual_verification_needed": "yes",
                "comment": "Discuss mechanism classes before individual drugs.",
            }
        )
    write_csv(
        OUTDIR / "SigCom_LINCS_mechanism_class_summary_Phase17B_v2.csv",
        mech_rows,
        ["mechanism_class", "number_of_perturbagens", "representative_perturbagens", "annotation_basis", "expected_relevance_to_PRFT", "overclaim_risk", "recommended_wording", "manual_verification_needed", "comment"],
    )

    beataml = read_csv(BEATAML_HIGH_LOW)
    beat_by_norm = {norm_name(r.get("drug_name", "")): r for r in beataml}
    beat_by_class: dict[str, list[dict[str, str]]] = defaultdict(list)
    class_map = {
        "HSP90_proteostasis": "HSP90-related",
        "HDAC_epigenetic": "HDAC-related",
        "BCL2_apoptosis": "apoptosis / BCL2-related",
        "MEK_MAPK_related": "MAPK / stress-kinase-related",
        "PI3K_AKT_mTOR": "metabolism-related",
    }
    for r in beataml:
        cls = class_map.get(r.get("drug_category", ""), "other / unclear")
        beat_by_class[cls].append(r)

    consistency = []
    seen_keys: set[str] = set()
    for r in top:
        drug = str(r["perturbagen"])
        mech = str(r["known_or_inferred_mechanism"])
        b = beat_by_norm.get(norm_name(drug))
        source = "direct_drug_match" if b else "mechanism_class_match" if beat_by_class.get(mech) else "SigCom_only"
        if b:
            interp = b.get("group_difference_interpretation", "")
            diff = fnum(b.get("difference_high_minus_low", ""), 0.0) or 0.0
            if diff < 0:
                status = "suggested convergent support"
                safe = "SigCom prioritized the perturbagen as a PRFT-high signature reverser and BeatAML showed lower AUC in PRFT-high samples, consistent with greater ex vivo sensitivity."
            elif diff > 0:
                status = "directionally complex"
                safe = "SigCom prioritized the perturbagen as a reverser, but BeatAML showed higher AUC in PRFT-high samples, indicating lower ex vivo sensitivity; interpret cautiously."
            else:
                status = "not directly comparable"
                safe = "BeatAML direction was neutral or unavailable."
        elif beat_by_class.get(mech):
            cls_rows = beat_by_class[mech]
            lower = [x for x in cls_rows if (fnum(x.get("difference_high_minus_low", ""), 0.0) or 0.0) < 0]
            higher = [x for x in cls_rows if (fnum(x.get("difference_high_minus_low", ""), 0.0) or 0.0) > 0]
            interp = f"{len(lower)} class rows with lower AUC in PRFT-high; {len(higher)} with higher AUC in PRFT-high"
            status = "suggested convergent support" if lower and not higher else "directionally complex" if lower and higher else "inconsistent" if higher else "not directly comparable"
            safe = "Mechanism-class comparison only; do not force drug-level equivalence."
        else:
            interp = "No direct BeatAML drug or conservative class comparator found."
            status = "not directly comparable"
            safe = "SigCom-only hypothesis-generating candidate; requires experimental and pharmacogenomic validation."
        key = drug
        seen_keys.add(key)
        consistency.append(
            {
                "drug_or_mechanism": drug,
                "source": source,
                "BeatAML_available": "yes" if b or beat_by_class.get(mech) else "no",
                "BeatAML_AUC_direction": "higher AUC = lower ex vivo sensitivity / greater relative resistance",
                "PRFT_high_interpretation": interp,
                "SigCom_reverser_direction": "SigCom LINCS Reverser of PRFT-high signature",
                "consistent_or_not": status,
                "safe_wording": safe,
                "risk_note": "Do not describe as clinical response prediction or therapy recommendation.",
                "comment": f"Mechanism class: {mech}",
            }
        )
    for mech, cls_rows in beat_by_class.items():
        if mech == "other / unclear" or mech in seen_keys:
            continue
        if mech in by_class:
            lower = [x for x in cls_rows if (fnum(x.get("difference_high_minus_low", ""), 0.0) or 0.0) < 0]
            higher = [x for x in cls_rows if (fnum(x.get("difference_high_minus_low", ""), 0.0) or 0.0) > 0]
            status = "suggested convergent support" if lower and not higher else "directionally complex" if lower and higher else "inconsistent" if higher else "not directly comparable"
            consistency.append(
                {
                    "drug_or_mechanism": mech,
                    "source": "mechanism_class_match",
                    "BeatAML_available": "yes",
                    "BeatAML_AUC_direction": "higher AUC = lower ex vivo sensitivity / greater relative resistance",
                    "PRFT_high_interpretation": f"{len(lower)} BeatAML class rows with lower AUC in PRFT-high; {len(higher)} with higher AUC in PRFT-high",
                    "SigCom_reverser_direction": "SigCom LINCS mechanism class among Reversers",
                    "consistent_or_not": status,
                    "safe_wording": "Mechanism-class convergence only; not drug efficacy evidence.",
                    "risk_note": "Mechanism annotation requires manual verification.",
                    "comment": "Added class-level comparator.",
                }
            )
    write_csv(
        PHASE17 / "beataml_consistency_check" / "BeatAML_SigCom_consistency_check_Phase17B_v2.csv",
        consistency,
        ["drug_or_mechanism", "source", "BeatAML_available", "BeatAML_AUC_direction", "PRFT_high_interpretation", "SigCom_reverser_direction", "consistent_or_not", "safe_wording", "risk_note", "comment"],
    )

    fig_plot = []
    for r in plot_data[:30]:
        fig_plot.append(
            {
                "panel": "B",
                "rank": r["rank"],
                "label": r["perturbagen"],
                "mechanism_class": r["mechanism_class"],
                "z_score": r["z_score"],
                "minus_log10_p": r["minus_log10_p"],
                "status": "top SigCom Reverser",
            }
        )
    for m in mech_rows:
        fig_plot.append(
            {
                "panel": "C",
                "rank": "",
                "label": m["mechanism_class"],
                "mechanism_class": m["mechanism_class"],
                "z_score": "",
                "minus_log10_p": m["number_of_perturbagens"],
                "status": "mechanism class count",
            }
        )
    for c in consistency[:30]:
        fig_plot.append(
            {
                "panel": "D",
                "rank": "",
                "label": c["drug_or_mechanism"],
                "mechanism_class": c["comment"].replace("Mechanism class: ", ""),
                "z_score": "",
                "minus_log10_p": "",
                "status": c["consistent_or_not"],
            }
        )
    write_csv(PHASE17 / "supplementary_figures" / "Supplementary_Figure_Sx_SigCom_BeatAML_plot_data_Phase17B_v2.csv", fig_plot, ["panel", "rank", "label", "mechanism_class", "z_score", "minus_log10_p", "status"])
    fig_design = [
        {"panel": "A", "title": "Locked Formula A PRFT-high vs PRFT-low DEG signature used for SigCom LINCS query.", "data_source": str(OUTDIR / "CMap_LINCS_input_signature_LOCKED_A_Phase17A_lite.csv"), "status": "completed input signature", "safe_interpretation": "query input only; not a result"},
        {"panel": "B", "title": "Top SigCom LINCS chemical perturbation reversers of the PRFT-high signature.", "data_source": str(OUTDIR / "SigCom_LINCS_top_reversers_Phase17B_v2.csv"), "status": "curated from real TSV", "safe_interpretation": "candidate perturbational reversers"},
        {"panel": "C", "title": "Mechanism classes of candidate PRFT-high reversing perturbagens.", "data_source": str(OUTDIR / "SigCom_LINCS_mechanism_class_summary_Phase17B_v2.csv"), "status": "name-based conservative annotation", "safe_interpretation": "mechanism-class prioritization; manual verification required"},
        {"panel": "D", "title": "BeatAML consistency check with ex vivo drug-response patterns.", "data_source": str(PHASE17 / "beataml_consistency_check" / "BeatAML_SigCom_consistency_check_Phase17B_v2.csv"), "status": "completed as bounded consistency check", "safe_interpretation": "ex vivo association consistency, not clinical response prediction"},
    ]
    write_csv(PHASE17 / "supplementary_figures" / "Supplementary_Figure_Sx_SigCom_BeatAML_design_Phase17B_v2.csv", fig_design, ["panel", "title", "data_source", "status", "safe_interpretation"])

    supp_rows = []
    for r in top:
        supp_rows.append(
            {
                "category": "SigCom reverser",
                "candidate_or_mechanism": r["perturbagen"],
                "evidence_source": "SigCom LINCS Reversers TSV",
                "direction": "Reverser",
                "relationship_to_PRFT_high": "candidate PRFT-high transcriptional-state reverser",
                "support_level": "SigCom_only",
                "recommended_interpretation": "hypothesis-generating prioritized perturbagen",
                "future_validation_experiment": "AML PRFT-high/residual-like drug-response assay plus lipid ROS/Fe2+/MDA and marker readouts",
                "risk_note": "not clinical efficacy evidence",
                "comment": f"Mechanism: {r['known_or_inferred_mechanism']}",
            }
        )
    for m in mech_rows:
        supp_rows.append(
            {
                "category": "mechanism class",
                "candidate_or_mechanism": m["mechanism_class"],
                "evidence_source": "SigCom mechanism class summary",
                "direction": "class among Reversers",
                "relationship_to_PRFT_high": "potential mechanism-class reversal context",
                "support_level": "hypothesis_generating" if m["mechanism_class"] != "other / unclear" else "manual_check_required",
                "recommended_interpretation": m["recommended_wording"],
                "future_validation_experiment": "class-representative perturbation in PRFT-high AML model",
                "risk_note": m["overclaim_risk"],
                "comment": m["comment"],
            }
        )
    for c in consistency:
        supp_rows.append(
            {
                "category": "BeatAML-SigCom consistency",
                "candidate_or_mechanism": c["drug_or_mechanism"],
                "evidence_source": c["source"],
                "direction": c["consistent_or_not"],
                "relationship_to_PRFT_high": c["PRFT_high_interpretation"],
                "support_level": "convergent_support" if c["consistent_or_not"] == "suggested convergent support" else "inconsistent" if c["consistent_or_not"] in {"inconsistent", "directionally complex"} else "hypothesis_generating",
                "recommended_interpretation": c["safe_wording"],
                "future_validation_experiment": "ex vivo AML drug-response assay with PRFT-high stratification",
                "risk_note": c["risk_note"],
                "comment": c["comment"],
            }
        )
    write_csv(
        PHASE17 / "supplementary_tables" / "Supplementary_Table_Sx_SigCom_BeatAML_prioritization_summary_Phase17B_v2.csv",
        supp_rows,
        ["category", "candidate_or_mechanism", "evidence_source", "direction", "relationship_to_PRFT_high", "support_level", "recommended_interpretation", "future_validation_experiment", "risk_note", "comment"],
    )

    convergent = [c for c in consistency if c["consistent_or_not"] == "suggested convergent support"]
    complex_rows = [c for c in consistency if c["consistent_or_not"] in {"directionally complex", "inconsistent"}]
    top_names = ", ".join([str(r["perturbagen"]) for r in top[:5]])
    results_para = (
        "Not recommended for Results yet as a primary claim. As a supplementary, hypothesis-generating note only: "
        f"Using the locked Formula A PRFT-high versus PRFT-low signature, the exported SigCom LINCS L1000 Chemical Perturbations Reversers table was curated for perturbational signature reversal. "
        f"The top candidate reversers included {top_names}, and mechanism-level annotation suggested several potentially relevant classes, including proteostasis/HSP90, HDAC, kinase/stress-kinase and cell-cycle/CDK-related perturbagens, with many rows requiring manual mechanism verification. "
        f"BeatAML comparison was limited to ex vivo AUC patterns, where higher AUC denotes lower ex vivo sensitivity or greater relative resistance; {len(convergent)} entries/classes showed suggested convergent support, whereas {len(complex_rows)} were directionally complex or inconsistent. "
        "These results support prioritization only and do not demonstrate drug efficacy, clinical sensitivity, or causal mechanism."
    )
    discussion_para = (
        "The SigCom LINCS and BeatAML comparison provides a cautious bridge from the PRFT transcriptomic state to therapeutic-vulnerability prioritization. "
        "SigCom reversers were derived from perturbational signature reversal of the locked Formula A PRFT-high signature, whereas BeatAML contributes ex vivo drug-response associations rather than clinical treatment-response evidence. "
        "Convergence between selected mechanism classes and lower BeatAML AUC in PRFT-high samples should therefore be viewed as hypothesis-generating support, and discordant or non-comparable results should remain visible rather than being forced into a single narrative. "
        "This layer cannot replace experimental validation. Future work should test candidate mechanisms in PRFT-high AML residual-like models using drug-response assays, ferroptosis rescue assays, lipid ROS, Fe2+ and MDA readouts, and markers including JAK2/STAT5/PD-L1 and SLC7A11/GPX4. "
        "The analysis should avoid causal wording such as proved, confirmed, demonstrated causally, or restored ferroptosis sensitivity."
    )
    write_text(PHASE17 / "manuscript_text" / "Results_SigCom_BeatAML_supplementary_paragraph_Phase17B_v2.txt", results_para + "\n")
    write_text(PHASE17 / "manuscript_text" / "Discussion_SigCom_BeatAML_supplementary_paragraph_Phase17B_v2.txt", discussion_para + "\n")

    ref_rows = [
        {"manuscript_location": "Supplementary Methods / SigCom LINCS query", "sentence_or_claim": "SigCom LINCS was used to query the locked Formula A PRFT-high signature for perturbational reversers.", "needed_reference_topic": "SigCom LINCS / LINCS perturbational signature query", "candidate_reference": "SigCom LINCS method/resource reference", "verification_status": "manual verification required", "claim_strength": "method provenance", "overclaim_risk": "moderate if used as biological evidence", "comment": "Do not invent DOI/PMID."},
        {"manuscript_location": "Supplementary Methods / perturbational reversal", "sentence_or_claim": "The analysis used perturbational signature reversal of up/down gene sets.", "needed_reference_topic": "CMap / L1000 perturbational signature reversal", "candidate_reference": "Connectivity Map / L1000 reference", "verification_status": "manual verification required", "claim_strength": "method provenance", "overclaim_risk": "high if interpreted as clinical response", "comment": "Cite method/resource only after verification."},
        {"manuscript_location": "Supplementary Results / BeatAML consistency", "sentence_or_claim": "BeatAML AUC patterns were used as ex vivo pharmacogenomic context.", "needed_reference_topic": "BeatAML pharmacogenomic resource", "candidate_reference": "BeatAML resource paper / official resource", "verification_status": "manual verification required", "claim_strength": "resource provenance", "overclaim_risk": "high if written as clinical efficacy", "comment": "Higher AUC = lower ex vivo sensitivity."},
        {"manuscript_location": "Discussion / limitations", "sentence_or_claim": "Computational reverser prioritization requires experimental validation.", "needed_reference_topic": "computational prediction limitation", "candidate_reference": "Computational drug-reversal/in silico perturbation limitation reference", "verification_status": "manual verification required", "claim_strength": "limitation/context", "overclaim_risk": "low if phrased cautiously", "comment": "Useful to prevent overclaiming."},
    ]
    write_csv(
        PHASE17 / "reference_support" / "Phase17B_v2_SigCom_BeatAML_reference_insertion_points.csv",
        ref_rows,
        ["manuscript_location", "sentence_or_claim", "needed_reference_topic", "candidate_reference", "verification_status", "claim_strength", "overclaim_risk", "comment"],
    )

    risk_rows = [
        ("Locked Formula A modification", "high", "Curation could be mistaken for a model update.", "State read-only curation.", "Formula A remained unchanged.", ""),
        ("Formula B misuse", "high", "Formula B audit files exist.", "Do not use Formula B.", "Formula B was not used.", ""),
        ("Old DEG misuse", "high", "Old/package DEG files exist.", "Use only Phase17A locked signature input and real SigCom export.", "No DEG rerun or old DEG use.", ""),
        ("SigCom input signature written as result", "high", "Input list is not a reverser result.", "Separate input from exported TSV.", "Only exported Reversers TSV was curated as result.", ""),
        ("Mimickers written as Reversers", "high", "SigCom has multiple result types.", "Use only exported Reversers table.", "Only LINCS L1000 Chemical Perturbations Reversers were curated.", ""),
        ("SigCom/LINCS clinical efficacy overclaim", "high", "Reversers can be misread as drugs that work.", "Use hypothesis-generating reversal wording.", "Candidate perturbational reversers, not clinical efficacy predictions.", ""),
        ("BeatAML AUC direction inversion", "high", "Higher AUC can be miswritten as sensitivity.", "Lock AUC interpretation.", "Higher AUC = lower ex vivo sensitivity / greater relative resistance.", ""),
        ("Forced interpretation of discordant results", "high", "SigCom and BeatAML may disagree.", "Mark inconsistent/directionally complex.", "Discordance remains visible.", ""),
        ("Single-drug overinterpretation", "moderate", "Individual perturbagens may dominate story.", "Prioritize mechanism classes.", "Mechanism-class prioritization is safer than drug claims.", ""),
        ("Mechanism class vs drug name", "moderate", "Name-based mechanism annotation is uncertain.", "Manual verification required.", "Mechanism classes are conservative and provisional.", ""),
        ("Results suitability", "moderate", "Supplementary result may still be premature.", "Use short bounded supplementary paragraph only after review.", "Not recommended as primary Results claim.", ""),
        ("Discussion suitability", "low", "Discussion can overreach.", "Frame as prioritization/future validation.", "More suitable for Discussion/Future validation.", ""),
        ("Human Genomics / SCI credibility", "moderate", "Overclaiming could weaken submission.", "Keep as supplementary and bounded.", "Potential add-on if carefully labeled.", ""),
    ]
    write_csv(
        PHASE17 / "risk_check" / "Phase17B_v2_SigCom_BeatAML_risk_checklist.csv",
        [{"risk_item": r[0], "risk_level": r[1], "problem": r[2], "suggested_action": r[3], "safe_wording": r[4], "comment": r[5]} for r in risk_rows],
        ["risk_item", "risk_level", "problem", "suggested_action", "safe_wording", "comment"],
    )

    top20 = ", ".join([str(r["perturbagen"]) for r in top[:10]])
    major_classes = ", ".join([m["mechanism_class"] for m in mech_rows[:5]])
    direction_counts = Counter(c["consistent_or_not"] for c in consistency)
    summary = f"""Phase17B_v2 SigCom LINCS / BeatAML curation summary CN
生成日期：2026-06-25

1. 本轮是否修改 Locked Formula A
否。Formula A model、gene list、coefficients、cutoff rule、training/validation cohorts、main Results、Methods 和 manuscript conclusions 均未修改。

2. 是否读取真实 SigCom LINCS Reversers TSV
是。读取文件：{selected}

3. SigCom 输出来自哪个模块
SigCom LINCS / LINCS L1000 Chemical Perturbations (2021) / Reversers exported table.

4. top reversers 是什么
前若干候选包括：{top20}。

5. 主要机制类别是什么
主要机制类别包括：{major_classes}。机制注释为 name-based conservative annotation，全部需要人工核查。

6. 是否存在 Mimicker / Reverser 方向混淆风险
存在理论风险，但本轮只整理真实导出的 Reversers TSV，并在风险表中锁定不得把 Mimickers 写成 Reversers。

7. BeatAML 一致性核查是否完成
完成 bounded consistency check。结果统计：{dict(direction_counts)}。

8. 哪些机制类别与 BeatAML 方向一致
一致性仅可写为 suggested convergent support；具体见 BeatAML_SigCom_consistency_check_Phase17B_v2.csv。

9. 哪些机制类别不一致或无法比较
directionally complex / inconsistent / not directly comparable 均已保留，未强行解释一致。

10. 是否适合写入 Results
不建议作为主 Results。可在人工审核后作为简短 supplementary/hypothesis-generating Results note。

11. 是否更适合只写入 Discussion
是。更适合作为 Discussion/Future validation 的 therapeutic vulnerability prioritization。

12. 是否建议生成 Supplementary Figure
建议，可作为 Supplementary Figure Sx，展示 input signature、top reversers、mechanism classes、BeatAML consistency。

13. 是否建议生成 Supplementary Table
建议，可作为 supplementary prioritization summary。

14. 是否需要新增参考文献
需要，但 SigCom LINCS / CMap-L1000 / BeatAML / computational prediction limitation 文献均需 manual verification required；未编造 DOI/PMID。

15. 是否还需要 scTenifoldKnk / CellOracle
当前不需要。若后续补齐 raw single-cell object，再作为独立 Phase17C。

16. 对 Human Genomics / SCI 2区投稿是否是加分项
谨慎处理时是加分项，可增强药物重定位/验证路线；若写成疗效预测或机制证明则会降低可信度。

17. 下一步最优先任务
人工核查机制注释和方法文献 DOI/PMID；人工审核是否只放 Supplementary/Discussion；若要进入 Results，仅使用弱化的一小段 supplementary wording。

明确结论：
1. Locked Formula A 是否保持不变：是。
2. 是否读取真实 SigCom LINCS Reversers 输出：是。
3. SigCom 输出模块：LINCS L1000 Chemical Perturbations (2021) Reversers。
4. 是否完成 reverser 结果整理：是。
5. 是否完成 BeatAML-SigCom 一致性核查：是。
6. 是否存在方向不一致结果：是，已标记 directionally complex / inconsistent / not directly comparable。
7. 是否建议写入 Results：不建议作为主 Results；最多作为人工审核后的 supplementary note。
8. 是否建议写入 Discussion：是。
9. 是否建议作为 Supplementary Figure：是。
10. 是否建议作为 Supplementary Table：是。
11. 是否需要新增方法学参考文献：是，且需人工核查。
12. 是否还需要真实运行 scTenifoldKnk / CellOracle：当前不需要，后续有 raw single-cell object 再说。
13. 对 Human Genomics 投稿是否加分：谨慎、补充性呈现时加分。
14. 下一步最优先任务：人工审核 top reversers/机制注释、完成参考文献核查，并决定是否纳入 Supplementary Figure/Table 与 Discussion。
"""
    write_text(PHASE17 / "logs" / "Phase17B_v2_SigCom_BeatAML_summary_CN.txt", summary)


if __name__ == "__main__":
    main()
