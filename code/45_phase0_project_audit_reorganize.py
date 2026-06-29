from __future__ import annotations

import csv
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "AML_PRFT_Human_Genomics_phase0_audit"

SUBDIRS = [
    "00_raw_data",
    "01_processed_data",
    "02_scripts",
    "03_results_tables",
    "04_figures",
    "05_logs",
    "06_manuscript_support",
]

PROJECT_DIRS = [
    "00_raw_data",
    "01_metadata",
    "02_processed_data",
    "03_gene_sets",
    "04_prft_score",
    "05_deg",
    "06_wgcna",
    "07_signature",
    "08_validation",
    "09_enrichment",
    "10_immune",
    "11_drug",
    "12_single_cell",
    "13_figures",
    "14_tables",
    "15_scripts",
    "16_logs",
    "Bioinformatics_audit_and_wet_validation_package",
    "Human_Genomics_PRFT_AML_OFFICIAL_submission_package",
    "Human_Genomics_PRFT_AML_INTERNAL_working_package",
    "Human_Genomics_PRFT_AML_submission_package_clean",
]

RAW_EXTS = {".gz", ".tsv", ".soft", ".mtx", ".h5", ".h5ad", ".loom", ".txt"}
PROCESSED_EXTS = {".rds", ".rdata"}
SCRIPT_EXTS = {".r", ".rmd", ".py"}
TABLE_EXTS = {".csv", ".tsv", ".xlsx"}
FIGURE_EXTS = {".png", ".pdf", ".tiff", ".tif", ".svg"}
MANUSCRIPT_EXTS = {".docx", ".md", ".txt"}

SKIP_DIR_PARTS = {
    ".git",
    ".codex",
    ".agents",
    "vendor",
    "outputs",
    "FQD_augmented_network_meta",
    "FQD_bioinformatics_work",
    "FQD_BMC_final_20260622",
    "FQD_original_reintegration",
    "JIMR_audit_extracted",
    "paper_rewriting_output",
    "results",
    "scripts",
    "17_tmp",
    OUT.name,
}


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace("\\", "/")


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for sub in SUBDIRS:
        (OUT / sub).mkdir(parents=True, exist_ok=True)


def log(message: str) -> None:
    log_path = OUT / "05_logs" / "phase0_audit_log.txt"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as f:
        f.write(f"[{now()}] {message}\n")


def should_scan(path: Path) -> bool:
    parts = set(path.relative_to(ROOT).parts)
    return not any(part in SKIP_DIR_PARTS for part in parts)


def scan_files() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for dirname in PROJECT_DIRS:
        base = ROOT / dirname
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if not path.is_file() or not should_scan(path):
                continue
            suffix = path.suffix.lower()
            category = categorize(path)
            records.append(
                {
                    "relative_path": rel(path),
                    "category": category,
                    "suffix": suffix,
                    "size_bytes": path.stat().st_size,
                    "last_modified": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
                }
            )
    return sorted(records, key=lambda x: (str(x["category"]), str(x["relative_path"])))


def categorize(path: Path) -> str:
    rp = rel(path).lower()
    suffix = path.suffix.lower()
    if rp.startswith("00_raw_data/"):
        if "geo" in rp:
            return "GEO raw data"
        if "tcga" in rp:
            return "TCGA raw data"
        return "raw data"
    if "beataml" in rp or "drug" in rp:
        if suffix in TABLE_EXTS:
            return "BeatAML/drug data or result"
    if "single_cell" in rp or "single-cell" in rp or "sc_" in path.name.lower():
        if suffix in TABLE_EXTS or suffix in PROCESSED_EXTS:
            return "single-cell data or result"
    if suffix in SCRIPT_EXTS:
        return "script"
    if suffix in PROCESSED_EXTS or rp.startswith(("01_metadata/", "02_processed_data/", "03_gene_sets/", "04_prft_score/")):
        return "processed expression/clinical/signature data"
    if suffix in TABLE_EXTS:
        return "result table"
    if suffix in FIGURE_EXTS:
        return "figure"
    if suffix == ".docx":
        return "manuscript/support document"
    if rp.startswith("16_logs/") or "sessioninfo" in rp or "log" in rp:
        return "log"
    return "other"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def copy_preserving_relative(src: Path, dest_root: Path) -> None:
    dest = dest_root / src.relative_to(ROOT)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


def copy_phase0_working_files(records: list[dict[str, object]]) -> dict[str, int]:
    copied = {sub: 0 for sub in SUBDIRS}
    for rec in records:
        src = ROOT / str(rec["relative_path"])
        category = str(rec["category"])
        suffix = src.suffix.lower()

        if category in {"TCGA raw data", "GEO raw data", "raw data"}:
            copy_preserving_relative(src, OUT / "00_raw_data")
            copied["00_raw_data"] += 1
        elif category == "processed expression/clinical/signature data":
            copy_preserving_relative(src, OUT / "01_processed_data")
            copied["01_processed_data"] += 1
        elif category == "script":
            copy_preserving_relative(src, OUT / "02_scripts")
            copied["02_scripts"] += 1
        elif category in {"result table", "BeatAML/drug data or result", "single-cell data or result"} and suffix in TABLE_EXTS | PROCESSED_EXTS:
            copy_preserving_relative(src, OUT / "03_results_tables")
            copied["03_results_tables"] += 1
        elif category == "figure":
            copy_preserving_relative(src, OUT / "04_figures")
            copied["04_figures"] += 1
        elif category == "log":
            copy_preserving_relative(src, OUT / "05_logs")
            copied["05_logs"] += 1
        elif category == "manuscript/support document":
            copy_preserving_relative(src, OUT / "06_manuscript_support")
            copied["06_manuscript_support"] += 1

    return copied


def run_r_environment_capture() -> None:
    session_path = OUT / "05_logs" / "sessionInfo_initial.txt"
    pkg_path = OUT / "05_logs" / "installed_R_packages_initial.csv"
    pkg_rel = "AML_PRFT_Human_Genomics_phase0_audit/05_logs/installed_R_packages_initial.csv"
    try:
        header = (
            "Phase 0 AML-PRFT Human Genomics project audit\n"
            f"Timestamp: {now()}\n\n"
        )
        session_path.write_text(header, encoding="utf-8")
        with session_path.open("a", encoding="utf-8") as out:
            subprocess.run(["Rscript", "-e", "sessionInfo()"], cwd=ROOT, check=True, stdout=out, stderr=out, text=True)
        pkg_code = (
            "pkgs <- as.data.frame(installed.packages()[, c('Package','Version','LibPath','Priority')]); "
            f"write.csv(pkgs, file={pkg_rel!r}, row.names=FALSE)"
        )
        subprocess.run(["Rscript", "-e", pkg_code], cwd=ROOT, check=True)
        log("Captured R sessionInfo_initial.txt and installed_R_packages_initial.csv.")
    except Exception as exc:
        session_path.write_text(f"R environment capture failed: {exc}\n", encoding="utf-8")
        log(f"R environment capture failed: {exc}")


def summarize(records: list[dict[str, object]], category: str, limit: int = 80) -> list[str]:
    rows = [r for r in records if r["category"] == category]
    out = []
    for r in rows[:limit]:
        out.append(f"- `{r['relative_path']}` ({r['size_bytes']} bytes)")
    if len(rows) > limit:
        out.append(f"- ... {len(rows) - limit} additional files listed in `05_logs/file_inventory_phase0.csv`")
    if not out:
        out.append("- Not detected in scanned project-relevant folders.")
    return out


def detect_missing_and_rerun_needs(records: list[dict[str, object]]) -> tuple[list[str], list[str]]:
    paths = {str(r["relative_path"]).lower() for r in records}
    names = {Path(str(r["relative_path"])).name.lower() for r in records}
    missing = []
    rerun = []

    expected = {
        "TCGA raw expression": any("tcga.laml.samplemap_hiseqv2.gz" in p for p in paths),
        "TCGA clinical matrix": any("clinicalmatrix" in p for p in paths),
        "GEO GPL570 validation matrices": any("gpl570_series_matrix" in p for p in paths),
        "GEO platform annotations GPL570/GPL96/GPL97": any("gpl570_family.soft.gz" in p for p in paths)
        and any("gpl96_family.soft.gz" in p for p in paths)
        and any("gpl97_family.soft.gz" in p for p in paths),
        "processed TCGA expression": "tcga_expr_hgnc_log2cpm.rds" in names,
        "processed clinical data": "tcga_clin_clean.rds" in names,
        "PRFT score table": "tcga_prft_score.csv" in names,
        "DEG full table": "deg_prft_high_vs_low_tcga_all.csv" in names or "03_all_tcga_prft_high_vs_low_deg_full.csv" in names,
        "WGCNA module table": "wgcna_module_assignments.csv" in names or "04_all_gene_module_membership.csv" in names,
        "candidate filtering trace": "05_candidate_gene_filtering_full_trace.csv" in names,
        "external validation results": "external_validation_summary.csv" in names or "07_external_validation_all_results.csv" in names,
        "BeatAML drug association table": "beataml_risk_score_drug_correlation_all.csv" in names or "09_beataml_all_drug_associations.csv" in names,
        "single-cell cell-type table": "sc_celltype_score_summary.csv" in names or "10_single_cell_all_celltype_scores.csv" in names,
        "main manuscript": any("main_manuscript" in n and n.endswith(".docx") for n in names),
        "Human Genomics official package": any("human_genomics_prft_aml_official_submission_package" in p for p in paths),
    }
    for label, ok in expected.items():
        if not ok:
            missing.append(label)

    if missing:
        rerun.append("Do not rerun yet. First recover or confirm missing files listed above.")
    if not missing:
        rerun.append("No mandatory module needs immediate rerun for Phase 0; rerun decisions should wait for Human Genomics manuscript-level gap review.")
    rerun.extend(
        [
            "Reference verification still needs manual or online bibliographic checking before journal upload.",
            "DOCX visual render QA needs Word/LibreOffice availability; local prior render attempts lacked LibreOffice/soffice.",
            "Layer 1 wet-lab validation is recommended for upgrade, but it is outside Phase 0 and should not be mixed with this audit.",
        ]
    )
    return missing, rerun


def build_report(records: list[dict[str, object]], copied: dict[str, int]) -> None:
    counts: dict[str, int] = {}
    sizes: dict[str, int] = {}
    for rec in records:
        cat = str(rec["category"])
        counts[cat] = counts.get(cat, 0) + 1
        sizes[cat] = sizes.get(cat, 0) + int(rec["size_bytes"])

    missing, rerun = detect_missing_and_rerun_needs(records)
    report = [
        "# AML-PRFT Human Genomics Phase 0 Project Audit",
        "",
        f"- Audit timestamp: {now()}",
        f"- Workspace root: `{ROOT}`",
        f"- New audit/reorganization directory: `{OUT}`",
        "- Scope: project audit and file reorganization only. No formal statistical analysis was started.",
        "- Safety: no original files were deleted or overwritten; relevant files were copied into the Phase 0 directory.",
        "",
        "## New Directory Structure",
        "",
    ]
    for sub in SUBDIRS:
        report.append(f"- `{sub}/` copied files: {copied.get(sub, 0)}")

    report.extend(["", "## File Inventory Summary", ""])
    for cat in sorted(counts):
        report.append(f"- {cat}: {counts[cat]} files, {sizes[cat]} bytes")

    report.extend(["", "## Existing Data Files", "", "### TCGA / raw data"])
    report.extend(summarize(records, "TCGA raw data"))
    report.extend(["", "### GEO data"])
    report.extend(summarize(records, "GEO raw data"))
    report.extend(["", "### Processed expression / clinical / PRFT data"])
    report.extend(summarize(records, "processed expression/clinical/signature data"))
    report.extend(["", "### BeatAML data/results"])
    report.extend(summarize(records, "BeatAML/drug data or result"))
    report.extend(["", "### Single-cell data/results"])
    report.extend(summarize(records, "single-cell data or result"))

    report.extend(["", "## Existing Scripts", ""])
    report.extend(summarize(records, "script", limit=120))
    report.extend(["", "## Existing Results", "", "### Result tables"])
    report.extend(summarize(records, "result table", limit=120))
    report.extend(["", "### Figures"])
    report.extend(summarize(records, "figure", limit=120))
    report.extend(["", "### Manuscript and support documents"])
    report.extend(summarize(records, "manuscript/support document", limit=120))

    report.extend(["", "## Missing Files", ""])
    if missing:
        report.extend(f"- {m}" for m in missing)
    else:
        report.append("- No mandatory Phase 0 project component was missing from the scanned project-relevant folders.")

    report.extend(["", "## Modules That May Need Rerun Or Follow-up", ""])
    report.extend(f"- {m}" for m in rerun)

    report.extend(
        [
            "",
            "## R Environment",
            "",
            "- R version and session information: `05_logs/sessionInfo_initial.txt`",
            "- Installed R packages: `05_logs/installed_R_packages_initial.csv`",
            "",
            "## Detailed Inventories",
            "",
            "- Full file inventory: `05_logs/file_inventory_phase0.csv`",
            "- Category counts: `05_logs/category_counts_phase0.json`",
            "- Operation log: `05_logs/phase0_audit_log.txt`",
        ]
    )
    (OUT / "project_audit_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    log("Wrote project_audit_report.md.")


def main() -> None:
    ensure_dirs()
    log("Started Phase 0 audit. No formal statistical analysis will be run.")
    log("Created/verified Phase 0 directory structure.")
    records = scan_files()
    log(f"Scanned project-relevant folders. Files detected: {len(records)}.")
    write_csv(OUT / "05_logs" / "file_inventory_phase0.csv", records)
    log("Wrote full file inventory CSV.")
    counts = {}
    for rec in records:
        counts[str(rec["category"])] = counts.get(str(rec["category"]), 0) + 1
    (OUT / "05_logs" / "category_counts_phase0.json").write_text(json.dumps(counts, indent=2, ensure_ascii=False), encoding="utf-8")
    log("Wrote category count JSON.")
    copied = copy_phase0_working_files(records)
    log(f"Copied relevant files into Phase 0 structure: {copied}.")
    run_r_environment_capture()
    build_report(records, copied)
    log("Completed Phase 0 audit.")
    print(OUT)
    print(OUT / "project_audit_report.md")
    print(OUT / "05_logs" / "sessionInfo_initial.txt")
    print(OUT / "05_logs" / "phase0_audit_log.txt")


if __name__ == "__main__":
    main()
