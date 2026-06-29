# GitHub Repository Final Replacement QA

Generated: 2026-06-30 04:30:42

## Scope

Repository figures, supplementary figures, and supplementary tables were replaced using the final formal upload package as the only source of truth. No analysis was rerun, no figure was redrawn, no manuscript scientific content was modified, Formula A was not changed, and statistical results were not recalculated. The repository was not made Public by this step.

## Replacement Summary

- Old GitHub figure/table folders were backed up locally under an ignored internal directory: `metadata/replaced_old_github_files_backup_do_not_upload/`.
- `figures/main/` was replaced with Figure 1-Figure 7 from the final upload package.
- `figures/final/` was synchronized with the same Figure 1-Figure 7 files for compatibility.
- `figures/supplementary_figures/` was replaced with Supplementary Figure S1-S8 from the final upload package.
- `supplementary_tables/` was replaced with Supplementary Tables S1-S15 from the final upload package.
- README and repository manifests were updated for the new paths.

## QA Answers

1. Main figures replaced from final upload package: yes.
2. Supplementary figures replaced from final upload package: yes.
3. Supplementary tables replaced from final upload package: yes.
4. Old Figure 3 wrong coefficient panel removed: yes.
5. Formula A correct in Figure 2: yes.
6. Formula B residue: no.
7. Old coefficient residue in public repository: no.
8. Old Figure 2 value residue: no.
9. Local path residue: no.
10. Figure 8/S9/S10 figure residue: no.
11. Fake DOI or Zenodo DOI residue: no.
12. Analysis rerun: no.
13. Ready to make repository public: yes.

## Visual/Provenance Check

The specified visual files were checked by source-of-truth replacement and SHA-256 equality against the final upload package:

- `figures/main/Figure 2.pdf`
- `figures/main/Figure 3.pdf`
- `figures/main/Figure 5.pdf`
- `figures/main/Figure 6.pdf`
- `figures/supplementary_figures/Supplementary Figure S2.pdf`
- `figures/supplementary_figures/Supplementary Figure S7.pdf`

Because these files now match the final upload package byte-for-byte, Figure 3 no longer uses the old public repository coefficient panel source, Figure 5/6 and Supplementary Figures S2/S7 are the final upload-package versions, and BeatAML remains assigned to Supplementary Figure S7.

Final conclusion: pass_ready_to_make_repository_public
