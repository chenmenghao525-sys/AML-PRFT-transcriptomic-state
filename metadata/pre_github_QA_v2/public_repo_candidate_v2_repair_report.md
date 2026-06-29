# Public Repo Candidate v2 Repair Report

Generated: 2026-06-29 15:21:54

Final conclusion: pass_ready_for_private_github_sync

## Scope

Created and repaired `AML_PRFT_Human_Genomics_public_repo_candidate_v2` as a new repository candidate. The source phase0 audit directory and v1 candidate were used as read-only inputs. No manuscript DOCX, figure content, supplementary table content, analysis, or statistic was modified.

## Repair Summary

- Excluded v1 code files flagged by the pre-GitHub author review scan.
- Retained mandatory final Figure 2 v3 and Figure 7 v2 files from v1.
- Added confirmed final Figure 3-Figure 6 files from the official submission package.
- Completed Supplementary Figures S1-S9 only.
- Completed Supplementary Tables S1-S15 and retained the legal tenth supplementary table.
- Regenerated repository-relative metadata files and v2 QA outputs.

## Postcheck

- Old Figure 2 final-value residuals clear: true
- ACOX2 obsolete coefficient residual clear: true
- Alternate formula label clear: true
- Prohibited S-figure reference clear: true
- Prohibited main-figure-eight reference clear: true
- 126 successfully fitted clear: true
- Local absolute paths clear: true
- Credential-like strings clear: true
- Figure 2 v3 present: true
- Figure 7 v2 present: true
- Figure 3-Figure 6 complete: true
- Supplementary Figures S1-S9 complete: true
- Prohibited S-figure file absent: true
- Supplementary Tables S1-S15 complete: true
- Supplementary Table S10 retained: true
- README and manifests use repository-relative paths: true
- phase0_audit modified: no
- v1 candidate modified: no

Supplementary table numeric advisory: legal supplementary tables contain ordinary numeric values that can match old-value substrings; these are recorded as nonblocking advisory rows in the risk scan because table content was not modified and no final Figure 2/legend context was detected.

## Required Answers

1. v2 candidate generated: yes.
2. phase0_audit modified: no.
3. v1 candidate modified: no.
4. Need rerun analysis: no.
5. Figure 8 code remnants cleared: yes.
6. Figure 3-Figure 6 completed: yes.
7. Supplementary Figures S1-S9 completed: yes.
8. Supplementary Tables S1-S15 completed and Supplementary Table S10 retained: yes.
9. Old values, prohibited S-figure, Figure 8, Formula B still found: no.
10. Can enter private GitHub sync: yes.
