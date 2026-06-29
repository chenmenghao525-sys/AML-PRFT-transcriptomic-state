# Private Git Sync Execution Report

Generated: 2026-06-29 16:15:46 +08:00

Final conclusion: fail_private_git_sync

## Scope

Created a true Git working copy, preserved its .git directory, and mirrored the QA-passed v2 candidate content from the read-only source candidate. No analysis was rerun and no manuscript, figure, or supplementary table content was edited.

## Execution Status

| Check | Status | Notes |
|---|---|---|
| True Git working copy with .git | pass | .git exists in the working copy. |
| v2 candidate content mirrored | pass | Root contains code, docs, figures, metadata, results, supplementary_tables, README.md, and LICENSE_note.md. |
| Local sync commit created | pass | 8b92f5984fc738500caff2ef7878757f82185062; subject: Sync clean v2 private QA candidate. |
| Local QA record commit created | pass | This report is included in the local QA record commit; subject: Record private Git sync network blocker QA. |
| Push to origin/main | fail | GitHub HTTPS/TCP 443 was unreachable; push failed before authentication/remote update. |
| Git working tree clean | pass_local | Final local status is clean after committing the QA record. |
| Remote private status | not_live_verified | The user had confirmed the repository is Private; no visibility-changing action was run. Live check was blocked by GitHub 443 connectivity failure. |

## Risk Scan Summary

Public-facing scan of the working copy found no blocking residuals for Formula B, Supplementary Figure S10 figure references, Figure 8, 126 successfully fitted, old Figure 2 values, local E:\ paths, or credential-like strings. Internal QA metadata contains historical risk strings only as audit records and is classified as advisory. Supplementary-table numeric coincidences are recorded separately as nonblocking advisory findings.

## Required Answers

1. 真正带 .git 的 Git working copy: yes.
2. v2 candidate 内容是否已同步进去: yes.
3. 是否已 push 到 origin/main: no; blocked by GitHub 443 network connectivity failure.
4. git status 是否 clean: yes locally after committing the QA record.
5. GitHub repo 是否仍为 Private: live verification blocked; no action was taken to change visibility, and user-provided precondition said Private.
6. 是否发现旧值、S10 figure、Figure 8、Formula B: no blocking public-facing residuals.
7. Figure 3-6 是否存在: yes.
8. Supplementary Figures S1-S9 是否齐全: yes.
9. Supplementary Tables S1-S15 是否齐全，S10 table 是否保留: yes; Supplementary Table S10 is retained.
10. 是否需要重跑分析: no.
11. 是否可以进入 DOI/URL decision: no, not until push succeeds and remote/private status is live-verified.
