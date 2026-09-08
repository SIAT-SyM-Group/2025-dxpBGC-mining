# Validation scope

## Completed locally

- **48 offline tests passed** on Python 3.12.2 with Biopython 1.88, requests 2.32.2 and pytest 7.4.4. The suite validates notebook schema/source synchronization and executes its setup/input/import cells without Colab or remote searches. Display tests cover empty results, HTML escaping, saved decisions and preview limits without truncating totals.
- The presentation revision preserves `workflow_support.py` byte for byte and preserves the Python AST of all parameter/analysis cells. Only the generated helper cell and the diagnosis display cell gain presentation code; the PfaA sequence and scientific defaults are unchanged.
- **Eight synthetic standard-annotation cases** were compared with the filter extracted from upstream commit `a67a4dc220b1090cc9e16356a1280bd4e3c94fa3`: gaps of 0, 24,999, 25,000 and 25,001 bp, each with one or two A/C domains and duplicate fullhmmer/clusterhmmer annotations. Original and revised pass/fail decisions agreed in all eight cases.
- **Actual cblaster 1.4.2 interfaces** were exercised with two synthetic GenBank inputs: native file-list input, batch size 1, SQLite generation, local-search output formats and downstream filtering. The DIAMOND executable was deliberately stubbed; this checks interfaces and data handoff, not homology search accuracy or performance. One of the two regions passed the independently specified gap rule.
- That interface check reproduced cblaster 1.4.2's multi-result `ndarray is not JSON serializable` plot failure. Separating search/plotting and normalizing rendering data resolved it. A regression test verifies that a plot failure does not repeat a completed search.
- **Actual clinker 0.0.32** generated HTML and session JSON from the synthetic GenBanks. Its registered `clinker.main:main` entry point was used; the package does not provide `python -m clinker`. The generated HTML loaded in a browser without console errors.
- The current antiSMASH 8.0.4 MITE implementation was checked against its upstream source: it reads the configured database directory. The previous notebook's destructive replacement of package-managed MITE data was removed; prerequisite failures now remain visible.

The included GitHub Actions workflow runs the offline suite on Linux with Python 3.10 and 3.12. Its actual result is reported by the PR checks; a local pass alone is not a claim of CI success.

## Not established by these checks

- No large or real NCBI BLAST/IPG/nuccore search was launched.
- The full Linux micromamba environment solve, antiSMASH database download and real antiSMASH execution have not been run as part of this revision.
- There is no verified reconstruction of the published study's complete candidate set; the repository does not include the separate cblaster marker FASTA or immutable database snapshots.
- There is no end-to-end speed or peak-memory benchmark, nor a measured equivalence comparison of `full` and `filter-required` annotation modes. Resource improvements are bounded scheduling, batched I/O, complete logs on disk and elimination of verified redundant work.

## Before scaling a scientific run

Use the intended Linux runtime and supply the exact anchor/marker FASTAs. Start with a small known set, review every stage's manifest, inspect predicted regions and diagnose both expected positives and negatives. Save the complete working directory, explicit environment export and reference database snapshot. Only then increase the search limits or change scientific thresholds.
