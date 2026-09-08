# Changes

## 2.0.0 — workflow reliability and Colab usability

### Correctness and recovery

- Implement environment-variable credentials and optional Colab Secrets; exclude credentials from provenance files.
- Validate every protein record and BLAST result row before writing output.
- Resume matching BLAST RIDs, enforce request/poll intervals, retry safe requests, distinguish no hits from transport failures, and retain ambiguous submissions.
- Replace append-only neighborhood downloads with atomic, verified per-segment files. Record returned coordinates, mapping decisions and failures; retain the original WP/taxid representative policy.
- Preserve antiSMASH annotation flags when database downloading is disabled. Check prerequisites for the chosen analysis mode without deleting package-managed data.
- Verify input/parameter/tool/database fingerprints before reusing antiSMASH or cblaster artifacts. Keep failed attempts separate from completed results.
- Save cblaster search results independently of plotting. Normalize the NumPy clustering array before HTML export to fix cblaster 1.4.2's multi-result JSON serialization failure; retry plotting without repeating a successful search.
- Use one GenBank-aware filter/diagnosis implementation. Keep the original thresholds; support wrapped qualifiers, compound cores and unversioned Pfam IDs, and prevent cross-record domain deduplication.
- Reject ambiguous region names and exports inside the input tree. Passing lists reference exported GBKs, with a separate source list.
- Remove clinker preview's process-wide working-directory change and close an earlier server on rerun.

### Usability and resource use

- Validate anchor and marker sets before expensive steps. Support an existing marker path, explicit upload, and importing BLAST TSVs.
- Keep the notebook self-contained while testing the same helper code as a Python module. Clear execution/widget metadata.
- Bound outstanding worker tasks, cap antiSMASH concurrency by CPUs and estimated memory, and stream complete logs to files.
- Skip dependency solving when tools exist; use cblaster's native file list and batched database writes; cache unchanged searches.
- Make clinker optional, with a display limit that leaves the full passing set intact.
- Record versions, environment exports, fingerprints and per-candidate diagnoses. Add offline regression tests and CI.

### Migration

- Use a new `WORK_DIR`. Older outputs lack verified completion manifests and are not automatically adopted as successful caches.
- Import an existing standard BLAST TSV through `BLAST_TSV_INPUT`. Keep old FASTAs, annotations and reports as historical records.
- The PfaA example, homology defaults, WP/taxid policy and filtering thresholds are unchanged. CPU scheduling and output organization change. Parser repairs can recover results previously missed by unsupported annotations.
- `DOWNLOAD_DB=False` means reuse prepared databases with the selected annotations; it no longer silently omits database-dependent analyses.
- `CLEAN_OLD` and unchecked `OVERWRITE` behavior are replaced by checkpoint recovery and isolated attempts/exports.
- Full Colab execution, the complete published candidate set and real-workload speed/memory improvements are not validated by the offline suite.
