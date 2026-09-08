# DxpBGC genome mining

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SIAT-SyM-Group/2025-dxpBGC-mining/blob/main/workflow.ipynb)

Find candidate biosynthetic gene clusters by combining **protein homology, genomic proximity and domain evidence**:

**BLASTP anchor → WP/IPG genome mapping → genomic neighborhoods → antiSMASH 8 → cblaster marker search → filtering → optional clinker view.**

The Colab notebook embeds its helper code and needs no separate download from a moving branch. The same code lives in [`workflow_support.py`](workflow_support.py) for inspection and offline testing.

**中文快速说明：**先配置 NCBI 联系邮箱，再准备两类输入——BLAST 锚定蛋白和 cblaster 多蛋白标记集。默认 `U=3` 要求至少三个不同标记蛋白。按步骤运行；重复运行会校验缓存并继续未完成的任务。通过与未通过的候选都有诊断记录，旧结果不会被覆盖。

## Start in Colab

1. Open the notebook in a **Linux x86_64 CPU runtime**. No GPU is needed. antiSMASH and its databases still require substantial setup time and disk space, even for the small example.
2. Run **0–1**. Set `WORK_DIR` and `NCBI_EMAIL`. Prefer **Colab Secrets** or environment variables named `NCBI_EMAIL` and `NCBI_API_KEY`; non-empty form fields take precedence. API keys are optional. Credentials are excluded from manifests; clear form values and execution outputs before sharing a notebook.
3. Run **2** to validate the BLAST anchor. The original **AAN54658.1-PfaA sequence (2,531 aa)** is unchanged. Raw amino-acid text or FASTA is accepted. Empty records, duplicate first-token IDs, gaps, stop symbols and invalid characters are reported before writing.
4. Run **3** to validate your **separate marker-protein FASTA**. Set `QUERY_FASTA` or enable `UPLOAD_MARKERS` and upload one file. The repository does **not** include a marker set. Default `U=3` requires at least three distinct marker IDs; one PfaA anchor is insufficient.
5. Run **4–8** in order, reviewing each step's configuration and result. Once inputs are configured, reruns can proceed without another upload. Set `BLAST_TSV_INPUT` to import a standard 12-column BLAST result and bypass the remote search.
6. Inspect **9** for explanations. Enable **10** only when a clinker view is wanted. Its default 20-region display limit leaves the full passing set intact.

The default **10 BLAST hits and one representative taxid** are a small demonstration, not a reconstruction of a published mining dataset. Study reproduction requires the exact queries, marker set, database snapshots, tool versions and parameters.

To survive a Colab runtime reset, use an already-mounted persistent `WORK_DIR` and save the complete directory. Files under `/content` disappear when the runtime resets. Tools and databases can remain on runtime-local disk, but may need reinstallation after a reset.

## Parameters

| Stage | Parameters | Effect |
| --- | --- | --- |
| BLAST | `DATABASE`, `HIT_TOP_N`, `EXPECT` | Defaults: `refseq_protein`, 10 hits, E-value `1e-5`. |
| BLAST waiting | `POLL_SEC`, `MAX_WAIT_H`, `FORCE_NEW_BLAST` | Polling is clamped to ≥60 s. Fractional waiting hours are supported. Rerun resumes a RID; force-new deliberately resubmits. |
| Neighborhoods | `FLANK`, `TOP_N_UNIQ_TAXID`, `PREFER_REFSEQ` | ±100,000 bp; one representative taxid; prefer `NC_`/`NZ_`. Taxid limit 0 means no limit. |
| Network | `THREADS`, `REFRESH_IPG` | Bounded download workers share an Entrez rate limiter. IPG XML is cached; refresh explicitly for updated mappings. |
| antiSMASH setup | `TOOL_ROOT`, `DBDIR`, `DOWNLOAD_DB` | Dedicated Linux environment and database paths. Downloading is independent of annotation modules. |
| Annotation | `ANNOTATION_MODE` | `full` preserves the original modules. Opt-in `filter-required` omits extra analyses while retaining full/cluster Pfam scans. |
| Resources | `MAX_JOBS`, `CPUS_PER_JOB`, `MEMORY_GB_PER_JOB` | CPU and estimated-memory caps. Defaults: 2 jobs, 1 CPU each, 6 GiB/job. This is a scheduling heuristic, not a memory guarantee. |
| cblaster | `U`, `MH`, `MI`, `MC`, `HS` | Original defaults: 3 unique queries, 3 hits, 20% identity, 50% coverage, 10,000 hit limit. Passed to cblaster 1.4.2, whose local parser uses strict `>` identity/coverage cutoffs. |
| Database build | `CBLASTER_CPUS`, `DB_BATCH_SIZE` | CPU cap and input files per parse/write batch; defaults 2 and 32. Native file lists avoid shell argument limits. |
| Filtering | `PROTO_GAP_MAX`, `COL`, `COPY_MODE`, `DRY_RUN` | Core-gap cutoff, one-based region-name column, export type and preview-only export. Defaults: 25,000 bp, column 1, complete report directory, false. |

The `full` mode retains `fullhmmer`, `clusterhmmer`, `tfbs`, `cb-general`, `cb-knownclusters`, `cb-subclusters`, `rre`, `cc-mibig`, `tigrfam`, `asf`, `pfam2go`, and `smcog-trees`. **`DOWNLOAD_DB=False` now preserves these options** and checks prerequisites for the chosen mode instead of silently dropping evidence needed downstream.

`filter-required` is an explicit alternative that omits comparative/extra annotations while retaining `fullhmmer` and `clusterhmmer`. It is not the default. Runtime improvement and candidate equivalence between modes have **not** been benchmarked on the published dataset.

## Filter definition

The biological thresholds are unchanged. Each matched region must meet all conditions:

1. GenBank text contains `PKS_AT` and `ketoacyl synthase`, case-insensitively.
2. Unique **`PFAM_domain` features** contain **≥2 PF00501 domains OR ≥2 AMP-binding domains**, and **≥2 condensation domains** (`PF00668`, `Condensation_LCL` or `Condensation_DCL`). A condensation feature matching both identifiers counts once.
3. NRPS and hglE-KS protocluster cores have a minimum gap **≤25,000 bp** by default. `core_location` uses zero-based half-open intervals; touching/overlapping cores have gap zero. Cores on separate sequence records are not paired.
4. An exact `NRPS` product annotation exists.

Identical Pfam annotations from fullhmmer and clusterhmmer are deduplicated by sequence record, location, Pfam IDs, domain name and protein coordinates. Biopython parses wrapped qualifiers, versionless Pfam IDs and compound core locations. These parser repairs may recover candidates that the former line-based parser missed; unchanged thresholds do **not** imply identical results for previously unsupported annotations.

`filter_diagnosis.tsv` records evidence and failed gates for **every matched region**. Selection and diagnosis use the same function and threshold. Read/parse errors and ambiguous/unmatched names stop export with an explicit report rather than being treated as biological negatives. Valid zero-hit or zero-passing results produce an empty passing list and a diagnosis file.

**`kept_region_files.abs.txt` is the authoritative passing list.** `COPY_MODE=dir` preserves complete antiSMASH folders, which can contain other, non-passing regions. Use `COPY_MODE=gbk` to export only passing GenBanks. The passing list references exported files; `kept_source_files.abs.txt` retains their sources. With `DRY_RUN=True`, files are not copied and the list refers to sources; clinker requires a real export.

Co-localization and domain evidence do not establish compound identity, bioactivity or experimentally demonstrated biosynthesis.

## Recovery and provenance

Stages isolate outputs by an input/parameter fingerprint. Downstream steps consume the **current manifest**, not a directory-wide collection of old outputs.

| Situation | Behavior |
| --- | --- |
| BLAST interrupted or waiting budget reached | Rerun with the same inputs to resume the RID, without another submission. |
| Submission timed out before a RID was confirmed | Preserve the ambiguous state. Inspect it before deliberately setting `FORCE_NEW_BLAST`; POST is never automatically retried. |
| RID failed or expired | Report the state. Force-new preserves the previous checkpoint/results before submitting again. |
| Temporary GET/Entrez failure | Bounded retries, request spacing and explicit timeouts; failures are recorded. |
| Repeated neighborhood download | Atomically written, hash-checked segments; no repeated FASTA appending. |
| Partial antiSMASH job | Preserve logs and partial output; retry in a fresh folder and reuse only verified successes. |
| Changed inputs, options, tools or database inventory | Invalidate affected caches and preserve earlier results. |
| Changed cblaster GBKs | Build a new database; check SQLite integrity and the presence of every source region. |
| cblaster plot export fails | Preserve the successful search and retry only plotting. NumPy clustering arrays are normalized before rendering to handle cblaster 1.4.2's multi-result serialization issue. |
| Changed filter parameters | Create a separate export and diagnosis. |

Inputs and cached outputs use **SHA-256**. Large antiSMASH reference databases use a fast **path/size/mtime inventory** to avoid repeatedly hashing many GB; this will not detect a same-size change with a deliberately preserved timestamp. Use an external immutable database snapshot for stricter reproducibility.

`runtime.json` records the helper-source hash, Python and package versions. `setup/tool-versions.json` and `setup/environment-explicit.txt` record the actual Linux environment. antiSMASH is constrained to major version 8: the initial solver selects its patch release, which is then recorded and reused. cblaster is constrained to 1.4.2. The original notebook did not record its full environment, so this does not claim to reconstruct it.

NCBI databases evolve. Keep the input FASTAs, saved BLAST results, cached IPG XML, selected accession versions, coordinates, environment exports and manifests together.

## Outputs

Under `WORK_DIR` (default `dxp-results/`):

```text
queries/                  anchor.faa and markers.faa
runtime.json              Python/packages and helper-source hash
blast/<id>/               job.json (RID), raw responses, validated hits.tsv
neighborhoods/ipg-cache/  accession-specific IPG XML
neighborhoods/<id>/       selection.tsv, manifest.json, download_errors.tsv, fasta/
setup/                    install/prerequisite logs, versions, explicit environment
antismash/                isolated attempts, reports, region GBKs, batch manifests
cblaster-db/<id>/         input list, SQLite/DIAMOND/FASTA, manifest and build log
cblaster-search/<id>/     session.json, summary.tsv, abspres.tsv, plot.html and log
final-target/<id>/        diagnosis, passing/source lists, export, optional clinker/
final-target/latest.json pointer to the most recently selected export
```

Neighborhood selection preserves the original policy: extract versioned `WP_` hits in BLAST row order, choose one CDS mapping per WP (prefer RefSeq, then a mapping with taxid), retain the first representative of each taxid. Missing taxids use a per-WP fallback and are marked explicitly. Other accession types are reported as unsupported. This samples representatives; it does not enumerate every genomic locus for every protein.

Downloads use **one-based inclusive genomic coordinates**. The left endpoint is clamped to 1 and actual returned length/end is recorded at contig boundaries. Sequences stay in forward genomic orientation regardless of the anchor strand; circular wraparound is not requested.

## Troubleshooting

| Symptom | Action |
| --- | --- |
| `U=3 ... query contains 1` | Supply the intended marker set, or deliberately change U for a different search. |
| No eligible WP/CDS mappings | Inspect `selection.tsv`, accession types and BLAST input. |
| IPG/download errors | Read the per-step report and rerun; verified completed items are retained. |
| No antiSMASH regions | Inspect the completed HTML reports and input sequence. This differs from a crashed job. |
| Prerequisite failure | Read `setup/prerequisites.log`; verify/prepare `DBDIR`. Package-managed MITE directories are not deleted or replaced. |
| Memory pressure | Lower concurrency/CPUs, raise the estimated per-job memory, reduce DB batch size, or use a larger runtime. |
| Region matching error | Inspect `matching_errors.tsv`; check COL and the current antiSMASH manifest. Duplicate basenames are rejected. |
| No passing candidates | Read the saved diagnosis and annotations; change thresholds only with a scientific reason. |

The [NCBI BLAST guidelines](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html) require ≥10 seconds between requests and ≥1 minute between polls of the same RID. The notebook enforces this across reruns within a workspace. An Entrez API key does **not** raise BLAST limits. For larger workloads, use a provisioned local/cloud BLAST workflow rather than increasing public-server polling or running concurrent copies of the notebook against the same service.

## Development and validation

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
python scripts/sync_notebook.py
python -m pytest -q
```

Edit helpers in `workflow_support.py`, then synchronize the embedded cell. Parameter/orchestration cells are edited in the notebook. Do not commit credentials, execution outputs, runtime folders or widget state.

Offline tests cover FASTA validation, RID recovery/rate limits, retries, idempotent downloads, partial failures, cache invalidation, stage handoff with synthetic GenBank/SQLite files, filter boundaries, duplicate domains, exports, notebook schema and preview behavior. Scientific tools are mocked in the stage-handoff test. **This is not full Colab/antiSMASH/database validation or a performance benchmark.** Validate a small dataset with known expected candidates in the intended Linux runtime before scaling up.

See [CHANGELOG.md](CHANGELOG.md) for migration notes. Older output folders without completion manifests are not automatically trusted as caches; old standard BLAST TSVs can be imported.

## Upstream interfaces

- [BLAST URL API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html) and [NCBI IPG XML schema](https://www.ncbi.nlm.nih.gov/data_specs/schema/other/seq_report/IPGReportSet.xsd)
- [antiSMASH installation](https://docs.antismash.secondarymetabolites.org/install/) and [release notes](https://github.com/antismash/antismash/releases)
- [cblaster](https://github.com/gamcil/cblaster) and [clinker](https://github.com/gamcil/clinker)

This repository is licensed under [Apache 2.0](LICENSE). External tools and databases retain their respective licenses.
