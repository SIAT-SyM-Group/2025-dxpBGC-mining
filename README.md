# Genome Mining Workflow (antiSMASH + cblaster + rule-based filtering)
<a href="https://colab.research.google.com/github/SIAT-SyM-Group/2025-dxpBGC-mining/blob/main/workflow.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

A Colab-ready notebook workflow for **BGC (biosynthetic gene cluster) genome mining**. Starting from protein queries, it automates **BLASTP hit retrieval → IPG/nuccore mapping and genomic neighborhood download (±FLANK) → batch antiSMASH annotation → local cblaster database build/search → rule-based filtering → (optional) clinker visualization**, enabling rapid triage of large candidate sets into higher-confidence cluster regions.

## Features
- **Run in Google Colab**: no local conda/antiSMASH required (dependencies/tools are installed by the notebook)
- **End-to-end automation**: BLASTP → neighborhood download → BGC annotation → similarity search → filtering/export
- **Scale control**: taxid-based deduplication, Top-N limits, tunable concurrency to avoid resource blowups
- **Interpretable filtering**: optional diagnostics table to explain why candidates failed and guide threshold tuning
- **Visualization-ready (optional)**: retained GBKs can be fed into clinker to produce interactive HTML views

## Pipeline Overview
1. **Dependencies & NCBI credentials**: set `NCBI_EMAIL` (required) and `NCBI_API_KEY` (recommended)
2. **Query FASTA**: provide protein sequences (raw or FASTA) and write a normalized FASTA
3. **BLASTP (URLAPI)**: run BLASTP against `refseq_protein` (default) and download TSV
4. **IPG/nuccore neighborhood**: map WP hits to genome records and download **±FLANK** neighborhoods (FASTA)
5. **antiSMASH batch annotation**: annotate neighborhoods to produce region GBKs and HTML reports
6. **Local cblaster search**: build a local database from antiSMASH region GBKs and search; export `summary/abspres/session`
7. **Rule-based filter & export**: apply rules and copy passing regions to a destination directory
8. **Diagnostics (optional)**: generate `filter_diagnosis.tsv` when nothing passes
9. **clinker visualization (optional)**: create interactive cluster alignment HTML for retained GBKs

## Requirements
- Recommended: **Google Colab**
- Python packages: `biopython`, `requests`, `tqdm`, etc. (auto-installed)
- Tools via micromamba:
  - **antiSMASH** for BGC prediction/annotation
  - **cblaster** for local gene-cluster similarity searches
  - **clinker (optional)** for visualization

## Quick Start (Colab)
1. Open `workflow.ipynb`
2. In the setup cell, provide:
   - `NCBI_EMAIL` (required)
   - `NCBI_API_KEY` (recommended for stability/throughput)
3. Paste/load protein query sequences and generate the query FASTA
4. Run subsequent cells (recommended: `Runtime → Run all`)
5. Inspect retained region GBKs and (optionally) clinker HTML outputs

## Key Parameters
- **BLAST**: `HIT_TOP_N`, `EXPECT`, `POLL_SEC`, `MAX_WAIT_H`
- **Neighborhood**: `FLANK` (bp), `THREADS`, `TOP_N_UNIQ_TAXID`
- **antiSMASH**: `MAX_JOBS`, `CPUS_PER_JOB`, `DOWNLOAD_DB`
- **cblaster**: `MI` (min identity), `MC` (min coverage)
- **Filtering**: `PROTO_GAP_MAX`, `COPY_MODE`, `COL`, `ROOT`, `DEST`

## Inputs & Outputs
### Inputs
- Protein query sequence(s) (raw/FASTA)
- NCBI Entrez credentials: `NCBI_EMAIL` (required), `NCBI_API_KEY` (recommended)
- BLAST TSV produced by the notebook and passed downstream

### Outputs (default folders/files)
- `wp_flanks_fna/`: downloaded neighborhoods (±FLANK)
- `antismash_out/`: antiSMASH outputs (including `index.html` and region GBKs)
- `final-target/`: retained GBKs/folders after filtering (depends on `COPY_MODE`)
- `final-target/filter_diagnosis.tsv`: (optional) failure diagnostics
- `summary-*.csv`, `abspres-*.csv`, `session-*.json`, `my_plot-*.html`: cblaster outputs

## Troubleshooting
- **No/slow BLAST hits**: reduce `HIT_TOP_N`, adjust `EXPECT`, increase `MAX_WAIT_H` if needed
- **Frequent Entrez/IPG failures**: ensure `NCBI_EMAIL`, strongly recommend `NCBI_API_KEY`, reduce `THREADS`
- **antiSMASH too slow/resource-limited**: lower `MAX_JOBS` or `CPUS_PER_JOB`; set `DOWNLOAD_DB=True` for first run
- **Nothing passes filters**: run Diagnostics to generate `filter_diagnosis.tsv`, then tune `PROTO_GAP_MAX` and domain thresholds
- **Region file matching fails**: confirm `ROOT` points to `antismash_out` (or GBK root) and `COL` matches the region-name column
