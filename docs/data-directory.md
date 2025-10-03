# Output data directory

EsViritu writes its primary tab-delimited outputs to the directory provided with `-o/--output_dir`.
Each file is prefixed with your `--sample` name.

Location: `<OUTPUT_DIR>/`

## Files

- `<SAMPLE>.detected_virus.info.tsv`
  - Per-contig summary table of detected viral contigs.
  - Contains detection and quantification metrics for each contig.

- `<SAMPLE>.detected_virus.assembly_summary.tsv`
  - Per-assembly summary table.
  - Aggregates contig-level results at the assembly level.

- `<SAMPLE>.tax_profile.tsv`
  -Taxonomic profile table.
  - Assigns taxonomy to records based on average nucleotide identity to reference
  - See `--species-threshold` (default 0.90)
  - See `--subspecies-threshold` (default 0.95)

- `<SAMPLE>.virus_coverage_windows.tsv`
  - Coverage in fixed windows across each reference contig.
  - Useful for visualizing coverage profiles and drops/gaps.

## Notes
- Temporary/intermediate files are written under `<OUTPUT_DIR>/<SAMPLE>_temp/` (unless you set `--temp`).
  These include aligned .bam files,read-sharing comparisons, clustering summaries, coverm-like tables, and read ANI per contig.
  They will be removed automatically unless you run with `--keep`.

## Column reference for main tables

### *SAMPLE*.detected_virus.info.tsv

| Column | Type | Description |
|---|---|---|
| sample_ID | string | Sample name provided via `--sample`. |
| Name | string | Virus name from database metadata. |
| description | string | Reference sequence description. |
| Length | int | Reference contig length (bp). |
| Segment | string | Segment identifier (if applicable). |
| Accession | string | Reference contig accession ID. |
| Assembly | string | Assembly name the contig belongs to. |
| Asm_length | int | Total assembly length (bp). |
| kingdom | string | Taxonomic rank. |
| phylum | string | Taxonomic rank. |
| tclass | string | Taxonomic class (named `tclass` to avoid keyword collision). |
| order | string | Taxonomic rank. |
| family | string | Taxonomic rank. |
| genus | string | Taxonomic rank. |
| species | string | Taxonomic rank. |
| subspecies | string | Taxonomic rank. |
| RPKMF | float | Reads per kilobase per million filtered reads: `(read_count/(Length/1000)) / (filtered_reads_in_sample/1e6)`. |
| read_count | int | Number of reads aligned to the contig. |
| covered_bases | int | Number of bases with coverage > 0 on the contig. |
| mean_coverage | float | Mean read depth across the contig. |
| avg_read_identity | float | Average per-read alignment identity for the contig. |
| Pi | float | Average nucleotide diversity across covered positions. |
| filtered_reads_in_sample | int | Total filtered reads used for normalization. |

### *SAMPLE*.detected_virus.assembly_summary.tsv

| Column | Type | Description |
|---|---|---|
| sample_ID | string | Sample name provided via `--sample`. |
| filtered_reads_in_sample | int | Total filtered reads used for normalization. |
| Assembly | string | Assembly name. |
| Asm_length | int | Total assembly length (bp). |
| kingdom | string | Taxonomic rank. |
| phylum | string | Taxonomic rank. |
| tclass | string | Taxonomic class (named `tclass` to avoid keyword collision). |
| order | string | Taxonomic rank. |
| family | string | Taxonomic rank. |
| genus | string | Taxonomic rank. |
| species | string | Taxonomic rank. |
| subspecies | string | Taxonomic rank. |
| read_count | int | Number of reads aligned across contigs in the assembly. |
| covered_bases | int | Number of bases with coverage > 0 across contigs in the assembly. |
| avg_read_identity | float | Mean read identity across contigs in the assembly. |
| Accession | string | Comma-separated list of accessions included. |
| Segment | string | Comma-separated list of segments included. |
| RPKMF | float | Assembly-level RPKMF using `Asm_length`. |

### *SAMPLE*.virus_coverage_windows.tsv

| Column | Type | Description |
|---|---|---|
| Accession | string | Reference accession ID. |
| window_index | int | Window index (0–99). |
| window_start | int | 0-based start coordinate (inclusive). |
| window_end | int | 0-based end coordinate (exclusive). |
| average_coverage | float | Mean read depth across the window. |

### *SAMPLE*.tax_profile.tsv

| Column | Type | Description |
|---|---|---|
| sample_ID | string | Sample name provided via `--sample`. |
| filtered_reads_in_sample | int | Total filtered reads used for normalization. |
| kingdom | string | Taxonomic rank. |
| phylum | string | Taxonomic rank. |
| tclass | string | Taxonomic class (named `tclass` to avoid keyword collision). |
| order | string | Taxonomic rank. |
| family | string | Taxonomic rank. |
| genus | string | Taxonomic rank. |
| species | string | Species classification; may be `s__unclassified <genus>` when `avg_read_identity` < `--species-threshold`. |
| subspecies | string | Subspecies classification; may be `t__unclassified <species>` when `avg_read_identity` < `--subspecies-threshold`. |
| read_count | int | Sum of reads aligned to assemblies contributing to this taxon. |
| RPKMF | float | Sum of RPKMF across assemblies in this taxon. |
| avg_read_identity | float | Mean read identity across assemblies in this taxon. |
| assembly_list | string | Comma-separated list of assemblies contributing to this row. |