# Making Custom Databases

EsViritu ships with a curated database of human, animal and plant virus genomes, but you can build and use a custom database if you need to target a different set of viruses (e.g. environmental viruses, phages, or a custom subset of NCBI sequences). One use case is improving consensus genome assembly by using lots of closely related genomes for a species of interest. See this example for [RSV](https://zenodo.org/records/18157379).

Two helper commands are provided:

- **`esv_create_taxonomy`** — Takes a table of accessions and generates a fully annotated EsViritu-style metadata TSV using NCBI taxonomy.
- **`esv_combine_tax`** — Merges two or more metadata TSVs into a single file, deduplicating records.

A custom database requires three files, both pointing to the same set of accessions:

1. A FASTA file of virus genome sequences (e.g. downloaded from NCBI).
 - `virus_pathogen_database.fna`
2. a minimap2 index (.mmi) with short-read settings.
 - `virus_pathogen_database.mmi`
 - via: `minimap2 -x sr -d virus_pathogen_database.mmi virus_pathogen_database.fna`
3. A metadata TSV produced by `esv_create_taxonomy` (and optionally combined with `esv_combine_tax`).
 - `virus_pathogen_database.all_metadata.tsv`

The directory with the three files are then passed to `EsViritu` via `--db`, respectively.

---

# Making an EsViritu-Style Database with NCBI Records

!!! note
    If you want to use custom records (not in NCBI), you'll have to generate a metadata file with the columns specified in the **Output columns** section.

---

### Install `pytaxonkit`

`esv_create_taxonomy` requires `pytaxonkit`, which is available on Bioconda:

```bash
conda install bioconda::pytaxonkit">=0.10"
```

### Download the NCBI accession-to-TaxID mapping

```bash
wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

zcat nucl_gb.accession2taxid.gz | cut -f 2,3 | gzip -c > acc2taxid.tsv.gz
```

Estimated size: ~900 MB (compressed).

### Download the NCBI taxonomy dump

```bash
mkdir taxdump && cd taxdump

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

tar -xvf taxdump.tar.gz
```

Estimated size: ~3 GB (uncompressed).

---

## Step 1: Prepare your accession table

Create a tab-delimited input TSV with at minimum these required columns:

| Column | Required | Description |
|---|---|---|
| `Accession` | **yes** | NCBI accession with version (e.g. `NC_001401.2`) |
| `Organism_Name` | **yes** | Virus name (becomes the `Name` field in output) |
| `Length` | **yes** | Sequence length in bp |
| `Assembly` | recommended | Assembly identifier for grouping multi-segment viruses (e.g. `GCF_000862125.1`). If absent, each accession gets its own assembly ID (`set:<Accession>`). |
| `Segment` | optional | Segment label (e.g. `S`, `M`, `L`). |
| `description` | optional | Free-text description. Defaults to `Organism_Name` if absent. |

This table is typically exported from the [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) web interface, but could be generated from [NCBI entrez](https://www.ncbi.nlm.nih.gov/books/NBK44864/).

!!! tip
    For segmented viruses, make sure all segments share the same `Assembly` value. `Asm_length` is computed as the sum of `Length` across all accessions with the same `Assembly`.

---

## Step 2: Run `esv_create_taxonomy`

```bash
esv_create_taxonomy \
    -i my_accessions.tsv \
    -o my_virus_metadata.tsv \
    -a /path/to/acc2taxid.tsv.gz \
    -t /path/to/taxdump
```

### Arguments

| Flag | Required | Description |
|---|---|---|
| `-i` / `--input` | **yes** | Input accession TSV (see Step 1). |
| `-o` / `--output` | **yes** | Output metadata TSV path. |
| `-a` / `--accfile` | **yes** | Path to `acc2taxid.tsv.gz` (see Prerequisites). |
| `-t` / `--taxonkit_dir` | **yes** | Path to the taxdump directory (see Prerequisites). |
| `-s` / `--subspecies-label` | no | How to populate the `subspecies` field. See below. |

### `--subspecies-label` options

By default (`-s` not specified), `esv_create_taxonomy` uses `pytaxonkit` to assign the subspecies/strain rank from NCBI taxonomy. This is the recommended option.

Alternatively, you can override the subspecies field using an existing column from your input:

| Value | Source column | Notes |
|---|---|---|
| *(default)* | NCBI taxonomy | Recommended; uses taxonkit `t__` format. |
| `organism` | `Organism_Name` | Useful when organism names encode strain info. |
| `genotype` | `Genotype` | Input must have a `Genotype` column. |
| `subspecies` | `subspecies` | Use a pre-existing subspecies column as-is. |
| `species` | `Species` | Input must have a `Species` column. |


### Output columns

The output TSV contains the following columns (plus any extra columns from the input):

| Column | Description |
|---|---|
| `Accession` | NCBI accession with version |
| `description` | Sequence description |
| `Name` | Virus name (from `Organism_Name`) |
| `Segment` | Segment label (null if not applicable) |
| `kingdom` | Taxonomic rank |
| `phylum` | Taxonomic rank |
| `tclass` | Taxonomic class |
| `order` | Taxonomic rank |
| `family` | Taxonomic rank |
| `genus` | Taxonomic rank |
| `species` | Taxonomic rank |
| `subspecies` | Subspecies/strain label |
| `Length` | Sequence length (bp) |
| `TaxID` | NCBI Taxonomy ID |
| `Assembly` | Assembly identifier |
| `Asm_length` | Total assembly length (sum of `Length` across segments) |

**NOTE:** There will be several extra columns from `pytaxonkit` that are not removed by default as they do not effect `EsViritu` processing.

---

## Step 3 (optional): Combine multiple metadata tables with `esv_combine_tax`

If you have built metadata tables from multiple input TSVs (e.g. different virus groups or NCBI datasets exports), combine them into a single file:

```bash
esv_combine_tax standard_DB/virus_pathogen_database.all_metadata.tsv table1.tsv table2.tsv -o custom_DB/virus_pathogen_database.all_metadata.tsv
```

- Duplicate rows (identical across all columns) are removed.
- Tables must all have the required columns (`Accession`, `description`, `Name`, `Segment`, `kingdom`, `phylum`, `tclass`, `order`, `family`, `genus`, `species`, `subspecies`, `Length`, `TaxID`, `Assembly`, `Asm_length`). Tables missing any of these will cause an error.
- Extra columns beyond the required set are preserved.

### Arguments

| Flag | Default | Description |
|---|---|---|
| `tables` (positional) | — | Two or more metadata TSV files to merge. |
| `-o` / `--output` | `virus_pathogen_database.all_metadata.tsv` | Output file path. |

---

## Step 4: Combine FASTA sequences and Index for minimap2

Combine all the FASTA files that correspond to records in your metadata table. 

```bash
cat standard_DB/virus_pathogen_database.fna custom_seqs.fna | seqkit rmdup | seqkit seq > custom_DB/virus_pathogen_database.fna

minimap2 -x sr -d custom_DB/virus_pathogen_database.mmi custom_DB/virus_pathogen_database.fna
```

!!! note
    The FASTA record IDs must match the `Accession` column in your metadata TSV exactly (including version suffix, e.g. `NC_001401.2`).

    A script to check this is not provided, so please do your due diligence.

---

## Step 5: Run EsViritu with the custom database

Pass your custom FASTA and metadata files via `--db`:

```bash
EsViritu \
    -r /path/to/reads1.fastq /path/to/reads2.fastq \
    -s my_sample \
    -o my_output_dir \
    --db /path/to/custom_DB
```

