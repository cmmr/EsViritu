# Running EsViritu


You might run this as part of a bash script, snakemake pipeline, do your own upstream read processing, etc, but these are the basic instructions.

*Required inputs:*

`-r reads file (.fastq)`

`-s sample name`

`-o output directory (may be shared with other samples)`

Activate the conda environment:

`conda activate EsViritu`

Individual samples can be run with the python script. E.g.:

**Basic run with 1 .fastq file:**

```bash
EsViritu -r /path/to/reads/myreads.fastq -s sample_ABC -o myproject_EsViritu1 -p unpaired
```

**Using paired end input .fastq files. Must be exactly 2 files.**

```bash
EsViritu -r /path/to/reads/myreads.R1.fastq /path/to/reads/myreads.R2.fastq -s sample_ABC -o myproject_EsViritu1 -p paired
```

**With pre-filtering steps:**

```bash
EsViritu -r /path/to/reads/myreads.fastq -s sample_ABC -o myproject_EsViritu1 -q True -f True -p unpaired
```

**Help menu**

```bash
EsViritu -h
```

## Making a summary report for a batch of samples

Run the batch summary script to collate reports from several sequencing libraries in a project:

Example:

Activate conda environment: `conda activate EsViritu`

Then run the `summarize_esv_runs` command with the relative path to the output directory as the only argument:

```bash
summarize_esv_runs myproject_EsViritu1
```

This command will generate the tables `myproject_EsViritu1.detected_virus.info.tsv`, `myproject_EsViritu1.detected_virus.assembly_summary.tsv` and the reactable `myproject_EsViritu1.batch_detected_viruses.html`, which summarize information about all the samples in the given directory.
