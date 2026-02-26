from ntpath import isfile
import polars as pl
import pysam
import logging
import subprocess
import json
import tempfile
from subprocess import Popen, PIPE, STDOUT
import sys, os
import shutil
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from collections import Counter
import statistics
from typing import Dict, Tuple
from collections import defaultdict, Counter
import shutil

import yaml

logger = logging.getLogger("esv_logger")



def str2bool(v):
    '''
    converts bool-like strings to True/False boolean objects
    '''
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def is_tool(name):
    """Check whether `name` is on PATH."""
    
    return shutil.which(name) is not None

def trim_filter(reads: list, outdir: str, tempdir: str, trim: bool, filter: bool, sample_name: str,
                filter_db: str = None, paired: str = "paired", threads: int = 4, mmk: str = "500M",
                mmp: str = "sr", dedup: bool = False) -> list:
    """
    Trim and/or filter .fastq reads using fastp (quality trim) and minimap2+pysam (host/spike-in filter).
    Args:
        reads: str, path to input fastq file(s). If paired, should be space-separated paths.
        outdir: str, output directory for summary files.
        tempdir: str, directory for intermediate files.
        trim: bool, whether to perform quality trimming with fastp.
        filter: bool, whether to filter reads mapping to filter_db with minimap2.
        filter_db: str, path to fasta file for minimap2 filtering (required if filter=True).
        paired: str, "paired" or "unpaired".
        threads: int, number of threads to use.
        mmk: str, minimap2 -K parameter
        mmp: str, minimap2 preset for sequencing tech
        dedup: bool, whether to remove PCR duplicates with fastp --dedup.
    Returns:
        str: path to output fastq file with processed reads.
    """

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(tempdir, exist_ok=True)

    # Set up file names
    trimmed_fastq = os.path.join(tempdir, f"{sample_name}.trimmed.fastq")
    filtered_fastq = os.path.join(tempdir, f"{sample_name}.filtered.fastq")
    fastp_html = os.path.join(tempdir, f"{sample_name}.fastp.html")
    fastp_json = os.path.join(tempdir, f"{sample_name}.fastp.json")


    input_fastq = list(reads)
    if int(threads) > 16:
        fastp_threads = 16
    else:
        fastp_threads = int(threads)
    # Step 1: Quality trim with fastp
    if trim:

        if paired == "paired":
            read1, read2 = reads
            trimmed1 = os.path.join(tempdir, f"{sample_name}_R1.trimmed.fastq")
            trimmed2 = os.path.join(tempdir, f"{sample_name}_R2.trimmed.fastq")
            fastp_cmd = [
                "fastp", "-i", read1, "-I", read2,
                "-o", trimmed1, "-O", trimmed2,
                "-w", str(fastp_threads), "--html", fastp_html, "--json", fastp_json
            ]
            if dedup:
                fastp_cmd.append("--dedup")
            try:
                subprocess.run(fastp_cmd, check=True)
            except Exception as e:
                logger.error(f"fastp paired trimming failed: {fastp_cmd}\nError: {e}")
                raise
            input_fastq = [trimmed1, trimmed2]
        else:
            fastp_cmd = [
                "fastp", "--in1", reads[0], "--out1", trimmed_fastq,
                "-w", str(fastp_threads), "--html", fastp_html, "--json", fastp_json
            ]
            if dedup:
                fastp_cmd.append("--dedup")
            try:
                subprocess.run(fastp_cmd, check=True)
            except Exception as e:
                logger.error(f"fastp single-end trimming failed: {fastp_cmd}\nError: {e}")
                raise
            input_fastq = [trimmed_fastq]
    else:
        logger.info("Not filtering for quality...")

    # Step 2: Filter with minimap2 + pysam
    if filter:
        if filter_db is None:
            raise ValueError("filter_db must be provided if filter=True")
        # Build commands for piping

        if not os.path.isfile(fastp_json):
            pre_filt_fastp_html = os.path.join(tempdir, f"{sample_name}.pre_filt.fastp.html")
            pre_filt_fastp_json = os.path.join(tempdir, f"{sample_name}.pre_filt.fastp.json")
            if paired == "paired":
                read1, read2 = input_fastq
                fastp_stat_paired_cmd = [
                    "fastp", "-i", read1, "-I", read2,
                    "-w", str(fastp_threads), "--html", pre_filt_fastp_html, 
                    "--json", pre_filt_fastp_json
                ]
                try:
                    subprocess.run(fastp_stat_paired_cmd, check=True)
                except Exception as e:
                    logger.error(f"fastp_stat_paired_cmd failed: {fastp_stat_paired_cmd}\nError: {e}")
                    raise
            else:
                fastp_stat_single_cmd = [
                    "fastp", "-i", input_fastq[0],
                    "-w", str(fastp_threads), "--html", pre_filt_fastp_html, 
                    "--json", pre_filt_fastp_json
                ]
                try:
                    subprocess.run(fastp_stat_single_cmd, check=True)
                except Exception as e:
                    logger.error(f"fastp_stat_single_cmd failed: {fastp_stat_single_cmd}\nError: {e}")
                    raise

        if paired == "paired":
            read1, read2 = input_fastq
            minimap2_cmd = [
                "minimap2", "-t", str(threads), "-ax", mmp, "-K", mmk, filter_db, read1, read2
            ]
        else:
            minimap2_cmd = [
                "minimap2", "-t", str(threads), "-ax", mmp, "-K", mmk, filter_db, input_fastq[0]
            ]
        samtools_view_cmd = ["samtools", "view", "-bS", "-"]
        samtools_sort_cmd = ["samtools", "sort", "-"]
        if paired == "paired":
            unmapped_fastq1 = filtered_fastq.replace(".fastq", "_R1.fastq")
            unmapped_fastq2 = filtered_fastq.replace(".fastq", "_R2.fastq")
            # the 12 flag is like 4 + 8 (keep if unmapped and pair is unmapped)
            samtools_fastq_cmd = ["samtools", "fastq", "-f", "12", "-1", unmapped_fastq1, "-2", unmapped_fastq2, "-"]
        else:
            samtools_fastq_cmd = ["samtools", "fastq", "-f", "12", "-o", filtered_fastq, "-"]

        # Pipe together: minimap2 | samtools view | samtools sort | samtools fastq
        p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE)
        p4 = subprocess.Popen(samtools_fastq_cmd, stdin=p3.stdout)
        p1.stdout.close()
        p2.stdout.close()
        p3.stdout.close()
        p4.communicate()

        if paired == "paired":
            output_fastq = [unmapped_fastq1, unmapped_fastq2]
        else: 
            output_fastq = [filtered_fastq]
        return output_fastq
    else:
        logger.info("Not filtering against host/spike-in reference...")
        return input_fastq

## fastp stats
def fastp_stats(reads: list, outdir: str, sample_name: str, trimarg: bool, filtarg: bool,
                tempdir: str, paired: str = "paired", threads: int = 4) -> str:
    """
    Getting stats for EsViritu denominator (total reads into main pipeline). 
    Running fastp if necessary.
    """
    pipeline_stats_fastp_html = os.path.join(outdir, f"{sample_name}.fastp.html")
    pipeline_stats_fastp_json = os.path.join(outdir, f"{sample_name}.fastp.json")

    if int(threads) > 16:
        fastp_threads = 16
    else:
        fastp_threads = int(threads)

    # Determine which existing fastp JSON to reuse, if any:
    # - trim=True, filter=False: reuse fastp.json from quality trimming
    # - trim=False, filter=True: reuse pre_filt.fastp.json from host filtering step
    # - Otherwise: run fastp on final reads to count them
    trim_json = os.path.join(tempdir, f"{sample_name}.fastp.json")
    pre_filt_json = os.path.join(tempdir, f"{sample_name}.pre_filt.fastp.json")
    trim_html = os.path.join(tempdir, f"{sample_name}.fastp.html")
    pre_filt_html = os.path.join(tempdir, f"{sample_name}.pre_filt.fastp.html")
    
    if trimarg and not filtarg and os.path.isfile(trim_json):
        # Reuse stats from quality trimming step
        shutil.copy(trim_json, pipeline_stats_fastp_json)
        if os.path.isfile(trim_html):
            shutil.copy(trim_html, pipeline_stats_fastp_html)
    elif filtarg and not trimarg and os.path.isfile(pre_filt_json):
        # Reuse stats from pre-filter step (avoids redundant fastp call)
        shutil.copy(pre_filt_json, pipeline_stats_fastp_json)
        if os.path.isfile(pre_filt_html):
            shutil.copy(pre_filt_html, pipeline_stats_fastp_html)
    elif paired == "paired":
        read1, read2 = reads
        fastp_stat_paired_cmd = [
            "fastp", "-i", read1, "-I", read2,
            "-w", str(fastp_threads), "--html", pipeline_stats_fastp_html, 
            "--json", pipeline_stats_fastp_json
        ]
        try:
            subprocess.run(fastp_stat_paired_cmd, check=True)
        except Exception as e:
            logger.error(f"fastp_stat_paired_cmd failed: {fastp_stat_paired_cmd}\nError: {e}")
            raise
    else:
        fastp_stat_single_cmd = [
            "fastp", "-i", reads[0],
            "-w", str(fastp_threads), "--html", pipeline_stats_fastp_html, 
            "--json", pipeline_stats_fastp_json
        ]
        try:
            subprocess.run(fastp_stat_single_cmd, check=True)
        except Exception as e:
            logger.error(f"fastp_stat_single_cmd failed: {fastp_stat_single_cmd}\nError: {e}")
            raise
    
    ## get total reads from fastp json
    # The pipeline_stats_fastp_json now contains stats for the final preprocessed reads.
    # - If we copied from trim fastp.json (trim=True, filter=False), use "after_filtering" 
    #   to get post-quality-filter count.
    # - Otherwise, we ran fastp on already-processed reads just to count them, 
    #   so use "before_filtering" (no filtering was applied in this counting run).
    total_reads = None
    copied_from_trim = (
        trimarg and 
        not filtarg and 
        os.path.isfile(os.path.join(tempdir, f"{sample_name}.fastp.json"))
    )
    which_reads = "after_filtering" if copied_from_trim else "before_filtering"
        
    try:
        with open(pipeline_stats_fastp_json, "r") as f:
            data = json.load(f)
            total_reads = data["summary"][which_reads]["total_reads"]
    except Exception as e:
        logger.warning(f"Could not parse total_reads from fastp JSON: {e}")
    return total_reads

def various_readstats(
    outdir: str, sample_name: str, trimarg: bool, 
    filtarg: bool,tempdir: str, filtered_reads: str
    ) -> str:
    """
    compile the stats from the various read filtering stages and return a yaml file.
    Args:
        reads: str, path to input fastq file(s). If paired, should be space-separated paths.
        outdir: str, output directory for summary files.
        tempdir: str, directory for intermediate files.
        trimarg: bool, whether to perform quality trimming with fastp.
        filtarg: bool, whether to filter reads mapping to filter_db with minimap2.
        filtered_reads: str, filtered read count from fastp_stats, used as the EsViritu denominator
    Returns:
        str: path to output yaml file with read stats.
    """
    pre_qual_json = os.path.join(tempdir, f"{sample_name}.fastp.json")
    pre_filt_json = os.path.join(tempdir, f"{sample_name}.pre_filt.fastp.json")
    yaml_readcount_path = os.path.join(outdir, f"{sample_name}_esviritu.readstats.yaml")

    readstats_dict = {}
    
    # Track reads at each filtering stage for full visibility into read attrition
    
    # 1. Pre-quality filter reads (only when trimming was performed)
    if trimarg and os.path.isfile(pre_qual_json):
        try:
            with open(pre_qual_json, "r") as f:
                data = json.load(f)
                readstats_dict["reads pre-quality filter"] = data["summary"]["before_filtering"]["total_reads"]
        except Exception as e:
            logger.warning(f"Could not parse pre-qual reads from fastp JSON: {e}")
    
    # 2. Pre-host removal reads (when host filtering was performed)
    if filtarg:
        if os.path.isfile(pre_filt_json):
            # trim=False, filter=True: pre_filt.fastp.json was created
            try:
                with open(pre_filt_json, "r") as f:
                    data = json.load(f)
                    readstats_dict["reads pre-host removal"] = data["summary"]["before_filtering"]["total_reads"]
            except Exception as e:
                logger.warning(f"Could not parse pre-host filter reads from fastp JSON: {e}")
        elif trimarg and os.path.isfile(pre_qual_json):
            # trim=True, filter=True: use post-trim count from fastp.json as pre-host count
            try:
                with open(pre_qual_json, "r") as f:
                    data = json.load(f)
                    readstats_dict["reads pre-host removal"] = data["summary"]["after_filtering"]["total_reads"]
            except Exception as e:
                logger.warning(f"Could not parse post-trim reads from fastp JSON: {e}")
    try:
        readstats_dict["reads for EsViritu denominator"] = filtered_reads
    except Exception as e:
        logger.warning(f"Could not parse filtered_reads from variable: {e}")

    try:
        with open(yaml_readcount_path, 'w') as yaml_file:
            yaml.dump(readstats_dict, yaml_file, default_flow_style=False)
    except Exception as e:
        logger.error(f"Failed to save readstats parameters YAML: {e}")

    return yaml_readcount_path


def minimap2_f(
    reference: str, reads: list, cpus: str, sorted_outbf, mmk: str = "500M", mmp: str = "sr"
    ) -> str:

    '''
    aligns read pairs to reference with minimap2, passes it to pysam,
    then filter low quality alignments.
    Args:
        reference: reference genome index.
        reads: list of fastq read file(s)
        cpus: number of CPUs to use for minimap2
        sorted_outbf: sorted bam file name
        mmk: -K parameter for minimap2
        mmp: str, minimap2 preset for sequencing tech
    Returns:
        sorted_outbf
    '''
    
    mini2_command = [
        'minimap2', '-t', cpus, 
        '-ax', mmp, 
        '--secondary=yes', 
        '--secondary-seq',
        '--sam-hit-only', "--MD",
        '-f', '10000',
        '-N' '100',
        '-p', '0.90',
        '-K', mmk,
        reference, 
        *reads
        ]

    # Launch process, to stdout and catch it with pysam
    with subprocess.Popen(
        mini2_command,
        stdout=subprocess.PIPE,
        #stderr=STDOUT
        ) as mini2_proc:

        fd_child = mini2_proc.stdout.fileno()
    
        # read minimap2 output with pysam

        with pysam.AlignmentFile(fd_child, 'r') as ministream:

            outbf = sorted_outbf.replace('.sorted.bam', '.bam')
            filtbam = pysam.AlignmentFile(outbf, "wb", template=ministream)

            for record in ministream:
                readLength: int = record.infer_read_length()
                alignLength: int = record.query_alignment_length

                # NM flag reports "edit distance" between read and ref
                try:
                    alignEdit: int=dict(record.get_tags(with_value_type=False)).get('NM')
                except:
                    alignEdit: int=100
                try:
                    alignProp: float=alignLength/readLength
                except:
                    alignProp: float=0
                # assume the edit distance is between aligned part of read and ref
                try:
                    alignAcc: float=(alignLength-alignEdit)/alignLength
                except:
                    alignAcc: float=0

                if alignLength >= 100 and alignProp >= 0.9 and alignAcc >= 0.8:

                    filtbam.write(record)

    filtbam.close()

    # make a sorted bam file

    pysam.samtools.sort("-@", str(cpus), "-o", sorted_outbf, outbf)

    pysam.index(sorted_outbf)

    # return an AlignmentFile
    return sorted_outbf


## take filtered .bam, make a coverm-like file
def _calculate_contig_stats_from_bam(bamfile: pysam.AlignmentFile, contig: str,
                                     include_secondary: bool = False) -> Dict:
    """Calculate statistics for a single contig using an open BAM handle."""
    length = bamfile.get_reference_length(contig)
    if length == 0:
        return None

    # Initialize coverage array
    coverage = np.zeros(length, dtype=np.uint32)
    pi_list = []

    # Get all reads for this contig
    reads = list(bamfile.fetch(contig))
    # only counting primary alignments
    primary_reads = [r for r in reads if not r.is_secondary and not r.is_unmapped]

    if not primary_reads:
        return None

    # Calculate coverage and base counts in one pass
    base_counts = [Counter() for _ in range(length)]

    for read in primary_reads:
        # Update coverage
        for block in read.get_blocks():
            coverage[block[0]:block[1]] += 1

        # Update base counts at each position
        if not read.is_unmapped and not read.is_secondary or include_secondary:
            ref_pos = read.reference_start
            for qpos, ref_pos, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                if qpos is not None and ref_pos is not None and ref_pos < length:
                    base = read.query_sequence[qpos].upper()
            if base in 'ACGTN':
                base_counts[ref_pos][base] += 1

    # Calculate Pi for each position with coverage
    for pos in range(length):
        if coverage[pos] > 0 and base_counts[pos]:
            counts = base_counts[pos]
            total = sum(counts.values())
            max_count = max(counts.values())
            pi = (total - max_count) / total if total > 0 else 0
            pi_list.append(pi)

    # Calculate statistics
    covered_bases = np.count_nonzero(coverage)
    mean_cov = np.mean(coverage) if length > 0 else 0
    read_count = len(primary_reads)
    avg_pi = round(statistics.mean(pi_list), 3) if pi_list else 0

    return {
        "Accession": contig,
        "contig_length": length,
        "covered_bases": int(covered_bases),
        "read_count": read_count,
        "mean_coverage": float(mean_cov),
        "Pi": avg_pi
    }


def calculate_contig_stats(bam_path: str, contig: str, include_secondary: bool = False) -> Dict:
    """Calculate statistics for a single contig."""
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        return _calculate_contig_stats_from_bam(bamfile, contig, include_secondary)


def _calculate_contig_stats_batch(args) -> list:
    """Process a batch of contigs with a single BAM handle for lower I/O overhead."""
    bam_path, contigs, include_secondary = args
    records = []
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        for contig in contigs:
            result = _calculate_contig_stats_from_bam(bamfile, contig, include_secondary)
            if result is not None:
                records.append(result)
    return records

def bam_to_coverm_table(bam_path: str, sample: str, max_workers: int, include_secondary: bool = False) -> pl.DataFrame:
    """
    Optimized version of bam_to_coverm_table that processes contigs in parallel.
    
    Args:
        bam_path: Path to the input BAM file
        sample: Sample name
        include_secondary: Whether to include secondary alignments
        max_workers: Maximum number of worker processes to use
        
    Returns:
        polars.DataFrame with coverage and diversity statistics
    """
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        # Get contigs with reads
        contigs_with_reads = set()
        for read in bamfile.fetch(until_eof=True):
            if not read.is_unmapped:
                contigs_with_reads.add(bamfile.get_reference_name(read.reference_id))

        contigs_with_reads = list(contigs_with_reads)
        if not contigs_with_reads:
            return pl.DataFrame([])

        # Process contigs in parallel with batched work to reduce per-contig overhead
        max_workers = max(1, int(max_workers))
        num_workers = min(max_workers, len(contigs_with_reads))
        batch_size = max(1, len(contigs_with_reads) // (num_workers * 4))
        batches = [
            contigs_with_reads[i:i + batch_size]
            for i in range(0, len(contigs_with_reads), batch_size)
        ]
        args_iter = [(bam_path, batch, include_secondary) for batch in batches]

        records = []
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            for batch_records in executor.map(_calculate_contig_stats_batch, args_iter):
                records.extend(batch_records)

    # Add sample name and convert to DataFrame
    for record in records:
        record["sample"] = sample
        
    return pl.DataFrame(records)

## take filtered .bam, make preliminary consensus .fastas
def bam_to_consensus_fasta(bam_path: str, output_fasta: str = None) -> str:
    """
    Calls samtools consensus on the given sorted BAM file and writes the consensus FASTA.
    Only outputs consensus for records with aligned reads.
    Appends "_consensus" to each FASTA record header.
    Returns the path to the consensus FASTA file.
    """

    if output_fasta is None:
        output_fasta = bam_path.replace('.bam', '.consensus.fasta')
    
    # Create a temporary output file for samtools consensus
    temp_output = output_fasta + '.temp'
    
    # samtools consensus command
    cmd = [
        'samtools', 'consensus',
        '-m', 'simple', 
        '-q', '-c', '0.51',
        '-o', temp_output,
        bam_path
    ]
    try:
        subprocess.run(cmd, check=True)
        
        # Process the FASTA file to append "_consensus" to each header
        with open(temp_output, 'r') as infile, open(output_fasta, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    # Append "_consensus" to the header line
                    outfile.write(f"{line.rstrip()}_consensus\n")
                else:
                    outfile.write(line)
                    
        # Remove temporary file
        os.remove(temp_output)
        
    except Exception as e:
        logger.error(f"samtools consensus failed: {cmd}\nError: {e}")
        raise
    return output_fasta

## select final records based on aniclust

def clust_record_getter(clust_tsv: str, field: int, vir_info_df: pl.DataFrame, db_fasta: str, output_fasta: str) -> str:
    """
    Select final records based on aniclust clusters.
    Args:
        clust_tsv: path to aniclust output tsv file (first field is assembly)
        field: which field has exemplar name
        vir_info_df: polars DataFrame with at least columns 'Assembly' and 'Accession'
        db_fasta: path to full fasta file
        output_fasta: path to write filtered fasta (optional)
    Returns:
        Path to new fasta file with selected accessions.
    """

    # Step 1: Get assemblies from aniclust tsv
    assemblies = set()
    with open(clust_tsv, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields:
                assemblies.add(fields[field])
    # Step 2: Get all rows from vir_info_df with those assemblies
    filtered_df = vir_info_df.filter(pl.col('Assembly').is_in(list(assemblies)))
    # Step 3: Get unique accessions
    unique_accessions = filtered_df['Accession'].unique().to_list()
    # Step 4: Write accessions to temp file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as acc_file:
        acc_file.write('\n'.join(unique_accessions) + '\n')
        acc_file.flush()
        acc_file_path = acc_file.name
    # Step 5: Run seqkit grep
    seqkit_cmd = ["seqkit", "grep", "-f", acc_file_path, db_fasta, "-o", output_fasta]
    try:
        subprocess.run(seqkit_cmd, check=True)
    except Exception as e:
        logger.error(f"seqkit_cmd failed: {seqkit_cmd}\nError: {e}")
        raise
    os.remove(acc_file_path)
    return output_fasta


## Make the main table and assembly-based output table
def assembly_table_maker(
    df1: pl.DataFrame, df2: pl.DataFrame, read_ani_df: pl.DataFrame,
    filtered_reads: int, sample: str) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """
    Merge three polars DataFrames on 'Accession' and add an 'RPKMF' column.
    RPKMF = (read_count / (contig_length / 1000)) / (filtered_reads / 1e6)
    Args:
        df1: polars DataFrame (must have 'Accession', 'read_count', 'contig_length')
        df2: polars DataFrame (must have 'Accession' and any columns to merge)
        read_ani_df: polars DataFrame (must have 'Accession' and any columns to merge)
        filtered_reads: total filtered reads (int)
        sample: sample_ID for run
    Returns:
        tuple of DataFrames, one summarized by accession, one by assembly.
    """

    # pare to necessary columns
    df1 = df1.select(["Accession", "covered_bases", "read_count", "mean_coverage", "Pi"])
    merged = df1.join(df2, on="Accession", how="right")
    # fill 0s
    merged = merged.with_columns(
        pl.col("read_count").fill_null(strategy="zero"),
        pl.col("covered_bases").fill_null(strategy="zero"),
        pl.col("mean_coverage").fill_null(strategy="zero"),
        pl.col("Pi").fill_null(strategy="zero")
        )
    # calculate/add more columns
    merged = merged.with_columns([
        (pl.col("read_count") / (pl.col("Length") / 1000) / (filtered_reads / 1e6)).alias("RPKMF"),
        pl.lit(filtered_reads).alias("filtered_reads_in_sample"),
        pl.lit(sample).alias("sample_ID")
    ])

    merged2 = merged.join(
        read_ani_df, on="Accession", how="left"
    ).select([
        "sample_ID", "Name", "description", "Length",
        "Segment", "Accession", "Assembly",
        "Asm_length", "kingdom", "phylum", "tclass", "order",
        "family", "genus", "species", "subspecies", "RPKMF",
        "read_count", "covered_bases", "mean_coverage", 
        "avg_read_identity", "Pi", "filtered_reads_in_sample"
    ]).sort([
        "family", "genus", "species", "Assembly", "Segment"
    ])
    # for the output "main" table we only want records with 1 or more reads
    merged_out = merged2.filter(pl.col("read_count") >= 1)

    # Group/summarize by assembly (for segmented viruses)
    assem_df = merged2.group_by(
        ["sample_ID", "filtered_reads_in_sample", "Assembly", 
        "Asm_length", "kingdom", "phylum", "tclass", "order", 
        "family", "genus", "species", "subspecies"] 
    ).agg(
        # read_count
        pl.col("read_count").sum().alias("read_count"),
        #covered_bases
        pl.col("covered_bases").sum().alias("covered_bases"),
        # read ANI
        pl.col("avg_read_identity").mean().alias("avg_read_identity"),
        #Accessions
        pl.col("Accession").flatten(),
        #Segments
        pl.col("Segment").flatten(),
    ).with_columns(
        (pl.col("read_count") / (pl.col("Asm_length") / 1000) / (filtered_reads / 1e6)).alias("RPKMF"),
        pl.col("Accession").list.eval(pl.element().cast(pl.String)),
        pl.col("Segment").list.eval(pl.element().cast(pl.String))
    ).with_columns(
        pl.col("Accession").list.join(","),
        pl.col("Segment").list.join(",")
    ).filter(pl.col("read_count") >= 1)


    return merged_out, assem_df

## recapitulate coverage table without bedtools
def _coverage_windows_batch(args) -> list:
    """Compute coverage windows for a batch of contigs with one BAM handle."""
    bam_path, contigs = args
    records = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for contig in contigs:
            contig_len = bam.get_reference_length(contig)
            window_size = contig_len // 100
            # If contig is shorter than 100bp, make window_size at least 1
            window_size = max(1, window_size)
            for i in range(100):
                start = i * window_size
                # Last window goes to end
                end = contig_len if i == 99 else min((i + 1) * window_size, contig_len)
                if start >= contig_len:
                    break
                # Get coverage for this window using pileup
                window_length = end - start
                coverage = [0] * window_length
                for pileupcolumn in bam.pileup(contig=contig, max_depth=100_000, start=start, end=end, stepper="all"):
                    pos_in_window = pileupcolumn.pos - start
                    if 0 <= pos_in_window < window_length:
                        coverage[pos_in_window] = pileupcolumn.nsegments
                covered_bases = sum(1 for d in coverage if d > 0)
                if covered_bases == 0:
                    avg_cov = 0.0
                else:
                    avg_cov = float(sum(coverage)) / window_length if window_length > 0 else 0.0
                records.append({
                    "Accession": contig,
                    "window_index": i,
                    "window_start": start,
                    "window_end": end,
                    "average_coverage": avg_cov
                })
    return records


def bam_coverage_windows(bam_path: str, max_workers: int = 1) -> pl.DataFrame:
    """
    Splits each reference contig in a BAM file into 100 equal windows and reports the average read coverage of each window.
    Returns a polars DataFrame with columns: contig, window_index, window_start, window_end, average_coverage
    """

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        contigs = list(bam.references)
    total_contigs = len(contigs)
    if total_contigs == 0:
        return pl.DataFrame([])

    max_workers = max(1, int(max_workers))
    if max_workers == 1 or total_contigs == 1:
        records = _coverage_windows_batch((bam_path, contigs))
        return pl.DataFrame(records)

    num_workers = min(max_workers, total_contigs)
    batch_size = max(1, total_contigs // (num_workers * 4))
    batches = [
        contigs[i:i + batch_size]
        for i in range(0, total_contigs, batch_size)
    ]
    args_iter = [(bam_path, batch) for batch in batches]

    records = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for batch, batch_records in zip(batches, executor.map(_coverage_windows_batch, args_iter)):
            records.extend(batch_records)

    return pl.DataFrame(records)

## read-based clustering function

def assembly_read_sharing_table(bam_path, meta_df: 'pl.DataFrame') -> pl.DataFrame:
    """
    For a given BAM file and meta_df (with 'Accession' and 'Assembly'),
    filter out contigs with no aligned reads, group contigs by Assembly, and compare each pair of contigs within each Assembly
    based on how many read IDs they share (primary or secondary alignments).
    Returns a polars DataFrame with Assembly, contig pairs, number of shared reads, share of total reads for each contig,
    and the average read identity for each assembly.
    """
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    # Step 1: Collect read IDs and read identities for each contig
    contig_reads = defaultdict(set)
    contig_read_identities = defaultdict(list)
    for read in bamfile.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        contig = bamfile.get_reference_name(read.reference_id)
        contig_reads[contig].add(read.query_name)
        # Calculate read identity if possible
        try:
            align_len = read.query_alignment_length
            nm = read.get_tag('NM') if read.has_tag('NM') else None
            if nm is not None and align_len > 0:
                identity = (align_len - nm) / align_len
                contig_read_identities[contig].append(identity)
        except Exception:
            pass
        # Also add for secondary alignments
        if read.has_tag("SA"):
            contig_reads[contig].add(read.query_name)
    bamfile.close()

    # Step 2: Filter out contigs with no reads
    contig_reads = {k: v for k, v in contig_reads.items() if v}
    contig_read_identities = {k: v for k, v in contig_read_identities.items() if v}
    if not contig_reads:
        return pl.DataFrame([])

    # Step 3: Map contigs to assemblies
    if not set(["Accession", "Assembly"]).issubset(set(meta_df.columns)):
        raise ValueError("meta_df must have columns 'Accession' and 'Assembly'")
    contig_to_assembly = dict(zip(meta_df["Accession"].to_list(), meta_df["Assembly"].to_list()))

    # Group contigs by assembly
    assembly_to_contigs = defaultdict(list)
    for contig in contig_reads:
        assembly = contig_to_assembly.get(contig, None)
        if assembly is not None:
            assembly_to_contigs[assembly].append(contig)

    # Step 4: Build assembly-level read sets and average identities
    assembly_reads = {}
    assembly_identities = {}
    for assembly, contigs in assembly_to_contigs.items():
        reads = set()
        identities = []
        for contig in contigs:
            reads.update(contig_reads[contig])
            identities.extend(contig_read_identities.get(contig, []))
        assembly_reads[assembly] = reads
        # Calculate average read identity for the assembly
        if identities:
            avg_identity = float(np.mean(identities))
        else:
            avg_identity = float('nan')
        assembly_identities[assembly] = avg_identity

    # Step 5: Build inverted index (read -> assemblies) and count shared reads per pair
    read_to_assemblies = defaultdict(list)
    for asm, reads in assembly_reads.items():
        for r in reads:
            read_to_assemblies[r].append(asm)

    # Count shared reads for each assembly pair (including self-self)
    pair_shared = Counter()
    for asms in read_to_assemblies.values():
        if len(asms) == 1:
            # self-self only
            pair_shared[(asms[0], asms[0])] += 1
        else:
            asms_sorted = sorted(asms)
            for i in range(len(asms_sorted)):
                for j in range(i, len(asms_sorted)):
                    pair_shared[(asms_sorted[i], asms_sorted[j])] += 1

    # Build output records from counted pairs
    records = []
    for (a, b), n_shared in pair_shared.items():
        tot_a = len(assembly_reads[a])
        tot_b = len(assembly_reads[b])
        share_a = n_shared / tot_a if tot_a else 0.0
        share_b = n_shared / tot_b if tot_b else 0.0
        avg_id_a = assembly_identities.get(a, float('nan'))
        avg_id_b = assembly_identities.get(b, float('nan'))
        records.append({
            "assembly_a": a,
            "assembly_b": b,
            "total_reads_a": tot_a,
            "total_reads_b": tot_b,
            "n_shared_reads": n_shared,
            "share_of_a": share_a,
            "share_of_b": share_b,
            "avg_read_identity_a": avg_id_a,
            "avg_read_identity_b": avg_id_b
        })
    return pl.DataFrame(records)


## cluster contigs on read sharing

def cluster_assemblies_by_read_sharing(df, threshold=0.33) -> pl.DataFrame:
    """
    Cluster assemblies that share at least `threshold` (default 33%) of their reads (based on the assembly with fewer reads).
    For each cluster, choose an exemplar assembly with the highest number of reads aligned.
    Only assemblies with direct connections to the exemplar are included in each cluster.
    Assemblies without direct connections to the exemplar are reclustered.
    Returns a polars DataFrame mapping each assembly to its cluster and exemplar.
    """

    # Build graph representation with direct connections
    direct_connections = defaultdict(set)
    
    # Map assemblies to their total reads
    total_reads = {}

    # Precompute per-assembly max avg_read_identity once (same semantics as prior get_avg_identity).
    has_avg = 'avg_read_identity_a' in df.columns
    max_avg_identity = {}
    def _update_max_identity(asm, val):
        if val is None:
            return
        cur = max_avg_identity.get(asm)
        if cur is None:
            max_avg_identity[asm] = val
        elif val > cur:
            max_avg_identity[asm] = val
    
    # Process edges and collect read counts
    for row in df.iter_rows(named=True):
        n_shared = row['n_shared_reads']
        min_reads = min(row['total_reads_a'], row['total_reads_b'])
        total_reads[row['assembly_a']] = row['total_reads_a']
        total_reads[row['assembly_b']] = row['total_reads_b']
        if has_avg:
            _update_max_identity(row['assembly_a'], row.get('avg_read_identity_a'))
        
        if min_reads == 0:
            continue
            
        share = n_shared / min_reads
        if share >= threshold:
            # Store direct connections in both directions
            direct_connections[row['assembly_a']].add(row['assembly_b'])
            direct_connections[row['assembly_b']].add(row['assembly_a'])

    # Second pass to match prior ordering: all assembly_a identities first, then assembly_b.
    if has_avg:
        for row in df.iter_rows(named=True):
            _update_max_identity(row['assembly_b'], row.get('avg_read_identity_b'))
    
    # Get all assemblies
    all_assemblies = set(df['assembly_a'].to_list()) | set(df['assembly_b'].to_list())
    
    # Track which assemblies have been assigned to clusters
    assigned = set()
    # Store final cluster assignments
    records = []
    # Generate unique cluster IDs
    next_cluster_id = 1
    
    # Keep trying to cluster until no assemblies remain or no more clusters can be formed
    remaining = list(all_assemblies)
    while remaining:
        # Sort remaining assemblies by read count (descending), breaking ties by max avg_read_identity
        # (descending), then by assembly name (ascending / forward alphabetical) for determinism.
        remaining.sort(
            key=lambda a: (-total_reads.get(a, 0), -max_avg_identity.get(a, float('-inf')), a)
        )
        
        # Try to create a new cluster with the assembly having most reads (break ties by avg_read_identity) as exemplar
        exemplar_candidate = remaining[0]
        
        # Find all assemblies directly connected to the exemplar
        directly_connected = direct_connections[exemplar_candidate]
        
        # If exemplar has connections, create a cluster
        if directly_connected:
            cluster_id = f"cluster_{next_cluster_id}"
            next_cluster_id += 1
            
            # Add exemplar to cluster
            records.append({
                'cluster_id': cluster_id,
                'assembly': exemplar_candidate,
                'exemplar': exemplar_candidate,
                'total_reads': total_reads.get(exemplar_candidate, 0),
                'is_exemplar': True
            })
            assigned.add(exemplar_candidate)
            
            # Add directly connected assemblies to cluster
            for assembly in directly_connected:
                if assembly not in assigned:  # Only add if not already in a cluster
                    records.append({
                        'cluster_id': cluster_id,
                        'assembly': assembly,
                        'exemplar': exemplar_candidate,
                        'total_reads': total_reads.get(assembly, 0),
                        'is_exemplar': False
                    })
                    assigned.add(assembly)
        else:
            # If no connections, this assembly forms a singleton cluster
            cluster_id = f"cluster_{next_cluster_id}"
            next_cluster_id += 1
            records.append({
                'cluster_id': cluster_id,
                'assembly': exemplar_candidate,
                'exemplar': exemplar_candidate,
                'total_reads': total_reads.get(exemplar_candidate, 0),
                'is_exemplar': True
            })
            assigned.add(exemplar_candidate)
        
        # Update list of remaining assemblies
        remaining = [a for a in all_assemblies if a not in assigned]
    
    return pl.DataFrame(records)

## read ANI function

def read_ani_from_bam(bam_path):
    """
    Calculate average read identity (ANI) and alignment length for each contig in a BAM file.
    Returns a polars DataFrame with columns: Accession, avg_read_identity, avg_align_len, n_reads
    """

    contig_identity_sum = defaultdict(float)
    contig_align_len_sum = defaultdict(int)
    contig_read_count = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            # only counting primary alignments
            if read.is_unmapped or read.is_secondary or read.reference_id is None:
                continue
            contig = bam.get_reference_name(read.reference_id)
            align_len = read.query_alignment_length
            nm = read.get_tag("NM") if read.has_tag("NM") else None
            if nm is not None and align_len > 0:
                identity = (align_len - nm) / align_len
            else:
                # fallback: identity unknown, skip
                continue
            contig_identity_sum[contig] += identity
            contig_align_len_sum[contig] += align_len
            contig_read_count[contig] += 1

    records = []
    for contig in contig_read_count:
        n = contig_read_count[contig]
        avg_identity = contig_identity_sum[contig] / n if n else 0.0
        avg_align_len = contig_align_len_sum[contig] / n if n else 0.0
        records.append({
            "Accession": contig,
            "avg_read_identity": avg_identity,
            "avg_align_len": avg_align_len,
            "n_reads": n
        })
    return pl.DataFrame(records)

## check for empty bam

def bam_has_alignments(bam_path):
    """
    Check if a BAM file has any aligned reads using pysam.
    Returns True if there are aligned reads, False otherwise.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for _ in bam.fetch(until_eof=True):
            return True
    return False


## get .fastq from aligned reads in bam

def bam_to_paired_fastq(bam_path, fastq_path) -> list:
    """
    Extracts paired-end reads from a BAM file and writes them to two FASTQ files (R1 and R2)
    using samtools fastq. Keeps all pairs where at least one mate is aligned.

    Args:
        bam_path (str): Path to input BAM file.
        fastq_r1_path (str): Path to output FASTQ file for read 1 (R1).
        fastq_r2_path (str): Path to output FASTQ file for read 2 (R2).
    Returns:
        list: [fastq_path]
    """

    cmd = [
        "samtools", "fastq",
        "-o", fastq_path,
        bam_path
    ]
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        logger.error(f"samtools fastq failed: {cmd}\nError: {e}")
        raise
    return [fastq_path]

## make subspecies tax profile
def tax_profile(assem_df: pl.DataFrame, sp_cutoff: float = 0.90, subsp_cutoff: float = 0.95) -> pl.DataFrame:
    """
    Create taxonomic profile with cutoffs at both subspecies and species levels from assembly dataframe.
    Records with avg_read_identity % < sp_cutoff are labeled as "unclassified" 
    at the subspecies level.
    Records with avg_read_identity % < subsp_cutoff are labeled as "unclassified" 
    at the species level.

    Args:
        assem_df: pl.DataFrame, Assembly dataframe from assembly_table_maker (second return value)
        sp_cutoff: float, average read ANI cutoff for species-level classification (default 90.0)
        subsp_cutoff: float, average read ANI cutoff for subspecies-level classification (default 95.0)
        
    Returns:
        pl.DataFrame: tax_df DataFrame
    """

    # Create taxonomical classification based on read identity threshold
    adjust_df = assem_df.with_columns(
        pl.when(pl.col("avg_read_identity") < sp_cutoff)
        .then(
            pl.lit("s__unclassified_") + 
            pl.col("genus").str.replace("^g__", "").str.replace("^unclassified_", "")
        )
        .otherwise(pl.col("species"))
        .alias("species")
    ).with_columns(
        pl.when(pl.col("avg_read_identity") < subsp_cutoff)
        .then(
            pl.lit("t__unclassified_") + 
            pl.col("species").str.replace("^s__", "").str.replace("^unclassified_", "")
        )
        .otherwise(pl.col("subspecies"))
        .alias("subspecies")
    )
    
    # Group by taxonomic levels including the new subspecies classification
    tax_df = adjust_df.group_by([
        "sample_ID", "filtered_reads_in_sample", "kingdom", "phylum", 
        "tclass", "order", "family", "genus", "species", "subspecies"
    ]).agg([
        pl.col("read_count").sum().alias("read_count"),
        pl.col("RPKMF").sum().alias("RPKMF"),
        pl.col("avg_read_identity").mean().alias("avg_read_identity"),
        pl.col("Assembly").unique().alias("assembly_list")
    ]).sort([
        "family", "genus", "species", "subspecies"
    ]).with_columns(
        pl.col("assembly_list").list.eval(pl.element().cast(pl.String))
    ).with_columns(
        pl.col("assembly_list").list.join(",")
    )
    
    return tax_df


def print_esviritu_banner():
    print(r"""
 _____   _____ __      __ _____  _____   _____  _______ _    _ 
| ____| / ____|\ \    / /|_   _||  __ \ |_   _||__   __| |  | |
| |__  | (_     \ \  / /   | |  | |__) |  | |     | |  | |  | |
|  __|  \___ \   \ \/ /    | |  |  _  /   | |     | |  | |  | |
| |___  ____) |   \  /    _| |_ | | \ \  _| |_    | |  | |__| |
|_____| |____/     \/    |_____||_|  \_\|_____|   |_|   \____/ 
                                                                 
""")

def ghost_banner():
    print(
        r"""

                                                     ░░▒▓▓▓▒░░▒▓▒░                                  
                                               ░▒▒▓▒▒░░░   ░     ░▓▒░                               
                                          ░▒▓▓▒░░░░░░░░  ▒    ░    ░▒▓░                             
                    ░░░░░░░░         ░▒▓▓▒░░░░░░░░░░░  ▒       ░      ▒▓▒                           
                    ░░░░░░░░░   ░▒▓▓▒░░░░░░░░░░░░░░ ░▒         ░░       ░▓▓░                        
                   ░░░░░░░░░░▒▓▓▒░░░░░░░░░░░░░░░░ ░▒            ▒         ░▒▒░                      
                    ░░░░░░▒▓░ ░░░░░░░░░░░░░░░░  ▒▒              ░▒           ▒▓▒                    
                     ░░░░▓░  ░░░░░░░░░░░░░░░  ▒▒                 ░░            ▒▓░                  
                        ▒▒░  ░░░░░░░░░░░░░  ▒░                    ▓              ░▓▒                
                        ▓░░░ ░░░░░░░░░░░  ▒░                      ░░               ░▓▒              
                       ░▓░░▒ ░░░░░░░░░  ▒░                         ▒                 ▒▓░            
                       ▒▒░░▒░░░░░░░   ▒░                            ▒                  ▒▓░ ░░░░░░░░ 
                       ▒▒░░░░ ░░░   ░░                              ░░                   ▓▒░░░░░░░░░
                       ▓░░░░░░░ ░░▒                                  ▒                    ░▓▒░░░░░░ 
                      ░▓░░░░░░░░░░░░                                  ░                     ▒▒░░░░░  
                      ░▒░░░░░░░░░░░░                                  ░                 ░░  ▓░░░░░  
                      ▒▒░░░░░░░░░░░░▒▒▒░░                              ░            ▒▒░    ░▓       
                      ▒░░░░░░░░░░░        ░░░▒▒▒▒▒░░                ░░░░░░░░    ░▒▒░       ▒▓       
                     ░▓░░░░░░░░░░                       ░░▒▒▒░     ░░░░░░░░░░▒▒        ▓▒  ▒▒       
                     ░▓░░░░░▒░░░░                                  ░░░░░░░░░░    ▒▒   ░▓▓ ░▓░       
                     ▒▒░░░░░░░░░░░                                   ░░░░░░░     ▓▓░  ▒▓▓ ▒▒        
                     ▓▒░░░░▒░░░░░▒                                 ░  ░░░░░     ░▓▓░  ░▓▓ ▒▒        
                    ░▓░░░░░▒ ░░░░▒░                               ▒        ▒    ░▓▓░  ░▓▒░▓░        
                    ░▓░░░░▒░ ░░░░░▒                     ░░░░░░  ░▒         ░     ▓▓░   ▒░▒▒         
                ░░  ▒▓░░░░▒ ░░░░░░▒░                  ░░░▒▒    ░░           ░    ▓▓      ▓▒         
             ░░░░░░░▒▒░░░▒  ░░░░░░░▒                ░░░▒▒     ▒             ▒░          ░▓░         
            ░░░░░░░░▓░░░▒░ ░░░░░░░░▒░          ░░░░░░▒▒      ▒░              ▒          ▒▒          
             ░░░░░░▓░░░░░ ░░░░░░░░░░▒             ░░       ░▒                ░░         ▒░          
              ░░░░▓▒░░░░ ░░░░░░░░░░░▒                     ░▒              ░   ▒░       ░▓           
                ░▓▒░░░░ ░░░░░░░░░░░░░░                   ░▒       ░       ░    ▒       ▒▒           
                ▓▒░░▒  ░░░░░░░░░▒░░░░░             ░    ▒░      ░░░      ░░░   ░░      ▓░           
               ▓▒░░▒  ░░░░░░░░░░░░░░░░░      ░░░▒░     ▓       ░░▒       ░░▒    ░░    ░▓            
             ░▓░░░░ ░░░░░░░░░░░▒ ░░░░░░░    ░▒▒       ░      ░░░░       ░░░░     ▒    ▒▒            
           ░▓▓░░░░░░░░░░░░░░░░▒ ░░░░░░░░░  ░▒       ░░      ░░░▒        ░░░       ░  ░▓░            
         ░▓▒░░░░░░░░░░░░░░░░░▒  ░░░░░░░░░  ░       ░       ░░░▒        ░░▒           ▒▒             
      ░▒▓░░░░░░░░░░░░░░░░░▒▒  ░░░░░░░░░░░                 ░░░▒░       ░░▒    ░░░▒▒░  ▒▒░░░░         
     ▓▒░░░░░░░░░░░░░░░░▒▒   ░░░░░░░░░░░░░░                ░░▒░    ░░░░░▒            ░▓░░░░░░        
   ▒▓░░░░░░░░░░░░░░░▒░   ░░░░░░░░░░░░░░░░                ░░▒░       ░░▒            ▒▒░░░░░░         
  ░▓░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▒         ░░░░    ░░▒░       ░░▒░           ▓▒░░░░░░          
  ░▓░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▒░                 ░░▒░       ░░▒░          ░▓▒  ░░░░           
   ▓▒░░░░░░░░░░░░░░░░░░░░░░░░░░░░▒▓░                  ░░▓        ░░░▒           ▓░                  
    ▓▒░░░░░░░░░▒▓▓▒░░░░░░░░░░▒▒░                     ░░▒        ░░░▒          ░▓░                   
      ▒▓▓▒▓▓▓▒░░░░░▓░                              ░░▒░        ░░░░░         ░▓░                    
         ░▓░░░░░░░▒▓                              ░░▒          ░░░░         ░▓░                     
          ▒▒░░░░░▒▓░                            ░░▒           ░░░▒          ▓░                      
           ░▓▓▒▒▓▒                           ░░░▒░            ░░░          ▓░                       
              ░▒▒                         ░░░░░░              ░░░         ▒▓                        
               ░▓         ░▒▓▓▓▓░                             ░░         ░▓░                        
                ░▓▒░   ░▒▓▒░   ░▓░                             ░         ▓░                         
                    ░░░░        ▓░                                      ▓▒                          
                               ▓▒                                      ▒▒                           
                              ░▓░            ░▓▓▓▓░                   ▒▒                            
                               ▒▒          ▒▓░░░░░░▓▓░               ▒▒                             
                                ░▓▒░    ░▒▓░░░░░░░░░▒▓              ▒▒                              
                                   ░▒▓▒▒░ ░▒░░░░░░░░░▓             ▒▒                               

"""
    )
