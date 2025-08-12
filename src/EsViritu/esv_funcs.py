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
                filter_db: str = None, paired: str = "paired", threads: int = 4) -> list:
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
            subprocess.run(fastp_cmd, check=True)
            input_fastq = [trimmed1, trimmed2]
        else:
            fastp_cmd = [
                "fastp", "--in1", reads[0], "--out1", trimmed_fastq,
                "-w", str(fastp_threads), "--html", fastp_html, "--json", fastp_json
            ]
            subprocess.run(fastp_cmd, check=True)
            input_fastq = [trimmed_fastq]

    # Step 2: Filter with minimap2 + pysam
    if filter:
        if filter_db is None:
            raise ValueError("filter_db must be provided if filter=True")
        # Build commands for piping
        if paired == "paired":
            read1, read2 = input_fastq
            minimap2_cmd = [
                "minimap2", "-t", str(threads), "-ax", "sr", filter_db, read1, read2
            ]
        else:
            minimap2_cmd = [
                "minimap2", "-t", str(threads), "-ax", "sr", filter_db, input_fastq[0]
            ]
        samtools_view_cmd = ["samtools", "view", "-bS", "-"]
        samtools_sort_cmd = ["samtools", "sort", "-"]
        if paired == "paired":
            unmapped_fastq1 = filtered_fastq.replace(".fastq", "_R1.fastq")
            unmapped_fastq2 = filtered_fastq.replace(".fastq", "_R2.fastq")
            samtools_fastq_cmd = ["samtools", "fastq", "-f", "4", "-1", unmapped_fastq1, "-2", unmapped_fastq2, "-"]
        else:
            samtools_fastq_cmd = ["samtools", "fastq", "-f", "4", "-o", filtered_fastq, "-"]

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
        return output_fastq
    else:
        return input_fastq

## fastp stats
def fastp_stats(reads: list, outdir: str, sample_name: str, trimarg: bool, filtarg: bool,
                tempdir: str, paired: str = "paired", threads: int = 4) -> str:

    pipeline_stats_fastp_html = os.path.join(outdir, f"{sample_name}.fastp.html")
    pipeline_stats_fastp_json = os.path.join(outdir, f"{sample_name}.fastp.json")

    if int(threads) > 16:
        fastp_threads = 16
    else:
        fastp_threads = int(threads)

    if trimarg and os.path.isfile(os.path.join(tempdir, f"{sample_name}.fastp.json")) and not filtarg:
        shutil.copy(
            os.path.join(tempdir, f"{sample_name}.fastp.json"),
            pipeline_stats_fastp_json
        )
        if os.path.isfile(os.path.join(tempdir, f"{sample_name}.fastp.html")):
            shutil.copy(
                os.path.join(tempdir, f"{sample_name}.fastp.html"),
                pipeline_stats_fastp_html
            )
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
    total_reads = None
    try:
        with open(pipeline_stats_fastp_json, "r") as f:
            data = json.load(f)
            total_reads = data["summary"]["before_filtering"]["total_reads"]
    except Exception as e:
        logging.warning(f"Could not parse total_reads from fastp JSON: {e}")
    return total_reads



def minimap2_f(reference: str, 
                   reads: list, 
                   cpus: str,
                   sorted_outbf) -> str:

    '''
    aligns read pairs to reference with minimap2, passes it to pysam,
    then filter low quality alignments
    '''
    
    mini2_command = [
        'minimap2', '-t', cpus, 
        '-ax', 'sr', 
        '--secondary=yes', 
        '--secondary-seq',
        '--sam-hit-only', "--MD",
        '-f', '10000',
        '-N' '100',
        '-p', '0.90',
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

    pysam.sort("-o", sorted_outbf, outbf)

    pysam.index(sorted_outbf)

    # return an AlignmentFile
    return sorted_outbf

## take filtered .bam, make a coverm-like file
def calculate_contig_stats(bam_path: str, contig: str, include_secondary: bool = False) -> Dict:
    """Calculate statistics for a single contig."""
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        length = bamfile.get_reference_length(contig)
        if length == 0:
            return None

        # Initialize coverage array
        coverage = np.zeros(length, dtype=np.uint32)
        pi_list = []
        
        # Get all reads for this contig
        reads = list(bamfile.fetch(contig))
        # only conunting primary alignments
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
        
        # Process contigs in parallel
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for contig in contigs_with_reads:
                futures.append(executor.submit(calculate_contig_stats, bam_path, contig, include_secondary))
            
            # Collect results
            records = []
            for future in futures:
                result = future.result()
                if result is not None:
                    records.append(result)

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
def bam_coverage_windows(bam_path: str) -> pl.DataFrame:
    """
    Splits each reference contig in a BAM file into 100 equal windows and reports the average read coverage of each window.
    Returns a polars DataFrame with columns: contig, window_index, window_start, window_end, average_coverage
    """

    bam = pysam.AlignmentFile(bam_path, "rb")
    records = []
    for contig in bam.references:
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
            for pileupcolumn in bam.pileup(contig=contig, max_depth = 100_000, start=start, end=end, stepper="all"):
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
    bam.close()
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

    # Step 5: Compare each pair of assemblies for shared reads
    assemblies = list(assembly_reads.keys())
    records = []
    for i, a in enumerate(assemblies):
        reads_a = assembly_reads[a]
        avg_id_a = assembly_identities.get(a, float('nan'))
        # this will do a self-self comparison so contigs without any sharing don't get missed
        for j in range(i, len(assemblies)):
            b = assemblies[j]
            reads_b = assembly_reads[b]
            avg_id_b = assembly_identities.get(b, float('nan'))
            shared = reads_a & reads_b
            if shared:
                n_shared = len(shared)
                share_a = n_shared / len(reads_a) if reads_a else 0.0
                share_b = n_shared / len(reads_b) if reads_b else 0.0
                tot_a = len(reads_a) if reads_a else 0.0
                tot_b = len(reads_b) if reads_b else 0.0
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
    
    # Process edges and collect read counts
    for row in df.iter_rows(named=True):
        n_shared = row['n_shared_reads']
        min_reads = min(row['total_reads_a'], row['total_reads_b'])
        total_reads[row['assembly_a']] = row['total_reads_a']
        total_reads[row['assembly_b']] = row['total_reads_b']
        
        if min_reads == 0:
            continue
            
        share = n_shared / min_reads
        if share >= threshold:
            # Store direct connections in both directions
            direct_connections[row['assembly_a']].add(row['assembly_b'])
            direct_connections[row['assembly_b']].add(row['assembly_a'])
    
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
        # Sort remaining assemblies by read count (descending), breaking ties by highest avg_read_identity
        def get_avg_identity(a):
            # Try both possible columns for avg_read_identity
            val = None
            if 'avg_read_identity_a' in df.columns:
                # Get max value for this assembly as either a or b
                vals = df.filter((pl.col('assembly_a') == a) | (pl.col('assembly_b') == a))
                vals_a = vals.filter(pl.col('assembly_a') == a)
                vals_b = vals.filter(pl.col('assembly_b') == a)
                id_a = vals_a['avg_read_identity_a'].to_list() if 'avg_read_identity_a' in vals_a.columns else []
                id_b = vals_b['avg_read_identity_b'].to_list() if 'avg_read_identity_b' in vals_b.columns else []
                all_ids = [x for x in id_a + id_b if x is not None]
                if all_ids:
                    val = max(all_ids)
            if val is None:
                return float('-inf')
            return val
        remaining.sort(key=lambda a: (total_reads.get(a, 0), get_avg_identity(a)), reverse=True)
        
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
    subprocess.run(cmd, check=True)
    return [fastq_path]


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
