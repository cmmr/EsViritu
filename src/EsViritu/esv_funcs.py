import polars as pl
import pysam
import logging
import subprocess
import json
import tempfile
from subprocess import Popen, PIPE, STDOUT
import sys, os
import statistics
from typing import Tuple
from collections import defaultdict, Counter
from distutils.spawn import find_executable
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
    
    return find_executable(name) is not None

def trim_filter(reads: list, outdir: str, tempdir: str, trim: bool, filter: bool, sample_name: str,
                filter_db: str = None, paired: bool = False, threads: int = 4) -> list:
    """
    Trim and/or filter .fastq reads using fastp (quality trim) and minimap2+pysam (host/spike-in filter).
    Args:
        reads: str, path to input fastq file(s). If paired, should be space-separated paths.
        outdir: str, output directory for summary files.
        tempdir: str, directory for intermediate files.
        trim: bool, whether to perform quality trimming with fastp.
        filter: bool, whether to filter reads mapping to filter_db with minimap2.
        filter_db: str, path to fasta file for minimap2 filtering (required if filter=True).
        paired: bool, whether reads are paired-end.
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
    # Step 1: Quality trim with fastp
    if trim:

        if paired:
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
        if paired:
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
        if paired:
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

        if paired:
            output_fastq = [unmapped_fastq1, unmapped_fastq2]
        return output_fastq
    else:
        return input_fastq

## fastp stats
def fastp_stats(reads: list, outdir: str, sample_name: str, 
                paired: bool = False, threads: int = 4) -> str:

    pipeline_stats_fastp_html = os.path.join(outdir, f"{sample_name}.fastp.html")
    pipeline_stats_fastp_json = os.path.join(outdir, f"{sample_name}.fastp.json")

    if int(threads) > 16:
        fastp_threads = 16

    if paired:
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
        '--sam-hit-only',
        '-f', '1000',
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
def bam_to_coverm_table(bam_path: str, sample: str, include_secondary: bool = True) -> pl.DataFrame:
    """
    Takes a sorted BAM file and sample name and generates a table with:
    contig length, covered bases, read count, mean coverage, and nucleotide diversity (π) for each contig.
    Efficiently skips contigs without alignments.
    Returns a polars DataFrame.
    """
    # Open the sorted BAM file
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    
    # First, identify contigs that have alignments
    contigs_with_reads = set()
    for read in bamfile.fetch():
        contigs_with_reads.add(bamfile.get_reference_name(read.reference_id))
    # Reset file pointer
    bamfile.reset()
    records = []
    stepper = "nofilter" if include_secondary else "all"

    for contig in bamfile.header.references:
        # Skip contigs without alignments
        if contig not in contigs_with_reads:
            continue
            
        length = bamfile.header.get_reference_length(contig)
        # Use pileup to get per-base coverage
        coverage = [0] * length
        pi_list = []
        for pileupcolumn in bamfile.pileup(contig=contig, start=0, end=length, stepper=stepper):
            coverage[pileupcolumn.pos] = pileupcolumn.nsegments
            # calculating pi 
            base_counts = Counter()
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    base_counts[base] += 1
            if base_counts: 
                tot_cov = sum(base_counts.values())
                maxN_cov = base_counts[max(base_counts, key=base_counts.get)]
                pi = (tot_cov - maxN_cov) / tot_cov
                pi_list.append(pi)
        
        avg_pi = round(statistics.mean(pi_list), 3)
        covered_bases = sum(1 for d in coverage if d > 0)
        mean_cov = sum(coverage) / length if length > 0 else 0
        read_count = bamfile.count(contig=contig)

        records.append({
            "sample": sample,
            "Accession": contig,
            "contig_length": length,
            "covered_bases": covered_bases,
            "read_count": read_count,
            "mean_coverage": mean_cov,
            "Pi": avg_pi
        })
    bamfile.close()

    return pl.DataFrame(records)

## take filtered .bam, make preliminary consensus .fastas
def bam_to_consensus_fasta(bam_path: str, output_fasta: str = None) -> str:
    """
    Calls samtools consensus on the given sorted BAM file and writes the consensus FASTA.
    Only outputs consensus for records with aligned reads.
    Returns the path to the consensus FASTA file.
    """

    if output_fasta is None:
        output_fasta = bam_path.replace('.bam', '.consensus.fasta')
    # samtools consensus command
    cmd = [
        'samtools', 'consensus',
        '-m', 'simple', 
        '-q', '-c', '0.51',
        '-o', output_fasta,
        bam_path
    ]
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        logger.error(f"samtools consensus failed: {cmd}\nError: {e}")
        raise
    return output_fasta

## (experimental) concatenate records from the same assembly
def concat_asm_accessions(fasta_path: str, df: pl.DataFrame, output_fasta: str = None) -> str:
    """
    Concatenate records from the same assembly in a fasta file based on a mapping DataFrame.
    - fasta_path: input fasta file
    - df: polars DataFrame with columns 'Accession' and 'Assembly'
    - output_fasta: path to write concatenated fasta (if None, auto-generate)
    Returns the path to the output fasta.
    """

    # Get accessions in fasta using seqkit
    with tempfile.NamedTemporaryFile(delete=False, mode='w+', suffix='.txt') as acc_file:
        seqkit_cmd = ["seqkit", "fx2tab", "-n", "-i", fasta_path]
        try:
            subprocess.run(seqkit_cmd, check=True, stdout=acc_file)
        except Exception as e:
            logger.error(f"seqkit_cmd failed: {seqkit_cmd}\nError: {e}")
            raise
        acc_file.flush()
        acc_file.seek(0)
        accessions = [line.strip() for line in acc_file.readlines() if line.strip()]
    os.remove(acc_file.name)

    # Filter DataFrame by accessions present in fasta
    df = df.filter(df['Accession'].is_in(accessions))

    # Build assembly: [Accession, ...] dictionary
    asm_dict = {}
    for row in df.iter_rows(named=True):
        asm = row['Assembly']
        acc = row['Accession']
        asm_dict.setdefault(asm, []).append(acc)

    # Prepare output
    if output_fasta is None:
        output_fasta = fasta_path.replace('.fasta', '.concat_asm.fasta')

    with open(output_fasta, 'w') as out_f:
        for asm, acc_list in asm_dict.items():
            # Use seqkit to pull out all accessions for this assembly
            with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as accs_f:
                accs_f.write('\n'.join(acc_list) + '\n')
                accs_f.flush()
                seqkit_cmd = ["seqkit", "grep", "-f", accs_f.name, fasta_path]
                try:
                    result = subprocess.run(seqkit_cmd, check=True, capture_output=True, text=True)
                except Exception as e:
                    logger.error(f"seqkit_cmd failed: {seqkit_cmd}\nError: {e}")
                    raise
            os.remove(accs_f.name)
            # Concatenate sequences (remove headers, join seqs)
            seqs = []
            for line in result.stdout.splitlines():
                if line.startswith('>'):
                    continue
                seqs.append(line.strip())
            concat_seq = ''.join(seqs)
            concat_seq = concat_seq.replace('NNN', '')
            # Write new header: assembly [accession1 accession2 ...]
            # only keep seqs over 100 nt
            if len(concat_seq) >= 100:
                out_f.write(f">{asm} [{' '.join(acc_list)}]\n{concat_seq}\n")
    return output_fasta


## use minimap2 to compare preliminary consensus .fastas
def minimap2_self_compare(fasta_path: str, threads: int = 2) -> pl.DataFrame:
    """
    Aligns a fasta file of long nucleotide sequences against itself using minimap2.
    Returns a polars DataFrame with alignment stats: query, reference, percent identity, percent query covered, percent reference covered.
    """

    # Best minimap2 settings for long nucleotide sequences: map-ont or asm20
    # Use PAF output for easy parsing
    paf_file = tempfile.NamedTemporaryFile(delete=False, suffix='.paf')
    cmd = [
        'minimap2',
        '-cx', 'asm20',  # or 'map-ont', but asm20 is for assembly/long reads
        '-t', str(threads),
        fasta_path, fasta_path
    ]
    with open(paf_file.name, 'w') as outpaf:
        try:
            subprocess.run(cmd, stdout=outpaf, check=True)
        except Exception as e:
            logger.error(f"minimap2 cmd failed: {cmd}\nError: {e}")
            raise

    # PAF format columns:
    # [0]query [1]qlen [2]qstart [3]qend [4]strand [5]target [6]tlen [7]tstart [8]tend [9]nmatch [10]alen [11]mapq ... tags

    # Aggregate all alignments for each (qname, tname) pair
    pair_records = defaultdict(list)
    with open(paf_file.name) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            qname = fields[0]
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            tname = fields[5]
            tlen = int(fields[6])
            tstart = int(fields[7])
            tend = int(fields[8])
            nmatch = int(fields[9])
            alen = int(fields[10])
            pair_records[(qname, tname)].append({
                'qlen': qlen,
                'qstart': qstart,
                'qend': qend,
                'tlen': tlen,
                'tstart': tstart,
                'tend': tend,
                'nmatch': nmatch,
                'alen': alen
            })

    # Summarize for each pair
    summary_records = []
    for (qname, tname), aligns in pair_records.items():
        total_alen = sum(a['alen'] for a in aligns)
        total_nmatch = sum(a['nmatch'] for a in aligns)
        qlen = aligns[0]['qlen']
        tlen = aligns[0]['tlen']
        # Collect all covered query and target positions
        q_covered = set()
        t_covered = set()
        for a in aligns:
            q_covered.update(range(a['qstart'], a['qend']))
            t_covered.update(range(a['tstart'], a['tend']))
        # Calculate summary stats
        pid = 100 * total_nmatch / total_alen if total_alen > 0 else 0
        qcov = 100 * len(q_covered) / qlen if qlen > 0 else 0
        tcov = 100 * len(t_covered) / tlen if tlen > 0 else 0
        num_alignments = len(aligns)
        summary_records.append({
            'qname': qname,
            'tname': tname,
            'num_alns': num_alignments,
            'ani': pid,
            'qcov': qcov,
            'tcov': tcov
        })
    os.remove(paf_file.name)
    df = pl.DataFrame(summary_records)
    # Filter out self-alignments where qname == tname
    df = df.filter(df['qname'] != df['tname'])
    return df


def blastn_self_compare(fasta_path: str, threads: int = 2) -> pl.DataFrame:
    """
    Aligns a fasta file of long nucleotide sequences against itself using BLASTN.
    Returns a polars DataFrame with alignment stats: query, subject, percent identity, percent query covered, percent subject covered, number of alignments.
    """

    # Create a temporary BLAST database
    db_dir = tempfile.TemporaryDirectory()
    db_prefix = os.path.join(db_dir.name, "blastdb")
    makeblastdb_cmd = [
        "makeblastdb",
        "-in", fasta_path,
        "-dbtype", "nucl",
        "-out", db_prefix
    ]
    try:
        subprocess.run(makeblastdb_cmd, check=True)
    except Exception as e:
        logger.error(f"makeblastdb_cmd failed: {makeblastdb_cmd}\nError: {e}")
        raise

    # Run BLASTN all-vs-all
    blast_out = tempfile.NamedTemporaryFile(delete=False, suffix='.blast.tsv')
    blastn_cmd = [
        "blastn",
        "-query", fasta_path,
        "-db", db_prefix,
        "-outfmt", "6 qseqid qlen qstart qend sseqid slen sstart send pident length",
        "-num_threads", str(threads),
        "-out", blast_out.name
    ]
    try:
        subprocess.run(blastn_cmd, check=True)
    except Exception as e:
        logger.error(f"blastn_cmd failed: {blastn_cmd}\nError: {e}")
        raise

    # Parse BLAST output
    pair_records = defaultdict(list)
    with open(blast_out.name) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            qname = fields[0]
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            sname = fields[4]
            slen = int(fields[5])
            sstart = int(fields[6])
            send = int(fields[7])
            pident = float(fields[8])
            alen = int(fields[9])
            pair_records[(qname, sname)].append({
                'qlen': qlen,
                'qstart': min(qstart, qend),
                'qend': max(qstart, qend),
                'slen': slen,
                'sstart': min(sstart, send),
                'send': max(sstart, send),
                'nmatch': int(round(pident * alen / 100)),
                'alen': alen,
                'pident': pident
            })

    # Summarize for each pair
    summary_records = []
    for (qname, sname), aligns in pair_records.items():
        total_alen = sum(a['alen'] for a in aligns)
        total_nmatch = sum(a['nmatch'] for a in aligns)
        qlen = aligns[0]['qlen']
        slen = aligns[0]['slen']
        # Collect all covered query and subject positions
        q_covered = set()
        s_covered = set()
        for a in aligns:
            q_covered.update(range(a['qstart'], a['qend']))
            s_covered.update(range(a['sstart'], a['send']))
        pid = 100 * total_nmatch / total_alen if total_alen > 0 else 0
        qcov = 100 * len(q_covered) / qlen if qlen > 0 else 0
        tcov = 100 * len(s_covered) / slen if slen > 0 else 0
        num_alignments = len(aligns)
        summary_records.append({
            'qname': qname,
            'tname': sname,
            'num_alns': num_alignments,
            'pid': float(f"{pid:.3g}"),
            'qcov': float(f"{qcov:.3g}"),
            'tcov': float(f"{tcov:.3g}")
        })
    # Cleanup
    blast_out.close()
    os.remove(blast_out.name)
    db_dir.cleanup()

    df = pl.DataFrame(summary_records)
    # Filter out self-alignments where qname == tname
    df = df.filter(df['qname'] != df['tname'])
    return df


## select final records based on aniclust

def final_record_getter(aniclust_tsv: str, field: int, vir_info_df: pl.DataFrame, db_fasta: str, output_fasta: str) -> str:
    """
    Select final records based on aniclust clusters.
    Args:
        aniclust_tsv: path to aniclust output tsv file (first field is assembly)
        field: which field has exemplar name
        vir_info_df: polars DataFrame with at least columns 'Assembly' and 'Accession'
        db_fasta: path to full fasta file
        output_fasta: path to write filtered fasta (optional)
    Returns:
        Path to new fasta file with selected accessions.
    """

    # Step 1: Get assemblies from aniclust tsv
    assemblies = set()
    with open(aniclust_tsv, 'r') as f:
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
            for pileupcolumn in bam.pileup(contig=contig, start=start, end=end, stepper="all"):
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

## Alternative, loose consensus func
def pileup_consensus(bam_path: str, output_fasta: str) -> str:
    """
    Generates a consensus FASTA sequence from a BAM file by examining the pileup at each position.
    For each position, the most frequent base is chosen.
    Writes the consensus sequences for each contig to the specified output FASTA file.
    Skips contigs with no aligned reads for efficiency.
    
    Args:
        bam_path: Path to the input BAM file.
        output_fasta: Path to the output FASTA file where consensus sequences will be written.
        
    Returns:
        The path to the output FASTA file containing the consensus sequences.
    """
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    
    # First, identify contigs that have alignments
    contigs_with_reads = set()
    for read in bamfile.fetch():
        contigs_with_reads.add(bamfile.get_reference_name(read.reference_id))
    
    # Reset file pointer
    bamfile.reset()
    
    with open(output_fasta, 'w') as outf:
        for contig in bamfile.header.references:
            # Skip contigs without alignments
            if contig not in contigs_with_reads:
                continue
                
            con_seq = []
            for pileupcolumn in bamfile.pileup(contig, stepper="nofilter"):
                base_counts = Counter()
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                        base_counts[base] += 1
                if base_counts: 
                    con_seq.append(max(base_counts, key=base_counts.get))
                else:
                    con_seq.append("N")
            
            if len(con_seq) > 100:
                print(f">{contig}", file = outf)
                seq = ''.join(con_seq)
                for i in range(0, len(seq), 60):
                    print(seq[i:i+60], file=outf)
    
    bamfile.close()

    return output_fasta

## read-based clustering function

def assembly_read_sharing_table(bam_path, meta_df: 'pl.DataFrame') -> pl.DataFrame:
    """
    For a given BAM file and meta_df (with 'Accession' and 'Assembly'),
    filter out contigs with no aligned reads, group contigs by Assembly, and compare each pair of contigs within each Assembly
    based on how many read IDs they share (primary or secondary alignments).
    Returns a polars DataFrame with Assembly, contig pairs, number of shared reads, and share of total reads for each contig.
    """
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    # Step 1: Collect read IDs for each contig
    contig_reads = defaultdict(set)
    for read in bamfile.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        contig = bamfile.get_reference_name(read.reference_id)
        contig_reads[contig].add(read.query_name)
        # Also add for secondary alignments
        if read.has_tag("SA"):
            contig_reads[contig].add(read.query_name)
    bamfile.close()

    # Step 2: Filter out contigs with no reads
    contig_reads = {k: v for k, v in contig_reads.items() if v}
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

    # Step 4: Build assembly-level read sets
    assembly_reads = {}
    for assembly, contigs in assembly_to_contigs.items():
        reads = set()
        for contig in contigs:
            reads.update(contig_reads[contig])
        assembly_reads[assembly] = reads

    # Step 5: Compare each pair of assemblies for shared reads
    assemblies = list(assembly_reads.keys())
    records = []
    for i, a in enumerate(assemblies):
        reads_a = assembly_reads[a]
        # this will do a self-self comparison so contigs without any sharing don't get missed
        for j in range(i, len(assemblies)):
            b = assemblies[j]
            reads_b = assembly_reads[b]
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
                    "share_of_b": share_b
                })
    return pl.DataFrame(records)


## cluster contigs on read sharing

def cluster_assemblies_by_read_sharing(df, threshold=0.33) -> pl.DataFrame:
    """
    Cluster assemblies that share at least `threshold` (default 33%) of their reads (based on the assembly with fewer reads).
    For each cluster, choose an exemplar assembly with the highest number of reads aligned.
    Returns a polars DataFrame mapping each assembly to its cluster and exemplar.
    """

    # Build undirected graph of assemblies to cluster
    edges = []
    for row in df.iter_rows(named=True):
        n_shared = row['n_shared_reads']
        min_reads = min(row['total_reads_a'], row['total_reads_b'])
        if min_reads == 0:
            continue
        share = n_shared / min_reads
        if share >= threshold:
            edges.append((row['assembly_a'], row['assembly_b']))

    # Union-find setup
    parent = {}
    def find(x):
        while parent.get(x, x) != x:
            x = parent[x]
        return x
    def union(x, y):
        xr, yr = find(x), find(y)
        if xr != yr:
            parent[yr] = xr

    # Union all connected assemblies
    for a, b in edges:
        union(a, b)

    # Assign cluster ids
    all_assemblies = set(df['assembly_a'].to_list()) | set(df['assembly_b'].to_list())
    clusters = defaultdict(list)
    for assembly in all_assemblies:
        cluster_id = find(assembly)
        clusters[cluster_id].append(assembly)

    # Get total reads for each assembly
    total_reads = {}
    for row in df.iter_rows(named=True):
        total_reads[row['assembly_a']] = row['total_reads_a']
        total_reads[row['assembly_b']] = row['total_reads_b']

    # Prepare output
    records = []
    for cluster_id, members in clusters.items():
        # Choose exemplar: assembly with max total reads
        exemplar = max(members, key=lambda c: total_reads.get(c, 0))
        for assembly in members:
            records.append({
                'cluster_id': cluster_id,
                'assembly': assembly,
                'exemplar': exemplar,
                'total_reads': total_reads.get(assembly, 0),
                'is_exemplar': assembly == exemplar
            })
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
            if read.is_unmapped or read.reference_id is None:
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
### declare ambiguity for high divergence?
### Calculate average nucleotide diversity per position
def calculate_nucleotide_diversity(bam_path: str, reference: str = None, min_mapq: int = 20, min_baseq: int = 20) -> pl.DataFrame:
    """
    Calculate nucleotide diversity (π) per contig from a BAM file using bcftools.
    
    Args:
        bam_path: Path to input BAM file (must be indexed)
        reference: Path to reference genome FASTA file (optional, but recommended)
        min_mapq: Minimum mapping quality
        min_baseq: Minimum base quality
        
    Returns:
        polars.DataFrame with columns: contig, length, n_sites, n_variant_sites, pi
        where pi is the average nucleotide diversity per site
    """

    # Check if bcftools is installed
    try:
        subprocess.run(['bcftools', '--version'], 
                      check=True, 
                      stdout=subprocess.PIPE, 
                      stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.error("bcftools is not installed or not in PATH. Please install bcftools.")
        raise RuntimeError("bcftools is required but not found. Please install bcftools.")
    
    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmpdir:
        # Generate VCF using bcftools mpileup
        vcf_path = os.path.join(tmpdir, 'variants.vcf.gz')
        
        # Build bcftools command
        cmd = [
            'bcftools', 'mpileup',
            '--output-type', 'z',
            '--min-MQ', str(min_mapq),
            '--min-BQ', str(min_baseq),
            '--output', vcf_path
        ]
        
        if reference:
            cmd.extend(['--fasta-ref', reference])
        
        cmd.append(bam_path)
        
        # Run bcftools mpileup
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"bcftools mpileup failed: {e.stderr}")
            raise
        
        # Index the VCF
        try:
            subprocess.run(['bcftools', 'index', vcf_path], check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"bcftools index failed: {e.stderr}")
            raise
        
        # Calculate π using bcftools stats
        try:
            result = subprocess.run(
                ['bcftools', 'stats', vcf_path],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse the output
            contig_data = []
            current_contig = None
            
            for line in result.stdout.split('\n'):
                if line.startswith('SN'):
                    fields = line.split('\t')
                    if len(fields) >= 4 and 'number of records' in fields[2]:
                        contig_data.append({
                            'contig': current_contig,
                            'n_variant_sites': int(fields[3])
                        })
                elif line.startswith('#'):
                    # Skip comment lines
                    continue
                elif line.startswith('ID'):
                    # This is a contig line
                    fields = line.split('\t')
                    if len(fields) >= 3:
                        current_contig = fields[1]
                        contig_data.append({
                            'contig': current_contig,
                            'length': int(fields[2])
                        })
            
            # Convert to a dictionary for easier processing
            contig_dict = {}
            for item in contig_data:
                contig = item['contig']
                if contig not in contig_dict:
                    contig_dict[contig] = {'contig': contig}
                contig_dict[contig].update(item)
            
            # Calculate π for each contig
            results = []
            for contig, data in contig_dict.items():
                if 'length' not in data:
                    continue
                    
                length = data['length']
                n_var_sites = data.get('n_variant_sites', 0)
                
                # π = (number of variant sites) / (total sites)
                # This is a simplified calculation - actual π would require allele frequencies
                pi = n_var_sites / length if length > 0 else 0.0
                
                results.append({
                    'contig': contig,
                    'length': length,
                    'n_sites': length,
                    'n_variant_sites': n_var_sites,
                    'pi': pi
                })
            
            return pl.DataFrame(results)
            
        except subprocess.CalledProcessError as e:
            logging.error(f"bcftools stats failed: {e.stderr}")
            raise