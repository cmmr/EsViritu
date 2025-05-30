import polars as pl
import pysam
import logging
import subprocess
from subprocess import Popen, PIPE, STDOUT
import sys, os
from collections import defaultdict
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

def trim_filter(reads: list, outdir: str, tempdir: str, trim: bool, filter: bool, 
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
    import os
    import subprocess
    import pysam

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(tempdir, exist_ok=True)

    # Set up file names
    sample_name = os.path.basename(reads[0]).split('.')[0]
    trimmed_fastq = os.path.join(tempdir, f"{sample_name}.trimmed.fastq")
    filtered_fastq = os.path.join(tempdir, f"{sample_name}.filtered.fastq")
    fastp_html = os.path.join(outdir, f"{sample_name}.fastp.html")
    fastp_json = os.path.join(outdir, f"{sample_name}.fastp.json")

    input_fastq = list(reads)
    # Step 1: Quality trim with fastp
    if trim:
        if threads > 16:
            fastp_threads = 16
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
        else:
            output_fastq = [filtered_fastq]
        return output_fastq
    else:
        return input_fastq

def minimap2_f(reference: str, 
                   reads: list, 
                   cpus: str,
                   sample: str,
                   tmpdir: str) -> str:
    '''
    aligns read pairs to reference with minimap2, passes it to pysam,
    then filter low quality alignments
    '''
    
    mini2_command = [
        'minimap2', '-t', cpus, 
        '-ax', 'sr', 
        '--secondary=yes', 
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

            outbf = os.path.join(tmpdir, f"{sample}.initial.filt.bam")
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
                ## assume the edit distance is between aligned part of read and ref
                try:
                    alignAcc: float=(alignLength-alignEdit)/alignLength
                except:
                    alignAcc: float=0

                if alignLength >= 100 and alignProp >= 0.9 and alignAcc >= 0.8:

                    filtbam.write(record)

    filtbam.close()

    # return an AlignmentFile
    return filtbam

# functions to add:
## take filtered .bam, make a coverm-like file
def bam_to_coverm_table(bam_path: str, sample: str) -> pl.DataFrame:
    """
    Takes a BAM file and sample name, sorts the BAM file, and generates a table with:
    contig length, covered bases, read count, and mean coverage for each contig.
    Returns a polars DataFrame.
    """
    # Sort the BAM file and index it
    sorted_bam_path = bam_path.replace('.bam', '.sorted.bam')
    pysam.sort('-o', sorted_bam_path, bam_path)
    pysam.index(sorted_bam_path)

    # Open the sorted BAM file
    bamfile = pysam.AlignmentFile(sorted_bam_path, "rb")
    
    records = []
    for contig in bamfile.header.references:
        length = bamfile.header.get_reference_length(contig)
        # Get coverage (depth) for each base in contig
        # count_coverage returns a tuple of 4 arrays (A,C,G,T)
        cov = bamfile.count_coverage(contig)
        # Sum across all bases
        total_cov = sum(cov)
        # Per-base depth
        per_base_coverage = [sum(bases) for bases in zip(*cov)]
        covered_bases = sum(1 for d in per_base_coverage if d > 0)
        mean_cov = sum(per_base_coverage) / length if length > 0 else 0
        # Read count for contig
        read_count = bamfile.count(contig=contig)
        records.append({
            "sample": sample,
            "contig": contig,
            "contig_length": length,
            "covered_bases": covered_bases,
            "read_count": read_count,
            "mean_coverage": mean_cov
        })
    bamfile.close()
    # Optionally remove temp sorted BAM files
    # os.remove(sorted_bam_path)
    # os.remove(sorted_bam_path + ".bai")
    return pl.DataFrame(records)

## take filtered .bam, make preliminary consensus .fastas
def bam_to_consensus_fasta(bam_path: str, output_fasta: str = None) -> str:
    """
    Calls samtools consensus on the given BAM file and writes the consensus FASTA.
    Only outputs consensus for records with aligned reads.
    Returns the path to the consensus FASTA file.
    """

    if output_fasta is None:
        output_fasta = bam_path.replace('.bam', '.consensus.fasta')
    # samtools consensus command
    cmd = [
        'samtools', 'consensus',
        '-o', output_fasta,
        bam_path
    ]
    subprocess.run(cmd, check=True)
    return output_fasta

## use minimap2 to compare preliminary consensus .fastas
def minimap2_self_compare(fasta_path: str, threads: int = 2) -> pl.DataFrame:
    """
    Aligns a fasta file of long nucleotide sequences against itself using minimap2.
    Returns a polars DataFrame with alignment stats: query, reference, percent identity, percent query covered, percent reference covered.
    """
    import subprocess
    import polars as pl
    from Bio import SeqIO
    import tempfile
    import os

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
        subprocess.run(cmd, stdout=outpaf, check=True)

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
        summary_records.append({
            'qname': qname,
            'tname': tname,
            'pid': pid,
            'qcov': qcov,
            'tcov': tcov
        })
    os.remove(paf_file.name)
    return pl.DataFrame(summary_records)

## use (updated?) anicalc and aniclust to cluster preliminary consensus .fastas
def run_aniclust(
    fna: str,
    ani: str,
    out: str,
    exclude: str = None,
    keep: str = None,
    min_ani: float = 95,
    min_qcov: float = 10,
    min_tcov: float = 70,
    min_length: int = 1
) -> None:
    """
    Runs the aniclust.py script as a subprocess with the given arguments.
    fna: path to nucleotide fasta
    ani: path to tab-delimited file with [qname, tname, num_alns, ani, qcov, tcov]
    out: output file path
    exclude: (optional) path to ids to exclude
    keep: (optional) path to ids to keep
    min_ani, min_qcov, min_tcov, min_length: clustering thresholds
    """
    import subprocess
    import sys
    script_path = os.path.join(os.path.dirname(__file__), 'utils', 'aniclust.py')
    cmd = [sys.executable, script_path, '--fna', fna, '--ani', ani, '--out', out,
           '--min_ani', str(min_ani), '--min_qcov', str(min_qcov), '--min_tcov', str(min_tcov), '--min_length', str(min_length)]
    if exclude:
        cmd += ['--exclude', exclude]
    if keep:
        cmd += ['--keep', keep]
    subprocess.run(cmd, check=True)

## re-align (mapped) reads to references based on preliminary consensus .fastas (fill N's??)
## take final .bam, make final consensus .fastas
## Make "main" output table
## recapitulate coverage table without bedtools
## Compare consensus .fastas to reference .fastas
## Make reactable (or python-based version?)