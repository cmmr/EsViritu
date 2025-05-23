import polars as pl
import pysam
import logging
import subprocess
from subprocess import Popen, PIPE, STDOUT
import sys, os

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
                   cpus: str) -> pl.DataFrame:
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

        filter_records = []


        with pysam.AlignmentFile(fd_child, 'r') as ministream:
            for record in ministream:
                readLength: int=record.infer_read_length()
                alignLength: int=record.query_alignment_length
                if record.is_paired:
                    if record.is_read1:
                        orient: str = "R1"
                    else:
                        orient: str = "R2"
                else:
                    orient: str = "N"

                # fixes cases when the reads have .1 or .2 at the end
                if record.query_name.endswith(".1")\
                  or record.query_name.endswith(".2"):
                    pairn: str = record.query_name[:-2]
                else:
                    pairn: str = record.query_name

                # NM flag reports "edit distance" between read and ref
                try:
                    alignEdit: int=dict(record.get_tags(with_value_type=False)).get('NM')
                except:
                    alignEdit: int=0
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

                    ### Instead of a polars dataframe, this should be a filtered .bam file
                    ### ALLOW MULTIPLE ALIGNMENTS
                    filter_records.append(
                        [
                            record.query_name,
                            record.reference_name,
                            readLength,
                            alignProp,
                            alignLength,
                            alignAcc,
                            orient,
                            pairn
                        ]
                            )
                    
    filter_df = pl.DataFrame(filter_records,
                            schema=['qname', 'rname', 
                                    'read_length', 'aln_prop', 
                                    'aln_len', 'aln_acc', 
                                    'orientation', 'pair_name'],
                            orient = "row")
    
    # return a polars dataframe
    return filter_df

# functions to add:
## take filtered .bam, make a coverm-like file
## take filtered .bam, make preliminary consensus .fastas
## use minimap2 to compare preliminary consensus .fastas
## use (updated?) anicalc and aniclust to cluster preliminary consensus .fastas
## re-align (mapped) reads to references based on preliminary consensus .fastas (fill N's??)
## take final .bam, make final consensus .fastas
## Make "main" output table
## recapitulate coverage table without bedtools
## Compare consensus .fastas to reference .fastas
## Make reactable (or python-based version?)