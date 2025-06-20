import os
import subprocess
import tempfile
from collections import defaultdict
from typing import Dict, List

import polars as pl
import pysam
from pysam import AlignmentFile

# Set up logger
import logging
logger = logging.getLogger("esv_logger")

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
