#!/usr/bin/env python

import argparse
import sys, os
import subprocess
import yaml
import logging
import polars as pl
try:
    from . import esv_funcs as esvf
except:
    import esv_funcs as esvf


def esviritu():

    pathname = os.path.dirname(__file__)  
    esviritu_script_path = os.path.abspath(pathname)      
    print(esviritu_script_path) 
    def_workdir = os.getcwd()

    __version__='0.9.0'

    parser = argparse.ArgumentParser(
        description='EsViritu is a read mapping pipeline for detection and measurement of \
            human and animal virus pathogens from short read metagenomic or clinical samples.\
            Version ' + str(__version__)
            )

    required_args = parser.add_argument_group(' REQUIRED ARGUMENTS for EsViritu ')

    required_args.add_argument(
        "-r", "--reads", nargs="+",
        dest="READS", required=True, 
        help='read file(s) in .fastq format. \
            You can specify more than one separated by a space'
            )
    required_args.add_argument(
        "-s", "--sample", 
        dest="SAMPLE", type=str, required=True, 
        help='Sample name. No space characters, please.'
        )

    required_args.add_argument(
        "-o", "--output_dir", 
        dest="OUTPUT_DIR", type=str, required=True, 
        help='Output directory name. Will be created if it does not exist. \
            Can be shared with other samples. No space characters, please. '
        )

    optional_args = parser.add_argument_group(' OPTIONAL ARGUMENTS for EsViritu.')

    optional_args.add_argument('--version', action='version', version=str(__version__))
    optional_args.add_argument(
        "-t", "--cpu", 
        dest="CPU", type=int, default = os.cpu_count(),
        help='Example: 32 -- Number of CPUs available for EsViritu.'
        )
    optional_args.add_argument(
        '-q', "--qual", dest="QUAL", 
        type=esvf.str2bool, default=False,
        help='True or False. Remove low-quality reads with fastp?'
        )
    optional_args.add_argument(
        '-f', "--filter_seqs", 
        dest="FILTER_SEQS", type=esvf.str2bool, default=False,
        help='True or False. Remove reads aligning to sequences at filter_seqs/filter_seqs.fna ?'
        )
    optional_args.add_argument(
        "--filter_dir", 
        dest="FILTER_DIR", type=str, default='default',
        help='path to directory of sequences to filter. If not set, EsViritu looks for environmental \
            variable ESVIRITU_FILTER. Then, if this variable is unset, it this is unset, \
            DB path is assumed to be ' + esviritu_script_path.replace("src/EsViritu", "filter_seqs")
            )
    optional_args.add_argument(
        '-i', "--compare", 
        dest="COMPARE", type=esvf.str2bool, default=True,
        help='True or False. Calculate percent identity between sample consensus and reference sequences?'
        )
    optional_args.add_argument(
        "-m", "--mode", 
        dest="MODE", type=str, default='general',
        help='Placeholder for future development. Do not change.'
        )
    optional_args.add_argument(
        "--temp", 
        dest="TEMP_DIR", type=str, default='default',
        help='path of temporary directory. Default is {OUTPUT_DIR}/{SAMPLE}_temp/'
        )
    optional_args.add_argument(
        "--keep", 
        dest="KEEP", type=esvf.str2bool, default=False,
        help='True of False. Keep the intermediate files, located in the temporary directory? \
            These can add up, so it is not recommended if space is a concern.'
        )
    optional_args.add_argument(
        "-p", "--read_format", 
        dest="READ_FMT", type=str, default='unpaired',
        help='unpaired or paired. Format of input reads. If paired, must provide 2 files \
            (R1, then R2) after -r argument.'
        )
    optional_args.add_argument(
        "--db", 
        dest="DB", type=str, default='default',
        help='path to sequence database. If not set, EsViritu looks for environmental \
            variable ESVIRITU_DB. Then, if this variable is unset, it this is unset, \
            DB path is assumed to be ' + esviritu_script_path.replace("src/EsViritu", "DBs/v2.0.2")
        )

    optional_args.add_argument(
        "-wd", "--working_directory", 
        dest="c_workdir", type=str, default=def_workdir, 
        help=f"Default: {def_workdir} -- Set working directory with absolute or relative path. \
            Run directory will be created within."
        )

    args = parser.parse_args()


    out_directory = os.path.join(str(args.c_workdir), str(args.OUTPUT_DIR))

    if not os.path.isdir(out_directory):
        os.makedirs(out_directory)

    #### define logger #####
    logger = logging.getLogger("esv_logger")
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Custom color formatter for info-level logs
    class ColorFormatter(logging.Formatter):
        BLUE = "\033[94m"
        BRIGHT_BLUE = "\033[94m"
        YELLOW = "\033[93m"
        RESET = "\033[0m"
        def format(self, record):
            msg = super().format(record)
            if record.levelno == logging.INFO:
                msg = f"{self.BRIGHT_BLUE}{msg}{self.RESET}"
            return msg

    # stream gets printed to terminal
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(ColorFormatter('%(asctime)s - %(levelname)s - %(message)s'))

    # file gets saved to a specified file
    log_file_path = os.path.join(
        out_directory, 
        f"esviritu_{args.SAMPLE}.log"
        )
    try:
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
    except Exception as e:
        logger.error(f"Failed to create log file handler: {e}")
        file_handler = None

    if not logger.hasHandlers():
        if file_handler:
            logger.addHandler(file_handler)
        logger.addHandler(stream_handler)

    logger.propagate = False
    #########################

    if str(args.TEMP_DIR) == "default":
        args.TEMP_DIR = os.path.join(
            args.OUTPUT_DIR,
            f"{args.SAMPLE}_temp"
        )
    
    logger.info(args.TEMP_DIR)

    READS = ' '.join(map(str,args.READS))

    if len(READS.split()) == 2 and str(args.READ_FMT).lower() == "paired":
        logger.info(len(READS.split()), str(args.READ_FMT).lower(), "read files")
    elif len(READS.split()) != 2 and str(args.READ_FMT).lower() == "paired":
        logger.info("if stating --read_format paired, must provide exactly 2 read files")
        sys.exit()

    for inr in args.READS:
        if not os.path.isfile(inr):
            logger.error(f'input read file not found at {inr}. exiting.')
            sys.exit()

    if args.DB == "default" and os.getenv('ESVIRITU_DB') != None:
        args.DB = os.getenv('ESVIRITU_DB')
    elif args.DB == "default":
        args.DB = esviritu_script_path.replace("src/EsViritu", "DBs/v2.0.2")

    db_index = os.path.join(args.DB, "virus_pathogen_database.mmi")
    if not os.path.isfile(db_index):
        logger.error(f'database file not found at {db_index}. exiting.')
        sys.exit()

    db_fasta = os.path.join(args.DB, "virus_pathogen_database.fna")
    if not os.path.isfile(db_fasta):
        logger.error(f'database file not found at {db_fasta}. exiting.')
        sys.exit()

    db_metadata = os.path.join(args.DB, "virus_pathogen_database.all_metadata.tsv")
    if not os.path.isfile(db_metadata):
        logger.error(f'database file not found at {db_metadata}. exiting.')
        sys.exit()

    if args.FILTER_DIR == "default" and os.getenv('ESVIRITU_FILTER') != None:
        args.FILTER_DIR = os.getenv('ESVIRITU_FILTER')
    elif args.FILTER_DIR == "default":
        args.FILTER_DIR = esviritu_script_path.replace("src/EsViritu", "filter_seqs")
    

    logger.info(f"DB: {str(args.DB)}")

    logger.info(f"filter seqs: {str(args.FILTER_DIR)}")

    logger.info(f"version {str(__version__)}")

    # check if R script with libraries returns good exit code
    completedProc = subprocess.run(['Rscript', str(esviritu_script_path) + '/check_R_libraries1.R'])

    #print(completedProc.returncode)
    if completedProc.returncode != 0 :
        logger.info("some required R packages are not found. Required:")
        logger.info("reactable, htmltools, dplyr, reactablefmtr, dataui, data.table, RColorBrewer, viridis, scales, knitr")
        logger.info("Did you activate the conda environment?")
        logger.info("see yml. Exiting")
        sys.exit()




    tool_dep_list = ['minimap2', 'fastp', 'seqkit', 'samtools']
    
    for tool in tool_dep_list:
        if not esvf.is_tool(tool):
            logger.warning(f"{tool} is not found. Exiting.")
            sys.exit()
    
    # parameters to yaml

    params_dict = vars(args).copy() 
    params_dict['version'] = str(__version__)
    params_dict['log_file'] = log_file_path

    # Save to YAML file in the same directory as the log file
    yaml_file_path = os.path.join(out_directory, f"esviritu_{args.SAMPLE}.params.yaml")
    try:
        with open(yaml_file_path, 'w') as yaml_file:
            yaml.dump(params_dict, yaml_file, default_flow_style=False)
        logger.info(f"Saved parameters to YAML: {yaml_file_path}")
    except Exception as e:
        logger.error(f"Failed to save parameters YAML: {e}")

    # trim and filter
    filter_db_fasta = os.path.join(
        args.FILTER_DIR,
        "filter_seqs.fna"
    )

    trim_filt_reads = esvf.trim_filter(
        args.READS,
        str(args.OUTPUT_DIR),
        str(args.TEMP_DIR),
        bool(args.QUAL),
        bool(args.FILTER_SEQS),
        str(args.SAMPLE),
        str(filter_db_fasta),
        bool(args.READ_FMT),
        str(args.CPU)
    )

    logger.info(trim_filt_reads)

    # get read stats for mapping pipeline
    filtered_reads = esvf.fastp_stats(
        trim_filt_reads,
        str(args.OUTPUT_DIR),
        str(args.SAMPLE),
        bool(args.READ_FMT),
        str(args.CPU)
    )

    # map reads to virus DB and filter for good alignments
    init_bam_f = os.path.join(str(args.TEMP_DIR), f"{str(args.SAMPLE)}.initial.filt.sorted.bam")
    initial_map_bam = esvf.minimap2_f(
        db_index,
        trim_filt_reads,
        str(args.CPU),
        init_bam_f
    )

    logger.info(initial_map_bam)

    # Make coverm-like table from initial bam
    init_coverm_like_dt = esvf.bam_to_coverm_table(
        initial_map_bam,
        str(args.SAMPLE)
    )

    init_coverm_like_dt.write_csv(
        file = os.path.join(
            args.TEMP_DIR,
            f"{str(args.SAMPLE)}_initial_coverm.tsv"
        ),
        separator = "\t"
    )

    # Load virus info metadata table in polars

    vir_meta_df = pl.read_csv(
        db_metadata, 
        separator='\t',
        schema_overrides={"TaxID": pl.String}
        )

    # make consensus .fasta from initial alignment
    initial_consensus = esvf.pileup_consensus(
        initial_map_bam,
        os.path.join(
            args.TEMP_DIR,
            f"{args.SAMPLE}.prelim.consensus.fasta"
        )
    )

    logger.info(initial_consensus)

    init_assembly_concat = esvf.concat_asm_accessions(
        initial_consensus,
        vir_meta_df
    )

    # compare initial consensus seqs to each other
    ## minimap2 method
    pairwise_mm_initial_dt = esvf.minimap2_self_compare(
        init_assembly_concat,
        str(args.CPU)
    )

    anicalc_mm_f = os.path.join(
            args.TEMP_DIR,
            f"{str(args.SAMPLE)}_mm_anicalc.tsv"
        )
    pairwise_mm_initial_dt.write_csv(
        file = anicalc_mm_f,
        separator = "\t"
    )
    logger.info(anicalc_mm_f)

    # Run aniclust.py
    aniclust_mm_f = os.path.join(
            args.TEMP_DIR,
            f"{str(args.SAMPLE)}_mm_aniclust.tsv"
        )
    aniclustpy_path = os.path.join(os.path.dirname(__file__), 'utils', 'aniclust.py')
    clustcmd = [
        sys.executable, aniclustpy_path, 
        '--fna', init_assembly_concat, 
        '--ani', anicalc_mm_f, 
        '--out', aniclust_mm_f,
        '--min_ani', str(98), 
        '--min_qcov', str(50), 
        '--min_tcov', str(0), 
        '--min_length', str(0)
        ]
    subprocess.run(clustcmd, check=True)

    logger.info(aniclust_mm_f)

    ## blastn method
    pairwise_bn_initial_dt = esvf.blastn_self_compare(
        init_assembly_concat,
        str(args.CPU)
    )

    anicalc_bn_f = os.path.join(
            args.TEMP_DIR,
            f"{str(args.SAMPLE)}_bn_anicalc.tsv"
        )
    pairwise_bn_initial_dt.write_csv(
        file = anicalc_bn_f,
        separator = "\t"
    )
    logger.info(anicalc_mm_f)

    # Run aniclust.py
    aniclust_bn_f = os.path.join(
            args.TEMP_DIR,
            f"{str(args.SAMPLE)}_bn_aniclust.tsv"
        )
    aniclustpy_path = os.path.join(os.path.dirname(__file__), 'utils', 'aniclust.py')
    clustcmd = [
        sys.executable, aniclustpy_path, 
        '--fna', init_assembly_concat, 
        '--ani', anicalc_mm_f, 
        '--out', aniclust_bn_f,
        '--min_ani', str(98), 
        '--min_qcov', str(50), 
        '--min_tcov', str(0), 
        '--min_length', str(0)
        ]
    subprocess.run(clustcmd, check=True)

    logger.info(aniclust_bn_f)

    ## make new fasta file from aniclust exemplars
    clust_db_fasta = esvf.final_record_getter(
        aniclust_bn_f,
        vir_meta_df,
        db_fasta,
        os.path.join(
            str(args.TEMP_DIR),
            f"{str(args.SAMPLE)}_clustered_refs.fasta"
        )
    )
    logger.info(clust_db_fasta)

    ## re-align (mapped) reads to references based on preliminary consensus .fastas (fill N's??)
    sec_bam_f = os.path.join(str(args.TEMP_DIR), f"{str(args.SAMPLE)}.second.filt.sorted.bam")

    second_map_bam = esvf.minimap2_f(
        clust_db_fasta,
        trim_filt_reads,
        str(args.CPU),
        sec_bam_f
    )

    ## take final .bam, make final consensus .fastas
    sec_con_f = os.path.join(
        str(args.OUTPUT_DIR),
        f"{str(args.SAMPLE)}_final_consensus.fasta"
    )
    second_consensus_fasta = esvf.bam_to_consensus_fasta(
        second_map_bam,
        sec_con_f
    )
    logger.info(second_consensus_fasta)

    # Make coverm-like table from final bam
    second_coverm_like_dt = esvf.bam_to_coverm_table(
        second_map_bam,
        str(args.SAMPLE)
    )

    logger.info(second_coverm_like_dt)

    second_coverm_like_dt.write_csv(
        file = os.path.join(
            args.TEMP_DIR,
            f"{str(args.SAMPLE)}_final_coverm.tsv"
        ),
        separator = "\t"
    )

    ## Make bedtools-like windows coverage
    windows_cov_df = esvf.bam_coverage_windows(
        second_map_bam
    )
    logger.info(windows_cov_df)

    ## Make "main" and "assembly" output table

    main_out_df, assem_out_df = esvf.assembly_table_maker(
        second_coverm_like_dt,
        vir_meta_df,
        filtered_reads,
        str(args.SAMPLE)
    )
    logger.info(main_out_df.schema)
    logger.info(main_out_df)

    logger.info(assem_out_df.schema)
    logger.info(assem_out_df)

    main_of = os.path.join(
        str(args.OUTPUT_DIR),
        f"{str(args.SAMPLE)}.detected_virus.info.tsv"
    )
    main_out_df.write_csv(
        file = main_of,
        separator = "\t"
    )

    logger.info(main_of)

    assem_of = os.path.join(
        str(args.OUTPUT_DIR),
        f"{str(args.SAMPLE)}.detected_virus.assembly_summary.json"
    )
    assem_out_df.write_ndjson(
        file = assem_of
    )

    logger.info(assem_of)


if __name__ == "__main__":
    esviritu()


