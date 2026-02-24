#!/usr/bin/env python

import argparse
import time
import sys, os
import subprocess
import yaml
import logging
import polars as pl
import shutil
try:
    from . import esv_funcs as esvf
    from .utils.timing import timed_function
except:
    import esv_funcs as esvf
    from utils.timing import timed_function


def esviritu():
    """EsViritu CLI entrypoint function"""
    pathname = os.path.dirname(__file__)  
    esviritu_script_path = os.path.abspath(pathname)      
    print(esviritu_script_path) 
    def_workdir = os.getcwd()

    __version__='1.2.0'

    esv_start_time = time.perf_counter()

    parser = argparse.ArgumentParser(
        description='EsViritu is a read mapping pipeline for detection and measurement of \
            human, animal, and plat virus pathogens from short read libraries. \
            It is useful for clinical and environmental samples. \
            Version ' + str(__version__)
            )

    required_args = parser.add_argument_group(' REQUIRED ARGUMENTS for EsViritu ')

    required_args.add_argument(
        "-r", "--reads", nargs="+",
        dest="READS", required=True, 
        help='read file(s) in .fastq format. \
            For unpaired reads of any kind, exactly one file path required. \
            For paired reads, exactly two file paths required, separated by a space. \
            Wildcards, e.g /path/to/reads/*fastq, will work but number of files must be correct.'
            )
    required_args.add_argument(
        "-s", "--sample", 
        dest="SAMPLE", type=str, required=True, 
        help='Sample name. No space characters, please.'
        )

    required_args.add_argument(
        "-o", "--output_dir", 
        dest="OUTPUT_DIR", type=str, required=True, 
        help='Output directory name (not a path). Will be created if it does not exist. \
            Can be shared with other samples. No space characters, please. \
            See also --working_directory to create at another path.'
        )

    optional_args = parser.add_argument_group(' OPTIONAL ARGUMENTS for EsViritu.')

    optional_args.add_argument('--version', action='version', version=str(__version__))
    optional_args.add_argument(
        "-t", "--cpu", 
        dest="CPU", type=int, default = os.cpu_count(),
        help='Example: 32 -- Number of CPUs available for EsViritu. Default = all available CPUs'
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
        "--temp", 
        dest="TEMP_DIR", type=str, default='default',
        help='path of temporary directory. Default is {OUTPUT_DIR}/{SAMPLE}_temp/'
        )
    optional_args.add_argument(
        "--keep", 
        dest="KEEP", type=esvf.str2bool, default=False,
        help='True of False. Keep the intermediate files located in the temporary directory? \
            These can add up, so it is not recommended if space is a concern.'
        )
    optional_args.add_argument(
        "-p", "--read_format", 
        dest="READ_FMT", type=str, default='paired',
        help='unpaired or paired. Format of input reads. If paired, must provide 2 files \
            (R1, then R2) after -r argument.'
        )
    optional_args.add_argument(
        "-mmP", "--minimap2-preset", dest="MM_SET", type=str, choices=['sr', 'map-hifi', 'lr:hq'], 
        default='sr',
        help=f'Default = sr -- minimap2 alignment preset. \
            sr: short read data (Illumina, other high-accuracy short reads). \
            map-hifi: HiFi PacBio reads. \
            lr:hq: Oxford Nanopore reads. \
            NOTE: You can only use -q True with "-mmP sr". fastp is not intended for long reads.'
        )
    optional_args.add_argument(
        "--db", 
        dest="DB", type=str, default='default',
        help='path to sequence database. If not set, EsViritu looks for environmental \
            variable ESVIRITU_DB. Then, if this variable is unset, it this is unset, \
            DB path is assumed to be ' + esviritu_script_path.replace("src/EsViritu", "DBs/v3.1.0. \
            As of EsViritu v1.0.0, DB v3.1.0 or higher is required.")
        )
    optional_args.add_argument(
        "-wd", "--working_directory", 
        dest="c_workdir", type=str, default=def_workdir, 
        help=f"Default: {def_workdir} -- Set working directory with absolute or relative path. \
            Run directory will be created within."
        )
    optional_args.add_argument(
        "-mmK", "--minimap2-K", 
        dest="MMK", type=str, default="50M", 
        help=f"Default: 500M -- minimap2 K parameter for Number of bases loaded \
            into memory to process in a mini-batch. Reducing this value lowers memory consumption"
        )
    optional_args.add_argument(
        "--species-threshold", 
        dest="spthresh", type=float, default=0.90, 
        help=f"Default: 0.90 -- minimum ANI of reads to reference to classify record at species level.\
            There is no perfect metric for this."
        )
    optional_args.add_argument(
        "--subspecies-threshold",
        dest="subspthresh", type=float, default=0.95,
        help=f"Default: 0.95 -- minimum ANI of reads to reference to classify record at subspecies level.\
            There is no perfect metric for this."
        )
    optional_args.add_argument(
        "--dedup",
        dest="DEDUP", type=esvf.str2bool, default=False,
        help='True or False. Remove PCR duplicates during fastp preprocessing? \
            This can reduce processing time and provide more accurate abundance estimates.'
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
        RED = "\033[31m"
        RESET = "\033[0m"
        def format(self, record):
            msg = super().format(record)
            if record.levelno == logging.INFO:
                msg = f"{self.YELLOW}{msg}{self.RESET}"
            elif record.levelno == logging.WARNING:
                msg = f"{self.BRIGHT_BLUE}{msg}{self.RESET}"
            elif record.levelno == logging.ERROR:
                msg = f"{self.RED}{msg}{self.RESET}"
            return msg

    # stream gets printed to terminal
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(ColorFormatter('%(asctime)s - %(levelname)s - %(message)s'))

    # file gets saved to a specified file
    log_file_path = os.path.join(
        out_directory, 
        f"{args.SAMPLE}_esviritu.log"
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

    esvf.print_esviritu_banner()
    esvf.ghost_banner()
    
    if str(args.TEMP_DIR) == "default":
        args.TEMP_DIR = os.path.join(
            out_directory,
            f"{args.SAMPLE}_temp"
        )
    

    if len(args.READS) != 2 and str(args.READ_FMT).lower() == "paired":
        logger.error("if stating --read_format paired, must provide exactly 2 read files")
        sys.exit()
    elif len(args.READS) != 1 and str(args.READ_FMT).lower() == "unpaired":
        logger.error("if stating --read_format unpaired, must provide exactly 1 read file")
        sys.exit()

    for inr in args.READS:
        if not os.path.isfile(inr):
            logger.error(f'input read file not found at {inr}. exiting.')
            sys.exit()

    if args.QUAL and str(args.MM_SET) != 'sr':
        logger.error(f'flag "-q True" is only compatible with "-mmP sr".')
        sys.exit()

    if args.DB == "default" and os.getenv('ESVIRITU_DB') != None:
        args.DB = os.getenv('ESVIRITU_DB')
    elif args.DB == "default":
        args.DB = esviritu_script_path.replace("src/EsViritu", "DBs/v3.1.0")

    if str(args.MM_SET) == 'sr':
        db_index = os.path.join(args.DB, "virus_pathogen_database.mmi")
    else:
        db_index = os.path.join(args.DB, "virus_pathogen_database.fna")
    if not os.path.isfile(db_index):
        logger.error(f'database file not found at {db_index}. Exiting. \
            As of EsViritu v1.0.0, DB v3.1.0 or higher is required.')
        sys.exit()

    db_fasta = os.path.join(args.DB, "virus_pathogen_database.fna")
    if not os.path.isfile(db_fasta):
        logger.error(f'database file not found at {db_fasta}. Exiting. \
            As of EsViritu v1.0.0, DB v3.1.0 or higher is required.')
        sys.exit()

    db_metadata = os.path.join(args.DB, "virus_pathogen_database.all_metadata.tsv")
    if not os.path.isfile(db_metadata):
        logger.error(f'database file not found at {db_metadata}. Exiting. \
            As of EsViritu v1.0.0, DB v3.1.0 or higher is required.')
        sys.exit()

    if args.FILTER_DIR == "default" and os.getenv('ESVIRITU_FILTER') != None:
        args.FILTER_DIR = os.getenv('ESVIRITU_FILTER')
    elif args.FILTER_DIR == "default":
        args.FILTER_DIR = esviritu_script_path.replace("src/EsViritu", "filter_seqs")
    
    logger.info(f"DB: {str(args.DB)}")

    logger.info(f"version {str(__version__)}")

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
    yaml_file_path = os.path.join(out_directory, f"{args.SAMPLE}_esviritu.params.yaml")
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



    ################################################
    ################################################
    ###### Beginning main alignment pipeline #######
    ################################################
    ################################################
    
    logger.info(f"main read processing starting for {args.SAMPLE}")

    trim_filter_fn = timed_function(logger=logger)(esvf.trim_filter)
    trim_filt_reads = trim_filter_fn(
        args.READS,
        str(out_directory),
        str(args.TEMP_DIR),
        bool(args.QUAL),
        bool(args.FILTER_SEQS),
        str(args.SAMPLE),
        str(filter_db_fasta),
        str(args.READ_FMT),
        str(args.CPU),
        str(args.MMK),
        str(args.MM_SET),
        bool(args.DEDUP)
    )

    logger.info(f"Main input reads: {trim_filt_reads}")

    # get read stats for mapping pipeline
    fastp_stats_fn = timed_function(logger=logger)(esvf.fastp_stats)
    filtered_reads = fastp_stats_fn(
        trim_filt_reads,
        str(out_directory),
        str(args.SAMPLE),
        bool(args.QUAL),
        bool(args.FILTER_SEQS),
        str(args.TEMP_DIR),
        str(args.READ_FMT),
        str(args.CPU)
    )

    # make a .yaml file with various readstats for later use
    out_readstats_yaml = esvf.various_readstats(
        str(out_directory),
        str(args.SAMPLE),
        bool(args.QUAL),
        bool(args.FILTER_SEQS),
        str(args.TEMP_DIR),
        filtered_reads
    )

    logger.info(f"read stats file: {out_readstats_yaml}")

    # map reads to virus DB and filter for good alignments
    init_bam_f = os.path.join(str(args.TEMP_DIR), f"{str(args.SAMPLE)}.initial.filt.sorted.bam")
    minimap2_f_fn = timed_function(logger=logger)(esvf.minimap2_f)
    initial_map_bam = minimap2_f_fn(
        db_index,
        trim_filt_reads,
        str(args.CPU),
        init_bam_f,
        str(args.MMK),
        str(args.MM_SET)
    )

    logger.info(f"initial bam: {initial_map_bam}")

    ## check if any reads aligned or quit
    if not esvf.bam_has_alignments(initial_map_bam):
        logger.error(
            f"No reads aligned to the EsViritu DB in {initial_map_bam}. Exiting..."
        )
        if args.KEEP:
            logger.info(f"keeping temp files in {args.TEMP_DIR}")
        else:
            try:
                if os.path.isdir(args.TEMP_DIR):
                    shutil.rmtree(args.TEMP_DIR)
                    logger.info(f"Removed temporary directory {args.TEMP_DIR}")
            except Exception as e:
                logger.warning(f"Failed to remove temp directory {args.TEMP_DIR}: {e}")
        sys.exit()

    # Load virus info metadata table in polars
    vir_meta_df = pl.read_csv(
        db_metadata, 
        separator='\t',
        schema_overrides={"TaxID": pl.String}
        )

    #########################
    ### 1 of 2 clustering ###
    #########################

    ## read-based clustering, attempt 1
    asm_read_sharing_table_fn = timed_function(logger=logger)(esvf.assembly_read_sharing_table)
    initial_read_comp_df = asm_read_sharing_table_fn(
        initial_map_bam,
        vir_meta_df
    )
    i_read_comp_of = os.path.join(
        str(args.TEMP_DIR),
        f"{str(args.SAMPLE)}.initial_read_comp_compare.tsv"
    )
    initial_read_comp_df.write_csv(
        file = i_read_comp_of,
        separator = "\t"
    )

    cluster_assemblies_by_read_sharinge_fn = timed_function(logger=logger)(esvf.cluster_assemblies_by_read_sharing)
    read_clust_init_df = cluster_assemblies_by_read_sharinge_fn(
        initial_read_comp_df
    )
    i_read_clust_of = os.path.join(
        str(args.TEMP_DIR),
        f"{str(args.SAMPLE)}.initial_read_comp_clust.tsv"
    )
    read_clust_init_df.write_csv(
        file = i_read_clust_of,
        separator = "\t"
    )

    ## make new fasta file database from aniclust exemplars, attempt 1
    clust_record_getter_fn = timed_function(logger=logger)(esvf.clust_record_getter)
    clust_db_fasta = clust_record_getter_fn(
        i_read_clust_of,
        2,
        vir_meta_df,
        db_fasta,
        os.path.join(
            str(args.TEMP_DIR),
            f"{str(args.SAMPLE)}_clustered_refs.fasta"
        )
    )

    ## re-align (mapped) reads to clustered references, attempt 1
    sec_bam_f = os.path.join(str(args.TEMP_DIR), f"{str(args.SAMPLE)}.second.filt.sorted.bam")

    minimap2_f_fn = timed_function(logger=logger)(esvf.minimap2_f)
    second_map_bam = minimap2_f_fn(
        clust_db_fasta,
        trim_filt_reads,
        str(args.CPU),
        sec_bam_f,
        str(args.MMK)
    )

    #########################
    ### 2 of 2 clustering ###
    #########################

    ## read-based clustering, attempt 2
    sec_read_comp_df = asm_read_sharing_table_fn(
        second_map_bam,
        vir_meta_df
    )
    sec_read_comp_of = os.path.join(
        str(args.TEMP_DIR),
        f"{str(args.SAMPLE)}.second_read_comp_compare.tsv"
    )
    sec_read_comp_df.write_csv(
        file = sec_read_comp_of,
        separator = "\t"
    )

    read_clust_sec_df = cluster_assemblies_by_read_sharinge_fn(
        sec_read_comp_df
    )
    sec_read_clust_of = os.path.join(
        str(args.TEMP_DIR),
        f"{str(args.SAMPLE)}.second_read_comp_clust.tsv"
    )
    read_clust_sec_df.write_csv(
        file = sec_read_clust_of,
        separator = "\t"
    )

    ## make new fasta file database from exemplars, attempt 2
    sec_clust_db_fasta = clust_record_getter_fn(
        sec_read_clust_of,
        2,
        vir_meta_df,
        clust_db_fasta,
        os.path.join(
            str(args.TEMP_DIR),
            f"{str(args.SAMPLE)}_second_clustered_refs.fasta"
        )
    )

    ## re-align (mapped) reads to clustered references, attempt 2
    third_bam_f = os.path.join(str(args.TEMP_DIR), f"{str(args.SAMPLE)}.third.filt.sorted.bam")

    third_map_bam = minimap2_f_fn(
        sec_clust_db_fasta,
        trim_filt_reads,
        str(args.CPU),
        third_bam_f,
        str(args.MMK)
    )

    ## take third .bam, make final consensus .fastas
    sec_con_f = os.path.join(
        str(out_directory),
        f"{str(args.SAMPLE)}_final_consensus.fasta"
    )
    bam_to_consensus_fasta_fn = timed_function(logger=logger)(esvf.bam_to_consensus_fasta)
    second_consensus_fasta = bam_to_consensus_fasta_fn(
        third_map_bam,
        sec_con_f
    )

    #########################

    # Make coverm-like table from final bam
    bam_to_coverm_table_fn = timed_function(logger=logger)(esvf.bam_to_coverm_table)
    third_coverm_like_dt = bam_to_coverm_table_fn(
        third_map_bam,
        str(args.SAMPLE),
        int(args.CPU),
        False
    )

    third_coverm_like_dt.write_csv(
        file = os.path.join(
            args.TEMP_DIR,
            f"{str(args.SAMPLE)}_final_coverm.tsv"
        ),
        separator = "\t"
    )

    ## Make bedtools-like windows coverage
    bam_coverage_windows_fn = timed_function(logger=logger)(esvf.bam_coverage_windows)
    windows_cov_df = bam_coverage_windows_fn(
        third_map_bam
    )

    windows_of = os.path.join(
        out_directory,
        f"{str(args.SAMPLE)}.virus_coverage_windows.tsv"
    )

    windows_cov_df.write_csv(
        file = windows_of,
        separator = "\t"
    )

    # get read ANI per contig
    read_ani_df = esvf.read_ani_from_bam(
        third_map_bam
    )
    rani_of = os.path.join(
        args.TEMP_DIR,
        f"{str(args.SAMPLE)}.read_ani.per_contig.tsv"
    )
    read_ani_df.write_csv(
        file = rani_of,
        separator = "\t"
    )

    ## Make "main" and "assembly" output table

    assembly_table_maker_fn = timed_function(logger=logger)(esvf.assembly_table_maker)
    main_out_df, assem_out_df = assembly_table_maker_fn(
        third_coverm_like_dt,
        vir_meta_df,
        read_ani_df,
        filtered_reads,
        str(args.SAMPLE)
    )

    main_of = os.path.join(
        str(out_directory),
        f"{str(args.SAMPLE)}.detected_virus.info.tsv"
    )
    main_out_df.write_csv(
        file = main_of,
        separator = "\t"
    )

    logger.info(f"summary table (per contig): {main_of}")

    assem_of = os.path.join(
        str(out_directory),
        f"{str(args.SAMPLE)}.detected_virus.assembly_summary.tsv"
    )
    assem_out_df.write_csv(
        file = assem_of,
        separator = "\t"
    )
    logger.info(f"summary table (per assembly): {assem_of}")

    tax_fn = timed_function(logger=logger)(esvf.tax_profile)
    tax_df = tax_fn(
        assem_out_df,
        float(args.spthresh),
        float(args.subspthresh)
    )

    tax_of = os.path.join(
        str(out_directory),
        f"{str(args.SAMPLE)}.tax_profile.tsv"
    )
    tax_df.write_csv(
        file = tax_of,
        separator = "\t"
    )

    # check if R script with libraries returns good exit code
    completedProc = subprocess.run(['Rscript', str(esviritu_script_path) + '/check_R_libraries1.R'])

    if completedProc.returncode == 0 :

        # make the html interactive table with sparklines
        logger.info(f"generating reactable report...")
        start_time = time.perf_counter()
        reactableR_path = os.path.join(
            os.path.dirname(__file__), 
            'EsViritu_general_reactable1.R'
            )
        reactablecmd = [
            'Rscript', reactableR_path,
            windows_of,
            main_of,
            str(out_directory),
            str(args.SAMPLE),
            str(filtered_reads)
            ]
        try:
            subprocess.run(reactablecmd, check=True)
        except Exception as e:
            logger.error(f"report generation failed: {reactablecmd}\nError: {e}")
            raise
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        logger.info(f"reactable report finished in {elapsed_time:.2f} seconds")

    else:
        logger.warning("some required R packages are not found. Required:")
        logger.warning("reactable, htmltools, reactablefmtr, scales, magrittr")
        logger.warning("HTML report not being created.")
        logger.warning("Did you activate the conda environment?")
        logger.warning("see EsViritu.yml.")

    # optionally removing temporary files
    if args.KEEP:
        logger.info(f"keeping temp files in {args.TEMP_DIR}")
    else:
        try:
            if os.path.isdir(args.TEMP_DIR):
                shutil.rmtree(args.TEMP_DIR)
                logger.info(f"Removed temporary directory {args.TEMP_DIR}")
        except Exception as e:
            logger.warning(f"Failed to remove temp directory {args.TEMP_DIR}: {e}")

    esv_end_time = time.perf_counter()
    esv_elapsed_time = esv_end_time - esv_start_time
    hours = int(esv_elapsed_time // 3600)
    minutes = int((esv_elapsed_time % 3600) // 60)
    seconds = int(esv_elapsed_time % 60)

    logger.info(f"EsViritu run for {args.SAMPLE} finished in {hours:02d}:{minutes:02d}:{seconds:02d}")

if __name__ == "__main__":
    esviritu()


