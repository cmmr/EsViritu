#!/usr/bin/env python

import argparse
from ntpath import isfile
import sys, os
import subprocess
try:
    from . import esv_funcs as esvf
except:
    import esv_funcs as esvf
import yaml
import logging

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
        type=esvf.str2bool, default='False',
        help='True or False. Remove low-quality reads with fastp?'
        )
    optional_args.add_argument(
        '-f', "--filter_seqs", 
        dest="FILTER_SEQS", type=esvf.str2bool, default='False',
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
        dest="COMPARE", type=esvf.str2bool, default='True',
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
        dest="KEEP", type=esvf.str2bool, default='False',
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
    #formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # stream gets printed to terminal
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    #stream_handler.setFormatter(formatter)

    # file gets saved to a specified file
    log_file_path = os.path.join(
        out_directory, 
        f"esviritu_{args.SAMPLE}.log"
        )
    try:
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.DEBUG)
        #file_handler.setFormatter(formatter)
    except Exception as e:
        logger.error(f"Failed to create log file handler: {e}")
        file_handler = None

    if not logger.hasHandlers():
        if file_handler:
            logger.addHandler(file_handler)
        logger.addHandler(stream_handler)

    logger.propagate = False
    #########################

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

    if args.FILTER_DIR == "default" and os.getenv('ESVIRITU_FILTER') != None:
        args.FILTER_DIR = os.getenv('ESVIRITU_FILTER')
    elif args.FILTER_DIR == "default":
        args.FILTER_DIR = esviritu_script_path.replace("src/EsViritu", "filter_seqs")
    

    logger.info("DB: ", str(args.DB))

    logger.info("filter seqs: ", str(args.FILTER_DIR))

    logger.info("version ", str(__version__))

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

    params_dict = vars(args).copy()  # Copy all argparse arguments
    #params_dict['DB'] = args.DB
    #params_dict['FILTER_DIR'] = args.FILTER_DIR
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
        args.OUTPUT_DIR,
        args.TEMP_DIR,
        args.QUAL,
        args.FILTER_SEQS,
        filter_db_fasta,
        args.READ_FMT,
        args.CPU,
    )

    logger.info(trim_filt_reads)

    # map reads to virus DB and filter for good alignments
    initial_map_bam = esvf.minimap2_f(
        db_index,
        trim_filt_reads,
        args.CPU,
        args.SAMPLE,
        args.TEMP_DIR
    )

    logger.info(initial_map_bam)

    initial_consensus = esvf.bam_to_consensus_fasta(
        initial_map_bam,
        os.path.join(
            args.TEMP_DIR,
            f"{args.SAMPLE}.prelim.consensus.fasta"
        )
    )

    logger.info(initial_consensus)


if __name__ == "__main__":
    esviritu()


