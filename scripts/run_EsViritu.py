#!/usr/bin/env python

import argparse
import sys, os
import subprocess

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        
pathname = os.path.dirname(sys.argv[0])  
esviritu_script_path = os.path.abspath(pathname)      
print(esviritu_script_path) 

parser = argparse.ArgumentParser(description='EsViritu is a read mapping pipeline for detection and measurement of human and animal virus pathogens from short read metagenomic or clinical samples. Version 0.1.1')

required_args = parser.add_argument_group(' REQUIRED ARGUMENTS for EsViritu ')

required_args.add_argument("-r", "--reads", nargs="+", action="append",
                           dest="READS", required=True, 
                           help='read file(s) in .fastq format. You can specify more than one separated by a space')
required_args.add_argument("-s", "--sample", 
                           dest="SAMPLE", type=str, required=True, 
                           help='Sample name. No space characters, please.')
required_args.add_argument("-t", "--cpu", 
                           dest="CPU", type=int, required=True, 
                           help='Example: 32 -- Number of CPUs available for EsViritu. ')
required_args.add_argument("-o", "--output_dir", 
                           dest="OUTPUT_DIR", type=str, required=True, 
                           help='Output directory name. Will be created if it does not exist. Can be shared with other samples. No space characters, please. ')



optional_args = parser.add_argument_group(' OPTIONAL ARGUMENTS for EsViritu.')

optional_args.add_argument('--version', action='version', version='v0.1.1')
optional_args.add_argument('-q', "--qual", dest="QUAL", type=str2bool, default='False',
                           help='True or False. Remove low-quality reads with fastp?')
optional_args.add_argument('-f', "--filter_seqs", dest="FILTER_SEQS", type=str2bool, default='False',
                           help='True or False. Remove reads aligning to sequences at filter_seqs/filter_seqs.fna ?')
optional_args.add_argument('-i', "--compare", dest="COMPARE", type=str2bool, default='False',
                           help='True or False. Calculate percent identity between sample consensus and reference sequences?')
optional_args.add_argument("-m", "--mode", 
                           dest="MODE", type=str, default='general',
                           help='general or curated. See README.md')
optional_args.add_argument("--temp", 
                           dest="TEMP_DIR", type=str, default='default',
                           help='path of temporary directory. Default is {OUTPUT_DIR}/{SAMPLE}_temp/')
optional_args.add_argument("--keep", 
                           dest="KEEP", type=str2bool, default='False',
                           help='True of False. Keep the intermediate files, located in the temporary directory? These can add up, so it is not recommended if space is a concern.')

args = parser.parse_args()

READS = [item for sublist in args.READS for item in sublist]

def is_tool(name):
	"""Check whether `name` is on PATH."""
	from distutils.spawn import find_executable
	return find_executable(name) is not None

if is_tool("coverm") :
	print ("coverm found")
else:
	print ("coverm is not found. Exiting.")
	quit()
if is_tool("minimap2") :
	print ("minimap2 found")
else:
	print ("minimap2 is not found. Exiting.")
	quit()
if is_tool("samtools") :
	print ("samtools found")
else:
	print ("samtools is not found. Exiting.")
	quit()
if is_tool("bioawk") :
	print ("bioawk found")
else:
	print ("bioawk is not found. Exiting.")
	quit()
if is_tool("seqkit") :
	print ("seqkit found")
else:
	print ("seqkit is not found. Exiting.")
	quit()
if is_tool("bedtools") :
	print ("bedtools found")
else:
	print ("bedtools is not found. Exiting.")
	quit()
if is_tool("fastp") :
	print ("fastp found")
else:
	print ("fastp is not found. Exiting.")
	quit()

if str(args.MODE) == "general" :
    subprocess.call(['bash', str(esviritu_script_path) + '/EsViritu_general.sh', 
		     str(READS), str(args.SAMPLE), str(args.CPU), str(args.OUTPUT_DIR), 
	         str(args.QUAL), str(args.FILTER_SEQS), str(args.COMPARE), 
	         str(args.TEMP_DIR), str(args.KEEP), str(esviritu_script_path)])
elif str(args.MODE) == "curated" :
    subprocess.call(['bash', str(esviritu_script_path) + '/EsViritu_curated.sh', 
		     str(READS), str(args.SAMPLE), str(args.CPU), str(args.OUTPUT_DIR), 
	        str(args.QUAL), str(args.FILTER_SEQS), str(args.COMPARE), 
            str(args.TEMP_DIR), str(args.KEEP), str(esviritu_script_path)])	
else:
	print("-m argument must be general or curated.")

