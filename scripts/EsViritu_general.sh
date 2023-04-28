#!/bin/bash

# scripts for EsViritu "general" mode. Run from "run_EsViritu.py" wrapper

# arguments

READS=$1
SAMPLE=$2
CPUS=$3
OUT_DIR=$4
QUAL=$5
FILTER_SEQS=$6
COMPARE=$7
TEMP_DIR=$8
KEEP=$9
ESVIRITU_DIR=${10}

MDYT=$( date +"%m-%d-%y---%T" )
echo "Time Update: Starting main bash script for EsViritu General Mode @ $MDYT"

#arguments check
if [ $# -ne 10 ] ; then 
	echo "expected 10 arguments passed on the command line:"
	echo "read file(s), sample ID, CPUs, output directory, trim by quality?, filter seqs?, compare consensus?, temp directory path, keep temp?, EsViritu script directory"
	echo "exiting"
	exit
fi

## temporary directory variable setup
if [ "$TEMP_DIR" == "default" ] ; then
	TEMP_DIR="${OUT_DIR}/${SAMPLE}_temp"
fi

# check filter_seqs
if [ "$FILTER_SEQS" == "True" ] && [ ! -s ${ESVIRITU_DIR%scripts}filter_seqs/filter_seqs.fna ]; then
	echo "-f True flag requires that this file exists and is not empty: "
	echo "${ESVIRITU_DIR%scripts}filter_seqs/filter_seqs.fna"
	echo "exiting"
	exit
fi

DB_DIR=$( find ${ESVIRITU_DIR%scripts}DBs/ -type d | tail -n1 )

if [ ! -d ${DB_DIR} ] ; then
	echo "can't find DB directory. should be a versioned directory here: ${ESVIRITU_DIR%scripts}DBs/ "
	echo "exiting"
	exit
fi

if [ ! -s ${DB_DIR}virus_pathogen_database.mmi ] ; then
	echo "can't find minimap2 index for pipeline. should be: ${DB_DIR}/virus_pathogen_database.mmi"
	echo "exiting"
	exit
fi


# check output directories
if [ ! -d ${OUT_DIR} ] ; then
	mkdir ${OUT_DIR}
fi
if [ ! -d ${OUT_DIR}/record/ ] ; then
	mkdir ${OUT_DIR}/record/
fi
if [ ! -d ${TEMP_DIR} ] ; then
	mkdir ${TEMP_DIR}
fi

# arguments file
date > ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "EsViritu general mode arguments" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Sample name:                $SAMPLE" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Read file(s):               $READS" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "CPUs:                       $CPUS" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Output directory:           $OUT_DIR" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Trim for quality:           $QUAL" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Remove host/spikein seqs:   $FILTER_SEQS" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Compare consensus to ref:   $COMPARE" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Temp directory path:        $TEMP_DIR" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Keep temp files:            $KEEP" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "EsViritu script directory:  $ESVIRITU_DIR" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt

cat ${OUT_DIR}/record/${SAMPLE}.arguments.txt

for READ_FILE in $READS ; do
	if [ -s $READ_FILE ] ; then
		echo $READ_FILE
	else
		echo "$READ_FILE not found."
		echo "exiting"
		exit
	fi
done

if [ "$QUAL" == "True" ] && [ "$FILTER_SEQS" == "True" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Trimming low quality reads with fastp THEN Aligning reads to filter_seqs.fna to remove host/spike-in @ $MDYT"

	##trim
	##filter
	cat "${READS}" | \
	fastp --stdin --stdout -w $CPUS -D 1 --html=${OUT_DIR}/record/${SAMPLE}.fastp.html --json=${OUT_DIR}/record/${SAMPLE}.fastp.html | \
	minimap2 -t $CPUS -ax sr ${ESVIRITU_DIR%scripts}filter_seqs/filter_seqs.fna - | \
	samtools fastq -n -f 4 - > ${TEMP_DIR}/${SAMPLE}.EV_input.fastq

elif [ "$QUAL" == "True" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Trimming low quality reads with fastp @ $MDYT"

	##trim
	cat "${READS}" | \
	fastp --stdin -o ${TEMP_DIR}/${SAMPLE}.EV_input.fastq -w $CPUS -D 1 --html=${OUT_DIR}/record/${SAMPLE}.fastp.html --json=${OUT_DIR}/record/${SAMPLE}.fastp.html

elif [ "$FILTER_SEQS" == "True" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Aligning reads to filter_seqs.fna to remove host/spike-in @ $MDYT"

	##filter
	cat "${READS}" | 
	minimap2 -t $CPUS -ax sr ${ESVIRITU_DIR%scripts}filter_seqs/filter_seqs.fna - | samtools fastq -n -f 4 - > ${TEMP_DIR}/${SAMPLE}.EV_input.fastq

else
	## cat
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Concatenating input reads @ $MDYT"

	cat "${READS}" > ${TEMP_DIR}/${SAMPLE}.EV_input.fastq
fi

if [ -s ${TEMP_DIR}/${SAMPLE}.EV_input.fastq ] ; then
	seqkit stats -T ${TEMP_DIR}/${SAMPLE}.EV_input.fastq > ${OUT_DIR}/${SAMPLE}.EV_input.seq_stats.tsv
	FILTERED_READS=$( tail -n1 ${OUT_DIR}/${SAMPLE}.EV_input.seq_stats.tsv | cut -f4  )

	# run initial coverm
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Running CoverM @ $MDYT"
	
	coverm contig --single ${TEMP_DIR}/${SAMPLE}.EV_input.fastq -r ${DB_DIR}virus_pathogen_database.mmi --minimap2-reference-is-index --min-read-percent-identity 90 --min-read-aligned-percent 90 -m length covered_bases count mean -o ${OUT_DIR}/${SAMPLE}.coverm.tsv -t $CPUS --bam-file-cache-directory ${TEMP_DIR} --discard-unmapped

	# Dereplicate hits from initial coverm
	# Checking for contigs with 1000bp of coverage or 50% breadth of coverage
	LINECOUNT=$( tail -n+2 ${OUT_DIR}/${SAMPLE}.coverm.tsv | awk '{OFS=FS="\t"}{ if (($3/$2) >= 0.5 || $3 >= 1000) {print} }'  | wc -l )
	if [ $LINECOUNT -ge 1 ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Time Update: Getting preliminary consensus sequences and all-vs-all comparison @ $MDYT"
		# list of accessions with >1000 nt or >50% coverage (passing threshold)
		tail -n+2 ${OUT_DIR}/${SAMPLE}.coverm.tsv | awk '{OFS=FS="\t"}{ if ($3>=1000 || ($3/$2)>=0.5) {print $1} }' > ${TEMP_DIR}/${SAMPLE}.accessions_prelim_threshold1.txt

		# sort bam files, keep alignments to accesstions passing threshold, get consensus fastas of those accessions
		## samtools view, consensus, sed NNs
		samtools sort -@ $CPUS ${OUT_DIR}/virus_pathogen_database.mmi.${SAMPLE}.fastq.bam -o ${TEMP_DIR}/virus_pathogen_database.mmi.${SAMPLE}.fastq.sort.bam
		samtools index ${TEMP_DIR}/virus_pathogen_database.mmi.${SAMPLE}.fastq.sort.bam
		REFS=$( cat ${TEMP_DIR}/${SAMPLE}.accessions_prelim_threshold1.txt | tr "\n" " " )
		samtools view -@ $CPUS -b ${TEMP_DIR}/virus_pathogen_database.mmi.${SAMPLE}.fastq.sort.bam ${REFS} > ${TEMP_DIR}/${SAMPLE}.accessions_prelim_threshold1.bam
		samtools consensus -@ $CPUS -f fasta ${TEMP_DIR}/${SAMPLE}.accessions_prelim_threshold1.bam | bioawk -c fastx '{print ">"$name ; print $seq}' |sed 's/NN//g' >  ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.fasta

		# pairwise comparison of consensus sequences
		## blastn -ungapped. anicalc, aniclust 98%ANI, 20% AF
		blastn -query ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.fasta -subject ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.fasta -outfmt '6 std qlen slen' -max_target_seqs 10000 -out ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.blastn.tsv -num_threads $CPUS -ungapped

		# calculate clusters of closely related sequences. Keep exemplars from each cluster
		python ${ESVIRITU_DIR}/anicalc.py -i ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.blastn.tsv -o ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.anicalc.tsv
		python ${ESVIRITU_DIR}/aniclust.py --fna ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.fasta --ani ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.anicalc.tsv --out ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.aniclust.98ANI_20AF.tsv --min_ani 98 --min_tcov 20 --min_qcov 0
		
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Time Update: Realigning reads to best references @ $MDYT"
		# Generate fasta of just the exemplars (dereplicated) sequences. Make minimap index
		cut -f1 ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.aniclust.98ANI_20AF.tsv > ${TEMP_DIR}/${SAMPLE}.prelim.derep_accessions1.txt
		seqkit grep -f ${TEMP_DIR}/${SAMPLE}.prelim.derep_accessions1.txt ${DB_DIR}/virus_pathogen_database.fna > ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.fna

		minimap2 -d ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.fna -t $CPUS

		# make fastq file from bam file of initial coverm alignment
		##	4) samtools fastq on mapped reads, coverm on these reads to new minimap index
		samtools fastq -@ $CPUS ${OUT_DIR}/virus_pathogen_database.mmi.${SAMPLE}.fastq.bam > ${TEMP_DIR}/${SAMPLE}.mapped_prelim.fastq

		# run coverm with just previously aligned reads and dereplicated reference sequences
		coverm contig --interleaved ${TEMP_DIR}/${SAMPLE}.mapped_prelim.fastq -r ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi --minimap2-reference-is-index --min-read-percent-identity 90 --min-read-aligned-percent 90 -m length covered_bases count mean -o ${TEMP_DIR}/${SAMPLE}.best.coverm.tsv -t $CPUS --bam-file-cache-directory ${TEMP_DIR} --discard-unmapped
		FINAL_HITS=$( tail -n+2 ${TEMP_DIR}/${SAMPLE}.best.coverm.tsv | awk '{OFS=FS="\t"}{ if ($3>=1000 || ($3/$2)>=0.5) {print} }' | wc -l )

		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Time Update: computing final consensus sequences @ $MDYT"
		# add another samtools consensus for final derep alignments
		samtools sort -@ $CPUS ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.bam -o ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.sort.bam
		samtools index ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.sort.bam
		F_CONS=$( tail -n+2 ${TEMP_DIR}/${SAMPLE}.best.coverm.tsv | awk '{OFS=FS="\t"}{ if ($3>=1000 || ($3/$2)>=0.5) {print $1} }' | tr "\n" " " )
		samtools view -@ $CPUS -b ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.sort.bam ${F_CONS} > ${TEMP_DIR}/accessions_final_threshold1.bam
		samtools consensus -@ $CPUS -f fasta ${TEMP_DIR}/accessions_final_threshold1.bam >  ${TEMP_DIR}/${OUT_DIR}.final.consensus.with_NNs.fasta

		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Time Update: Generating detection and abundance table @ $MDYT"

		# make an info table on all the genomes detected in the sample based on the coverm stats, number of filtered reads, and the metadata table Comp_viral_panel_w_curated_plus_extra_v2.0.2.all_metadata.tsv
		echo -e "accession\treference_length\tcovered_bases\treads_aligned\tmean_coverage\tRPKM\tsequence_name\ttaxid\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tsubspecies\tstrain\tSAMPLE\ttotal_filtered_reads_in_sample" > ${OUT_DIR}/${SAMPLE}.detected_virus.info.tsv
		tail -n+2 ${TEMP_DIR}/${SAMPLE}.best.coverm.tsv | awk -v readq="$FILTERED_READS" '{OFS=FS="\t"}{ if (($3/$2) >= 0.5 || $3 >= 1000) {print $0, (($4/($2/1000)/readq)*1000000)} }' | while read LINE ; do 
			ACC=$( echo "$LINE" | cut -f1 )
			### update this file
			TAXINFO=$( grep "^$ACC" ${DB_DIR}/virus_pathogen_database.all_metadata.tsv | cut -f 3- )
			echo -e "${LINE}\t${TAXINFO}\t${SAMPLE}\t${FILTERED_READS}"
		done >> ${OUT_DIR}/${SAMPLE}.detected_virus.info.tsv
	else
		echo "no reference genome had more than 1000 positions covered or >=50% length covered. Skipping report."
	fi



else
	##
	echo "${TEMP_DIR}/${SAMPLE}.EV_input.fastq"
	echo "cannot generate reports"

fi

if [ -s ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.fna ] && [ -s ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.sort.bam ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Caclulating average coverage across detected viruses @ $MDYT"

	bioawk -c fastx '{print $name, length($seq)}' ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.fna > ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.genome.tsv

	bedtools makewindows -g ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.genome.tsv -n 100 > ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.100windows.bed

	bedtools coverage -mean -a ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.100windows.bed -b ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.sort.bam > ${OUT_DIR}/${SAMPLE}.100windows.mean_cov.tsv
else
	echo "Can't find .bam or .fna reference for coverage calculation"
fi

if [ -s ${OUT_DIR}/${SAMPLE}.100windows.mean_cov.tsv ] && [ -s ${OUT_DIR}/${SAMPLE}.detected_virus.info.tsv ] ; then
	# runs the R script that makes a reactable for each sequence passing the threshold

	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Running reactablefmtr to make html table @ $MDYT"

	Rscript ${ESVIRITU_DIR}/EsViritu_general_reactable1.R ${OUT_DIR}/${SAMPLE}.100windows.mean_cov.tsv ${OUT_DIR}/${SAMPLE}.detected_virus.info.tsv ${OUT_DIR} ${SAMPLE_ID}

else
	echo "either ${OUT_DIR}/${SAMPLE}.100windows.mean_cov.tsv or ${OUT_DIR}/${SAMPLE}.detected_virus.info.tsv were not found. Could not make reactable."
fi

### ! insert compare consensus to reference

# delete temp
if [ "$KEEP" == "True" ] ; then
	echo "Keeping temporary files in: ${TEMP_DIR}"
else
	rm ${TEMP_DIR}/*
fi

MDYT=$( date +"%m-%d-%y---%T" )
echo "EsViritu general mode is finishing @ $MDYT"

