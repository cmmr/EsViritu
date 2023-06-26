#!/bin/bash

# scripts for EsViritu "general" mode. Run from "run_EsViritu.py" wrapper

# arguments

READS=$1
SAMPLE=$2
CPUS=$3
OUT_DIR=$4
QUAL=$5
FILTER_SEQS=$6
FILTER_LOC=$7
COMPARE=$8
TEMP_DIR=$9
KEEP=${10}
READ_FMT=${11}
ES_VERSION=${12}
DB_DIR=${13}
ESVIRITU_DIR=${14}

MDYT=$( date +"%m-%d-%y---%T" )
echo "Time Update: Starting main bash script for EsViritu General Mode @ $MDYT"

#arguments check
if [ $# -ne 14 ] ; then 
	echo "expected 14 arguments passed on the command line:"
	echo "read file(s), sample ID, CPUs, output directory, trim by quality?, filter seqs?, filter seqs directory, "
	echo "compare consensus?, temp directory path, keep temp?, read format, version, DB directory, EsViritu script directory"
	echo "exiting"
	exit
fi

## temporary directory variable setup
if [ "$TEMP_DIR" == "default" ] ; then
	TEMP_DIR="${OUT_DIR}/${SAMPLE}_temp"

elif [ -d $TEMP_DIR ] ; then
	MDYT_g=$( date +"%m-%d-%y---%T" | sed 's/:/_/g' )
	TEMP_DIR="${TEMP_DIR}/${SAMPLE}_temp_${MDYT_g}"

fi

# check filter_seqs
if [ "$FILTER_SEQS" == "True" ] && [ ! -s ${FILTER_LOC}/filter_seqs.fna ]; then
	echo "-f True flag requires that this file exists and is not empty: "
	echo "${FILTER_LOC}/filter_seqs.fna"
	echo "exiting"
	exit
fi


if [ ! -d ${DB_DIR} ] ; then
	echo "can't find DB directory. should be a versioned directory here: ${DB_DIR}/ "
	echo "exiting"
	exit
fi

if [ ! -s ${DB_DIR}/virus_pathogen_database.mmi ] ; then
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
echo "filter seqs directory:      $FILTER_LOC" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Compare consensus to ref:   $COMPARE" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Temp directory path:        $TEMP_DIR" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Keep temp files:            $KEEP" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "Read format:                $READ_FMT" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "EsViritu Version:           $ES_VERSION" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "EsViritu script directory:  $ESVIRITU_DIR" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt
echo "EsViritu DB version:        ${DB_DIR}/virus_pathogen_database.mmi" >> ${OUT_DIR}/record/${SAMPLE}.arguments.txt


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


### parsing paired/qual/pre-filt arguments and process reads accordingly

if [ "$QUAL" == "True" ] && [ "$FILTER_SEQS" == "True" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Trimming low quality reads with fastp THEN Aligning reads to filter_seqs.fna to remove host/spike-in @ $MDYT"

	##trim
	##filter
	if [ "$READ_FMT" == "paired" ] ; then
		READ1=$( echo $READS | cut -d " " -f1 )
		READ2=$( echo $READS | cut -d " " -f2 )

		fastp -i $READ1 -I $READ2 --stdout -w $CPUS -D 1 --html=${OUT_DIR}/record/${SAMPLE}.fastp.html --json=${OUT_DIR}/record/${SAMPLE}.fastp.json | \
		minimap2 -t $CPUS -ax sr ${FILTER_LOC}/filter_seqs.fna - | \
		samtools collate -u -O - | \
		samtools fastq -n -f 4 - > ${TEMP_DIR}/${SAMPLE}.EV_input.fastq			
	else
		cat ${READS} | \
		fastp --stdin --stdout -w $CPUS -D 1 --html=${OUT_DIR}/record/${SAMPLE}.fastp.html --json=${OUT_DIR}/record/${SAMPLE}.fastp.json | \
		minimap2 -t $CPUS -ax sr ${FILTER_LOC}/filter_seqs.fna - | \
		samtools fastq -n -f 4 - > ${TEMP_DIR}/${SAMPLE}.EV_input.fastq
	fi

elif [ "$QUAL" == "True" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Trimming low quality reads with fastp @ $MDYT"

	##trim
	if [ "$READ_FMT" == "paired" ] ; then
		READ1=$( echo $READS | cut -d " " -f1 )
		READ2=$( echo $READS | cut -d " " -f2 )
		fastp -i $READ1 -I $READ2 --stdout -w $CPUS -D 1 --html=${OUT_DIR}/record/${SAMPLE}.fastp.html --json=${OUT_DIR}/record/${SAMPLE}.fastp.html > ${TEMP_DIR}/${SAMPLE}.EV_input.fastq		
	else
		cat ${READS} | \
		fastp --stdin --stdout -w $CPUS -D 1 --html=${OUT_DIR}/record/${SAMPLE}.fastp.html --json=${OUT_DIR}/record/${SAMPLE}.fastp.html > ${TEMP_DIR}/${SAMPLE}.EV_input.fastq
	fi

elif [ "$FILTER_SEQS" == "True" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Aligning reads to filter_seqs.fna to remove host/spike-in @ $MDYT"

	##filter
	if [ "$READ_FMT" == "paired" ] ; then
		minimap2 -t $CPUS -ax sr ${FILTER_LOC}/filter_seqs.fna ${READS} | \
		samtools collate -u -O - | \
		samtools fastq -n -f 4 - > ${TEMP_DIR}/${SAMPLE}.EV_input.fastq
	else
		minimap2 -t $CPUS -ax sr ${FILTER_LOC}/filter_seqs.fna ${READS} | \
		samtools fastq -n -f 4 - > ${TEMP_DIR}/${SAMPLE}.EV_input.fastq
	fi

else
	## cat
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Concatenating input reads @ $MDYT"

	if [ "$READ_FMT" == "paired" ] ; then
		READ1=$( echo $READS | cut -d " " -f1 )
		READ2=$( echo $READS | cut -d " " -f2 )
		seqfu interleave -1 $READ1 -2 $READ2 -o ${TEMP_DIR}/${SAMPLE}.EV_input.fastq
	else
		cat ${READS} > ${TEMP_DIR}/${SAMPLE}.EV_input.fastq
	fi
fi

### main mapping/dereplicating/remapping part

if [ -s ${TEMP_DIR}/${SAMPLE}.EV_input.fastq ] ; then
	seqkit stats -T ${TEMP_DIR}/${SAMPLE}.EV_input.fastq > ${OUT_DIR}/${SAMPLE}.EV_input.seq_stats.tsv
	FILTERED_READS=$( tail -n1 ${OUT_DIR}/${SAMPLE}.EV_input.seq_stats.tsv | cut -f4  )

	# run initial coverm
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Running CoverM @ $MDYT"
	
	if [ "$READ_FMT" == "paired" ] ; then 
		coverm contig --interleaved ${TEMP_DIR}/${SAMPLE}.EV_input.fastq -r ${DB_DIR}/virus_pathogen_database.mmi --minimap2-reference-is-index --min-read-percent-identity 90 --min-read-aligned-percent 90 -m length covered_bases count mean -o ${TEMP_DIR}/${SAMPLE}.coverm.tsv -t $CPUS --bam-file-cache-directory ${TEMP_DIR} --discard-unmapped
	else
		coverm contig --single ${TEMP_DIR}/${SAMPLE}.EV_input.fastq -r ${DB_DIR}/virus_pathogen_database.mmi --minimap2-reference-is-index --min-read-percent-identity 90 --min-read-aligned-percent 90 -m length covered_bases count mean -o ${TEMP_DIR}/${SAMPLE}.coverm.tsv -t $CPUS --bam-file-cache-directory ${TEMP_DIR} --discard-unmapped
	fi

	# Dereplicate hits from initial coverm
	# Checking for contigs with 1000bp of coverage or 50% breadth of coverage
	LINECOUNT=$( tail -n+2 ${TEMP_DIR}/${SAMPLE}.coverm.tsv | awk '{OFS=FS="\t"}{ if (($3/$2) >= 0.5 || $3 >= 1000) {print} }'  | wc -l )
	if [ $LINECOUNT -ge 1 ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Time Update: Getting preliminary consensus sequences and all-vs-all comparison @ $MDYT"
		# list of accessions with >1000 nt or >50% coverage (passing threshold)
		tail -n+2 ${TEMP_DIR}/${SAMPLE}.coverm.tsv | awk '{OFS=FS="\t"}{ if ($3>=1000 || ($3/$2)>=0.5) {print $1} }' > ${TEMP_DIR}/${SAMPLE}.accessions_prelim_threshold1.txt

		# sort bam files, keep alignments to accesstions passing threshold, get consensus fastas of those accessions
		## samtools view, consensus, sed NNs
		samtools sort -@ $CPUS ${TEMP_DIR}/virus_pathogen_database.mmi.${SAMPLE}.EV_input.fastq.bam -o ${TEMP_DIR}/virus_pathogen_database.mmi.${SAMPLE}.EV_input.fastq.sort.bam
		samtools index ${TEMP_DIR}/virus_pathogen_database.mmi.${SAMPLE}.EV_input.fastq.sort.bam
		REFS=$( cat ${TEMP_DIR}/${SAMPLE}.accessions_prelim_threshold1.txt | tr "\n" " " )
		samtools view -@ $CPUS -b ${TEMP_DIR}/virus_pathogen_database.mmi.${SAMPLE}.EV_input.fastq.sort.bam ${REFS} > ${TEMP_DIR}/${SAMPLE}.accessions_prelim_threshold1.bam
		samtools consensus -@ $CPUS -f fasta ${TEMP_DIR}/${SAMPLE}.accessions_prelim_threshold1.bam | bioawk -c fastx '{print ">"$name ; print $seq}' |sed 's/NN//g' >  ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.fasta

		# pairwise comparison of consensus sequences
		## blastn -ungapped. anicalc, aniclust 98%ANI, 20% AF
		blastn -query ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.fasta -subject ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.fasta -outfmt '6 std qlen slen' -max_target_seqs 10000 -out ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.blastn.tsv -ungapped

		# calculate clusters of closely related sequences. Keep exemplars from each cluster
		python ${ESVIRITU_DIR}/utils/anicalc.py -i ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.blastn.tsv -o ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.anicalc.tsv
		python ${ESVIRITU_DIR}/utils/aniclust.py --fna ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.fasta --ani ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.anicalc.tsv --out ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.aniclust.98ANI_20AF.tsv --min_ani 98 --min_tcov 20 --min_qcov 0
		
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Time Update: Realigning reads to best references @ $MDYT"
		# Generate fasta of just the exemplars (dereplicated) sequences. Make minimap index
		cut -f1 ${TEMP_DIR}/${SAMPLE}.prelim.consensus.noNNs.aniclust.98ANI_20AF.tsv > ${TEMP_DIR}/${SAMPLE}.prelim.derep_accessions1.txt
		seqkit grep -f ${TEMP_DIR}/${SAMPLE}.prelim.derep_accessions1.txt ${DB_DIR}/virus_pathogen_database.fna > ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.fna

		minimap2 -d ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.fna -t $CPUS

		# make fastq file from bam file of initial coverm alignment
		##	4) samtools fastq on mapped reads, coverm on these reads to new minimap index
		samtools fastq -@ $CPUS ${TEMP_DIR}/virus_pathogen_database.mmi.${SAMPLE}.EV_input.fastq.bam > ${TEMP_DIR}/${SAMPLE}.mapped_prelim.fastq

		# run coverm with just previously aligned reads and dereplicated reference sequences
		if [ "$READ_FMT" == "paired" ] ; then 
			coverm contig --interleaved ${TEMP_DIR}/${SAMPLE}.mapped_prelim.fastq -r ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi --minimap2-reference-is-index --min-read-percent-identity 90 --min-read-aligned-percent 90 -m length covered_bases count mean -o ${TEMP_DIR}/${SAMPLE}.best.coverm.tsv -t $CPUS --bam-file-cache-directory ${TEMP_DIR} --discard-unmapped
		else
			coverm contig --single ${TEMP_DIR}/${SAMPLE}.mapped_prelim.fastq -r ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi --minimap2-reference-is-index --min-read-percent-identity 90 --min-read-aligned-percent 90 -m length covered_bases count mean -o ${TEMP_DIR}/${SAMPLE}.best.coverm.tsv -t $CPUS --bam-file-cache-directory ${TEMP_DIR} --discard-unmapped
		fi
		
		FINAL_HITS=$( tail -n+2 ${TEMP_DIR}/${SAMPLE}.best.coverm.tsv | awk '{OFS=FS="\t"}{ if ($3>=1000 || ($3/$2)>=0.5) {print} }' | wc -l )

		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Time Update: computing final consensus sequences @ $MDYT"
		# add another samtools consensus for final derep alignments
		samtools sort -@ $CPUS ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.bam -o ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.sort.bam
		samtools index ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.sort.bam
		F_CONS=$( tail -n+2 ${TEMP_DIR}/${SAMPLE}.best.coverm.tsv | awk '{OFS=FS="\t"}{ if ($3>=1000 || ($3/$2)>=0.5) {print $1} }' | tr "\n" " " )
		samtools view -@ $CPUS -b ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.sort.bam ${F_CONS} > ${TEMP_DIR}/accessions_final_threshold1.bam
		samtools consensus -@ $CPUS -f fasta ${TEMP_DIR}/accessions_final_threshold1.bam | bioawk -c fastx '{print ">"$name"__consensus" ; print $seq}' >  ${OUT_DIR}/${SAMPLE}.final.consensus.with_NNs.fasta

		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Time Update: Generating detection and abundance table @ $MDYT"

		# make an info table on all the genomes detected in the sample based on the coverm stats, number of filtered reads, and the metadata table Comp_viral_panel_w_curated_plus_extra_v2.0.2.all_metadata.tsv
		echo -e "accession\treference_length\tcovered_bases\treads_aligned\tmean_coverage\tRPKMF\tsequence_name\ttaxid\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tsubspecies\tstrain\tsample_ID\ttotal_filtered_reads_in_sample" > ${OUT_DIR}/${SAMPLE}.detected_virus.info.tsv
		tail -n+2 ${TEMP_DIR}/${SAMPLE}.best.coverm.tsv | awk -v readq="$FILTERED_READS" '{OFS=FS="\t"}{ if (($3/$2) >= 0.5 || $3 >= 1000) {print $0, (($4/($2/1000)/readq)*1000000)} }' | while read LINE ; do 
			ACC=$( echo "$LINE" | cut -f1 )

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

### calculate coverage across genome windows
if [ -s ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.fna ] && [ -s ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.sort.bam ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Caclulating average coverage across detected viruses @ $MDYT"

	bioawk -c fastx '{print $name, length($seq)}' ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.fna > ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.genome.tsv

	bedtools makewindows -g ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.genome.tsv -n 100 > ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.100windows.bed

	bedtools coverage -mean -a ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.100windows.bed -b ${TEMP_DIR}/${SAMPLE}.best_ref_seqs.mmi.${SAMPLE}.mapped_prelim.fastq.sort.bam > ${OUT_DIR}/${SAMPLE}.100windows.mean_cov.tsv
else
	echo "Can't find .bam or .fna reference for coverage calculation"
fi

### run the R script that makes a reactable for each sequence passing the threshold
if [ -s ${OUT_DIR}/${SAMPLE}.100windows.mean_cov.tsv ] && [ -s ${OUT_DIR}/${SAMPLE}.detected_virus.info.tsv ] ; then
	

	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Time Update: Running reactablefmtr to make html table @ $MDYT"

	Rscript ${ESVIRITU_DIR}/EsViritu_general_reactable1.R ${OUT_DIR}/${SAMPLE}.100windows.mean_cov.tsv ${OUT_DIR}/${SAMPLE}.detected_virus.info.tsv ${OUT_DIR} ${SAMPLE}

else
	echo "either ${OUT_DIR}/${SAMPLE}.100windows.mean_cov.tsv or ${OUT_DIR}/${SAMPLE}.detected_virus.info.tsv were not found. Could not make reactable."
fi

### compare consensus to reference
if [ "$COMPARE" == "True" ] ; then
	##
	if [ -s ${OUT_DIR}/${SAMPLE}.final.consensus.with_NNs.fasta ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Time Update: comparing consensus seqs to references @ $MDYT"

		sed 's/NN//g' ${OUT_DIR}/${SAMPLE}.final.consensus.with_NNs.fasta | \
		blastn -query - -subject ${DB_DIR}/virus_pathogen_database.fna -outfmt '6 std qlen slen' -max_target_seqs 10000 -out ${TEMP_DIR}/${SAMPLE}.consensus.remove_NNs.VS.refs.blastn.tsv

		python ${ESVIRITU_DIR}/utils/anicalc.py -i ${TEMP_DIR}/${SAMPLE}.consensus.remove_NNs.VS.refs.blastn.tsv -o ${TEMP_DIR}/${SAMPLE}.consensus.remove_NNs.VS.refs.anicalc.tsv

		if [ -s ${TEMP_DIR}/${SAMPLE}.consensus.remove_NNs.VS.refs.anicalc.tsv ] ; then
			Rscript ${ESVIRITU_DIR}/Consensuses_vs_Refs_comparison.R ${DB_DIR}/virus_pathogen_database.all_metadata.tsv ${TEMP_DIR}/${SAMPLE}.consensus.remove_NNs.VS.refs.anicalc.tsv ${OUT_DIR} ${SAMPLE}
		else
			echo "${TEMP_DIR}/${SAMPLE}.consensus.remove_NNs.VS.refs.anicalc.tsv NOT FOUND"
			echo "Cannot generate Consensuses_vs_Refs_comparison table"
		fi

	else
		echo "${OUT_DIR}/${SAMPLE}.final.consensus.with_NNs.fasta NOT FOUND"
		
	fi
fi

### delete temp
if [ "$KEEP" == "True" ] ; then
	echo "Keeping temporary files in: ${TEMP_DIR}"
else
	echo "Removing temp files"
	rm ${TEMP_DIR}/*
fi


echo "Main Output Files Generated in ${OUT_DIR}/"

find ${OUT_DIR}/ -type f -name "${SAMPLE}.detected_virus.info.tsv"
find ${OUT_DIR}/ -type f -name "${SAMPLE}.final.consensus.with_NNs.fasta"
find ${OUT_DIR}/ -type f -name "${SAMPLE}_EsViritu_reactable.html"
find ${OUT_DIR}/ -type f -name "${SAMPLE}_consensus_seqs_vs_ref_seqs.tsv"


MDYT=$( date +"%m-%d-%y---%T" )
echo "EsViritu general mode is finishing @ $MDYT"

