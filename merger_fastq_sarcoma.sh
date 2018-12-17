#!/bin/sh

## this script merges illumina NGS fasq files into one fastq file for Sarcomas

echo $(date)

FASTQ_FOLDER="/ext/path/DRL/Molecular/NGS/NextSeq_Fastq"
OUTPUT_FASTQ_FOLDER="/ext/path/DRL/Molecular/NGS/Sarcoma_Merged_Fastq"

folders=`find $FASTQ_FOLDER -type d -name "*_?????[ABC]??_*"`


for fold in $folders; do
	out_folder=`echo $fold |rev| cut -d'/' -f 1 | rev` 
	file_prefix=`echo $fold |rev| cut -d'/' -f 1 | rev|cut -d'_' -f-2`
	OUT_DIR=$OUTPUT_FASTQ_FOLDER/$out_folder
	#FINAL_OUT_DIR=$OUTPUT_FASTQ_FOLDER/$out_folder

	if [ -d "$OUT_DIR" ]; then
		continue
	fi

	mkdir $OUT_DIR
    STATUS_FILE=$OUT_DIR/"status.txt"
	R1_FILE=$OUT_DIR/$file_prefix"_R1.fastq.gz"
	R2_FILE=$OUT_DIR/$file_prefix"_R2.fastq.gz"

    echo "Merging Job Started: $(date)" >$STATUS_FILE
	#echo $fold
	#echo $out_folder
	#echo $file_prefix

	R1_files=`ls $fold/* | grep "_R1_" | sort -n`
	R2_files=`ls $fold/* | grep "_R2_" | sort -n`
	#echo $R1_files
	#echo $R2_files

	cat $R1_files >$R1_FILE
	cat $R2_files >$R2_FILE
	R1_COUNT=`zcat $R1_FILE  |wc -l`
	R2_COUNT=`zcat $R2_FILE  |wc -l`

	if [ $R1_COUNT -ne $R2_COUNT ]
	then
		echo "R1 and R2 files do not have same number of reads" 
		rm -r $OUT_DIR
	fi

    echo "Merging Job Successfully Done: $(date)" >>$STATUS_FILE

done

