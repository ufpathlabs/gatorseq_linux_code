#!/bin/bash
IFS=$'\n'
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG_FILE="$SCRIPT_DIR/linux_gatorseq.config.yaml"


HPC_SFTP=$(cat $CONFIG_FILE |grep "^HPC_SFTP:" |cut -d"'" -f2)
LINUX_PATHOLOGY_FASTQ_FOLDER=$(cat $CONFIG_FILE |grep "^LINUX_PATHOLOGY_FASTQ_FOLDER:" |cut -d"'" -f2)
HPC_FASTQ_FOLDER=$(cat $CONFIG_FILE |grep "^HPC_FASTQ_FOLDER:" |cut -d"'" -f2)

#echo $HPC_SFTP
#echo $LINUX_PATHOLOGY_FASTQ_FOLDER
#echo $HPC_FASTQ_FOLDER


echo $(date)

if [ ! -d $LINUX_PATHOLOGY_FASTQ_FOLDER ]; then
	echo "ERROR: $LINUX_PATHOLOGY_FASTQ_FOLDER folder does not exist."
	exit 1
fi

#if [ ! -d $LINUX_HPC_FASTQ_FOLDER ]; then
#	echo "ERROR: $LINUX_HPC_FASTQ_FOLDER folder does not exit."
#	exit 1
#fi

#SYNC_FLAG=$(eval rsync -ltgoDvzrn --progress --stats $LINUX_PATHOLOGY_FASTQ_FOLDER/ $LINUX_HPC_FASTQ_FOLDER/ |grep "Total transferred file size: 0 bytes")
#rsync_dry="rsync -ltgoDvzrn --progress --stats $LINUX_PATHOLOGY_FASTQ_FOLDER/ $LINUX_HPC_FASTQ_FOLDER/ "
#SYNC_FLAG=$(eval $rsync_dry |grep "Total transferred file size: 0 bytes")
#echo $rsync_dry
#eval $rsync_dry
#echo "" 
#echo "" 
#SYNC_FLAG=$(eval $rsync_dry | grep "Total transferred file size: 0 bytes")

rsync_dry="rsync -ltgoDvzrn --progress --stats $LINUX_PATHOLOGY_FASTQ_FOLDER/ $HPC_SFTP:$HPC_FASTQ_FOLDER/ "
SYNC_FLAG=$(eval $rsync_dry |grep "Total transferred file size: 0 bytes")
echo $rsync_dry
#eval $rsync_dry
echo "" 


echo "SYNC_FLAG"
echo $SYNC_FLAG

if [ "$SYNC_FLAG" == "" ]; then
    #echo "Need to sync folders $LINUX_PATHOLOGY_FASTQ_FOLDER $LINUX_HPC_FASTQ_FOLDER"
    #echo "rsync -ltgoDvzr --progress --stats $LINUX_PATHOLOGY_FASTQ_FOLDER/ $HPC_SFTP:$HPC_FASTQ_FOLDER/ "
    echo "Need to sync folders $LINUX_PATHOLOGY_FASTQ_FOLDER $HPC_SFTP:$HPC_FASTQ_FOLDER"
    rsync_cmd="rsync -ltgoDvzr --progress --stats $LINUX_PATHOLOGY_FASTQ_FOLDER/ $HPC_SFTP:$HPC_FASTQ_FOLDER/ "
    echo $rsync_cmd
    eval $rsync_cmd
else 
    echo "Folders are already in sync $LINUX_PATHOLOGY_FASTQ_FOLDER $HPC_SFTP:$HPC_FASTQ_FOLDER"
fi

echo "" 

