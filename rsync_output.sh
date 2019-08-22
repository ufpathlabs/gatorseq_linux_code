#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG_FILE="$SCRIPT_DIR/linux_gatorseq.config.yaml"

LINUX_HPC_ANALYSIS_FOLDER=$(cat $CONFIG_FILE |grep "^LINUX_HPC_ANALYSIS_FOLDER:" |cut -d"'" -f2)

HPC_ANALYSIS_FOLDER=$(cat $CONFIG_FILE |grep "^HPC_ANALYSIS_FOLDER:" |cut -d"'" -f2)
LINUX_ANALYSIS_OUT_FOLDER=$(cat $CONFIG_FILE |grep "^LINUX_ANALYSIS_OUT_FOLDER:" |cut -d"'" -f2)
HPC_SFTP=$(cat $CONFIG_FILE |grep "^HPC_SFTP:" |cut -d"'" -f2)


echo $(date)

if [ ! -d $LINUX_HPC_ANALYSIS_FOLDER ];then
	echo "ERROR: $LINUX_HPC_ANALYSIS_FOLDER folder does not exit."
	exit 1
fi;

if [ ! -d $LINUX_ANALYSIS_OUT_FOLDER ];then
	echo "ERROR: $LINUX_ANALYSIS_OUT_FOLDER folder does not exit."
	exit 1
fi;

IFS=$'\n'

for d in $LINUX_HPC_ANALYSIS_FOLDER/*; do
    if [[ -d "$d" && ! -L "$d" ]]; then
        echo $d
        f_s=$d\/SUCCESS.txt
        f_f=$d\/FAILED.txt
        if [ -f "$f_s" ] || [ -f "$f_f" ]; then
            #rm -f -r $d/.snakemake
            RUN_NAME=$(echo $d  | rev | cut -d '/' -f1 |rev | cut -d '_' -f1 )
            SAMPLE_FOLDER=$(echo $d  | rev | cut -d '/' -f1 |rev )
            RUN_GO_FOLDER=$LINUX_ANALYSIS_OUT_FOLDER/$RUN_NAME
            HPC_SAMPLE_FOLDER=$HPC_ANALYSIS_FOLDER/$SAMPLE_FOLDER

            mkdir -p $RUN_GO_FOLDER

            #syncronize output folder, if already not synced 
            #rsync_dry="rsync -ltgoDvzrn --progress --stats $d/ $RUN_GO_FOLDER/$SAMPLE_FOLDER/"
            rsync_dry="rsync -ltgoDvzrn --progress --stats $HPC_SFTP:$HPC_SAMPLE_FOLDER/ $RUN_GO_FOLDER/$SAMPLE_FOLDER/ "
            SYNC_FLAG=$(eval $rsync_dry |grep "Total transferred file size: 0 bytes")
            echo $rsync_dry
            echo $SYNC_FLAG

            if [ "$SYNC_FLAG" == "" ]; then
                mv $d\.cronjob.sh  $d/
                mv $d\.cronjob.sh.slurm*  $d/
                #mv $d\.slurm.* $d/
                rm $d\.cronjob.sh.log
                rsync_cmd="rsync -ltgoDvzr --progress --stats $HPC_SFTP:$HPC_SAMPLE_FOLDER/ $RUN_GO_FOLDER/$SAMPLE_FOLDER/ "
                echo $rsync_cmd
                eval $rsync_cmd
            else
                echo "rm -r $d " 
                rm -r $d 
            fi

        fi
    fi
done

echo "" 
