#!/bin/bash

CODE_PROD_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo $(date)
echo "" 

/usr/bin/flock -n $CODE_PROD_FOLDER/LOCKS/gatorseq_input_excel_to_db.py.lock /home/path-svc-mol/Software/miniconda3/bin/python $CODE_PROD_FOLDER/gatorseq_input_excel_to_db.py &>> $CODE_PROD_FOLDER/LOGS/gatorseq_input_excel_to_db.py.log.txt 

/usr/bin/flock -n $CODE_PROD_FOLDER/LOCKS/gatorseq_fastq_to_analyze.py.lock /home/path-svc-mol/Software/miniconda3/bin/python $CODE_PROD_FOLDER/gatorseq_fastq_to_analyze.py >> $CODE_PROD_FOLDER/LOGS/gatorseq_fastq_to_analyze.py.log.txt 2>&1

/usr/bin/flock -n $CODE_PROD_FOLDER/LOCKS/QCI_upload.py.lock /home/path-svc-mol/Software/miniconda3/bin/python $CODE_PROD_FOLDER/QCI_upload.py &>> $CODE_PROD_FOLDER/LOGS/QCI_upload.py.log.txt

/usr/bin/flock -n $CODE_PROD_FOLDER/LOCKS/QCI_download.py.lock /home/path-svc-mol/Software/miniconda3/bin/python $CODE_PROD_FOLDER/QCI_download.py &>> $CODE_PROD_FOLDER/LOGS/QCI_download.py.log.txt

/usr/bin/flock -n $CODE_PROD_FOLDER/LOCKS/EPIC_upload.py.lock /home/path-svc-mol/Software/miniconda3/bin/python $CODE_PROD_FOLDER/EPIC_upload.py &>> $CODE_PROD_FOLDER/LOGS/EPIC_upload.py.log.txt 

/usr/bin/flock -n $CODE_PROD_FOLDER/LOCKS/gatorseq_db_to_excel.py.lock /home/path-svc-mol/Software/miniconda3/bin/python $CODE_PROD_FOLDER/gatorseq_db_to_excel.py &>> $CODE_PROD_FOLDER/LOGS/gatorseq_db_to_excel.py.log.txt

