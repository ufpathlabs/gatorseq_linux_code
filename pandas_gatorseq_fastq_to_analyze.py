import openpyxl as pyxl
import sys
import os
import datetime
import glob
import subprocess
import shlex
import yaml
import pandas as pd

print(str(datetime.datetime.now()) + "\n")

script_path = os.path.dirname(os.path.abspath( __file__ ))
CONFIG_FILE=script_path+"/linux_gatorseq.config.yaml"

#CONFIG_FILE=os.path.expanduser("~/gatorseq.config.1.1.yaml")


#FASTQ_FILES_DIR=HPC_FASTQ_FOLDER+""
CODE_ENV="DevEnv"

config_dict=dict()
with open(CONFIG_FILE, 'r') as stream:
    try:
        config_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()


# CODE_ENV=script_path.split('/')[-2]
USER_NAME=os.environ['USER']

def replace_env(strname):
    strname=strname.replace("USER_NAME",USER_NAME).replace("CODE_ENV",CODE_ENV)
    return strname
    

LINUX_PATHOLOGY_FASTQ_FOLDER = replace_env(config_dict['LINUX_PATHOLOGY_FASTQ_FOLDER'])
LINUX_ANALYSIS_OUT_FOLDER = replace_env(config_dict['LINUX_ANALYSIS_OUT_FOLDER'])
#ToDo: remove this
GATOR_SEQ_SAMPLE_INPUT_FILE = '/ext/path/DRL/Molecular/NGS/GatorSeq/DevEnv/3_GS_Fastq_to_Analyze_V1_0.xlsx'#replace_env(config_dict['GATOR_SEQ_SAMPLE_INPUT_FILE'])
LINUX_HPC_ANALYSIS_FOLDER = replace_env(config_dict['LINUX_HPC_ANALYSIS_FOLDER'])
HPC_ANALYSIS_FOLDER = replace_env(config_dict['HPC_ANALYSIS_FOLDER'])
LINUX_HPC_FASTQ_FOLDER = replace_env(config_dict['LINUX_HPC_FASTQ_FOLDER'])
HPC_FASTQ_FOLDER = replace_env(config_dict['HPC_FASTQ_FOLDER'])
GATOR_SEQ_HPC_CODE_DIR = replace_env(config_dict['GATOR_SEQ_HPC_CODE_DIR'])
HPC_NEXTFLOW_PROGRAM = replace_env(config_dict['HPC_NEXTFLOW_PROGRAM'])
HPC_SFTP = replace_env(config_dict['HPC_SFTP'])
LOCAL_PATHOLOGY_MNT = replace_env(config_dict['LOCAL_PATHOLOGY_MNT'])
LOCAL_HPC_MNT = replace_env(config_dict['LOCAL_HPC_MNT'])
GSBW_VERSION = replace_env(config_dict['GSBW_VERSION'])
NEXTFLOW_GIT_REPO = replace_env(config_dict['NEXTFLOW_GIT_REPO'])


#LINUX_HPC_ANALYSIS_FOLDER="/Users/path-svc-mol/Documents/mnt/HPC/GatorSeq/GatorSeq_V1_1/GatorSeq_Analysis"
#HPC_ANALYSIS_FOLDER="/ufrc/chamala/path-svc-mol/GatorSeq/GatorSeq_V1_1/GatorSeq_Analysis"

#LINUX_HPC_FASTQ_FOLDER="/Users/path-svc-mol/Documents/mnt/HPC/GatorSeq/GatorSeq_Fastq"
#HPC_FASTQ_FOLDER="/ufrc/chamala/path-svc-mol/GatorSeq/GatorSeq_Fastq"

#LINUX_ANALYSIS_OUT_FOLDER="/Users/path-svc-mol/Documents/mnt/PATHOLOGY/DRL/Molecular/NGS/GenomOncology/NextSeq"

#GATOR_SEQ_HPC_CODE_DIR="/ufrc/chamala/path-svc-mol/GatorSeq/GatorSeq_V1_1/gatorseq_hpc_code_v1_1"

#LINUX_PATHOLOGY_FASTQ_FOLDER="/Users/path-svc-mol/Documents/mnt/PATHOLOGY/DRL/Molecular/NGS/NextSeq_Fastq"
#GATOR_SEQ_SAMPLE_INPUT_FILE="/Users/path-svc-mol/Documents/mnt/PATHOLOGY/DRL/Molecular/NGS/NextSeq/GS_Fastq_to_Analyze_V1_1.xlsx"
#HPC_SNAKEMAKE_PROGRAM="/ufrc/chamala/path-svc-mol/GatorSeq/GatorSeq_V1_1/GatorSeq_Software/miniconda3/bin/snakemake"



def check_folders_exist():
    if not os.path.isfile(GATOR_SEQ_SAMPLE_INPUT_FILE):
        sys.exit("ERROR: Does not have access to following folder: " + GATOR_SEQ_SAMPLE_INPUT_FILE + "\n") 

    if not os.path.isdir(LINUX_HPC_ANALYSIS_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + LINUX_HPC_ANALYSIS_FOLDER + "\n")

    if not os.path.isdir(LINUX_HPC_FASTQ_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + LINUX_HPC_FASTQ_FOLDER + "\n")

    if not os.path.isdir(LINUX_ANALYSIS_OUT_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + LINUX_ANALYSIS_OUT_FOLDER + "\n")

    if not os.path.isdir(LINUX_PATHOLOGY_FASTQ_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + LINUX_PATHOLOGY_FASTQ_FOLDER + "\n") 


check_folders_exist()



def save_workbook(df):
    try:
        df.to_excel(GATOR_SEQ_SAMPLE_INPUT_FILE)
    except:
        print("could not save excel")
        sys.exit()

    
    
def check_success(log_file):
    datafile = open(log_file, 'rU')
    for line in datafile:
        if 'Ran Successfully' in line:
            return True
    return False

def check_failure(log_file):
    datafile = open(log_file, 'rU')
    for line in datafile:
        if 'ERROR' in line:
            return True
    return False

def populate_error_msg(row, error_message, xldf):
    row["STATUS"] = "FAILED"
    row["MESSAGE"] = error_message
    save_workbook(xldf)

# This checks if workbook is open; if so it will exit out

if __name__ == "__main__":
    
    try: 
        excel_file = open(GATOR_SEQ_SAMPLE_INPUT_FILE, "r+")
    except:
        print(" Could not open file! Please close Excel!")
        sys.exit()

    try:
            xldf_full = pd.read_excel(GATOR_SEQ_SAMPLE_INPUT_FILE)
            xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
            excel_file.close()
    except:
        print("Problem Reading Excel")
        sys.exit()
    


    for index, row in xldf.iterrows():
        
        #print(row)
        sample_path=row['SAMPLE_DIR_PATH']
        status=row["STATUS"]
        #print(sample_path)
        #print(status)
        #print(list(status))

        ## Check if the sample directory exist
        # Strips from both sites for following characaters ' ', '\t', '\r', '\/', '.', '*', '+'
        try:
            sample_path = sample_path.strip(' \t\n\r\/.\*\+')
        except Exception:
            populate_error_msg(row, "ERROR: Empty directory", xldf)
            continue

        if sample_path == '':
            row["STATUS"] = "FAILED"
            row["MESSAGE"] = "ERROR: Invalid directory"
            save_workbook(xldf)
            continue

        run_prefix = sample_path.split('/')[0]
        sample_prefix = sample_path.split('/')[1]
        
        sample_path = LINUX_PATHOLOGY_FASTQ_FOLDER + "/" + run_prefix + "-*/" +sample_prefix + "-*"
        listing = glob.glob(sample_path)
        directories=[d for d in listing if os.path.isdir(d)]
        run_name = ''
        sample_name = ''

        if len(directories) == 0:
            #ToDO: replace error messages with definition
            row["STATUS"] = "FAILED"
            row["MESSAGE"] = "ERROR: Following directory extention not found:" + sample_path
            save_workbook(xldf)
            continue
        elif len(directories) > 1:
            row["STATUS"] = "FAILED"
            row["MESSAGE"] = "ERROR: More than one directory extention found:" + ';'.join(directories)
            save_workbook(xldf)
            continue
        else:
            sample_path = directories[0]
            #run_name = sample_path.split('/')[-2]
            #sample_name = sample_path.split('/')[-1]

        linux_hpc_sample_path =  LINUX_HPC_FASTQ_FOLDER + "/" + run_prefix + "-*/" +sample_prefix + "-*"
        hpc_sample_path = ""
        listing = glob.glob(linux_hpc_sample_path)
        directories=[d for d in listing if os.path.isdir(d)]
        if len(directories) == 0:
            row["STATUS"] = "RUN"
            row["MESSAGE"] = "WARNING: Directory is not yet uploaded to HPC:" + linux_hpc_sample_path
            save_workbook(xldf)
            continue
        elif len(directories) > 1:
            row["STATUS"] = "FAILED"
            row["MESSAGE"] = "ERROR: More than one directory extention found on HPC:" + ';'.join(directories)
            save_workbook(xldf)
            continue
        else:
            linux_hpc_sample_path = directories[0]
            sample_name = sample_path.split('/')[-1]
            run_name = sample_path.split('/')[-2]
            hpc_sample_path =  HPC_FASTQ_FOLDER + "/" + run_name + "/" +sample_name 




        ## Check if the sample directory exist
        if status == 'RUN':
            print("running from pandas")
            ### Check sample folders are synced properly 
            rsync_check_cmd =  'rsync -ltgoDvzrn --progress --stats ' +\
                sample_path + "/ " +\
                linux_hpc_sample_path + '/ |grep "Total transferred file size: 0 bytes"'
            proc = subprocess.Popen(rsync_check_cmd, shell=True, stdout=subprocess.PIPE)
            #need to decode from bytes to str
            stdout_list = proc.communicate()[0].decode()
            if len(stdout_list) == 0:
                row["STATUS"] = "RUN"
                row["MESSAGE"] = "WARNING: Directory is not completely in sync with HPC folder:" + linux_hpc_sample_path
                save_workbook(xldf)

                rsync_hpc_check_cmd =  'rsync -ltgoDvzrn --progress --stats ' +\
                    sample_path + "/ " + HPC_SFTP + ":" +\
                    hpc_sample_path + '/ |grep "Total transferred file size: 0 bytes"'
                proc = subprocess.Popen(rsync_hpc_check_cmd, shell=True, stdout=subprocess.PIPE)
                #need to decode from bytes to str
                stdout_list = proc.communicate()[0].decode()
                if len(stdout_list) != 0:
                    umount_cmd = '/sbin/umount '+ LOCAL_HPC_MNT + ';'
                    proc = subprocess.Popen(umount_cmd, shell=True, stdout=subprocess.PIPE)
                    sys.exit("ERROR: mount snyc problem; executed "+ umount_cmd + "\n")

                continue

            time_stamp = str(datetime.datetime.now()).replace('-', '').replace(' ','').replace(':', '').replace('.', '')
            time_stamp= time_stamp +CODE_ENV+GSBW_VERSION

            run_log = sample_path + "/" + time_stamp + ".log"

            row['STATUS'] = "SUBMITTED"
            row['TIME_STAMP'] = time_stamp 
            row['MESSAGE'] = "Currently job is running."
            save_workbook(xldf)

            HPC_JOB_PREFIX= sample_prefix+ "_"+ time_stamp
            HPC_RUN_DIR=HPC_ANALYSIS_FOLDER+"/"+ HPC_JOB_PREFIX
            LINUX_HPC_RUN_DIR_CRONJOB=LINUX_HPC_ANALYSIS_FOLDER+"/"+ HPC_JOB_PREFIX + ".testJob.sh"
            LINUX_HPC_RUN_DIR_CRONJOB_LOG = LINUX_HPC_ANALYSIS_FOLDER + "/" + HPC_JOB_PREFIX + ".testJob.sh.log"
            #FASTQ_FILES_DIR=HPC_FASTQ_FOLDER+'/'+run_prefix+'*/'+sample_prefix+'*'
            FASTQ_ROOT_DIR_SUFFIX=run_prefix+'*/'+sample_prefix+'*'
            #print(FASTQ_FILES_DIR)

            HPC_RUN_CMD =   "#!/bin/bash\n" 

            HPC_RUN_CMD = HPC_RUN_CMD + "SAMPLE_RUN_DIR=" + HPC_RUN_DIR +  ";\n"
            HPC_RUN_CMD = HPC_RUN_CMD + "rm -r -f $SAMPLE_RUN_DIR" +  ";\n"
            HPC_RUN_CMD = HPC_RUN_CMD + "mkdir $SAMPLE_RUN_DIR" + ";\n"
            HPC_RUN_CMD = HPC_RUN_CMD + "chmod -R 750 $SAMPLE_RUN_DIR" + ";\n"
            HPC_RUN_CMD = HPC_RUN_CMD + "cd $SAMPLE_RUN_DIR" ";\n"
            HPC_RUN_CMD = HPC_RUN_CMD + "\n\n"


            HPC_RUN_CMD = HPC_RUN_CMD + "export NXF_TEMP=$SAMPLE_RUN_DIR"+"\n"
            HPC_RUN_CMD = HPC_RUN_CMD + "export NXF_ASSETS=$SAMPLE_RUN_DIR"+"\n"
            HPC_RUN_CMD = HPC_RUN_CMD + HPC_NEXTFLOW_PROGRAM + " run \\\n" +\
            GATOR_SEQ_HPC_CODE_DIR +" -r " + GSBW_VERSION + " \\\n"+\
            "  -with-report -with-trace -with-timeline -with-dag flowchart.html \\\n" +\
            "  -resume \\\n" +\
            " --FASTQ_ROOT_DIR_SUFFIX=" + FASTQ_ROOT_DIR_SUFFIX + " \\\n" +\
            " --SAMPLE_DIR=$SAMPLE_RUN_DIR" + " \\\n" +\
            " --RUN_NAME=" + run_prefix + " \\\n" +\
            " --SAMPLE_NAME=" + sample_prefix + " \\\n" +\
            " --TIME_STAMP=" + time_stamp + " \\\n" +\
            " --CODE_ENV=" + CODE_ENV + " ;\n"



            HPC_RUN_CMD = HPC_RUN_CMD + "rm -r -f " + HPC_RUN_DIR + "/.nextflow* ;\n"
            HPC_RUN_CMD = HPC_RUN_CMD + "rm -r -f " + HPC_RUN_DIR + "/work ;\n"
            HPC_RUN_CMD = HPC_RUN_CMD + "rm -r -f " + HPC_RUN_DIR + "/" + NEXTFLOW_GIT_REPO + ";\n"

            print(HPC_RUN_CMD)
            cron_fw = open(LINUX_HPC_RUN_DIR_CRONJOB,'w')
            cron_fw.write(HPC_RUN_CMD)
            cron_fw.close()
            #print(HPC_RUN_CMD)
            print("Submitted: "+ LINUX_HPC_RUN_DIR_CRONJOB + "\n")

        elif status == 'SUBMITTED':
            time_stamp = row["TIME_STAMP"]
            linux_analysis_out_sample_dir = LINUX_ANALYSIS_OUT_FOLDER + "/" + run_prefix + "/" + sample_prefix + "_" + time_stamp
            hpc_analysis_out_sample_dir = LINUX_HPC_ANALYSIS_FOLDER + "/" + run_prefix + "/" + sample_prefix + "_" + time_stamp

            #run_log = sample_path + "/" + time_stamp + ".log"
            if os.path.isdir(linux_analysis_out_sample_dir):
                if os.path.isdir(hpc_analysis_out_sample_dir):
                    pass
                elif os.path.isfile(linux_analysis_out_sample_dir + "/SUCCESS.txt"):
                    row["STATUS"] = "DONE"
                    row["MESSAGE"] = "Successfully processed"
                    save_workbook(xldf)
                elif os.path.isfile(linux_analysis_out_sample_dir + "/FAILED.txt"):
                    row["STATUS"] = "FAILED"
                    row["MESSAGE"] = "ERROR: Failed while running on HPC"
                    save_workbook(xldf)

        elif status == 'FAILED':
            pass

        elif status == 'DONE':
            pass

        else:
            row["MESSAGE"] = "ERROR: Invalid Status." 
            save_workbook(xldf)
        

