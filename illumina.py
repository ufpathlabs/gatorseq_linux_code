import openpyxl as pyxl
import sys
import os
import datetime
import glob
import subprocess
import numpy as np
import pandas as pd
from filelock import FileLock
import traceback
import sqlite3
import yaml
import mysql.connector
import json

PROJECT_ID_MAP = {
    "CLIN_WGS": "152973821"
}
APPLICATION_ID_MAP = {
    "CLIN_WGS": "9650641"
}

print(str(datetime.datetime.now()) + "\n")

script_path = os.path.dirname(os.path.abspath( __file__ ))
CONFIG_FILE=script_path+"/linux_gatorseq.config.yaml"

config_dict=dict()
with open(CONFIG_FILE, 'r') as stream:
    try:
        config_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()


CODE_ENV=script_path.split('/')[-2]
USER_NAME=os.environ['USER']

def replace_env(strname):
    strname=strname.replace("USER_NAME",USER_NAME).replace("CODE_ENV",CODE_ENV)
    return strname

ILLUMINA_SAMPLE_FILE = replace_env(config_dict['ILLUMINA_SAMPLE_FILE'])
ILLUMINA_TABLE_NAME = replace_env(config_dict['ILLUMINA_TABLE_NAME'])

MYSQL_HOST = config_dict['MYSQL_HOST']
MYSQL_USERNAME = config_dict['MYSQL_USERNAME']
# MYSQL_PASSWAORD = config_dict['MYSQL_PASSWAORD']
MYSQL_DATABASE = config_dict['MYSQL_DATABASE']

CONFIG_TOKENS_FILE = script_path + "/" + config_dict['CONFIG_TOKENS_FILE']
config_token_dict=dict()
with open(CONFIG_TOKENS_FILE, 'r') as stream:
    try:
        config_token_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

MYSQL_PASSWAORD = config_token_dict['MYSQL_PASSWAORD']

if CODE_ENV == "ProdEnv":
    MYSQL_HOST = config_dict['PROD_MYSQL_HOST']
    MYSQL_USERNAME = config_dict['PROD_MYSQL_USERNAME']
    MYSQL_PASSWAORD = config_token_dict['PROD_MYSQL_PASSWAORD']
    MYSQL_DATABASE = config_dict['PROD_MYSQL_DATABASE']


file_to_lock = ILLUMINA_SAMPLE_FILE + '.lock'
lock = FileLock(file_to_lock)
try:
    lock.acquire(timeout=1)
except:
    print("some other script is using it")
    sys.exit()

try: 
    excel_file = open(ILLUMINA_SAMPLE_FILE, "r+")
except:
    print(" Could not open file! Please close Excel!")
    sys.exit()


def create_connection():
    conn = None
    try:
        conn = mysql.connector.connect(
            host=MYSQL_HOST,
            user=MYSQL_USERNAME,
            passwd=MYSQL_PASSWAORD,
            database=MYSQL_DATABASE
        )
    except:
        print(traceback.format_exc())
 
    return conn


def read_excel_and_upsert(conn):
    xldf_full = pd.read_excel(ILLUMINA_SAMPLE_FILE)
    xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    xldf = xldf.replace(np.nan, '', regex=True)
    for index, row in xldf.iterrows():
        # print(row['SAMPLE_DIR_PATH'], row['STATUS'], row['TIME_STAMP'], row['MESSAGE'], row['PLMO_Number'], row['Test_Product_Profile'], row['Test_Product_Code'], row['Diagnosis'], row['Primary_Tumor_Site'], row['Pre_Filter'], row['Report_Template'], row['QCIType'], row['Treatments_Policy'], row['Reporting_Method'])
        # update_task(conn, (row['SAMPLE_DIR_PATH'], row['STATUS'], row['PLMO_Number'], row['Test_Product_Profile'], row['Test_Product_Code'], row['Diagnosis'], row['Primary_Tumor_Site'], row['Pre_Filter'], row['Report_Template'], row['QCIType'], row['Treatments_Policy'], row['Reporting_Method']))
        cur = conn.cursor()
        cur.execute("SELECT * FROM "+ILLUMINA_TABLE_NAME+" where SAMPLE_NAME = '" + row['SAMPLE_NAME'] + "'")
        rows = cur.fetchall()
        if len(rows) > 0:
            sql = '''
            UPDATE '''+ILLUMINA_TABLE_NAME+'''
            SET SAMPLE_NAME = %s,
                PROJECT_NAME = %s,
                GENDER = %s
                
            WHERE SAMPLE_NAME = %s;''' 
            #print(sql)
            cur2 = conn.cursor()
            cur2.execute(sql, (row['SAMPLE_NAME'], row['PROJECT_NAME'], row['GENDER'], row['SAMPLE_NAME'] ))
            conn.commit()
            cur2.close()
        else:
            sql = ''' INSERT into '''+ILLUMINA_TABLE_NAME+'''(SAMPLE_NAME, PROJECT_NAME, GENDER) values(%s,%s,%s); '''
            
            cur2 = conn.cursor()
            #print(sql)
            cur2.execute(sql, (row['SAMPLE_NAME'], row['PROJECT_NAME'], row['GENDER'] ))
            conn.commit()
            cur2.close()
    
    
def getJSONFromBashCommand(cmd):
    process = subprocess.Popen(cmd.split(), 
                           stdout=subprocess.PIPE)
                        #    universal_newlines=True)
    jsonBin, errors = process.communicate()
   
    if not errors:
        return json.loads(jsonBin)
    else:
        print("error while executing the command-----------> ", cmd)
        return None

def findStatus(appSessionId):
    statusJson = getJSONFromBashCommand("/home/path-svc-mol/Illumina_Binary/bin/bs --config UFMOL_ENTERPRISE appsession get --id=219690534 --format=json")
    print(statusJson["Id"])

def updateRowWithStatus(sampleName, status, appSessionLabel, conn):
    cur = conn.cursor()
    updateSql = "update "+ ILLUMINA_TABLE_NAME +" set APP_SESSION_STATUS = %s, APP_SESSION_ID = %s where SAMPLE_NAME = %s;"
    cur.execute(updateSql, (status, appSessionLabel, sampleName))
    conn.commit()
    cur.close()

def submitJob(applicationId, projectId, appSessionName, appSessionLabel):
    gender = "auto"
    # if genderProvided:
    #     gender = genderProvided
    bashCommand = """/home/path-svc-mol/Illumina_Binary/bin/bs \
        --config UFMOL_ENTERPRISE 
        application launch 
        --id=""" + applicationId + """ 
        --option=project-id:""" + projectId + """ \
        --option=ht-ref:hg19-altaware-cnv-anchor.v8 \
        --option=cnv_checkbox:1 \
        --option=cnv_ref:1 \
        --option=cnv_segmentation_mode:slm \
        --option=sv_checkbox:1 \
        --option=eh_checkbox:1 \
        --option=output_format:BAM \
        --option=vcf_or_gvcf:GVCF \
        --option=dupmark_checkbox:1 \
        --option=vc_enable_bqd_checkbox:1 \
        --option=metrics_checkbox:1 \
        --option=md5_all:1 \
        --option=automation_checkbox:1 \
        --option=automation-sample-id:216890698 \
        --option=automation-sex:""" + gender + """ \
        --option=app-session-name: """ + appSessionName + """ \
        --appsession-label=""" + appSessionLabel + """ \
        --format=json 
    """
    #print(bashCommand)
    submitJobResponseJson = getJSONFromBashCommand(bashCommand)
    return submitJobResponseJson
         



#read the database for any jobs with APP_SESSION_STATUS == "" 
def getJobsToSubmit(conn):
    cur = conn.cursor()
    cur.execute("SELECT * FROM "+ILLUMINA_TABLE_NAME+" where APP_SESSION_STATUS = '';")
    rows = cur.fetchall()
    cur.close()
    if len(rows):
        for row in rows:
            sampleId = None
            sampleName = row[0]
            getSampleIdResponseJson = getJSONFromBashCommand("/home/path-svc-mol/Illumina_Binary/bin/bs --config UFMOL_ENTERPRISE  biosample get --name=" + sampleName + " --format=json")
            if getSampleIdResponseJson:
                sampleId = getSampleIdResponseJson["Id"]
                #submitting the job
                projectId = PROJECT_ID_MAP.get(row[1]) 
                applicationId = APPLICATION_ID_MAP.get(row[1])
                if projectId and applicationId and sampleId:
                    time_stamp = str(datetime.datetime.now()).replace('-', '').replace(' ','').replace(':', '').replace('.', '')
                    time_stamp = time_stamp + CODE_ENV
                    appSessionName = str(row[0]) + "_AppSessionName_" + time_stamp
                    appSessionLabel = str(row[0]) + "_AppSessionLabel_" + time_stamp

                    jobSubmitted = submitJob(applicationId, projectId, appSessionName, appSessionLabel)
                    if jobSubmitted:
                        print(jobSubmitted)
                        updateRowWithStatus(row[0], "Submitted", appSessionLabel, connection)
                    else:
                        print("error while submitting the job")
                    
            else:
                print("unable to get sample Id for sample name: " + row[0])
                continue
                           
    


if __name__ == "__main__":
    connection = create_connection()
    #read_excel_and_upsert(connection)

    toSubmitJobs = getJobsToSubmit(connection)
    # findStatus(1234)
    connection.close()



