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
parent_path = os.path.abspath(os.path.join(script_path, '..'))

CONFIG_FILE = parent_path + "/linux_gatorseq.config.yaml"

config_dict=dict()
with open(CONFIG_FILE, 'r') as stream:
    try:
        config_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()


CODE_ENV=script_path.split('/')[-3]
USER_NAME=os.environ['USER']

def replace_env(strname):
    strname=strname.replace("USER_NAME",USER_NAME).replace("CODE_ENV",CODE_ENV)
    return strname

ILLUMINA_SAMPLE_FILE = replace_env(config_dict['ILLUMINA_BASESPACE_APP_SUBMISSION_INPUT_FILE'])
ILLUMINA_TABLE_NAME = replace_env(config_dict['ILLUMINA_TABLE_NAME'])

MYSQL_HOST = config_dict['DEV_MYSQL_HOST']
MYSQL_USERNAME = config_dict['DEV_MYSQL_USERNAME']
# MYSQL_PASSWAORD = config_dict['MYSQL_PASSWAORD']
MYSQL_DATABASE = config_dict['DEV_MYSQL_DATABASE']

CONFIG_TOKENS_FILE = parent_path + "/" + config_dict['CONFIG_TOKENS_FILE']
config_token_dict=dict()
with open(CONFIG_TOKENS_FILE, 'r') as stream:
    try:
        config_token_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

MYSQL_PASSWAORD = config_token_dict['DEV_MYSQL_PASSWORD']

if CODE_ENV == "ProdEnv":
    MYSQL_HOST = config_dict['PROD_MYSQL_HOST']
    MYSQL_USERNAME = config_dict['PROD_MYSQL_USERNAME']
    MYSQL_PASSWAORD = config_token_dict['PROD_MYSQL_PASSWORD']
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

# read the excel for all rows and add any new samples to database or update the rows in database.
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
        cur.close()
        if len(rows) > 0:
            sql = '''
            UPDATE '''+ILLUMINA_TABLE_NAME+'''
            SET SAMPLE_NAME = %s,
                PROJECT_NAME = %s
                
            WHERE SAMPLE_NAME = %s;''' 
            #print(sql)
            cur2 = conn.cursor()
            cur2.execute(sql, (row['SAMPLE_NAME'], row['PROJECT_NAME'], row['SAMPLE_NAME'] ))
            conn.commit()
            cur2.close()
        else:
            sql = ''' INSERT into '''+ILLUMINA_TABLE_NAME+'''(SAMPLE_NAME, PROJECT_NAME) values(%s,%s); '''
            
            cur2 = conn.cursor()
            #print(sql)
            cur2.execute(sql, (row['SAMPLE_NAME'], row['PROJECT_NAME'] ))
            conn.commit()
            cur2.close()
    return xldf
    
    
def getJSONFromBashCommand(cmd):
    process = subprocess.Popen(cmd.split(), 
                           stdout=subprocess.PIPE)
                        #    universal_newlines=True)
    jsonBin, errors = process.communicate()
    print(jsonBin)
    print(errors) 
    if not errors and len(jsonBin):
        return json.loads(jsonBin)
    else:
        print("error while executing the command-----------> ", cmd)
        return None

def findStatus(appSessionId):
    bashCommand = """/home/path-svc-mol/Illumina_Binary/bin/bs \
        --config UFMOL_ENTERPRISE \
        appsession get \
        --name=""" + appSessionId + """ \
        --format=json"""
    statusJson = getJSONFromBashCommand(bashCommand)
    if statusJson:
        return statusJson["ExecutionStatus"]
    else:
        return None

def updateRowWithStatus(sampleName, status, appSessionLabel, conn):
    cur = conn.cursor()
    updateSql = "update "+ ILLUMINA_TABLE_NAME +" set APP_SESSION_STATUS = %s, APP_SESSION_ID = %s where SAMPLE_NAME = %s;"
    cur.execute(updateSql, (status, appSessionLabel, sampleName))
    conn.commit()
    cur.close()

def updateRowWithError(sampleName, status, conn):
    cur = conn.cursor()
    updateSql = "update "+ ILLUMINA_TABLE_NAME +" set APP_SESSION_STATUS = %s where SAMPLE_NAME = %s;"
    cur.execute(updateSql, (status, sampleName))
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
         



#read the database for any jobs with APP_SESSION_STATUS == "" and tries to submit them. 
def submitNewJobs(conn):
    cur = conn.cursor()
    cur.execute("SELECT * FROM "+ILLUMINA_TABLE_NAME+" where APP_SESSION_STATUS is null;")
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
                        updateRowWithStatus(row[0], "In Progress", appSessionLabel, connection)
                    else:
                        print("-----------------error while submitting the job---------------------")
                        updateRowWithError(row[0], "Error while submitting, please check with tech team", connection)
                    
            else:
                print("unable to get sample Id for sample name: " + row[0])
                           
    
def checkAndUpdateStatus(conn):
    cur = conn.cursor()
    cur.execute("SELECT * FROM "+ILLUMINA_TABLE_NAME+" where APP_SESSION_STATUS = 'In Progress';")
    rows = cur.fetchall()
    cur.close()
    if len(rows):
        for row in rows:
            status = findStatus(row[6])
            #ToDo: check for other statuses
            if status and status == "Complete":
                updateRowWithStatus(row[0], status, row[6], conn)
            
def populateStatusInExcel(conn, df):
    for index, row in df.iterrows():
        cur = conn.cursor()
        cur.execute("SELECT * FROM "+ILLUMINA_TABLE_NAME+" where SAMPLE_NAME = '" + row['SAMPLE_NAME'] + "'")
        rows = cur.fetchall()
        if len(rows) > 0:
            row = rows[0]
            df.at[index, "APP_SESSION_STATUS"] = row[7]
            df.at[index, "APP_SESSION_LABEL_ID"] = row[6]
    df.to_excel(ILLUMINA_SAMPLE_FILE, index=False)

if __name__ == "__main__":
    connection = create_connection()

    df = read_excel_and_upsert(connection)

    submitNewJobs(connection)

    checkAndUpdateStatus(connection)

    populateStatusInExcel(connection, df)

    connection.close()



