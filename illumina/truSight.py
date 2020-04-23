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
import requests



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

TRUSIGHT_APP_SUBMISSION_INPUT_FILE = replace_env(config_dict['TRUSIGHT_APP_SUBMISSION_INPUT_FILE'])
TRUSIGHT_TABLE_NAME = replace_env(config_dict['TRUSIGHT_TABLE_NAME'])

TRUSIGHT_NEW_CASE_URL = "https://ufl-tss.trusight.illumina.com/crs/api/v1/cases"

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
TRUSIGHT_API_KEY = config_token_dict['TRUSIGHT_API_KEY']

if CODE_ENV == "ProdEnv":
    MYSQL_HOST = config_dict['PROD_MYSQL_HOST']
    MYSQL_USERNAME = config_dict['PROD_MYSQL_USERNAME']
    MYSQL_PASSWAORD = config_token_dict['PROD_MYSQL_PASSWORD']
    MYSQL_DATABASE = config_dict['PROD_MYSQL_DATABASE']
    TRUSIGHT_TABLE_NAME = replace_env(config_dict['TRUSIGHT_TABLE_NAME_PROD'])


file_to_lock = TRUSIGHT_APP_SUBMISSION_INPUT_FILE + '.lock'
lock = FileLock(file_to_lock)
try:
    lock.acquire(timeout=1)
except:
    print("some other script is using it")
    sys.exit()

try: 
    excel_file = open(TRUSIGHT_APP_SUBMISSION_INPUT_FILE, "r+")
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
    xldf_full = pd.read_excel(TRUSIGHT_APP_SUBMISSION_INPUT_FILE)
    xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    xldf = xldf.replace(np.nan, '', regex=True)
    for index, row in xldf.iterrows():
        # print(row['SAMPLE_DIR_PATH'], row['STATUS'], row['TIME_STAMP'], row['MESSAGE'], row['PLMO_Number'], row['Test_Product_Profile'], row['Test_Product_Code'], row['Diagnosis'], row['Primary_Tumor_Site'], row['Pre_Filter'], row['Report_Template'], row['QCIType'], row['Treatments_Policy'], row['Reporting_Method'])
        # update_task(conn, (row['SAMPLE_DIR_PATH'], row['STATUS'], row['PLMO_Number'], row['Test_Product_Profile'], row['Test_Product_Code'], row['Diagnosis'], row['Primary_Tumor_Site'], row['Pre_Filter'], row['Report_Template'], row['QCIType'], row['Treatments_Policy'], row['Reporting_Method']))
        cur = conn.cursor()
        cur.execute("SELECT * FROM "+TRUSIGHT_TABLE_NAME+" where SAMPLE_NAME = '" + row['SAMPLE_NAME'] + "'")
        rows = cur.fetchall()
        cur.close()
        if len(rows) > 0:
            sql = '''
            UPDATE '''+TRUSIGHT_TABLE_NAME+'''
            SET SAMPLE_NAME = %s,
                GENDER = %s
                
            WHERE SAMPLE_NAME = %s;''' 
            #print(sql)
            cur2 = conn.cursor()
            cur2.execute(sql, (row['SAMPLE_NAME'], row['GENDER'], row['SAMPLE_NAME'] ))
            conn.commit()
            cur2.close()
        else:
            sql = ''' INSERT into '''+TRUSIGHT_TABLE_NAME+'''(SAMPLE_NAME, GENDER, STATUS) values(%s,%s, ""); '''
            
            cur2 = conn.cursor()
            #print(sql)
            cur2.execute(sql, (row['SAMPLE_NAME'], row['GENDER'] ))
            conn.commit()
            cur2.close()
    return xldf
    
    
def runBashCommand(cmd):
    process = subprocess.Popen(cmd.split(), 
                           stdout=subprocess.PIPE)
                        #    universal_newlines=True)
    responseBin, errors = process.communicate()
    #print(responseBin)
    #print(errors) 
    if not errors and len(responseBin):
        return responseBin
    else:
        print("error while executing the command-----------> ", cmd)
        return None

def updateRowWithStatus(sampleName, status, conn):
    cur = conn.cursor()
    updateSql = "update "+ TRUSIGHT_TABLE_NAME +" set STATUS = %s where SAMPLE_NAME = %s;"
    cur.execute(updateSql, (status, sampleName))
    conn.commit()
    cur.close()


            
def populateStatusInExcel(conn, df):
    for index, row in df.iterrows():
        cur = conn.cursor()
        cur.execute("SELECT * FROM "+TRUSIGHT_TABLE_NAME+" where SAMPLE_NAME = '" + row['SAMPLE_NAME'] + "'")
        rows = cur.fetchall()
        if len(rows) > 0:
            row = rows[0]
            df.at[index, "STATUS"] = row[2]
    df.to_excel(TRUSIGHT_APP_SUBMISSION_INPUT_FILE, index=False)

# create a sample using REST API of trusight 
# use upload command of trusight to upload the sample
def createSample(conn, baseMountDir):
    cur = conn.cursor()
    cur.execute("SELECT * FROM "+TRUSIGHT_TABLE_NAME+" where STATUS = 'FASTQ_UPLOADED';")
    rows = cur.fetchall()
    cur.close()
    for row in rows:
        sampleName = row[0]
        gender = row[1]
        headers = {
            "x-auth-token": "ApiKey " + TRUSIGHT_API_KEY,
            "X-ILMN-Domain": "ufl-tss",
            "X-ILMN-Workgroup": "51b925ca-56ed-37c9-89d7-85b83a1f7e55",
            "Content-Type": "application/json",
            "Accept": "application/json"
        }
        fastq_list = os.listdir(baseMountDir + "/Projects/WGS/Samples/"+ sampleName +"/Files")
        print("fastq files list : ", fastq_list)
        fastqs = []
        for f in fastq_list:
            fastqs.append({
                "fastqLink": sampleName + "/" + f,
                "userProvidedDirectory": sampleName
            })

        data = {
            "displayId": ("auto_" + sampleName)[:12],
            "testDefinitionId": "1cb2c841-9a25-4cb3-8b55-c3ea0d639086",
            "subjects": [{
                "isAffected": "AFFECTED",
                "relationshipToProband": "PROBAND",
                "previousTestHistory": "",
                "reportTypes": ["10443391-216f-4ccc-9513-35d2777fb17f"],
                "gender": gender,
                "samples": [
                    {
                        "externalSampleId": sampleName,
                        "molecularData": fastqs
                    }]
            }]
        }
        response = requests.post(TRUSIGHT_NEW_CASE_URL, headers=headers, data=json.dumps(data))
        
  #ToDO: add success logic
        if response.status_code in [200, 201]:
            caseId = response.json()["id"]
            updateRowWithStatus(sampleName, "SUBMITTED", conn)
            processCase(caseId, sampleName, connection)
        else:
            updateRowWithStatus(sampleName, "ERROR_WHILE CREATING_CASE", conn)
            print("error while creating case:", response.json())

def processCase(caseId, sampleName, conn):
    print("going to process casewith (caseId, sampleName):", caseId, sampleName)
    headers = {
            "X-Auth-Token": "ApiKey " + TRUSIGHT_API_KEY,
            "X-ILMN-Domain": "ufl-tss",
            "X-ILMN-Workgroup": "51b925ca-56ed-37c9-89d7-85b83a1f7e55",
            "Content-Type": "application/json"
        }
    data = {
        
    }
    response = requests.post(TRUSIGHT_NEW_CASE_URL+"/"+caseId+"/process", headers=headers, data=data)
    if response.status_code in [200,201]:
        updateRowWithStatus(sampleName, "DONE", conn)
    else:
        print(response.json())
        updateRowWithStatus(sampleName, "ERROR_WHILE_PROCESSING_CASE", conn)

# mounts the basespace folder in the 'basDir' folder
def mountBaseSpace(basDir):
    bashResponse = runBashCommand("basemount --unmount " + script_path + "/" + basDir)
    if bashResponse:
        bashCommand = """basemount \
            --config UFMOL_ENTERPRISE """ + basDir
        bashResponse = runBashCommand(bashCommand)
        return True if bashResponse else False
    else:
        return False

# gets all sample names from DB with no status
# checks in basemount dir if fastq files exists
# returns list of samplenames which have status as "RUN" and FastQ files in basemount
def checkFastqExists(conn, baseMountDir):
    cur = conn.cursor()
    cur.execute("SELECT * FROM "+TRUSIGHT_TABLE_NAME+" where STATUS = '';")
    rows = cur.fetchall()
    cur.close()
    for row in rows:
        sampleName = row[0]
        if os.path.isdir(baseMountDir + "/Projects/WGS/Samples/"+sampleName+"/Files"):
            updateRowWithStatus(sampleName, "UPLOAD_FASTQ_PENDING", conn)
        else:
            updateRowWithStatus(sampleName, "FASTQ_NOT_FOUND", conn)
        

if __name__ == "__main__":
    connection = create_connection()

    df = read_excel_and_upsert(connection)

    baseMountDir = "BaseMount" #+ datetime.datetime.now.strftime("%d/%m/%Y %H:%M:%S")

    if mountBaseSpace(baseMountDir):

        checkFastqExists(connection, baseMountDir)

        createSample(connection, baseMountDir)

        

        populateStatusInExcel(connection, df)
    #processCase("40e0e508-0052-4c59-bd5e-6d6749337f8e", "NQ-20-02_BC712507_Z8-5x", connection)
    connection.close()

