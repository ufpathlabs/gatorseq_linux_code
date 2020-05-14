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
from truSight import runBashCommand
from truSight import mountBaseSpace
from truSight import updateRowWithStatus
from truSight import populateStatusInExcel
from truSight import updateRowWithStatusAndMessage
import time
import multiprocessing as mp




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
TRUSIGHT_CLI = replace_env(config_dict['TRUSIGHT_CLI'])

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
TRUSIGHT_API_KEY = config_token_dict['DEV_MYSQL_PASSWORD']

if CODE_ENV == "ProdEnv":
    MYSQL_HOST = config_dict['PROD_MYSQL_HOST']
    MYSQL_USERNAME = config_dict['PROD_MYSQL_USERNAME']
    MYSQL_PASSWAORD = config_token_dict['PROD_MYSQL_PASSWORD']
    MYSQL_DATABASE = config_dict['PROD_MYSQL_DATABASE']
    TRUSIGHT_TABLE_NAME = replace_env(config_dict['TRUSIGHT_TABLE_NAME_PROD'])



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

def processCode(sampleName, directoryName, filep, index, conn):
    startTime = time.time()
    cmd = "java -jar "+ TRUSIGHT_CLI + " stage --stageDirectory=" + directoryName + "/ --localDirectory=" +  filep + "/" + sampleName + "/"
    print("running the following command: ", cmd)
    #updateRowWithStatusAndMessage(sampleName, "STARTED_UPLOAD", "uopload started at " + time.ctime(), directoryName, conn)
    print('update DB that I am starting:', sampleName)
    statusJson = runBashCommand(cmd, index)
    endTime = time.time()
    #updateRowWithStatusAndMessage(sampleName, "STARTED_UPLOAD", "uopload started at " + time.ctime(), conn)
    timeForExecution = round((endTime - startTime)/60)
    print("---------------------------------", sampleName) 
    print("-----status received from bash command", statusJson)
    new_conn = create_connection()
    if statusJson:
        # updateRowWithStatus(sampleName, "FASTQ_UPLOADED", conn)
        updateRowWithStatusAndMessage(sampleName, "FASTQ_UPLOADED", "file uploaded in: " + str(timeForExecution) + " minutes", directoryName, new_conn)
    else:
        print("status is None")
        updateRowWithStatusAndMessage(sampleName, "ERROR_UPLOADING", "Error uploading the file to trusight", directoryName, new_conn)
    new_conn.close()
    return sampleName
    

def uploadFastQ(conn):
    cur = conn.cursor()
    cur.execute("SELECT * FROM "+TRUSIGHT_TABLE_NAME+" where  STATUS = 'UPLOAD_FASTQ_PENDING' and FILE_PATH is not null;")
    rows = cur.fetchall()
    cur.close()

    while len(rows) > 0:
        print("uploading fastq first 3")
        cur = rows[:3]
        pool = mp.Pool(processes=len(cur))

        results = [pool.apply_async(processCode, args=(row[0], row[4], row[5], i+1, conn,)) for i, row in enumerate(cur)]
        
        output = [p.get() for p in results]
        print("output for current set of cur:", output)
        populateStatusInExcel(connection, df)
        
        rows = rows[len(cur):]
        
# gets all sample names from DB with no status
# checks in basemount dir if fastq files exists
# returns list of samplenames which have status as "RUN" and FastQ files in basemount
def checkFastqExists(conn):
    cur = conn.cursor()
    cur.execute("SELECT * FROM "+TRUSIGHT_TABLE_NAME+" where STATUS = 'RUN' and FILE_PATH is not null;")
    rows = cur.fetchall()
    cur.close()
    for row in rows:
        sampleName = row[0]
        filepath = row[5]
        if os.path.isdir(filepath + "/" + sampleName + "/"):
            updateRowWithStatus(sampleName, "UPLOAD_FASTQ_PENDING", row[4], conn)
        else:
            updateRowWithStatus(sampleName, "FASTQ_NOT_FOUND", row[4], conn)        

def read_excel_and_upsert(conn):
    xldf_full = pd.read_excel(TRUSIGHT_APP_SUBMISSION_INPUT_FILE)
    xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    xldf = xldf.replace(np.nan, '', regex=True)
    for index, row in xldf.iterrows():
        print("row === ")
        print(row)
        cur = conn.cursor()
        selectSql = "SELECT * FROM "+TRUSIGHT_TABLE_NAME+" where SAMPLE_NAME = %s and DIRECTORY_NAME = %s;"

        cur.execute(selectSql, (row['SAMPLE_NAME'], row['DIRECTORY_NAME'] ))
        rows = cur.fetchall()
        cur.close()
        if len(rows) > 0:
            sql = '''
            UPDATE '''+TRUSIGHT_TABLE_NAME+'''
            SET SAMPLE_NAME = %s,
                GENDER = %s,
                FILE_PATH = %s

            WHERE SAMPLE_NAME = %s and DIRECTORY_NAME = %s ;''' 
            #print(sql)
            cur2 = conn.cursor()
            cur2.execute(sql, (row['SAMPLE_NAME'], row['GENDER'], row['FILE_PATH'], row['SAMPLE_NAME'], row['DIRECTORY_NAME'] ))
            conn.commit()
            cur2.close()
        else:
            sql = "INSERT into "+TRUSIGHT_TABLE_NAME+"(SAMPLE_NAME, GENDER, STATUS, DIRECTORY_NAME, FILE_PATH) values(%s,%s, %s, %s, %s); "
            cur2 = conn.cursor()
            #print(sql)
            cur2.execute(sql, (row['SAMPLE_NAME'], row['GENDER'], row['STATUS'], row['DIRECTORY_NAME'], row['FILE_PATH'] ))
            conn.commit()
            cur2.close()
    return xldf

if __name__ == "__main__":
    connection = create_connection()
    #baseMountDir = "BaseMount"
    print("Reading the excel and upsert")
    df = read_excel_and_upsert(connection)
    print("reading excel completed")
    print("checking for fastq exists")
    checkFastqExists(connection)
    print("uploading fastq starts")
    uploadFastQ(connection)
    print("starting to populate status to excel")
    populateStatusInExcel(connection, df)   

    '''if  mountBaseSpace(baseMountDir):
        print("checking for fastq exists")
        checkFastqExists(connection, baseMountDir)
        print("uploading fastq starts")
        uploadFastQ(connection, baseMountDir)
        print("starting to populate status to excel")
        populateStatusInExcel(connection, df)'''

    connection.close()
