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
from truSight import updateRowWithStatusAndMessage
from truSight import updateRowWithStatus
import time



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

TRUSIGHT_APP_SUBMISSION_INPUT_FILE = replace_env(config_dict['ILLUMINA_BASESPACE_APP_SUBMISSION_INPUT_FILE'])
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

def uploadFastQ(conn, baseMountDir):
    cur = conn.cursor()
    cur.execute("SELECT * FROM "+TRUSIGHT_TABLE_NAME+" where STATUS = 'UPLOAD_FASTQ_PENDING';")
    rows = cur.fetchall()
    cur.close()
    for row in rows:
        sampleName = row[0]
        directoryName = row[4]
        startTime = time.time()
        cmd = "java -jar "+ TRUSIGHT_CLI + " stage --stageDirectory=" + directoryName + "/ --localDirectory=" +  baseMountDir + "/Projects/WGS/Samples/" + sampleName + "/Files/"
        print("running the following command: ", cmd)
        statusJson = runBashCommand(cmd, False)
        endTime = time.time()
        timeForExecution = round((endTime - startTime)/60)
        if statusJson:
            # updateRowWithStatus(sampleName, "FASTQ_UPLOADED", conn)
            updateRowWithStatusAndMessage(sampleName, "FASTQ_UPLOADED", "file uploaded in: " + str(timeForExecution) + " minutes", directoryName, conn)
        else:
            print("status is None")
            updateRowWithStatusAndMessage(sampleName, "ERROR_UPLOADING", "Error uploading the file to trusight", directoryName, conn)

if __name__ == "__main__":
    connection = create_connection()
    baseMountDir = "BaseMount"
    uploadFastQ(connection, baseMountDir)

    connection.close()



