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
from illumina.truSight import runBashCommand



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


def updateRowWithStatus(sampleName, status, conn):
    cur = conn.cursor()
    updateSql = "update "+ TRUSIGHT_TABLE_NAME +" set STATUS = %s where SAMPLE_NAME = %s;"
    cur.execute(updateSql, (status, sampleName))
    conn.commit()
    cur.close()



# create a sample using REST API of trusight 
# use upload command of trusight to upload the sample
def createSampleAndUpload(conn, baseMountDir):
    cur = conn.cursor()
    cur.execute("SELECT * FROM "+TRUSIGHT_TABLE_NAME+" where STATUS = 'FASTQ_UPLOADED';")
    rows = cur.fetchall()
    cur.close()
    for row in rows:
        sampleName = row[0]
        headers = {
            "Authorization": "ApiKey " + TRUSIGHT_API_KEY,
            "X-ILMN-Domain": "ufl-tss",
            "X-ILMN-Workgroup": "51b925ca-56ed-37c9-89d7-85b83a1f7e55",
            "Content-Type": "application/json"
        }
        data = {
            "displayID": "api_created_case_" + sampleName,
            "testDefinitionId": "1cb2c841-9a25-4cb3-8b55-c3ea0d639086",
            "subjects": [
                {
                "gender": "Male",
                "samples": [
                    {
                    "externalSampleId": sampleName
                    }
                ]
                }
            ]
        }
        response = requests.post(TRUSIGHT_NEW_CASE_URL, headers=headers, params=data)
        #ToDO: add success logic
        if response.json() is not None:
            updateRowWithStatus(sampleName, "UPLOAD_PENDING", conn)
            print("added a sample successfully")
        return response.json()["status"]

def uploadFastQ(conn, baseMountDir):
    cur = conn.cursor()
    cur.execute("SELECT * FROM "+TRUSIGHT_TABLE_NAME+" where STATUS = 'UPLOAD_FASTQ_PENDING';")
    rows = cur.fetchall()
    cur.close()
    for row in rows:
        sampleName = row[0]
        cmd = "java -jar "+ TRUSIGHT_CLI + " stage --stageDirectory=" + sampleName + "/ --localDirectory=" +  baseMountDir + "/Projects/WGS/Samples/" + sampleName + "/Files/"
        statusJson = runBashCommand(cmd)
        if statusJson:
            updateRowWithStatus(sampleName, "FASTQ_UPLOADED", conn)
        else:
            print("status is None")
            # updateRowWithStatus(sampleName, "ERROR_UPLOADING", conn)

if __name__ == "__main__":
    connection = create_connection()
    baseMountDir = "BaseMount"
    uploadFastQ(connection, baseMountDir)

    connection.close()



