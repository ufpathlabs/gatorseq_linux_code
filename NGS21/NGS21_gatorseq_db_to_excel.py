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
import NGS21_database_connection

print(str(datetime.datetime.now()) + "\n")

script_path = os.path.dirname(os.path.abspath( __file__ ))
script_path = os.path.abspath(os.path.join(script_path, '..'))
CONFIG_FILE=script_path+"/linux_gatorseq.config.yaml"
config_dict=dict()
with open(CONFIG_FILE, 'r') as stream:
    try:
        config_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

CONFIG_TOKENS_FILE = script_path + "/" + config_dict['CONFIG_TOKENS_FILE']
CODE_ENV=script_path.split('/')[-2]
USER_NAME=os.environ['USER']

def replace_env(strname):
    strname=strname.replace("USER_NAME",USER_NAME).replace("CODE_ENV",CODE_ENV)
    return strname

GATOR_SEQ_SAMPLE_INPUT_FILE = replace_env(config_dict['GATOR_SEQ_SAMPLE_INPUT_FILE_NGS21'])
TABLE_NAME = replace_env(config_dict['TABLE_NAME'])
if CODE_ENV == "ProdEnv":
    TABLE_NAME = replace_env(config_dict['TABLE_NAME_PROD'])

file_to_lock = GATOR_SEQ_SAMPLE_INPUT_FILE + '.lock'
lock = FileLock(file_to_lock)
try:
    lock.acquire(timeout=1)
except:
    print("some other script is using it")
    sys.exit()

try: 
    excel_file = open(GATOR_SEQ_SAMPLE_INPUT_FILE, "r+")
except:
    print(" Could not open file! Please close Excel!")
    sys.exit()


def create_connection():
    return NGS21_database_connection.getSQLConnection(CONFIG_FILE, CONFIG_TOKENS_FILE, CODE_ENV)

def read_excel_and_upsert(conn):
    xldf_full = pd.read_excel(GATOR_SEQ_SAMPLE_INPUT_FILE,  engine='openpyxl')
    xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    xldf = xldf.replace(np.nan, '', regex=True)

    for index, row in xldf.iterrows():
        cur = conn.cursor()
        cur.execute("SELECT * FROM "+TABLE_NAME+" where SAMPLE_DIR_PATH = '" + row['SAMPLE_DIR_PATH'] + "'")
        rows = cur.fetchall()
        # print(rows)
        # print("-----------------")
        if len(rows) > 0:
            row = rows[0]
            xldf.at[index, "STATUS"] = row[2]
            xldf.at[index, "TIME_STAMP"] = row[3]
            xldf.at[index, "MESSAGE"] = row[4]
            # xldf.at[index, "PLMO_Number"] = row[5]
            # xldf.at[index, "Test_Product_Profile"] = row[6]
            # xldf.at[index, "Test_Product_Code"] = row[7]
            # xldf.at[index, "Diagnosis"] = row[8]
            # xldf.at[index, "Primary_Tumor_Site"] = row[9]
            # xldf.at[index, "Pre_Filter"] = row[10]
            # xldf.at[index, "Report_Template"] = row[11]
            # xldf.at[index, "QCIType"] = row[12]
            # xldf.at[index, "Treatments_Policy"] = row[13]
            # xldf.at[index, "Reporting_Method"] = row[14]

            xldf.at[index, "QCI_Upload_Message"] = row[15]
            xldf.at[index, "QCI_Download_Message"] = row[16]
            xldf.at[index, "EPIC_Upload_Message"] = row[17]
            xldf.at[index, "QCI_Re_Run"] = row[18]
            xldf.at[index, "EPIC_Re_Run"] = row[19]

        cur.close()
            
    try:
        xldf.to_excel(GATOR_SEQ_SAMPLE_INPUT_FILE, index=False)
    except:
        print("could not save excel")
        sys.exit()
        
connection = create_connection()
read_excel_and_upsert(connection)
connection.close()

