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

GATOR_SEQ_SAMPLE_INPUT_FILE = replace_env(config_dict['GATOR_SEQ_SAMPLE_INPUT_FILE'])
TABLE_NAME = replace_env(config_dict['TABLE_NAME'])
SQLITE_DB = replace_env(config_dict['SQLITE_DB'])

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


def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except:
        print(traceback.format_exc())
 
    return conn

def read_excel_and_upsert(conn):
    xldf_full = pd.read_excel(GATOR_SEQ_SAMPLE_INPUT_FILE)
    xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    xldf = xldf.replace(np.nan, '', regex=True)

    for index, row in xldf.iterrows():
        cur = conn.cursor()
        cur.execute("SELECT * FROM "+TABLE_NAME+" where SAMPLE_DIR_PATH = '" + row['SAMPLE_DIR_PATH'] + "'")
        rows = cur.fetchall()
        if len(rows) > 0:
            row = rows[0]
            xldf.at[index, "STATUS"] = row[1]
            xldf.at[index, "TIME_STAMP"] = row[2]
            xldf.at[index, "MESSAGE"] = row[3]
            # xldf.at[index, "PLMO_Number"] = row[4]
            # xldf.at[index, "Test_Product_Profile"] = row[5]
            # xldf.at[index, "Test_Product_Code"] = row[6]
            # xldf.at[index, "Diagnosis"] = row[7]
            # xldf.at[index, "Primary_Tumor_Site"] = row[8]
            # xldf.at[index, "Pre_Filter"] = row[9]
            # xldf.at[index, "Report_Template"] = row[10]
            # xldf.at[index, "QCIType"] = row[11]
            # xldf.at[index, "Treatments_Policy"] = row[12]
            # xldf.at[index, "Reporting_Method"] = row[13]

            xldf.at[index, "QCI_Upload_Message"] = row[14]
            xldf.at[index, "QCI_Download_Message"] = row[15]
            xldf.at[index, "EPIC_Upload_Message"] = row[16]
            xldf.at[index, "QCI_Re_Run"] = row[19]
            xldf.at[index, "EPIC_Re_Run"] = row[20]

        cur.close()
            
    try:
        xldf.to_excel(GATOR_SEQ_SAMPLE_INPUT_FILE, index=False)
    except:
        print("could not save excel")
        sys.exit()
        
connection = create_connection(SQLITE_DB)
read_excel_and_upsert(connection)
connection.close()

