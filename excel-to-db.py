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
        # print(row['SAMPLE_DIR_PATH'], row['STATUS'], row['TIME_STAMP'], row['MESSAGE'], row['PLMO_Number'], row['Test_Product_Profile'], row['Test_Product_Code'], row['Diagnosis'], row['Primary_Tumor_Site'], row['Pre_Filter'], row['Report_Template'], row['QCIType'], row['Treatments_Policy'], row['Reporting_Method'])
        # update_task(conn, (row['SAMPLE_DIR_PATH'], row['STATUS'], row['PLMO_Number'], row['Test_Product_Profile'], row['Test_Product_Code'], row['Diagnosis'], row['Primary_Tumor_Site'], row['Pre_Filter'], row['Report_Template'], row['QCIType'], row['Treatments_Policy'], row['Reporting_Method']))
        cur = conn.cursor()
        cur.execute("SELECT * FROM "+TABLE_NAME+" where SAMPLE_DIR_PATH = '" + task[0] + "'")
        rows = cur.fetchall()
        if len(rows) > 0:
            sql = '''
            UPDATE '''+TABLE_NAME+'''
            SET SAMPLE_DIR_PATH = ?,
                STATUS = ?,
                PLMO_Number = ?,
                Test_Product_Profile = ?,
                Test_Product_Code = ?,
                Diagnosis = ?,
                Primary_Tumor_Site = ?,
                Pre_Filter = ?,
                Report_Template = ?,
                QCIType = ?,
                Treatments_Policy = ?,
                Reporting_Method = ?
            WHERE SAMPLE_DIR_PATH = ''' + "'" + row['SAMPLE_DIR_PATH'] + "';"
            cur2 = conn.cursor()
            cur2.execute(sql, (row['SAMPLE_DIR_PATH'], row['STATUS'], row['PLMO_Number'], row['Test_Product_Profile'], row['Test_Product_Code'], row['Diagnosis'], row['Primary_Tumor_Site'], row['Pre_Filter'], row['Report_Template'], row['QCIType'], row['Treatments_Policy'], row['Reporting_Method']))
            conn.commit()
            cur2.close()
        else:
            sql = ''' INSERT into '''+TABLE_NAME+'''(SAMPLE_DIR_PATH, STATUS, TIME_STAMP, MESSAGE, PLMO_Number, 
        Test_Product_Profile , Test_Product_Code, Diagnosis, Primary_Tumor_Site, Pre_Filter, 
                    Report_Template, QCIType, Treatments_Policy, Reporting_Method, QCI_Upload_Message, QCI_Download_Message, EPIC_Upload_Message) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?); '''
            #     ON CONFLICT(SAMPLE_DIR_PATH)
            #     DO UPDATE
            #   SET STATUS = excluded.STATUS
            #   ;'''
            cur2 = conn.cursor()
            cur2.execute(sql, (row['SAMPLE_DIR_PATH'], row['STATUS'], row['TIME_STAMP'], row['MESSAGE'], row['PLMO_Number'],  row['Test_Product_Profile'], row['Test_Product_Code'], row['Diagnosis'], row['Primary_Tumor_Site'], row['Pre_Filter'], row['Report_Template'], row['QCIType'], row['Treatments_Policy'], row['Reporting_Method'], "", "", ""))
            conn.commit()
            cur2.close()
    
    conn.close()

connection = create_connection(SQLITE_DB)
read_excel_and_upsert(connection)

