import requests 
from requests.exceptions import HTTPError 
import time 
import ast
import sys
import yaml
import numpy
import pandas as pd
import hl7
import os
import datetime
from shutil import copyfile
from shutil import move
import csv
import traceback
import mysql.connector
from filelock import FileLock
from pathlib import Path
import os.path
import re
import math
print("Run start time: ", str(datetime.datetime.now()) + "\n")

script_path = os.path.dirname(os.path.abspath( __file__ ))
parent_path = os.path.abspath(os.path.join(script_path, '..'))

CONFIG_FILE = parent_path +"/linux_gatorseq.config.yaml"
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

COVID_19_EPIC_UPLOAD_TABLE = replace_env(config_dict['COVID_19_EPIC_UPLOAD_TABLE'])
COVID_19_TEST_INPUT_FOLDER = replace_env(config_dict['POOL_COVID_19_TEST_INPUT_FOLDER']) #'G:\DRL\Molecular\Assays\PGX\PGX_Beaker_Interface' 
CONFIG_TOKENS_FILE = parent_path  + "/" + config_dict['CONFIG_TOKENS_FILE']
MIRTH_GATORSEQ = config_dict['MIRTH_GATORSEQ']
if CODE_ENV=='ProdEnv':
    MIRTH_GATORSEQ += '/PROD'
    COVID_19_EPIC_UPLOAD_TABLE = replace_env(config_dict['COVID_19_EPIC_UPLOAD_TABLE_PROD'])
else:
    MIRTH_GATORSEQ += '/TEST'



config_token_dict=dict()
with open(CONFIG_TOKENS_FILE, 'r') as stream:
    try:
        config_token_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

if CODE_ENV == "ProdEnv":
    MYSQL_HOST = config_dict['PROD_MYSQL_HOST']
    MYSQL_USERNAME = config_dict['PROD_MYSQL_USERNAME']
    MYSQL_PASSWORD = config_token_dict['PROD_MYSQL_PASSWORD']
    MYSQL_DATABASE = config_dict['PROD_MYSQL_DATABASE']

if CODE_ENV == "DevEnv":
    MYSQL_HOST = config_dict['DEV_MYSQL_HOST']
    MYSQL_USERNAME = config_dict['DEV_MYSQL_USERNAME']
    MYSQL_PASSWORD = config_token_dict['DEV_MYSQL_PASSWORD']
    MYSQL_DATABASE = config_dict['DEV_MYSQL_DATABASE']


def check_folders_exist():
    #if not os.path.isfile(GATOR_SEQ_SAMPLE_INPUT_FILE):
        #sys.exit("ERROR: Does not have access to following folder: " + GATOR_SEQ_SAMPLE_INPUT_FILE + "\n") 

    if not os.path.isdir(MIRTH_GATORSEQ):
        sys.exit("ERROR: Does not have access to following folder: " + MIRTH_GATORSEQ + "\n") 

    if not os.path.isdir(COVID_19_TEST_INPUT_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + COVID_19_TEST_INPUT_FOLDER + "\n") 


check_folders_exist()
def format_for_unity(x):
    if(x<10):
        return "0"+str(x)
    else:
        return str(x)

def get_current_formatted_date():
    currentDT = datetime.datetime.now()
    if currentDT:
        data = format_for_unity(currentDT.year) +format_for_unity(currentDT.month) +format_for_unity(currentDT.month)+format_for_unity(currentDT.hour)+format_for_unity(currentDT.minute)+format_for_unity(currentDT.second)
        
        return data
    return str(currentDT)

def get_ticks(dt):
    return (dt - datetime.datetime(1, 1, 1)).total_seconds() * 10000000

class hl7update:
    
    def __init__(self, h):
        self.h = h

    def update_msh_segment(self):
        if self.h and self.h['MSH']:
            for msh_segment in self.h['MSH']:
                if msh_segment:
                    msh_segment[7] = get_current_formatted_date()
                    msh_segment[8] = ''
                    msh_segment[9][0][0] = 'ORU'
                    msh_segment[9][0][1] = 'R01'
                    msh_segment[10] = get_ticks(datetime.datetime.now())

    def update_orc_segment(self):
        if self.h and self.h['ORC']:
            for orc_segment in self.h['ORC']:
                orc_segment[1] = 'RE'

    def update_obr_segment(self):
        if self.h and self.h['OBR']:
            for obr_segment in self.h['OBR']:
                obr_segment[22] = get_current_formatted_date()
                obr_segment[25] = 'P'
                obr_segment[27] = '^^^^^R^^'

    def update_obx_segment(self):
        if self.h and self.h['OBX']:
            for obx_segment in self.h['OBX']:
                obx_segment[2] = 'ST'
                obx_segment[11] = 'P'
                obx_segment[14] = get_current_formatted_date()
                if(len(obx_segment)==19):
                    obx_segment.append(obx_segment[14])
                elif(len(obx_segment)>=19):
                    obx_segment[19] = obx_segment[14]
   
    def update_obx_seg_containing_gene(self, result):
        updates = 0
        temp_obx = self.h[:]
        l = len(self.h)
        for i in range(l):
            del temp_obx[l-i-1]
        new_obx_index = 1
        for obxSegment in self.h['OBX']:
            if obxSegment[3][0][1][0] == "SARS-COV-2, NAA":
                obxSegment[5][0] = result
                obxSegment[1] = new_obx_index
                new_obx_index +=1 
                temp_obx.append(obxSegment) 
            
        h_t = self.h[:]
        l = len(self.h)
        for i in range(l):
            del h_t[l-i-1]
        for i in range(len(self.h)):
            if(self.h[i][0][0]!="OBX"):
                h_t.append(self.h[i])
        h_t.extend(temp_obx)
        return h_t

    def get_first_obx_index(self):
        idx = 0
        for seg in self.h:
            if seg[0][0] == 'OBX':
                return idx
            idx += 1
        return -1
    # Assuming insertion is just above first OBX segment
    def update_comments(self, comments):
        comments_arr = comments.split("\n")
        obx_idx = self.get_first_obx_index()
        if (obx_idx == -1):
            print("OBX segment not found, so appending it")
            obx_idx = len(self.h) - 1 
        i=1
        for comment in comments_arr:
            self.h.append('NTE|{}|L|{}'.format(i,comment))
            obx_idx += 1
            i += 1


def create_connection():
    conn = None
    try:
        conn = mysql.connector.connect(
            host=MYSQL_HOST,
            user=MYSQL_USERNAME,
            passwd=MYSQL_PASSWORD,
            database=MYSQL_DATABASE
        )
    except:
        print(traceback.format_exc())
    return conn
SQL_CONNECTION = create_connection()


def updateRowInDatabase(containerId, PLMO, MRN, ptName, ptSex, ptAge, ordDept, excelFileName):
    cur = SQL_CONNECTION.cursor()

    updateSql = "UPDATE "+ COVID_19_EPIC_UPLOAD_TABLE +" set PLMO_Number= %s, MRN = %s, PATIENT_NAME = %s, PATIENT_SEX = %s, PATIENT_AGE = %s, ORDERING_DEPARTMENT = %s, EPIC_UPLOAD_TIMESTAMP = %s WHERE CONTAINER_ID = %s and SOURCE_EXCEL_FILE = %s;"
    #print(sample.completeSampleName, type(sample.name), "!!!!!!",(sample.name, PLMO, get_current_formatted_date(), sample.nCoV_N1, sample.nCoV_N2, sample.nCoV_N3, sample.RP, sample.result))
    cur.execute(updateSql, (PLMO[0], MRN, ptName, ptSex, ptAge, ordDept, str(datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S")), containerId, excelFileName))
    SQL_CONNECTION.commit()
    cur.close()


# arguments are <class Sample> sample, str PLMO
def addRowInDatabase(containerId, result, PLMO, MRN, ptName, ptSex, ptAge, ordDept, excelFileName):
    cur = SQL_CONNECTION.cursor()

    findsql = "SELECT * from " + COVID_19_EPIC_UPLOAD_TABLE + " where CONTAINER_ID = %s and SOURCE_EXCEL_FILE = %s;"
    cur.execute(findsql, (containerId, excelFileName))
    rows = cur.fetchall()
    # print(rows)
    # print("-----------------")
    if len(rows) > 0:
        updateSql = "UPDATE " + COVID_19_EPIC_UPLOAD_TABLE + " set QUANTSTUDIO_SPECIMEN_ID = %s, RESULT = %s where CONTAINER_ID = %s and SOURCE_EXCEL_FILE = %s;" 
        cur.execute(updateSql, ("pooled", result, containerId, excelFileName))
        SQL_CONNECTION.commit()
    else:
        insertSql = "INSERT INTO "+ COVID_19_EPIC_UPLOAD_TABLE +" VALUES(%s, %s,%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);"
        cur.execute(insertSql, ("pooled", containerId, PLMO, MRN, ptName, ptSex, ptAge, ordDept, excelFileName, "", "", "", "", "", result))
        SQL_CONNECTION.commit()
    cur.close()

# method to write the entire database table to excel
def writeDataToExcel(excelName, containerToResult):
    xldf = pd.read_sql_query('select * from '+ COVID_19_EPIC_UPLOAD_TABLE +' where SOURCE_EXCEL_FILE = "'+ excelName +'" ;', SQL_CONNECTION)
    xldf = xldf.drop("SOURCE_EXCEL_FILE", 1)
    xldf = xldf.drop("ORDERING_DEPARTMENT", 1)
    xldf = xldf.drop("2019nCoV_N3", 1)
    xldf = xldf.drop("2019nCoV_N2", 1)
    xldf = xldf.drop("2019nCoV_N1", 1)
    xldf = xldf.drop("RP", 1)

    cols = xldf.columns.tolist()
    upload_col = cols.pop(cols.index("EPIC_UPLOAD_TIMESTAMP"))
    cols.insert(len(cols), upload_col)
    xldf = xldf.reindex(columns= cols)
    RESULT_LOG = COVID_19_TEST_INPUT_FOLDER + "/" + \
        excelName.split("/")[-1].replace("_SAMPLE_RESULTS_UPDATED_ID","_SAMPLE_EPIC_UPLOAD_LOG")
    SAMPLE_MAP_FILE = excelName.replace("_SAMPLE_RESULTS_UPDATED_ID", "_SAMPLE_MAP")

    sampleMapDf = pd.read_excel(SAMPLE_MAP_FILE)
    sampleMapDf["UPLOADED_TO_EPIC"] = "No"
    for index, row in sampleMapDf.iterrows():
        if xldf['QUANTSTUDIO_SPECIMEN_ID'].str.contains(str(row["Internal_Sample_ID"])).any() :
            sampleMapDf.at[index, "UPLOADED_TO_EPIC"] = "Yes"
    
    try:
        # xldf.to_excel(RESULT_LOG , index=False)
        with pd.ExcelWriter(RESULT_LOG) as writer:
            xldf.to_excel(writer, index=False, sheet_name='Sheet 1')
            #sampleMapDf.to_excel(writer, index=False, sheet_name='Sheet 2')
        print("done writeToExcel method and writing done to -->", RESULT_LOG )
    except:
        print("unable to save status excel, please close it")

def checkIncomingHl7(sampleDict, excelFile):
    UPLOAD_PATH = MIRTH_GATORSEQ + '/RESULTS'
    ORDERS_ARCHIVE_DIR = MIRTH_GATORSEQ + '/ORDERS_ARCHIVE/'
    ORDERS_DIR = MIRTH_GATORSEQ + '/ORDERS/'
    allhl7filenames = []
    for (dirpath, dirnames, filenames) in os.walk(ORDERS_DIR):
        allhl7filenames.extend(filenames)
        break
    for hl7_file_name in allhl7filenames:
        try:
            hl7file = open(ORDERS_DIR + hl7_file_name, mode="r").read()
        except:
            continue
        arr = hl7file.split("\n\n") #split by blank lines
        #Iterate all HL7 messages in file
        c=0
        for hl7msg in arr: 
            if hl7msg: #check if message not empty
                msg_unix_fmt = hl7msg.replace("\n","\r")
                h = hl7.parse(msg_unix_fmt)
                newHl7 = hl7update(h)
                #Read message id from HL7 message
                try:
                    messageId = str(h['OBR'][0][3]).replace("^", "")
                except:
                    continue
                if (not messageId):
                    continue
                plm = None
                try:
                    plm = h['ORC'][0][2]
                except:
                    print("---------could not find PLMO in hl7 file: ", hl7_file_name)
                    continue
                
                mrn = ""
                try:
                    mrn = h['PID'][0][3][0][0]
                except:
                    print("---------could not find MRN in hl7 file: ", hl7_file_name)
        
                ptName = ""
                try:
                    ptName = h['PID'][0][5][0]
                except:
                    print("---------could not find PATIENT_NAME in hl7 file: ", hl7_file_name) 

                ptSex = ""
                try:
                    ptSex = h['PID'][0][8][0]
                except:
                    print("---------could not find PATIENT_SEX in hl7 file: ", hl7_file_name) 

                ptAge = -1
                try:
                    ptAge = 2020 - int(h['PID'][0][7][0][:4]) 
                except:
                    print("---------could not find PATIENT_AGE in hl7 file: ", hl7_file_name) 

                ordDept = ""
                try:
                    ordDept = h['OBR'][0][15][0]
                except:
                    print("---------could not find Ordering_DEPT in hl7 file: ", hl7_file_name) 


                # search for messageId in the sampleDict
                #if messageId == "100047187": #100047166  100047187
                if messageId in sampleDict or plm[0] in sampleDict:
                   # print("--------found----------")
                    if sampleDict.get(messageId) is not None:
                        givenSample =  sampleDict.get(messageId)
                    else:
                        givenSample = sampleDict.get(plm[0])
                    #print("processing hl7 input file: ", hl7_file_name)                   
                    newHl7.update_msh_segment()
                    newHl7.update_orc_segment()
                    newHl7.update_obr_segment()
                    newHl7.update_obx_segment()
                    h = newHl7.update_obx_seg_containing_gene( givenSample.result )
                    
                    out_file_path = UPLOAD_PATH + '/hl7-pooled-COVID_19-{}-output.txt'.format(messageId)
                    if h:
                        with open(out_file_path, 'w' ,  encoding='utf-8') as f:
                            f.write(str(h))
                        print("Out file available at :",out_file_path)
                        move(ORDERS_DIR + hl7_file_name, ORDERS_ARCHIVE_DIR + 'POOLED_COVID_19_processed_' + get_current_formatted_date() + "-" + hl7_file_name) 
                        if plm:
                            updateRowInDatabase(messageId, plm, str(mrn), str(ptName), str(ptSex), str(ptAge), str(ordDept), excelFile )
                    

class Sample:
    def __init__(self, sample_name, completeSampleName):
        self.name = str(sample_name)
        self.nCoV_N1 = None
        self.nCoV_N2 = None
        self.nCoV_N3 = None
        self.RP = None
        self.result = None
        self.completeSampleName = completeSampleName

    def __str__(self):
        return str( self.completeSampleName) + ": " + str(self.nCoV_N1) + " & " + str(self.nCoV_N2) + " & "  + str(self.nCoV_N3) + " & " + str(self.RP) + " & " + str(self.result)

    def __repr__(self):
        return str( self.completeSampleName) + ": " + str(self.nCoV_N1) + " & " + str(self.nCoV_N2) + " & "  + str(self.nCoV_N3) + " & " + str(self.RP) + " & " + str(self.result) + " | "

def isFloatValue(value, maxThreshold):
    try:
        val = float(value)
        #ToDo: check RP condition, if RP > 35 then is it considered to be Undetermined
        if maxThreshold and val > maxThreshold:
            return False
        return True

    except:
        if value == "Undetermined":
            return False
        else:
            print("-----------ERROR: unable to identify the value-------------->", value)
            return False
    
def addSampleDictToDatabase(sampleResult, excelName):
    for containerId in sampleResult.keys():
        result = sampleResult[containerId]
        addRowInDatabase(containerId, result, "", "", "", "", "", "", excelName)
    
    print("-> all samples added to database <-")


if __name__ == "__main__":
    os.chdir(COVID_19_TEST_INPUT_FOLDER)
    _files = filter(os.path.isfile, os.listdir(COVID_19_TEST_INPUT_FOLDER))
    excel_files = [os.path.join(COVID_19_TEST_INPUT_FOLDER, f) for f in _files if "$" not in f] # add path to each file

    fileNames = []
    toProcess = []
    for f in excel_files:
        if "_SAMPLE_POOL" in f:
            sampleGroupName = f[:f.index("_SAMPLE_POOL")]
            if sampleGroupName + "_FINAL.xlsx" not in excel_files:
                fileNames.append(sampleGroupName)
    
    for f in fileNames:
        sampleResult = f + "_SAMPLE_POOL_RESULTS.xlsx"
        sampleMap = f + "_SAMPLE_POOL.xlsx"
        sampleToPool = {}
        df = pd.read_excel(sampleMap)
        for i, row in df.iterrows():
            sampleToPool[row["Container_ID"]] = row["Pooled_ID"]
        
        results_df = pd.read_excel(sampleResult)#, skiprows=range(0,35))
        pool_results = {}
        for i, row in results_df.iterrows():
            pool_results[row["Pooled_ID"]] = row["Result"]
        
        #contains results corresponding to each containerId
        sampleToResult = {}
        for i, row in df.iterrows():
            sampleToResult[row["Container_ID"]] = pool_results[row["Pooled_ID"]]

        containerToResult = {}
        for containerId in sampleToResult:
            if sampleToResult[containerId] == "Negative":
                containerToResult[containerId.replace("\\", "")  ] = "Not Detected"
        
        print(containerToResult)

        #add all samples to database
        addSampleDictToDatabase(containerToResult, f)

        #logic to add the corresponding hl7 file to RESULTS folder
        checkIncomingHl7(containerToResult, f)

        print(f)
        writeDataToExcel(f, containerToResult)


        
    #writeDataToExcel("/ext/path/DRL/Molecular/COVID19/COVID_19_QuantStudio/ProdEnv/Results/2020-03-20 203810_QuantStudio_export_UPDATED_CONTAINER_ID.xlsx")
    SQL_CONNECTION.close()
