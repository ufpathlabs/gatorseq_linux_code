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
COVID_19_TEST_INPUT_FOLDER = replace_env(config_dict['COVID_19_TEST_INPUT_FOLDER']) #'G:\DRL\Molecular\Assays\PGX\PGX_Beaker_Interface' 
COVID_19_TEST_SAMPLE_LOG = replace_env(config_dict['COVID_19_TEST_SAMPLE_LOG'])
CONFIG_TOKENS_FILE = parent_path  + "/" + config_dict['CONFIG_TOKENS_FILE']
MIRTH_GATORSEQ = config_dict['MIRTH_GATORSEQ']
if CODE_ENV=='ProdEnv':
    MIRTH_GATORSEQ += '/PROD'
    COVID_19_EPIC_UPLOAD_TABLE = replace_env(config_dict['COVID_19_EPIC_UPLOAD_TABLE_PROD'])
    COVID_19_TEST_INPUT_FOLDER = replace_env(config_dict['COVID_19_TEST_INPUT_FOLDER_PROD']) 
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
   
    def update_obx_seg_containing_gene(self, sample):
        updates = 0
        temp_obx = self.h[:]
        l = len(self.h)
        for i in range(l):
            del temp_obx[l-i-1]
        new_obx_index = 1
        for obxSegment in self.h['OBX']:
            if obxSegment[3][0][1][0] == "SARS-COV-2, NAA":
                obxSegment[5][0] = sample.result
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
# arguments are <class Sample> sample, str PLMO
def addRowInDatabase(sample, PLMO, MRN, ptName, excelFileName):
    cur = SQL_CONNECTION.cursor()
    updateSql = "INSERT INTO "+ COVID_19_EPIC_UPLOAD_TABLE +" VALUES(%s,%s, %s, %s, %s, %s, %s, %s, %s, %s, %s);"
    #print(sample.completeSampleName, type(sample.name), "!!!!!!",(sample.name, PLMO, get_current_formatted_date(), sample.nCoV_N1, sample.nCoV_N2, sample.nCoV_N3, sample.RP, sample.result))
    cur.execute(updateSql, (sample.completeSampleName, PLMO[0], MRN, ptName, excelFileName, str(datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S")), sample.nCoV_N1, sample.nCoV_N2, sample.nCoV_N3, sample.RP, sample.result))
    SQL_CONNECTION.commit()
    cur.close()

# method to write the entire database table to excel
def writeDataToExcel(excelName):
    xldf = pd.read_sql_query('select * from '+ COVID_19_EPIC_UPLOAD_TABLE +' where SOURCE_EXCEL_FILE = "'+ excelName +'" ;', SQL_CONNECTION)
    xldf = xldf.drop("SOURCE_EXCEL_FILE", 1)
    try:
        xldf.to_excel(COVID_19_TEST_SAMPLE_LOG + "_1" + excelName.split("/")[-1], index=False)
    except:
        print("unable to save status excel, please close it")
    print("--------done writeToExcel method------------")
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
                # search for messageId in the sampleDict
                #if messageId == "100047187": #100047166  100047187
                if messageId in sampleDict or plm[0] in sampleDict:
                   # print("--------found----------")
                    if sampleDict.get(messageId) is not None:
                        givenSample =  sampleDict.get(messageId)
                    else:
                        givenSample = sampleDict.get(plm[0])
                    print("processing hl7 input file: ", hl7_file_name)                   
                    newHl7.update_msh_segment()
                    newHl7.update_orc_segment()
                    newHl7.update_obr_segment()
                    newHl7.update_obx_segment()
                    h = newHl7.update_obx_seg_containing_gene( givenSample )
                    
                    out_file_path = UPLOAD_PATH + '/hl7-{}-output.txt'.format(messageId)
                    if h:
                        with open(out_file_path, 'w' ,  encoding='utf-8') as f:
                            f.write(str(h))
                        print("---> Out file available at :",out_file_path, "<---")
                        move(ORDERS_DIR + hl7_file_name, ORDERS_ARCHIVE_DIR + 'COVID_19_processed_' + get_current_formatted_date() + "-" + hl7_file_name) 
                        if plm:
                            addRowInDatabase(givenSample, plm, str(mrn), str(ptName), excelFile )
                    

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
    
def addSampleDictToExcel(sampleDict, excelName, writeFlag):
    #print(status_df.head())
    #print(sampleDict)
    status_df = pd.DataFrame(columns=['SAMPLE_NAME', 'PLMO', 'MRN', 'TIMESTAMP', 'N1', 'N2', 'RP', 'RESULT'])
    for key in sampleDict.keys():
        curSample = sampleDict[key]
        if curSample.completeSampleName in status_df["SAMPLE_NAME"]:
            continue
        else:
            #addRowInDatabase(curSample, "", "")
            status_df.loc[len(status_df)] = [curSample.completeSampleName, "", "", "", curSample.nCoV_N1, curSample.nCoV_N2, curSample.RP, curSample.result]
    #print(len(status_df))
    if writeFlag:
        status_df.to_excel(COVID_19_TEST_SAMPLE_LOG + "_" + excelName.split("/")[-1], index=False)
        print("-> status file written successfully <-")


if __name__ == "__main__":
    os.chdir(COVID_19_TEST_INPUT_FOLDER)
    _files = filter(os.path.isfile, os.listdir(COVID_19_TEST_INPUT_FOLDER))
    excel_files = [os.path.join(COVID_19_TEST_INPUT_FOLDER, f) for f in _files if "$" not in f] # add path to each file

    fileNames = []
    toProcess = []
    for f in excel_files:
        if "_SAMPLE_MAP" in f:
            sampleGroupName = f[:f.index("_SAMPLE_MAP")]
            if sampleGroupName + "_RESULTED.xlsx" not in excel_files:
                fileNames.append(sampleGroupName)
    
    for f in fileNames:
        sampleResult = f + "_SAMPLE_RESULTS.xlsx"
        sampleMap = f + "_SAMPLE_MAP.xlsx"
        sampleToContainer = {}
        df = pd.read_excel(sampleMap)
        for i, row in df.iterrows():
            sampleToContainer[row["Internal_Sample_ID"]] = row["Container_ID"]
        
        results_df = pd.read_excel(sampleResult, skiprows=range(0,35))
        results_df["CONTAINER_ID"] = None
        
        
        for index, row in results_df.iterrows():
            results_df.at[index, "CONTAINER_ID"] = sampleToContainer[row["Sample Name"]]
        
        print(results_df.head())
        toProcess.append(f + "_UPDATED_CONTAINER_ID.xlsx")
        results_df.to_excel(f + "_UPDATED_CONTAINER_ID.xlsx", index=False)
        


    for f in toProcess:
        sampleDict = {}
        plmoDict = {}
        eachExcel = os.path.join(COVID_19_TEST_INPUT_FOLDER, f)
        xldf = pd.read_excel(eachExcel)
        #xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
        # generate sample dictionary in below format:
        # {<sample name>: <Class Sample>}
        #print(xldf.head())  
        for index, row in xldf.iterrows():
            sampleName = str(row["CONTAINER_ID"])
            if 'PLMO' in sampleName:
                plm = re.findall(r"PLMO\d+-\d+" , sampleName)
                if len(plm):
                    sampleName = plm[0]

            sampleName = sampleName.replace("\\", "")            

            targetName = row["Target Name"]
            value = row["CT"]
            
            #ToDo: sampleName = <PLMO of a sample> or <id number of a sample>
            if sampleDict.get(sampleName) is None:
                sampleDict[sampleName] = Sample(sampleName, str(row["Sample Name"]))
            if targetName == "RP":# and not math.isnan(value):
                setattr(sampleDict[sampleName], "%s" % targetName, value)
            else:# not math.isnan(value):
                #ToDO: is there any better of deriving 'nCoV_N1' from '2019nCoV_N1'? (python does not allow variable names to start with number)
                setattr(sampleDict[sampleName], "%s" % targetName[4:], value)

        for sampleName in sampleDict.keys():
            sample = sampleDict[sampleName]

            if all([isFloatValue(sample.nCoV_N1, 40), isFloatValue(sample.nCoV_N2, 40)]):#, isFloatValue(sample.nCoV_N3, None)]):
                sample.result = "Detected"
                continue
            elif any([not isFloatValue(sample.nCoV_N1, 40), not isFloatValue(sample.nCoV_N2, 40)]) and not all([not isFloatValue(sample.nCoV_N1, 40), not isFloatValue(sample.nCoV_N2, 40)]):
                sample.result = "Indeterminate"
                continue
            
            if (sample.RP and not isFloatValue(sample.RP, 40.0)) or ( type(sample.nCoV_N1) == "float" and sample.nCoV_N1 > 40 and type(sample.nCoV_N2) == "float" and sample.nCoV_N2 > 40  ):
                #INVALID Case is to be handled by pathologists separately and hence just continuing
                sample.result = "Invalid"
                continue
            if all([not isFloatValue(sample.nCoV_N1, None), not isFloatValue(sample.nCoV_N2, None)]):
                sample.result = "Not Detected"
            
            if sample.result is None:
                print("------unable to identify result for the sample-----", sample)
                #del sampleDict[sampleName]
                sample.result = "Invalid"
        #print("below is the dictionary of all samples:")
        #print(sampleDict["PLMO20-000129"])
        #print(sampleDict)
        checkIncomingHl7(sampleDict, f)
        writeDataToExcel(f)  
        #addSampleDictToExcel(sampleDict, f, True) 
        
    #writeDataToExcel("/ext/path/DRL/Molecular/COVID19/COVID_19_QuantStudio/ProdEnv/Results/2020-03-20 203810_QuantStudio_export_UPDATED_CONTAINER_ID.xlsx")
    SQL_CONNECTION.close()
