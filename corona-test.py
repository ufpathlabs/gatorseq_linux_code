import requests 
from requests.exceptions import HTTPError 
import time 
import ast
import sys
import yaml
import numpy
import pandas as pd
import hl7update
import hl7
import os
import datetime
from shutil import copyfile
from shutil import move
import csv

print(str(datetime.datetime.now()) + "\n")

MIRTH_GATORSEQ = 'Z:/MIRTH_GATORSEQ'
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


LINUX_ANALYSIS_OUT_FOLDER = replace_env(config_dict['LINUX_ANALYSIS_OUT_FOLDER'])
#GATOR_SEQ_SAMPLE_INPUT_FILE = replace_env(config_dict['GATOR_SEQ_SAMPLE_INPUT_FILE'])
GATOR_PGX_SAMPLE_INPUT_FOLDER = replace_env(config_dict['GATOR_PGX_SAMPLE_INPUT_FOLDER']) #'G:\DRL\Molecular\Assays\PGX\PGX_Beaker_Interface' 
CONFIG_TOKENS_FILE = script_path + "/" + config_dict['CONFIG_TOKENS_FILE']
MIRTH_GATORSEQ = config_dict['MIRTH_GATORSEQ']
if CODE_ENV=='ProdEnv':
    MIRTH_GATORSEQ += '/PROD'
else:
    MIRTH_GATORSEQ += '/TEST'

def check_folders_exist():
    #if not os.path.isfile(GATOR_SEQ_SAMPLE_INPUT_FILE):
        #sys.exit("ERROR: Does not have access to following folder: " + GATOR_SEQ_SAMPLE_INPUT_FILE + "\n") 

    if not os.path.isdir(LINUX_ANALYSIS_OUT_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + LINUX_ANALYSIS_OUT_FOLDER + "\n")


    if not os.path.isdir(MIRTH_GATORSEQ):
        sys.exit("ERROR: Does not have access to following folder: " + MIRTH_GATORSEQ + "\n") 

    if not os.path.isdir(GATOR_PGX_SAMPLE_INPUT_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + GATOR_PGX_SAMPLE_INPUT_FOLDER + "\n") 


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

            if obxSegment[3][0][1][0] == "SARS-CoV-2, NAA":
                #print("tumor type found")
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


def checkIncomingHl7(sampleDict):
    UPLOAD_PATH = MIRTH_GATORSEQ + '/RESULTS1'
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

                # search for messageId in the sampleDict
                #if messageId == "100047187": #100047166  100047187
                foundMessageId = False
                if messageId in sampleDict:
                    foundMessageId = True
                    newHl7.update_msh_segment()
                    newHl7.update_orc_segment()
                    newHl7.update_obr_segment()
                    newHl7.update_obx_segment()
                    h = newHl7.update_obx_seg_containing_gene( sampleDict[messageId] )
                    
                    out_file_path = UPLOAD_PATH + '/hl7-{}-output.txt'.format(messageId)
                    if h:
                        with open(out_file_path, 'w' ,  encoding='utf-8') as f:
                            f.write(str(h))
                        print("Out file available at :",out_file_path)
                        # move(ORDERS_DIR + hl7_file_name, ORDERS_ARCHIVE_DIR + 'COVID-processed-' + get_current_formatted_date() + "-" + hl7_file_name) 

                    

class Sample:
    def __init__(self, sample_name):
        self.name = sample_name
        self.nCoV_N1 = None
        self.nCoV_N2 = None
        self.nCoV_N3 = None
        self.RP = None
        self.result = None

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
            print("------------unable to identify the value-------------->", value)
            return False
    
    
if __name__ == "__main__":
    CORONAVIRUS_TEST_INPUT_FOLDER = replace_env(config_dict['CORONAVIRUS_TEST_INPUT_FOLDER'])
    _files = filter(os.path.isfile, os.listdir(CORONAVIRUS_TEST_INPUT_FOLDER))
    excel_files = [os.path.join(CORONAVIRUS_TEST_INPUT_FOLDER, f) for f in _files] # add path to each file
    sampleDict = {}
    for eachExcel in excel_files:
        xldf_full = pd.read_excel(eachExcel, skiprows=range(0,34))
        xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
        # generate sample dictionary in below format:
        # {<sample name>: <Class Sample>}

        for index, row in xldf.iterrows():
            sampleName = row["Sample Name"]
            targetName = row["Target Name"]
            value = row["CT"]
            if sampleDict.get(sampleName) is None:
                sampleDict[sampleName] = Sample(sampleName)

            if targetName == "RP":
                sampleDict[sampleName][targetName] = value
            else:
                #ToDO: is there any better of deriving 'nCoV_N1' from '2019nCoV_N1'? (python does not allow variable names to start with number)
                sampleDict[sampleName][targetName[4:]] = value

    for sampleName in sampleDict.keys():
        sample = sampleDict[sampleName]

        if all([isFloatValue(sample.nCoV_N1, None), isFloatValue(sample.nCoV_N2, None), isFloatValue(sample.nCoV_N3, None)]):
            sample.result = "Positive"
            continue
        elif any([not isFloatValue(sample.nCoV_N1, None), not isFloatValue(sample.nCoV_N2, None), not isFloatValue(sample.nCoV_N3, None)]):
            sample.result = "Indeterminate"
            continue
        
        if sample.RP and isFloatValue(sample.RP, 35.0):
            sample.result = "Invalid"
            continue
        if all([not isFloatValue(sample.nCoV_N1, None), not isFloatValue(sample.nCoV_N2, None), not isFloatValue(sample.nCoV_N3, None)]):
            sample.result = "Negative"



        # if sample.RP and isFloatValue(sample.RP, 35.0):
        #     sample.result = "Invalid"
        #     continue

        # if sample.nCoV_N1 and sample.nCoV_N2 and sample.nCoV_N3 and sample.RP:
        #     if all([not isFloatValue(sample.nCoV_N1, None), not isFloatValue(sample.nCoV_N2, None), not isFloatValue(sample.nCoV_N3, None)]):
        #         sample.result = "Negative"
            
        #     elif any([not isFloatValue(sample.nCoV_N1, None), not isFloatValue(sample.nCoV_N2, None), not isFloatValue(sample.nCoV_N3, None)]):
        #         sample.result = "Indeterminate"

        #     elif all([isFloatValue(sample.nCoV_N1, None), isFloatValue(sample.nCoV_N2, None), isFloatValue(sample.nCoV_N3, None)]):
        #         sample.result = "Positive"
        if sample.result is None:
            print("------unable to identify result for the sample-----", sample)
            del sampleDict[sampleName]

    checkIncomingHl7(sampleDict)