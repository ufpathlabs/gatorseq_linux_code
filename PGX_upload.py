import requests 
from requests.exceptions import HTTPError 
import xmlschema 
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
GATOR_SEQ_SAMPLE_INPUT_FILE = replace_env(config_dict['GATOR_SEQ_SAMPLE_INPUT_FILE'])
GATOR_PGX_SAMPLE_INPUT_FOLDER = replace_env(config_dict['GATOR_SEQ_SAMPLE_INPUT_FILE']) #'G:\DRL\Molecular\Assays\PGX\PGX_Beaker_Interface' 
CONFIG_TOKENS_FILE = script_path + "/" + config_dict['CONFIG_TOKENS_FILE']
MIRTH_GATORSEQ = config_dict['MIRTH_GATORSEQ']
if CODE_ENV=='DevEnv':
    MIRTH_GATORSEQ += '/TEST'
else:
    MIRTH_GATORSEQ += '/PROD'

def check_folders_exist():
    if not os.path.isfile(GATOR_SEQ_SAMPLE_INPUT_FILE):
        sys.exit("ERROR: Does not have access to following folder: " + GATOR_SEQ_SAMPLE_INPUT_FILE + "\n") 

    if not os.path.isdir(LINUX_ANALYSIS_OUT_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + LINUX_ANALYSIS_OUT_FOLDER + "\n")

    if not os.path.isfile(CONFIG_TOKENS_FILE):
        sys.exit("ERROR: Does not have access to following folder: " + CONFIG_TOKENS_FILE + "\n") 

    if not os.path.isdir(MIRTH_GATORSEQ):
        sys.exit("ERROR: Does not have access to following folder: " + MIRTH_GATORSEQ + "\n") 

    if not os.path.isdir(GATOR_PGX_SAMPLE_INPUT_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + GATOR_PGX_SAMPLE_INPUT_FOLDER + "\n") 


#check_folders_exist()
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
   
    def update_obx_seg_containing_gene(self, row):
        updates = 0
        temp_obx = self.h[:]
        l = len(self.h)
        for i in range(l):
            del temp_obx[l-i-1]
        new_obx_index = 1
        report_genotype_only_arr = ["CYP2C_CLUSTER", "CYP4F2", "VKORC1"]
        for obxSegment in self.h['OBX']:
            isGenoType = "GENOTYPE" in obxSegment[3][0][1][0]
            gene = obxSegment[3][0][1][0].replace(" GENOTYPE", "").replace(" PHENOTYPE", "")
            gene = gene.replace(' ', '_')
            if "CYP2C_" in gene:
                print(gene)
            #if we find the value for that gene in csv, add it to the OBX segment, else add it as is.
            if row.get(gene) is not None and (gene in report_genotype_only_arr and not isGenoType):
                print("totest")
            if row.get(gene) is not None and not (gene in report_genotype_only_arr and not isGenoType):
                try:
                    obxSegment[5][0] = row[gene].iloc[0].split("; ")[0] if isGenoType else row[gene].iloc[0].split("; ")[1]
                except:
                    print('------------------------------', row[gene], isGenoType, obxSegment[3][0][1][0])
            obxSegment[1] = new_obx_index
            new_obx_index +=1
            temp_obx.append(obxSegment)
            updates += 1

            
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

def main():
    UPLOAD_PATH = MIRTH_GATORSEQ + '/RESULTS'
    ORDERS_ARCHIVE_DIR = MIRTH_GATORSEQ + '/ORDERS_ARCHIVE/'
    ORDERS_DIR = MIRTH_GATORSEQ + '/ORDERS/'
    try: 
        #ToDo: check if they are sorted correctly

        os.chdir(GATOR_PGX_SAMPLE_INPUT_FOLDER)
        csv_files = filter(os.path.isfile, os.listdir(GATOR_PGX_SAMPLE_INPUT_FOLDER))
        csv_files = [os.path.join(GATOR_PGX_SAMPLE_INPUT_FOLDER, f) for f in csv_files] # add path to each file
        csv_files.sort(key=lambda x: os.path.getctime(x))

    except:
        print(" Could not read folder")
        sys.exit()

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
                #Read PLM id from HL7 message
                try:
                    messageId = str(h['OBR'][0][3]).replace("^", "")
                except:
                    continue
                if (not messageId):
                    continue

                # search for messageId in the csv_files of GATOR_PGX_SAMPLE_INPUT_FOLDER
                #if messageId == "100047187": #100047166  100047187
                print("Im here")
                foundMessageId = False
                for f in csv_files:
                    #First, find the header in the csv (using the constraint that header starts with 'sample ID')
                    # And then create a data frame using pandas
                    if foundMessageId:
                        continue
                    with open(f) as csvFile:
                        readCSV = csv.reader(csvFile, delimiter=',')
                        rowCount = -1
                        for row in readCSV:
                            rowCount += 1
                            if len(row) > 0 and row[0] == 'sample ID':
                                break 
                    df = pd.read_csv(f, header=rowCount)
                    df.columns = df.columns.str.strip().str.replace(' ', '_')
                    # not found in this file and hence moving on to next file
                    if df[df['sample_ID'] == messageId].empty:
                        print("sample_ID " + messageId +" not found in this file: ", f)
                    else:
                        #ToDO: check for duplicate entries of sample Id

                        # found messageId and generating a new hl7 file
                        print("processing ", len(df[df['sample_ID'] == messageId]), df[df['sample_ID'] == messageId])
                        foundMessageId = True
                        newHl7.update_msh_segment()
                        newHl7.update_orc_segment()
                        newHl7.update_obr_segment()
                        #newHl7.update_comments(open( r"C:\Users\s.majety\Desktop\PGX project\test.txt", mode="r",  encoding='utf-8').read())
                        newHl7.update_obx_segment()
                        h = newHl7.update_obx_seg_containing_gene( df[df['sample_ID'] == messageId])
                        
                        out_file_path = UPLOAD_PATH + '/hl7-{}-output.txt'.format(messageId)
                        if h:
                            with open(out_file_path, 'w' ,  encoding='utf-8') as f:
                                f.write(str(h))
                            print("Out file available at :",out_file_path)
                            move(ORDERS_DIR + hl7_file_name, ORDERS_ARCHIVE_DIR + 'processed-' + hl7_file_name) 
                            
main()
