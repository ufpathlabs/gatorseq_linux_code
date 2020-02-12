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
import xmltodict

import datetime
import traceback
import sqlite3
import mysql.connector

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


LINUX_ANALYSIS_OUT_FOLDER = replace_env(config_dict['LINUX_ANALYSIS_OUT_FOLDER'])
GATOR_SEQ_SAMPLE_INPUT_FILE = replace_env(config_dict['GATOR_SEQ_SAMPLE_INPUT_FILE'])
CONFIG_TOKENS_FILE = script_path + "/" + config_dict['CONFIG_TOKENS_FILE']
MIRTH_GATORSEQ = config_dict['MIRTH_GATORSEQ']

TABLE_NAME = replace_env(config_dict['TABLE_NAME'])
SQLITE_DB = replace_env(config_dict['SQLITE_DB'])

MYSQL_HOST = config_dict['MYSQL_HOST']
MYSQL_USERNAME = config_dict['MYSQL_USERNAME']
# MYSQL_PASSWAORD = config_dict['MYSQL_PASSWAORD']
MYSQL_DATABASE = config_dict['MYSQL_DATABASE']

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


check_folders_exist()

config_token_dict=dict()
with open(CONFIG_TOKENS_FILE, 'r') as stream:
    try:
        config_token_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

QCI_CLIENT_ID = config_token_dict['QCI_CLIENT_ID']
QCI_CLIENT_ID_KEY = config_token_dict['QCI_CLIENT_ID_KEY']
MYSQL_PASSWAORD = config_token_dict['MYSQL_PASSWAORD']


#populates a map with each drug name. it is required as we need to show the treatments grouped by drugnames
def getDrugMaps(treatmentsList):
    drugMap = {}
    if type(treatmentsList) != list:
        treatmentsList = [treatmentsList]
    for treatment in treatmentsList:
        if drugMap.get(treatment.get("drug").get("drugname")) is None:
            isMatch = treatment["diagnosismatch"] == 'true' or treatment["tissuematch"]  == 'true' or False
            drugMap[treatment.get("drug").get("drugname")] = {'isMatch': isMatch, 'listTreatments':[]}
        else:
            isMatch = treatment["diagnosismatch"] == 'true' or treatment["tissuematch"]  == 'true' or False
        drugMap[treatment.get("drug").get("drugname")]['isMatch'] = drugMap[treatment.get("drug").get("drugname")]['isMatch'] or isMatch
        drugMap[treatment.get("drug").get("drugname")]['listTreatments'].append(treatment)
    return drugMap

def create_connection(db_file):
    conn = None
    # try:
    #     conn = sqlite3.connect(db_file)
    # except:
    #     print(traceback.format_exc())

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

def updateStatus(SAMPLE_DIR_PATH, message, con):
    cursor = con.cursor()
    sql_update_query = 'Update '+ TABLE_NAME +'  set EPIC_Upload_Message = %s, EPIC_Re_Run = "" where SAMPLE_DIR_PATH = %s ;'
    cursor.execute(sql_update_query, (message, SAMPLE_DIR_PATH))
    con.commit()
    cursor.close()

# 1. Tries to open excel and exits if already open
# 2. Iterates over the folder and tries to read PLMO number from each hl7 file
# 3. Checks if there is an entry in excel for that PLMO and retrieves the corresponding accession id.
# 4. queries the Qiagen to see if the accessionid is in final state
# 5. downloads the xml, parses it and generates comments
# 5. processes the genes and moves the new hl7 file to a results folder where it gets pushed to EPIC
# 6. archives the initial file 
def main():
        
    UPLOAD_PATH = MIRTH_GATORSEQ + '/RESULTS'
    ORDERS_ARCHIVE_DIR = MIRTH_GATORSEQ + '/ORDERS_ARCHIVE/'
    ORDERS_DIR = MIRTH_GATORSEQ + '/ORDERS/'
    #Check if excel file is opened by any other user
    # try: 
    #     excel_file = open(GATOR_SEQ_SAMPLE_INPUT_FILE, "r+")
    # except:
    #     print(" Could not open file! Please close Excel!")
    #     sys.exit()
    # try:
    #     xldf_full = pd.read_excel(GATOR_SEQ_SAMPLE_INPUT_FILE)
    #     xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    # except:
    #     print("Problem Reading Excel")
    #     sys.exit()

    xldf = pd.read_sql_query('select * from '+ TABLE_NAME +' where status =  "DONE" and PLMO_Number != "" and (EPIC_Upload_Message = "" or EPIC_Re_Run = "yes") ;', create_connection(SQLITE_DB))
    conn = create_connection(SQLITE_DB)

    #re_run logic
    # if a row is to be rerun, move the input hl7 file from ORDERS_ARCHIVE folder to ORDERS folder
    # the rows to re-run is a rare case and hence, each time, we are querying all files in archive (will not be a time complexity issue)
    rerun_barcodes = pd.read_sql_query('select * from '+ TABLE_NAME +' where status =  "DONE" and PLMO_Number != "" and (EPIC_Re_Run = "yes") ;', create_connection(SQLITE_DB))
    for index, row in rerun_barcodes.iterrows():
        excel_row_plmo = row['PLMO_Number']
        for (dirpath, dirnames, filenames) in os.walk(ORDERS_ARCHIVE_DIR):
            for fileName in filenames:
                try:
                    hl7file = open(ORDERS_ARCHIVE_DIR + fileName, mode="r").read()
                except:
                    continue
                arr = hl7file.split("\n\n") #split by blank lines
                c=0
                for hl7msg in arr: 
                    if hl7msg: 
                        msg_unix_fmt = hl7msg.replace("\n","\r")
                        h = hl7.parse(msg_unix_fmt)
                        # Read PLM id from HL7 message and then compare it with the plmo in reRun_excel
                        try:
                            plm = h['ORC'][0][2]
                        except:
                            continue
                        if (plm) and str(plm) == excel_row_plmo:
                            try:
                                os.rename(ORDERS_ARCHIVE_DIR + fileName, ORDERS_DIR + fileName)
                            except:
                                print("cannot move file from archive to orders during re-run case")

    #---------- additional step removing initially created hl7.txt file is below ------------
            

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
        #move('/ext/mirth_gatorseq/MIRTH_GATORSEQ/TEST/ORDERS/61-4caadebf-812c-4b31-b30e-a95673c85e5f.txt', '/ext/mirth_gatorseq/MIRTH_GATORSEQ/TEST/ORDERS_ARCHIVE/61-4caadebf-812c-4b31-b30e-a95673c85e5f.txt')
        c=0
        for hl7msg in arr: 
            if hl7msg: #check if message not empty
                msg_unix_fmt = hl7msg.replace("\n","\r")
                h = hl7.parse(msg_unix_fmt)
                #Read PLM id from HL7 message
                try:
                    plm = h['ORC'][0][2]
                except:
                    continue
                if (not plm):
                    print('PLM not found in HL7 message!')
                    continue
                if xldf[xldf['PLMO_Number'] == str(plm)].empty:
                    # print('PLMO '+ str(plm) +' is not found in excel')
                    continue
                if len(xldf[xldf['PLMO_Number'] == str(plm)]) > 1:
                    print(str(plm), " has duplicate entries in excel")
                    continue
                #print((plm))
                sample_dir_path = xldf[xldf['PLMO_Number'] == str(plm)]['SAMPLE_DIR_PATH'].item()
                if (not sample_dir_path):
                    continue
                #if xldf[xldf['PLMO_Number'] == str(plm)]['downloadedXML'].item() == 0:
                accessionId = sample_dir_path.split("/")[1] + "_" + xldf[xldf['PLMO_Number'] == str(plm)]['TIME_STAMP'].item()
                vcfFolder = LINUX_ANALYSIS_OUT_FOLDER + "/" +  xldf[xldf['PLMO_Number'] == str(plm)]['SAMPLE_DIR_PATH'].item() + '_' + xldf[xldf['PLMO_Number'] == str(plm)]['TIME_STAMP'].item()  + "/"
                Perc_Target_Cells =  xldf[xldf['PLMO_Number'] == str(plm)]['Perc_Target_Cells'].item()
                if Perc_Target_Cells == "":
                    Perc_Target_Cells = None

                Perc_Tumor =  xldf[xldf['PLMO_Number'] == str(plm)]['Perc_Tumor'].item()
                if Perc_Tumor == "":
                    Perc_Tumor = None
                
                force_re_run = True if xldf[xldf['PLMO_Number'] == str(plm)]['EPIC_Re_Run'].item() == "yes" else False
                if force_re_run:
                    try:
                        os.remove(vcfFolder+accessionId+".hl7.txt")
                    except:
                        pass
                    try:
                        os.remove(vcfFolder+accessionId+".QCIXml.xml")
                    except:
                        pass
                    print(sample_dir_path)
                    cursor = conn.cursor()
                    sql_update_query = 'Update '+ TABLE_NAME +'  set QCI_Download_Message = "", EPIC_Re_Run = "FETCHING_REPORT"  where SAMPLE_DIR_PATH ="'+sample_dir_path+'" ;'
                    cursor.execute(sql_update_query)
                    print(sql_update_query)
                    conn.commit()
                    cursor.close()
                
                
                if not os.path.isfile(vcfFolder+accessionId+".hl7.txt") and os.path.isfile(vcfFolder+accessionId+".QCIXml.xml"):  #accessionIdStatusMap.get(accessionId) is not None:
                    if os.path.isfile(vcfFolder+accessionId+".QCIXml.xml") and os.path.isfile(vcfFolder+accessionId+".QCIreport.txt"):
                        genes_list, diagnosis = hl7update.find_genes_from_XML(vcfFolder+accessionId+".QCIXml.xml")
                        gene_map={}
                        if(genes_list):
                            for gene in genes_list:
                                x = gene.split(" ", 1)
                                if x[0] in gene_map.keys():
                                    gene_map[x[0]].append(x[1])
                                else:
                                    gene_map[x[0]] = [x[1]]
                            #gene_map = dict(gene.split(" ", 1) for gene in genes_list)
                        hl7update.update_msh_segment(h)
                        hl7update.update_orc_segment(h)
                        hl7update.update_obr_segment(h)
                        hl7update.update_comments(h, open( vcfFolder+accessionId+".QCIreport.txt", mode="r",  encoding='utf-8').read())
                        hl7update.update_obx_segment(h)
                        h = hl7update.update_obx_seg_containing_gene(h, gene_map, accessionId, diagnosis, Perc_Target_Cells, Perc_Tumor)
                        
                        out_file_path = UPLOAD_PATH + '/hl7-{}-output.txt'.format(plm)
                        if h:
                            with open(out_file_path, 'w' ,  encoding='utf-8') as f:
                                f.write(str(h))
                            print("Out file available at :",out_file_path)
                            move(ORDERS_DIR + hl7_file_name, ORDERS_ARCHIVE_DIR + 'processed-' + hl7_file_name) 
                            copyfile(out_file_path, vcfFolder+accessionId+".hl7.txt")
                            updateStatus(sample_dir_path, "Successfully added hl7 file at " + time.ctime(), conn) 
                        else:
                            print("Couldn't replace '-' in hl7. Check logs for more details!")
                            updateStatus(sample_dir_path, "Error in processing hl7 file", conn) 
                    else:
                        print("XML was not yet generated for the  " + accessionId)

    conn.close()
    #logging.debug('=======================Execution ends===========================')
    # excel_file.close()

#for handler in logging.root.handlers[:]:
#   logging.root.removeHandler(handler)
#logging.basicConfig(format='%(asctime)s - %(message)s', filename='loggingForGetXML.log', level=logging.DEBUG)

main()


#accessionId = "NS-19-13_BC710502_53_20190416154205825653DevEnv3.0"
#accessToken = generateAccessToken()
#y, status = getXMLForAccessionId(accessionId, accessToken["access_token"])
#if y and status == 200:
#    print(parseXML(y, accessionId, 'fdadfs'))

