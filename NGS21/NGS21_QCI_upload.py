import requests 
from requests.exceptions import HTTPError 
import ast
from pprint import pprint
import xml.etree.cElementTree as ET
import datetime
import numpy as np
import pandas as pd
import zipfile
import os
import logging
import logging.config
import sys
import yaml
import math
from shutil import copyfile
from shutil import move
import time
import datetime
import traceback
import sqlite3
import NGS21_database_connection

print("\n", str(datetime.datetime.now()) + "\n")
#GATOR_SEQ_SAMPLE_INPUT_FILE = r'C:\Users\s.majety\Desktop\Copy of Sheet1.xlsx'
#LINUX_ANALYSIS_OUT_FOLDER = r'G:/DRL/Molecular/NGS/GenomOncology/NextSeq/'

script_path = os.path.dirname(os.path.abspath( __file__ ))
script_path = os.path.abspath(os.path.join(script_path, '..'))
CONFIG_FILE=script_path+"/linux_gatorseq.config.yaml"
QCI_ACCESS_TOKEN_URL = "https://api.ingenuity.com/v1/oauth/access_token"
QCI_UPLOAD_URL = "https://api.ingenuity.com/v1/datapackages"
QCI_QUERY_URL = "https://api.ingenuity.com/datastream/api/v1/clinical?startReceivedDate=2018-04-28"
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
    


LINUX_ANALYSIS_OUT_FOLDER = replace_env(config_dict['LINUX_ANALYSIS_OUT_FOLDER_NGS21'])
GATOR_SEQ_SAMPLE_INPUT_FILE = replace_env(config_dict['GATOR_SEQ_SAMPLE_INPUT_FILE_NGS21'])
CONFIG_TOKENS_FILE = script_path + "/" + config_dict['CONFIG_TOKENS_FILE'] 

TABLE_NAME = replace_env(config_dict['TABLE_NAME'])
if CODE_ENV == "ProdEnv":
    TABLE_NAME = replace_env(config_dict['TABLE_NAME_PROD'])
def check_folders_exist():
    if not os.path.isfile(GATOR_SEQ_SAMPLE_INPUT_FILE):
        sys.exit("ERROR: Does not have access to following folder: " + GATOR_SEQ_SAMPLE_INPUT_FILE + "\n") 

    if not os.path.isdir(LINUX_ANALYSIS_OUT_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + LINUX_ANALYSIS_OUT_FOLDER + "\n")

    if not os.path.isfile(CONFIG_TOKENS_FILE):
        sys.exit("ERROR: Does not have access to following folder: " + CONFIG_TOKENS_FILE + "\n") 

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

PROCESSING_STATUS_URLS = []

# Helper function to make a HTTP call
def makePostCall(url, files, headers, SAMPLE_DIR_PATH):
    response = requests.post(url, files=files, headers=headers)
    if response.json() and response.json().get("status-url") is not None:
        PROCESSING_STATUS_URLS.append((response.json()["status-url"], SAMPLE_DIR_PATH))
    
    return response.json()["status"]

# Function to generate a new AccessToken. AccessToken is expired every 300 sec from Qiagen's end
def generateAccessToken():
    response = requests.post(QCI_ACCESS_TOKEN_URL, 
        {"grant_type": "client_credentials", 
        "client_id": QCI_CLIENT_ID, 
        "client_secret": QCI_CLIENT_ID_KEY}, 
        {"content-type": "text/plain;charset=UTF-8"} )

    return ast.literal_eval(response.content.decode("utf-8"))
    
# uploads the zip to Qiagen API
def uploadToQiagen(zipFile, SAMPLE_DIR_PATH):
    accessToken = generateAccessToken()
    url = QCI_UPLOAD_URL
    files = {'file': open(zipFile,'rb')}
    return makePostCall(url, files, {"Accept": "application/json", "Authorization": accessToken["access_token"]}, SAMPLE_DIR_PATH)
   
# Generates xml from a dictionary 
def createXMLFromDict(map, xmlFileName):
    if map["QCIType"] == "QCISomaticTest":
        xsiSchemaLocation = "http://qci.qiagen.com/xsd/interpret ../schema/QCISomaticTest.xsd "
    elif map["QCIType"] == "QCIHereditaryTest":
        xsiSchemaLocation = "http://qci.qiagen.com/xsd/interpret ../schema/QCIHereditaryTest.xsd "
    
    d = datetime.datetime.today()
    root = ET.Element(map.get("QCIType"), version="1.12.0", xmlns="http://qci.qiagen.com/xsd/interpret")
    root.attrib['xmlns:xsi'] = "http://www.w3.org/2001/XMLSchema-instance"
    root.attrib['xsi:schemaLocation'] = "http://qci.qiagen.com/xsd/interpret ../schema/QCISomaticTest.xsd "
    
    testProduct = ET.SubElement(root, "TestProduct")
    ET.SubElement(testProduct, "Code").text = map.get("Test_Product_Code")
    ET.SubElement(testProduct, "Profile").text = map.get("Test_Product_Profile ")
    ET.SubElement(testProduct, "ReportTemplate").text = map.get("Report_Template")
    
    
    if map.get("ReportingMethod"):
        ET.SubElement(testProduct, "ReportingMethod").text = map.get("ReportingMethod")
    
    if map.get("TreatmentsPolicy"):
        ET.SubElement(testProduct, "TreatmentsPolicy").text = map.get("TreatmentsPolicy")
    if map.get("Pre_Filter"):
        filter = ET.SubElement(testProduct, "Filters")
        ET.SubElement(filter, "Prefilter").text = map.get("Pre_Filter")
    
    test = ET.SubElement(root, "Test")
    ET.SubElement(test, "AccessionId").text = map.get("accessionId")
    ET.SubElement(test, "VariantsFilename").text = map.get("vcfFileName")
    ET.SubElement(test, "TestDate").text = map.get("TestDate") or d.strftime('%Y-%m-%d')
    ET.SubElement(test, "Diagnosis").text = map.get("Diagnosis")
    ET.SubElement(test, "PrimarySourceTissue").text = map.get("Primary_Tumor_Site")
       
    reviewers = ET.SubElement(root, "Reviewers")
    ET.SubElement(reviewers, "Reviewer").text = "University_of_Florida_TPP"

    tree = ET.ElementTree(root)
    tree.write(xmlFileName)
    return xmlFileName

# takes in a row from the excel sheet and generates a XML file filling in all the details
def generateXMLFileFromRow(row, xmlFilename):
    map = {}
    map["Diagnosis"] = row["Diagnosis"]
    map["Test_Product_Code"] = row['Test_Product_Code']
    map["Primary_Tumor_Site"] = row['Primary_Tumor_Site']
    map["Report_Template"] = row['Report_Template']
    map["accessionId"] =  row['SAMPLE_DIR_PATH'].split("/")[1].strip() + '_' + row['TIME_STAMP']
    map["vcfFileName"] = row['SAMPLE_DIR_PATH'].split("/")[1].strip() + '_' + row['TIME_STAMP'] + '.vcf'
    if row.get("QCIType") is not None:
        map["QCIType"] = (row['QCIType'])
    else:
        map["QCIType"] = "QCISomaticTest"
    map["Pre_Filter"] = row["Pre_Filter"]
    map["TreatmentsPolicy"] = row["Treatments_Policy"]
    map["ReportingMethod"] = row["Reporting_Method"]
    map["Test_Product_Profile "] = row["Test_Product_Profile"]
    xmlFileName = createXMLFromDict(map, xmlFilename)
    return xmlFileName
    
# Gets the status of all the accessionIds. It is used because checkStatus() is taking awfully lomg time to query by accessionIds
def getAllStatus():
    url = QCI_QUERY_URL
    headers = {"Accept": "application/json", "Authorization": generateAccessToken()["access_token"]}
    response = requests.get(url, headers = headers)
    return ast.literal_eval(response.content.decode("utf-8"))

def populateStatusMap():
    array = getAllStatus()
    map = {}
    for x in array:
        map[x.get("accessionID")] = x.get("state")
    return map

def checkStatus(url):
    headers = {"Accept": "application/json", "Authorization": generateAccessToken()["access_token"]}
    response = requests.get(url, headers = headers) 
    return (response.json()["status"], response.json())

def create_connection():
    return NGS21_database_connection.getSQLConnection(CONFIG_FILE, CONFIG_TOKENS_FILE, CODE_ENV)

def updateStatus(SAMPLE_DIR_PATH, message, con):
    cursor = con.cursor()
    sql_update_query = 'Update '+ TABLE_NAME +'  set QCI_Upload_Message = %s, QCI_Re_Run = "" where SAMPLE_DIR_PATH = %s ;'
    cursor.execute(sql_update_query, (message, SAMPLE_DIR_PATH))
    con.commit()
    cursor.close()

# main method
# 1. checks if excel is open and exits if it is open
# 2. reads from excel
# 3. gets all the statuses of accesion Id, used to upload only those accessionId which are not in Qiagen
# 4. for each row in excel
#   4.1 check if PLMO and other entries exists and then check if its not yet uploaded to Qiagen
#   4.2 upload to Qiagen by generating a XML file from row and a VCF file at the specified location in row.
# 5. closes the excel in the end 
if __name__ == "__main__":
    # try: 
    #     excel_file = open(GATOR_SEQ_SAMPLE_INPUT_FILE, "r+")
    # except:
    #     print("Could not open file! Please close Excel!")
    #     sys.exit(0)
    # try:
    #     df_full = pd.read_excel(GATOR_SEQ_SAMPLE_INPUT_FILE)
    #     df = df_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    # except:
    #     print("could not read excel")
    #     sys.exit()

    df = pd.read_sql_query('select * from '+ TABLE_NAME +' where status =  "DONE" and PLMO_Number != "" and (QCI_Upload_Message = "" or QCI_Re_Run = "yes");', create_connection())
    conn = create_connection()
    statusChanged = False
    accessionIdMap = populateStatusMap()

    for index, row in df.iterrows():
        if row["STATUS"] == "DONE" and type(row.get("PLMO_Number")) == str:#  math.isnan(float(row.get("PLMO_Number"))):
            assay_folder = row['ASSAY_DIR'].strip().split('-')[0]
            vcfFolder = LINUX_ANALYSIS_OUT_FOLDER + "/" +  assay_folder + "/" + row['ASSAY_DIR'].strip() + "/" + row['SAMPLE_DIR_PATH'].strip() + '_' + row['TIME_STAMP']  + "/"
            accessionId = row['SAMPLE_DIR_PATH'].split("/")[1].strip() + '_' + row['TIME_STAMP']
            vcfFileName = accessionId + ".vcf"
            if row['QCI_Re_Run'].lower() != "yes" and accessionIdMap.get(accessionId) is not None:
                print(accessionId, " is already uploaded and hence not uploading again")
                updateStatus(row['SAMPLE_DIR_PATH'], "Successfully Uploaded", conn)
                continue
            
            if row['QCI_Re_Run'].lower() == "yes":
                try:
                    os.remove(vcfFolder + accessionId + ".QCIUpload.xml")
                    os.remove(vcfFolder + accessionId + ".QCIUpload.zip")
                except OSError:
                    pass
                print("files removed, if exists as rerunning QCI_upload")


            if os.path.isfile(vcfFolder + vcfFileName):
                
                xmlFileName = generateXMLFileFromRow(row,  vcfFolder + accessionId + ".QCIUpload.xml")
                try:
                    with zipfile.ZipFile( vcfFolder + accessionId + ".QCIUpload.zip", 'w') as zip: 
                        zip.write(xmlFileName, accessionId + ".xml")
                        zip.write(vcfFolder + vcfFileName, accessionId + ".vcf")
                    try:
                        status_code = uploadToQiagen(vcfFolder + accessionId + ".QCIUpload.zip", row['SAMPLE_DIR_PATH'].strip())
                        print("UPLOADED ", accessionId, " to QCI with status: ", status_code)
                        statusChanged = True
                    except:
                        print("error while uploading to Qiagen for accessionId: ", accessionId , " with exception: ", sys.exc_info()[0] )

                except:
                    print("could not generate a zip file for accessionId: ", accessionId , " with exception: ", sys.exc_info()[0] )

            else:
                print("COULD NOT find vcf file: ",  vcfFolder + vcfFileName)
    

    while len(PROCESSING_STATUS_URLS) > 0:
        time.sleep(30)
        for url, SAMPLE_DIR_PATH in PROCESSING_STATUS_URLS: 
            status, response = checkStatus(url)
            if status == "DONE" or status == "FAILED":
                print("final status of Data packet with url: ", url, " is ", status)
                PROCESSING_STATUS_URLS.remove((url, SAMPLE_DIR_PATH))
                if status == "DONE":
                    updateStatus(SAMPLE_DIR_PATH, "Successfully Uploaded", conn)
                else:
                    print(str(response.get("errors")))
                    updateStatus(SAMPLE_DIR_PATH, "Error while Uploading: " + str(response.get("errors"))[:200], conn)
            else:
                print("retrying as status is not yet done or failed and url: ", url, " is ", status)

    conn.close()
    # excel_file.close()
   # logging.basicConfig(filename='uploadToQiagen.log',level=logging.ERROR)
