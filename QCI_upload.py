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
print("\n", str(datetime.datetime.now()) + "\n")
#GATOR_SEQ_SAMPLE_INPUT_FILE = r'C:\Users\s.majety\Desktop\Copy of Sheet1.xlsx'
#LINUX_ANALYSIS_OUT_FOLDER = r'G:/DRL/Molecular/NGS/GenomOncology/NextSeq/'

script_path = os.path.dirname(os.path.abspath( __file__ ))
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
    


LINUX_ANALYSIS_OUT_FOLDER = replace_env(config_dict['LINUX_ANALYSIS_OUT_FOLDER'])
GATOR_SEQ_SAMPLE_INPUT_FILE = replace_env(config_dict['GATOR_SEQ_SAMPLE_INPUT_FILE'])
CONFIG_TOKENS_FILE = script_path + "/" + config_dict['CONFIG_TOKENS_FILE'] 

def check_folders_exist():
    if not os.path.isfile(GATOR_SEQ_SAMPLE_INPUT_FILE):
        sys.exit("ERROR: Does not have access to following folder: " + GATOR_SEQ_SAMPLE_INPUT_FILE + "\n") 

    if not os.path.isdir(LINUX_ANALYSIS_OUT_FOLDER):
        sys.exit("ERROR: Does not have access to following folder: " + LINUX_ANALYSIS_OUT_FOLDER + "\n")

    if not os.path.isfile(CONFIG_TOKENS_FILE):
        sys.exit("ERROR: Does not have access to following folder: " + CONFIG_TOKENS_FILE + "\n") 

check_folders_exist()

def save_workbook(df):
    try:
        df.to_excel(GATOR_SEQ_SAMPLE_INPUT_FILE, index=False)
    except:
        print("could not save excel")
        sys.exit()

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
def makePostCall(url, files, headers):
    response = requests.post(url, files=files, headers=headers)
    if response.json() and response.json().get("status-url") is not None:
        PROCESSING_STATUS_URLS.append(response.json()["status-url"])
    
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
def uploadToQiagen(zipFile):
    accessToken = generateAccessToken()
    url = QCI_UPLOAD_URL
    files = {'file': open(zipFile,'rb')}
    return makePostCall(url, files, {"Accept": "application/json", "Authorization": accessToken["access_token"]})
   
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
    map["Test_Product_Profile "] = row["Test_Product_Profile "]
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
    return response.json()["status"]

def populateQCIMessage(df, index, msg):
    df.at[index, "QCI_Message"] = msg

# main method
# 1. checks if excel is open and exits if it is open
# 2. reads from excel
# 3. gets all the statuses of accesion Id, used to upload only those accessionId which are not in Qiagen
# 4. for each row in excel
#   4.1 check if PLMO and other entries exists and then check if its not yet uploaded to Qiagen
#   4.2 upload to Qiagen by generating a XML file from row and a VCF file at the specified location in row.
# 5. closes the excel in the end 
if __name__ == "__main__":
    try: 
        excel_file = open(GATOR_SEQ_SAMPLE_INPUT_FILE, "r+")
    except:
        print("Could not open file! Please close Excel!")
        sys.exit(0)
    try:
        df_full = pd.read_excel(GATOR_SEQ_SAMPLE_INPUT_FILE)
        df = df_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    except:
        print("could not read excel")
        sys.exit()

    statusChanged = False
    accessionIdMap = populateStatusMap()


    for index, row in df.iterrows():
        if row["STATUS"] == "DONE" and type(row.get("PLMO_Number")) == str:#  math.isnan(float(row.get("PLMO_Number"))):
            vcfFolder = LINUX_ANALYSIS_OUT_FOLDER + "/" +  row['SAMPLE_DIR_PATH'].strip() + '_' + row['TIME_STAMP']  + "/"
            accessionId = row['SAMPLE_DIR_PATH'].split("/")[1].strip() + '_' + row['TIME_STAMP']
            vcfFileName = accessionId + ".vcf"
            if accessionIdMap.get(accessionId) is not None:
               # print(accessionId, " is already uploaded and hence not uploading again")
                continue

            if os.path.isfile(vcfFolder + vcfFileName):
                try: 
                    xmlFileName = generateXMLFileFromRow(row,  vcfFolder + accessionId + ".QCIUpload.xml")
                    with zipfile.ZipFile( vcfFolder + accessionId + ".QCIUpload.zip", 'w') as zip: 
                        zip.write(xmlFileName, accessionId + ".xml")
                        zip.write(vcfFolder + vcfFileName, accessionId + ".vcf")
                    try:
                        status_code = uploadToQiagen(vcfFolder + accessionId + ".QCIUpload.zip")
                        print("UPLOADED ", accessionId, " to QCI with status: ", status_code)
                        populateQCIMessage(df, index, "UPLOADED ", accessionId, " to QCI")
                        statusChanged = True
                    except:
                        print("error while uploading to Qiagen for accessionId: ", accessionId , " with exception: ", sys.exc_info()[0] )
                        populateQCIMessage(df, index, "error while uploading to Qiagen for accessionId: ", accessionId , " with exception: ", sys.exc_info()[0])
                except:
                    print("could not generate a zip file for accessionId: ", accessionId , " with exception: ", sys.exc_info()[0] )
                    populateQCIMessage(df, index, "could not generate a zip file for accessionId: ", accessionId , " with exception: ", sys.exc_info()[0])
            else:
                print("COULD NOT find vcf file: ",  vcfFolder + vcfFileName)
                populateQCIMessage(df, index, "COULD NOT find vcf file: ",  vcfFolder + vcfFileName)
    while len(PROCESSING_STATUS_URLS) > 0:
        time.sleep(60)
        for url in PROCESSING_STATUS_URLS: 
            status = checkStatus(url)
            if status != "PREPROCESSING":
                print("final status of Data packet with url: ", url, " is ", status)
                PROCESSING_STATUS_URLS.remove(url)

    excel_file.close()
    save_workbook(df)
   # logging.basicConfig(filename='uploadToQiagen.log',level=logging.ERROR)
