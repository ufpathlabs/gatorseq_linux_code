import requests 
from requests.exceptions import HTTPError 
import ast
import sys
import yaml
import numpy
import pandas as pd
import os
import datetime
from shutil import copyfile
from shutil import move
import xmltodict
import time
import datetime
print(str(datetime.datetime.now()) + "\n")

script_path = os.path.dirname(os.path.abspath( __file__ ))
CONFIG_FILE=script_path+"/linux_gatorseq.config.yaml"
config_dict=dict()
QCI_ACCESS_TOKEN_URL = "https://api.ingenuity.com/v1/oauth/access_token"
QCI_DOWNLOAD_URL =  "https://api.ingenuity.com/v1/export/"
QCI_QUERY_FINAL_URL =  "https://api.ingenuity.com/datastream/api/v1/clinical?state=final&startReceivedDate=2018-04-28"
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

config_token_dict=dict()
with open(CONFIG_TOKENS_FILE, 'r') as stream:
    try:
        config_token_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

QCI_CLIENT_ID = config_token_dict['QCI_CLIENT_ID']
QCI_CLIENT_ID_KEY = config_token_dict['QCI_CLIENT_ID_KEY']


# Gets all the accessionIds with 'final' status. Used to check if a accessionId is ready to be pulled from Qiagen
def getAllStatus():
    url = QCI_QUERY_FINAL_URL
    headers = {"Accept": "application/json", "Authorization": generateAccessToken()["access_token"]}
    response = requests.get(url, headers = headers)
    return ast.literal_eval(response.content.decode("utf-8"))

def populateStatusMap():
    array = getAllStatus()
    map = {}
    for x in array:
        map[x.get("accessionID")] = x.get("state")
    return map

# helper call to make a GET HTTP call
def makeGetCall(url, params, headers):
    response = requests.get(url, params=params, headers=headers)
    return response.content.decode("utf-8"), response.status_code 

def generateAccessToken():
    response = requests.post(QCI_ACCESS_TOKEN_URL, {"grant_type": "client_credentials", "client_id": QCI_CLIENT_ID, "client_secret": QCI_CLIENT_ID_KEY}, {"content-type": "text/plain;charset=UTF-8"} )
    return ast.literal_eval(response.content.decode("utf-8"))

# downloads the xml for a given accessionId
def getXMLForAccessionId(accessionId, accessToken):
    url = QCI_DOWNLOAD_URL
    url = url + accessionId
    return makeGetCall(url, {"view" : "reportXml"}, {"Accept": "application/xml", "Authorization": accessToken})

# compares the xml with the xsd file and genreates a dictionary
def parseXML(xml, accessionId, plm, accessionIdPath):
    filename = accessionIdPath +  ".QCIXml.xml"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w", encoding='utf-8') as f:
        f.write(xml)
        f.close()
    with open(filename, encoding='utf-8') as fd:
        map = xmltodict.parse(fd.read())
    return generateTxtFileAndSave(map['report'], accessionId, plm, accessionIdPath)

def newLine():
    return "  \n" 

# helper function to generate comments
def populateIndividualText2(text, list, gene_map):
    foundAtleastOne = False
    for iterator in range(0, len(list)):
        variant = list[iterator]
        
        foundAtleastOne = True
        text += "    " +  str(variant["gene"]) + " " 
        if variant.get("transcriptchange"):
            text += variant.get("transcriptchange").get("change") + " "
        if variant.get("proteinchange"):
            text += str(variant.get("proteinchange").get("change")) + " "
        if variant.get("allelefraction"):
            text += "VAF: " + str(variant.get("allelefraction")) + "%"
        text += "    "
        
        if variant.get("chromosome"):
            text += "chr"+variant.get("chromosome") + ":"
        if variant.get("genomicchange"):
            text += str(variant.get("genomicchange").get("change")) + "  "
        if variant.get("actionability") is not None and variant["assessment"] != "Uncertain Significance":
            text += "(Tier " + variant.get("actionability") + ")"
        text += " \n"
         
        gene_map[str(variant["gene"])] = variant
    if not foundAtleastOne:
        text += "    None" + " \n"
    return text, gene_map



def addHeading(text, heading):
    text += "------------------------------------------------------------" + "\n"
    text += heading + " \n"
    text += "------------------------------------------------------------" + "\n"
    return text

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

#this method is used to generate comments for matched diagnosis and as well as not matched diagnosis
def generateCommentsForDrug(text, drugMap, gene_map, isMatch):
    foundAtleastOne = False
    for drugname in drugMap:
        drug = drugMap[drugname]
        if drug['isMatch'] == isMatch:
            foundAtleastOne = True
            for treatment in drug['listTreatments']:
                text += "    " + treatment.get("gene") + " "
                variant = gene_map.get(treatment.get("gene"))
                if treatment.get("proteinchange"):
                    text += str(treatment.get("proteinchange").get("change")) + " "
                if variant.get("actionability"):
                    text += "Tier " + variant.get("actionability") + " "
                if variant.get("assessment"):
                    text += variant.get("assessment") + " "
                text += " \n"
            treatment = drug['listTreatments'][0]
            if treatment.get("type") is not None:
                text += "    Response=" 
                if  treatment.get("type") == "candidate":
                    text += "Sensitive" + "\n"
                elif treatment.get("type") == "caution":
                    text += "Resistant" + "\n"
                elif treatment.get("type") == "predictive":
                    text += "predictive" + "\n"
                elif treatment.get("type") == "not recommended":
                    text += "not recommended" + "\n"
            if treatment.get("drug") is not None:
                text += "    Drug: " + treatment.get("drug").get("drugname") + "\n"
            if treatment.get("rcomment"):# and treatment.get("rcomment")[0]:
                # text += "    Therapies description: " + treatment.get("rcomment")[0].get("text") + "\n"
                text += "    Therapies description: " + treatment.get("rcomment").get("text") + "\n"
            text += " \n"
    if not foundAtleastOne:
        text += "    None" + "\n"
    return text

# generate comments from the dictionary created from downloaded xml   
def generateTxtFileAndSave(map, accessionId, plm, accessionIdPath):
    text = ''
    gene_map = {}
    custom_order = ["Pathogenic", "Likely Pathogenic", "Uncertain Significance"]
    order = {key: i for i, key in enumerate(custom_order)}
    variants = None
    if map.get("variant"):
        if type(map.get("variant")) != list:
            map["variant"] = [map.get("variant")]
        variants = sorted(map.get("variant"), key=lambda d: order[d['assessment']])
    
    text += "UF HEALTH PATHLABS GATORSEQ177 NGS SCREENING PANEL" + newLine()
    text += "Sample: " + plm + "\n"
    if map.get("accession"):
        text += "Accession Id: " + map.get("accession") + "\n"
    if map.get("diagnosis"):
        text += "Diagnosis: " + map.get("diagnosis") + "\n"
        
    
    pVar, likVar, uncVar = [], [], []
    if variants is not None:
        for variant in variants:
            if variant["assessment"] == "Pathogenic":
                pVar.append(variant)
            elif variant["assessment"] == "Likely Pathogenic":
                likVar.append(variant)
            elif variant["assessment"] == "Uncertain Significance":
                uncVar.append(variant)
	
    text += newLine()
    text = addHeading(text, "1. Pathogenic variants")
    if variants is not None:
        text, gene_map = populateIndividualText2(text, pVar, gene_map)
    else:
        text += "    None" + " \n"
    text += " \n"
    text = addHeading(text, "2. Likely pathogenic variants")
    if variants is not None:
        text, gene_map = populateIndividualText2(text, likVar, gene_map)
    else:
        text += "    None" + " \n"
    text += " \n"
    text = addHeading(text, "3. Variants of unknown significance (Tier 3)")
    if variants is not None:
        text, gene_map = populateIndividualText2(text, uncVar, gene_map)
    else:
        text += "    None" + " \n"

    text += newLine()
    drugMap = {}
    if map.get("treatment") is not None:
        drugMap = getDrugMaps(map.get("treatment"))
    text = addHeading(text, "4. Therapeutic implications for " + map.get("diagnosis")) 
    try:

        text = generateCommentsForDrug(text, drugMap, gene_map, 1)
        text += newLine()

        text = addHeading(text, "5. Therapeutic implications for other indications") 
        text = generateCommentsForDrug(text, drugMap, gene_map, 0)
    except:
        print("error while getting all treatments")
    text += newLine()

    text = addHeading(text, "6. Individual Variant Interpretations")            
    if variants is not None:
        for variant in variants:
            if variant["assessment"] == "Uncertain Significance":
                continue
            text += "Gene: " + str(variant["gene"]) + "\n"
            if variant.get("transcriptchange"):
                text += "    " + "Exon: " + str(variant.get("transcriptchange").get("exonNumber")) +  "; " + "Nucleotide: " + str(variant.get("transcriptchange").get("transcript")) + ":" + variant.get("transcriptchange").get("change") + ";  "
            
            if variant.get("genomicchange"):
                text += " " + "Genomic Location: " + str(variant.get("genomicchange").get("change")) + ";"
                
            if variant.get("proteinchange"):
                text += " " + "Amino acid: " + str(variant.get("proteinchange").get("change")) + ";"
                
            if variant.get("function"):
                text += " " + "Function: " + str(variant.get("function")) + ";"
                
            if variant.get("assessment"):
                text += " " + "Assessment: " + str(variant.get("assessment")) + ";"
            if variant.get("actionability"):
                text += " " + "Classification: Tier " + str(variant.get("actionability")) + ";"
                
            if variant.get("allelefraction"):
                text += " " + "Allele Fraction: " + str(variant.get("allelefraction")) + "%(of "+ str(variant.get("readdepth")) +" reads)" + ";"
                   
            if variant.get("variation"):
                text += " " + "Variation: " + str(variant.get("variation")) + "\n"
            if variant.get("rcomment"):# and variant.get("rcomment")[0]:
                #text += "    " + "Interpretation: " + str(variant.get("rcomment")[0].get("text")) + "\n"
                text += "   " + "Interpretation: " + str(variant.get("rcomment").get("text")) + "\n"
            text += " \n"
    else:
        text += "    None"  + " \n"  
    text += newLine()

    text = addHeading(text, "7. Methods")       
    text += """The area of tumor is localized on an H/E slide and microdissected. Genomic DNA extracted from the tissue is amplified using the GatorSeq NGS Panel and sequenced on the Illumina NextSeq to high uniform depth (targeting 500x coverage by non-PCR duplicate read pairs with >99% of exons at coverage >100x). Sequence data is processed using a customized analysis pipeline (GatorSeq v3.0) designed to accurately detect base substitutions and insertions/deletions. Human genome version hg19 is used as the reference and was downloaded from http://ftp.broadinstitute.org/bundle/2.8/hg19. The mutation nomenclature is based on the convention recommended by the Human Genome Variation Society (http://www.hgvs.org/mutnomen/). This assay will only detect alterations in the exons of genes listed in Genes assayed. The analytical sensitivity limit of the assay is approximately 5% at 500x coverage. With lower coverage, the sensitivity worsens especially in tissues with low tumor cell percentage, and mutations present below the level of detection of this assay will not be identified. This assay is intended to detect somatic alterations present in tumor cells and is not to be used to determine or report germline alterations. 

Annotated reports are generated using QIAGEN Clinical Insight (QCI) software and databases*. QIAGEN Clinical Insight (QCI) is a variant analysis, interpretation and decision support tool for research and clinical labs analyzing human genetics data and is not intended to be used for diagnostic purposes. QCI uses rules-based approaches to automatically compute pathogenicity classifications and actionability classifications (Tier 1 to 3) for each alteration according to the 2015 professional guidelines from the American College of Medical Genetics and Association for Molecular Pathology (ACMG/AMP) and the 2017 guideline from the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists (AMP/ASCO/CAP). 

Pathogenicity/Actionability Classification
Tier 1 Variants: Variants with Strong Clinical Significance
Tier 2 Variants: Variants with Potential Clinical Significance
Tier 3 Variants: Variants of Unknown Significance.

The QIAGEN treatment reporting policy will report approved treatments associated with 1A/B or 2C/D actionability tiers based on evidence from FDA, NCCN, and primary literature. For LOF alterations that fall into oncogenes only exact variant matches to treatments associated with 1A/B or 2C/D will be reported. Treatments associated with 2C ("off-label treatments") are only reported if 1A/B evidence in the patient's diagnosis is not available. If contradictory evidence exists for a drug (sensitive and resistant information) the treatment with the highest level of evidence (the highest tier) will be reported out (e.g. 1A > 1B). If the same level of evidence (e.g. 1B resistance and 1B sensitivity) for the same drug is present, sensitivity trumps resistance.
"""
    text += map.get("version") + " \n"
    text += newLine()

    text = addHeading(text, "8. Genes Tested")    
    text += str(map.get("labTestedGenes")).replace(";", ", ") + " \n"
    
    text += newLine()

    text = addHeading(text, "9. Appendix") 
    text += """For all non-FDA approved tests listed above unless specified above, UF Health Pathology Laboratories (UFHPL) test have been developed and their performance characteristics determined by UFHPL. These tests have not been cleared or approved by the U.S. Food and Drug Administration. The FDA has determined that such clearance or approval is not necessary. Tests are used for clinical purposes andshould not be regarded as investigational or for research purposes. UFHPL is certified under the Clinical Laboratory Improvement Amendments of 1988 (CLIA) as qualified to perform high complexity clinical testing. Decisions on patient care and treatment must be based on the independent medical judgment of the treating physician, taking into considerations all applicable information concerning the patient's condition including family history, physical examination, patient preferences and other diagnostic tests in accordance with the applicable standard of care. Treatment decisions should not be made based solely on the information contained in this report. The presence of a genomic alteration that has been reported to show response to a particular therapy is no guarantee that the agent will show efficacy in any given patient. Any reference to available clinical trials or treatment options is based on the currently available public knowledge aggregated by publicly available websites which we do not provide or update content nor are we responsible for monitoring the update of content. Treating physician certifies by receiving this report that the patient has been adequately counseled about the results and implications of this genetic test in accordance with all applicable laws."""
    
    # text += newLine()
    # text += "Selected Citations" + "\n"
    # if map.get("article") is not None:
    #     for article in map.get("article"):
    #         text += "    " +  str(article.get("@author")) + " " + str(article.get("@title")) + " " + str(article.get("@citation")) + " " + str(article.get("@date")) + " " + str(article.get("@url")) + "\n\n" 

    filename = accessionIdPath + ".QCIreport.txt"#os.path.join(os.path.dirname(__file__) + "data/txt/", "txt-" + accessionId +".txt")
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w+", encoding='utf-8') as f:
        f.write(text)
        f.close()
    #filename = r"E:\Sid Workspace\projects\hl7-parse-update\HL7_Python_parsing_script\txt\txt-" + accessionId +".txt"
    #file = open(filename, 'w+', encoding='utf-8')
    #file.write(text)
    #file.close()
    return filename

# move the processed file to archive
def moveHL7ToArchive(src, dest):
    # os.remove(dest)
    move(src, dest)

# calls the Qiagen API for the report
def callQCIApi(accessionId, plm, accessionIdPath):
    accessToken = generateAccessToken()
    accessionId = accessionId or "NS-19-11_BC707507_92_20190319110904776453DevEnv3.0"
    y, status = getXMLForAccessionId(accessionId, accessToken["access_token"])
    if y and status == 200:
        return parseXML(y, accessionId, plm, accessionIdPath)
    return False

# 1. Tries to open excel and exits if already open
# 2. Iterates over the folder and tries to read PLMO number from each hl7 file
# 3. Checks if there is an entry in excel for that PLMO and retrieves the corresponding accession id.
# 4. queries the Qiagen to see if the accessionid is in final state
# 5. downloads the xml, parses it and generates comments
# 5. processes the genes and moves the new hl7 file to a results folder where it gets pushed to EPIC
# 6. archives the initial file 
def main():
    #Check if excel file is opened by any other user
    try: 
        excel_file = open(GATOR_SEQ_SAMPLE_INPUT_FILE, "r+")
    except:
        print(" Could not open file! Please close Excel!")
        sys.exit()
    try:
        xldf_full = pd.read_excel(GATOR_SEQ_SAMPLE_INPUT_FILE)
        xldf = xldf_full.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    except:
        print("Problem Reading Excel")
        sys.exit()


    accessionIdStatusMap = populateStatusMap()

    for index, row in xldf.iterrows():
        if row["STATUS"] == "DONE" and type(row.get("PLMO_Number")) == str:#  math.isnan(float(row.get("PLMO_Number"))):
            vcfFolder = LINUX_ANALYSIS_OUT_FOLDER + "/" +  row['SAMPLE_DIR_PATH'].strip() + '_' + row['TIME_STAMP']  + "/"
            accessionId = row['SAMPLE_DIR_PATH'].split("/")[1].strip() + '_' + row['TIME_STAMP']
            if accessionIdStatusMap.get(accessionId) is not None and not os.path.isfile(vcfFolder+accessionId+".QCIXml.xml"):
                text_file = callQCIApi(accessionId, row.get("PLMO_Number"), vcfFolder + accessionId)
                if not text_file:
                    print("could not pull XML from QCI")
    
    #time.sleep(600)    
    #logging.debug('=======================Execution ends===========================')
    excel_file.close()

#for handler in logging.root.handlers[:]:
#   logging.root.removeHandler(handler)
#logging.basicConfig(format='%(asctime)s - %(message)s', filename='loggingForGetXML.log', level=logging.DEBUG)

main()


#accessionId = "NS-19-13_BC710502_53_20190416154205825653DevEnv3.0"
#accessToken = generateAccessToken()
#y, status = getXMLForAccessionId(accessionId, accessToken["access_token"])
#if y and status == 200:
#    print(parseXML(y, accessionId, 'fdadfs'))

