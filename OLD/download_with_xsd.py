import requests 
from requests.exceptions import HTTPError 
import xmlschema 
import time 
import logging
import logging.config
import ast
import sys
import yaml
import numpy
import pandas as pd
import hl7
import hl7update
import os
import datetime
from shutil import copyfile
from shutil import move
import untangle
import xmltodict

UPLOAD_DOWNLOAD_PATH = 'Z:\MIRTH_GATORSEQ\TEST'
EXCEL_FILE_PATH = r'C:\Users\s.majety\Desktop\Copy of Sheet1.xlsx'
VCF_DIRECTORY_PATH = r'G:/DRL/Molecular/NGS/GenomOncology/NextSeq/'

schema_file_name = os.path.join(os.path.dirname(__file__), 'reportXML.xsd')
qci_schema = xmlschema.XMLSchema(schema_file_name)

# Gets all the accessionIds with 'final' status. Used to check if a accessionId is ready to be pulled from Qiagen
def getAllStatus():
    url = "https://api.ingenuity.com/datastream/api/v1/clinical?state=final&startReceivedDate=2018-04-28"
    headers = {"Accept": "application/json", "Authorization": generateAccessToken()["access_token"]}
    response = requests.get(url, headers = headers)
    #print(ast.literal_eval(response.content.decode("utf-8")))
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
    response = requests.post("https://api.ingenuity.com/v1/oauth/access_token", {"grant_type": "client_credentials", "client_id": "bd7a1714af065a9219d6d31ea6a1da22", "client_secret": "83c821283fb21969494bb205863e5211"}, {"content-type": "text/plain;charset=UTF-8"} )
    return ast.literal_eval(response.content.decode("utf-8"))

# downloads the xml for a given accessionId
def getXMLForAccessionId(accessionId, accessToken):
    url = "https://api.ingenuity.com/v1/export/"
    accessionId = accessionId or "NS-19-11_BC712508_27T_20190317175414382540ProdEnv3.0"
    url = url + accessionId
    return makeGetCall(url, {"view" : "reportXml"}, {"Accept": "application/xml", "Authorization": accessToken})

# compares the xml with the xsd file and genreates a dictionary
def parseXML(xml, accessionId):
    filename = os.path.join(os.path.dirname(__file__) + "\\data\\xml", "xml-" + accessionId + ".xml")
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w", encoding='utf-8') as f:
        f.write(xml)
        f.close()
    map = qci_schema.to_dict(filename)
    return generateTxtFileAndSave(map, accessionId), filename

def newLine():
    return " \n \n" 

# helper function to generate comments
def populateIndividualText(text, i, list, nextSeperator, gene_map):
    nextSeperator = nextSeperator.split(",")
    foundAtleastOne = False
    for iterator in range(i, len(list)):
        variant = list[iterator]
        
        if variant["assessment"] in nextSeperator:
            break
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
    return text, iterator, gene_map

def addHeading(text, heading):
    text += "------------------------------------------------------------" + "\n"
    text += heading + " \n"
    text += "------------------------------------------------------------" + "\n"
    return text

#populates a map with each drug name. it is required as we need to show the treatments grouped by drugnames
def getDrugMaps(treatmentsList):
    drugMap = {}
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
                if variant.get("proteinchange"):
                    text += str(variant.get("proteinchange").get("change")) + " "
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
            if treatment.get("rcomment") and treatment.get("rcomment")[0]:
                text += "    Therapies description: " + treatment.get("rcomment")[0].get("text") + "\n"
            text += " \n"
    if not foundAtleastOne:
        text += "    None" + "\n"
    return text

# generate comments from the dictionary created from downloaded xml   
def generateTxtFileAndSave(map, accessionId):
    text = ''
    gene_map = {}
    custom_order = ["Pathogenic", "Likely Pathogenic", "Uncertain Significance"]
    order = {key: i for i, key in enumerate(custom_order)}
    variants = None
    if map.get("variant"):
        variants = sorted(map.get("variant"), key=lambda d: order[d['assessment']])
    
    text += "UF HEALTH PATHLABS GATORSEQ177 NGS SCREENING PANEL" + newLine()
    if map.get("accession"):
        text += "Accession Id: " + map.get("accession") + "\n"
    if map.get("diagnosis"):
        text += "Diagnosis: " + map.get("diagnosis") + "\n"
        
    
    text += newLine()
    text = addHeading(text, "1. Pathogenic variants")
    if variants is not None:
        text, iterator, gene_map = populateIndividualText(text, 0, variants, "Likely Pathogenic,Uncertain Significance", gene_map)
    else:
        text += "    None" + " \n"
    text += " \n"
    text = addHeading(text, "2. Likely pathogenic variants")
    if variants is not None:
        text, iterator, gene_map = populateIndividualText(text, iterator, variants, "Uncertain Significance", gene_map)
    else:
        text += "    None" + " \n"
    text += " \n"
    text = addHeading(text, "3. Variants of unknown significance (Tier 3)")
    if variants is not None:
        text, iterator, gene_map = populateIndividualText(text, iterator, variants, "", gene_map)
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
            text += "Gene: " + str(variant["gene"]) + "\n"
            if variant.get("transcriptchange"):
                text += "    " + "Exon: " + str(variant.get("transcriptchange").get("exonNumber")) + "\n"
                text += "    " + "Nucleotide: " + str(variant.get("transcriptchange").get("transcript")) + ":" + variant.get("transcriptchange").get("change") + "\n"
            
            if variant.get("genomicchange"):
                text += "    " + "Genomic Location: " + str(variant.get("genomicchange").get("change")) + "\n"
                
            if variant.get("proteinchange"):
                text += "    " + "Amino acid: " + str(variant.get("proteinchange").get("change")) + "\n"
                
            if variant.get("function"):
                text += "    " + "Function: " + str(variant.get("function")) + "\n"
                
            if variant.get("assessment"):
                text += "    " + "Assessment: " + str(variant.get("assessment")) + "\n"
            if variant.get("actionability"):
                text += "    " + "Classification: Tier " + str(variant.get("actionability")) + "\n"
                
            if variant.get("allelefraction"):
                text += "    " + "Allele Fraction: " + str(variant.get("allelefraction")) + "(of "+ str(variant.get("readDepth")) +" reads)" + "\n"
                   
            if variant.get("variation"):
                text += "    " + "Variation: " + str(variant.get("variation")) + "\n"
            if variant.get("rcomment") and variant.get("rcomment")[0]:
                text += "    " + "Interpretation: " + str(variant.get("rcomment")[0].get("text")) + "\n"
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
    text += str(map.get("labTestedGenes")) + " \n"
    
    text += newLine()

    text = addHeading(text, "9. Appendix") 
    text += """For all non-FDA approved tests listed above unless specified above, UF Health Pathology Laboratories (UFHPL) test have been developed and their performance characteristics determined by UFHPL. These tests have not been cleared or approved by the U.S. Food and Drug Administration. The FDA has determined that such clearance or approval is not necessary. Tests are used for clinical purposes andshould not be regarded as investigational or for research purposes. UFHPL is certified under the Clinical Laboratory Improvement Amendments of 1988 (CLIA) as qualified to perform high complexity clinical testing. Decisions on patient care and treatment must be based on the independent medical judgment of the treating physician, taking into considerations all applicable information concerning the patient's condition including family history, physical examination, patient preferences and other diagnostic tests in accordance with the applicable standard of care. Treatment decisions should not be made based solely on the information contained in this report. The presence of a genomic alteration that has been reported to show response to a particular therapy is no guarantee that the agent will show efficacy in any given patient. Any reference to available clinical trials or treatment options is based on the currently available public knowledge aggregated by publicly available websites which we do not provide or update content nor are we responsible for monitoring the update of content. Treating physician certifies by receiving this report that the patient has been adequately counseled about the results and implications of this genetic test in accordance with all applicable laws."""
    
    # text += newLine()
    # text += "Selected Citations" + "\n"
    # if map.get("article") is not None:
    #     for article in map.get("article"):
    #         text += "    " +  str(article.get("@author")) + " " + str(article.get("@title")) + " " + str(article.get("@citation")) + " " + str(article.get("@date")) + " " + str(article.get("@url")) + "\n\n" 

    filename = os.path.join(os.path.dirname(__file__) + "\\data\\txt", "txt-" + accessionId +".txt")
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
def callQCIApi(accessionId):
    accessToken = generateAccessToken()
    accessionId = accessionId or "NS-19-11_BC707507_92_20190319110904776453DevEnv3.0"
    y, status = getXMLForAccessionId(accessionId, accessToken["access_token"])
    if y and status == 200:
        return parseXML(y, accessionId)
    return False

# 1. Tries to open excel and exits if already open
# 2. Iterates over the folder and tries to read PLMO number from each hl7 file
# 3. Checks if there is an entry in excel for that PLMO and retrieves the corresponding accession id.
# 4. queries the Qiagen to see if the accessionid is in final state
# 5. downloads the xml, parses it and generates comments
# 5. processes the genes and moves the new hl7 file to a results folder where it gets pushed to EPIC
# 6. archives the initial file 
def main():
    # if( len(sys.argv)>=3):
        # UPLOAD_DOWNLOAD_PATH =  sys.argv[1]
        # EXCEL_FILE_PATH = sys.argv[2]
        
    UPLOAD_PATH = UPLOAD_DOWNLOAD_PATH + '\RESULTS'
    ORDERS_ARCHIVE_DIR = UPLOAD_DOWNLOAD_PATH + '\ORDERS_ARCHIVE\\'
    ORDERS_DIR = UPLOAD_DOWNLOAD_PATH + '\ORDERS\\'
    #Check if excel file is opened by any other user
    try: 
        excel_file = open(EXCEL_FILE_PATH, "r+")
    except:
        logging.error(" Could not open file! Please close Excel!")
        sys.exit()
    try:
        xldf = pd.read_excel(EXCEL_FILE_PATH)
    except:
        logging.error("Problem Reading Excel")
        sys.exit()

    accessionIdStatusMap = populateStatusMap()
    allhl7filenames = []
    for (dirpath, dirnames, filenames) in os.walk(ORDERS_DIR):
        allhl7filenames.extend(filenames)
        break
    for hl7_file_name in allhl7filenames:
        try:
            hl7file = open(ORDERS_DIR + hl7_file_name, mode="r").read()
        except:
            logging.error(hl7_file_name + " file not found")
            continue
        arr = hl7file.split("\n\n") #split by blank lines
        #Iterate all HL7 messages in file

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
                    logging.error('PLM not found in HL7 message!')
                    continue
                if not(xldf['PLMO_Number'].str.contains( str(plm)).any()):
                    logging.error('PLM is not found in excel')
                    continue
                #print((plm))
                sample_dir_path = xldf[xldf['PLMO_Number'] == str(plm)]['SAMPLE_DIR_PATH'].item()
                if (not sample_dir_path):
                    continue
                #if xldf[xldf['PLMO_Number'] == str(plm)]['downloadedXML'].item() == 0:
                accessionId = sample_dir_path.split("/")[1] + "_" + xldf[xldf['PLMO_Number'] == str(plm)]['TIME_STAMP'].item()
                if accessionIdStatusMap.get(accessionId) is not None:
                    text_file, xml_fileName = callQCIApi(accessionId)
                    if(text_file):
                        vcfFolder = VCF_DIRECTORY_PATH + xldf[xldf['PLMO_Number'] == str(plm)]['SAMPLE_DIR_PATH'].item() + '_' + xldf[xldf['PLMO_Number'] == str(plm)]['TIME_STAMP'].item()  + "/"
                        copyfile(text_file, vcfFolder+accessionId+".QCIreport.txt")
                        copyfile(xml_fileName, vcfFolder+accessionId+".QCIXml.xml")
                        genes_list = hl7update.find_genes_from_XML(xml_fileName)
                        #Map structure {"KIT": "c.1739_1740insTGACCCAACACAACTTCCTTATGATCA p.D572_H580dup"}
                        gene_map={}
                        if(genes_list):
                            gene_map = dict(gene.split(" ", 1) for gene in genes_list)
                        
                            
                        hl7update.update_msh_segment(h)
                        hl7update.update_orc_segment(h)
                        hl7update.update_obr_segment(h)
                        hl7update.update_comments(h, open( text_file, mode="r",  encoding='utf-8').read())
                        hl7update.update_obx_segment(h)
                        h = hl7update.update_obx_seg_containing_gene(h, gene_map)
                        
                        out_file_path = UPLOAD_PATH + '\hl7-{}-output.txt'.format(plm)
                        if h:
                            with open(out_file_path, 'w' ,  encoding='utf-8') as f:
                                f.write(str(h))
                            print("Out file available at :",out_file_path)
                            moveHL7ToArchive(ORDERS_DIR + hl7_file_name, ORDERS_ARCHIVE_DIR + hl7_file_name)  
                        else:
                            logging.error("Couldn't replace '-' in hl7. Check logs for more details!")
                    else:
                        logging.error("could not generate comment files for accessionId: " + accessionId)

    logging.debug('=======================Execution ends===========================')
    excel_file.close()

for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
logging.basicConfig(format='%(asctime)s - %(message)s', filename='loggingForGetXML.log', level=logging.DEBUG)
main()
# generateTxtFileAndSave(qci_schema.to_dict(r'C:\Users\s.majety\Desktop\hl7 project\Scripts to download\data\xml\xml-NS-19-20_BC704504_880_20190512190730780905ProdEnv3.0.xml'), "NS-19-20_BC704504_880_20190512190730780905ProdEnv3.0")
# listAccessions = ["NS-19-21_BC710506_203_20190520074304544964ProdEnv3.0",
# "NS-19-21_BC709505_199_20190520074304282257ProdEnv3.0",
# "NS-19-21_BC711507_219_20190520074305012015ProdEnv3.0",
# "NS-19-21_BC712508_278_20190520074305918928ProdEnv3.0", 
# "NS-19-21_BC701501_286_20190520074306468020ProdEnv3.0", 
# "NS-19-21_BC703503_290_20190520074307102419ProdEnv3.0",
# "NS-19-21_BC704504_433_20190520074307527555ProdEnv3.0",
# "NS-19-21_BC705505_442_20190520074307890863ProdEnv3.0", 
# "NS-19-21_BC707503_097T_20190520120403576779ProdEnv3.0",
# "NS-19-21_BC702502_287_20190520120404805194ProdEnv3.0"] 

# for accessionId in listAccessions:
#     accessToken = generateAccessToken()
#     accessionId = accessionId or "NS-19-11_BC707507_92_20190319110904776453DevEnv3.0"
#     y, status = getXMLForAccessionId(accessionId, accessToken["access_token"])
#     if y and status == 200:
#         print(parseXML(y, accessionId))
        