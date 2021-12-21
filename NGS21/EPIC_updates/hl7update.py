import hl7
import argparse
import xml.etree.ElementTree as ET
import glob, os
import yaml
from html.parser import HTMLParser
import datetime
import sys
import xmltodict


def format_for_unity(x):
    if (x < 10):
        return "0" + str(x)
    else:
        return str(x)


def get_current_formatted_date():
    currentDT = datetime.datetime.now()
    if currentDT:
        data = format_for_unity(currentDT.year) + format_for_unity(currentDT.month) + format_for_unity(
            currentDT.month) + format_for_unity(currentDT.hour) + format_for_unity(currentDT.minute) + format_for_unity(
            currentDT.second)

        return data
    return str(currentDT)


def get_ticks(dt):
    return (dt - datetime.datetime(1, 1, 1)).total_seconds() * 10000000


# Return example: KIT c.1739_1740insTGACCCAACACAACTTCCTTATGATCA p.D572_H580dup
def find_genes_from_XML(xmlFile):
    with open(xmlFile, encoding='utf-8') as fd:
        map = xmltodict.parse(fd.read())
        map = map['report']
    diagnosis = ""
    if map.get("diagnosis") is not None:
        diagnosis = map.get("diagnosis")
    genes_list = []
    if map.get("variant") is not None:
        if type(map.get("variant")) != list:
            map["variant"] = [map.get("variant")]
        for variant in map.get("variant"):
            if variant["assessment"] == "Uncertain Significance":
                if variant.get("gene") is not None:
                    print("ignoring uncertain imortance genes: ", str(variant["gene"]))
                continue
            if variant.get("gene") is not None:
                gene = ''
                transcriptChange = ""
                proteinChange = ""
                if variant.get("transcriptchange"):
                    transcriptChange = variant["transcriptchange"]["change"]
                if variant.get("proteinchange"):
                    proteinChange = variant["proteinchange"]["change"]
                gene += str(variant["gene"]) + " " + transcriptChange + " " + proteinChange
                genes_list.append(gene)
    return genes_list, diagnosis


def get_first_obx_index(h):
    idx = 0
    for seg in h:
        if seg[0][0] == 'OBX':
            return idx
        idx += 1
    return -1


# Assuming insertion is just above first OBX segment
def update_comments(h, comments):
    comments_arr = comments.split("\n")
    obx_idx = get_first_obx_index(h)
    if (obx_idx == -1):
        print("OBX segment not found, so appending it")
        obx_idx = len(h) - 1
    i = 1
    for comment in comments_arr:
        h.append('NTE|{}|L|{}'.format(i, comment))
        obx_idx += 1
        i += 1


def update_msh_segment(h):
    if h and h['MSH']:
        for msh_segment in h['MSH']:
            if msh_segment:
                msh_segment[7] = get_current_formatted_date()
                msh_segment[8] = ''
                msh_segment[9][0][0] = 'ORU'
                msh_segment[9][0][1] = 'R01'
                msh_segment[10] = get_ticks(datetime.datetime.now())


def update_orc_segment(h):
    if h and h['ORC']:
        for orc_segment in h['ORC']:
            orc_segment[1] = 'RE'


def update_obr_segment(h):
    if h and h['OBR']:
        for obr_segment in h['OBR']:
            obr_segment[22] = get_current_formatted_date()
            obr_segment[25] = 'P'
            obr_segment[27] = '^^^^^R^^'


def update_obx_segment(h):
    if h and h['OBX']:
        for obx_segment in h['OBX']:
            obx_segment[2] = 'ST'
            obx_segment[11] = 'P'
            obx_segment[14] = get_current_formatted_date()
            if (len(obx_segment) == 19):
                obx_segment.append(obx_segment[14])
            elif (len(obx_segment) >= 19):
                obx_segment[19] = obx_segment[14]


def update_obx_seg_containing_gene(h, gene_map, accessionId, diagnosis, Perc_Target_Cells, Perc_Tumor):
    updates = 0
    temp_obx = h[:]
    l = len(h)
    for i in range(l):
        del temp_obx[l - i - 1]
    new_obx_index = 1
    toReplace = {"ERBB2": "ERBB2/HER2", "NSD2": "WHSC1", "MRTFA": "MKL1"}
    for gene in gene_map.keys():
        if gene in toReplace.keys():
            gene_map[toReplace[gene]] = gene_map[gene]
            del gene_map[gene]
    for obxSegment in h['OBX']:
        # print(obxSegment)
        # print(obxSegment[3])
        # print(obxSegment[3][0][1][0])
        # print(obxSegment[7][0])
        if (obxSegment[3][0][1][0] in gene_map.keys() and obxSegment[7][0] == "-"):
            obxSegment[5][0] = ", ".join(gene_map[obxSegment[3][0][1][0]])
            obxSegment[1] = new_obx_index
            new_obx_index += 1
            temp_obx.append(obxSegment)
            updates += 1
            # print(obxSegment[3][0][1][0])
        elif obxSegment[3][0][1][0] == "BARCODE":
            # print("barcode found")
            obxSegment[5][0] = accessionId
            obxSegment[1] = new_obx_index
            new_obx_index += 1
            temp_obx.append(obxSegment)
        elif obxSegment[3][0][1][0] == "TUMOR TYPE":
            # print("tumor type found")
            obxSegment[5][0] = diagnosis
            obxSegment[1] = new_obx_index
            new_obx_index += 1
            temp_obx.append(obxSegment)
        elif obxSegment[3][0][1][0] == "% TARGET CELLS":
            # print("tumor type found")
            if Perc_Target_Cells:
                obxSegment[5][0] = Perc_Target_Cells
                obxSegment[1] = new_obx_index
                new_obx_index += 1
                temp_obx.append(obxSegment)
        elif obxSegment[3][0][1][0] == "ATYPICAL LYMPHOCYTES IN MICRODISSECTED TISSUE":
            # print("tumor type found")
            if Perc_Tumor:
                obxSegment[5][0] = Perc_Tumor
                obxSegment[1] = new_obx_index
                new_obx_index += 1
                temp_obx.append(obxSegment)

    h_t = h[:]
    l = len(h)
    for i in range(l):
        del h_t[l - i - 1]
    for i in range(len(h)):
        if (h[i][0][0] != "OBX"):
            h_t.append(h[i])
    h_t.extend(temp_obx)
    # ToDo: remove the comment for False
    if (updates < len(gene_map.keys())):
        print('Some genes were not found in OBX segment! Please check logs.')
        return False
    else:
        return h_t


def search_file(file_name_query, path):
    os.chdir(path)
    search_result = glob.glob(file_name_query)
    if (not search_result):
        print("Couldn't find a file similar to {}".format(file_name_query))
        return None
    return search_result[0]

