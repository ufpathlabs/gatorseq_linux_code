import os
import shutil
import hashlib

GATORSEQ_CANCER_BASEMOUNT_FOLDER_NAME = "/home/path-svc-mol/Replica_BaseSpace/Projects/Gatorseq_NGS/Samples/"
#LINUX_PATHOLOGY_FASTQ_FOLDER = "../test/test2"
LINUX_PATHOLOGY_FASTQ_FOLDER = "/ext/path/DRL/Molecular/NGS/NextSeq_Fastq"

#converting NQ-20-02_BC702506_239 into NQ-20-02-000000/NQ-20-02_BC702506_239-000000
def getNameInFastQFormat(name):
    TAIL = "-000000"
    split = name.split("_")
    batch = split[0]

    sample = "_".join(split[1:])

    #print("received ->" + name + "<- returning -->" + batch + TAIL + "/" + name + TAIL)
    
    return batch + TAIL + "/" + name + TAIL

# get all sample names from GATORSEQ_CANCER_BASEMOUNT_FOLDER_NAME
def getAllSampleNames():
    outList = []
    files = os.listdir(GATORSEQ_CANCER_BASEMOUNT_FOLDER_NAME)
    for f in files:
        if os.path.isdir(GATORSEQ_CANCER_BASEMOUNT_FOLDER_NAME + f):
            outList.append(f)
    return outList

#check if the sampleNames exist in the LINUX_PATHOLOGY_FASTQ_FOLDER Folder and if not present, then copy the contents creating new folders
# folderName eg: NQ-20-02_BC702506_239
def checkIfFolderExists(folderName):
    fileNameNGS = getNameInFastQFormat(folderName)

    if os.path.isfile(LINUX_PATHOLOGY_FASTQ_FOLDER +"/" + fileNameNGS + "/success.txt"):
        print(folderName + " already exists as " + LINUX_PATHOLOGY_FASTQ_FOLDER + fileNameNGS)
    else:
        copyFastqFiles(GATORSEQ_CANCER_BASEMOUNT_FOLDER_NAME + folderName + "/Files/", LINUX_PATHOLOGY_FASTQ_FOLDER + "/"+ fileNameNGS + "/" )

# copy all the fastQ files from src to dest
def copyFastqFiles(src, dest):
    sourceFiles = os.listdir(src)
    if not os.path.exists(dest):
        os.makedirs(dest)
    issueDetected = False
    print("copying-->" + "&".join(sourceFiles))
    for fileName in sourceFiles:
        srcFileHash = hashlib.md5(open(src + fileName, 'rb').read()).hexdigest()
        if os.path.isfile(src + fileName):
            shutil.copy(src + fileName, dest)
        destFileHash = hashlib.md5(open(dest + fileName, 'rb').read()).hexdigest()
        if srcFileHash != destFileHash:
            print("not copied correctly")
            issueDetected = True
    if not issueDetected:
        f = open(dest + "success.txt", "w")
        f.write("FASTQ FILES COPIED successfully")
        f.close()

if __name__ == "__main__":
   # print("im here")
    sampleNames = getAllSampleNames()
   # print(sampleNames)
   # print("sfg")
   # if len(sampleNames)> 0:
   #     checkIfFolderExists(sampleNames[1])
    for sample in sampleNames:
        checkIfFolderExists(sample)
    pass
    



