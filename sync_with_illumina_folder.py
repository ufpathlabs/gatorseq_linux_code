import os
import shutil

ILLUMINA_FOLDER_NAME = "/home/path-svc-mol/Replica_BaseSpace/Projects/Gatorseq_NGS/Samples/"
PATH_TO_COPY = "/ext/path/DRL/Molecular/NGS/NextSeq_Fastq/"

#converting NQ-20-02_BC702506_239 into NQ-20-02-000000/NQ-20-02_BC702506_239-000000
def getNameInFastQFormat(name):
    TAIL = "-000000"
    split = name.split("_")
    batch = split[0]

    sample = "_".join(split[1:])

    print("received ->" + name + "<- returning -->" + batch + TAIL + "/" + name + TAIL)
    
    return batch + TAIL + "/" + name + TAIL

# get all sample names from ILLUMINA_FOLDER_NAME
def getAllSampleNames():
    outList = []
    files = os.listdir(ILLUMINA_FOLDER_NAME)
    for f in files:
        if os.path.isdir(f):
            outList.append(f)
    return outList

#check if the sampleNames exist in the PATH_TO_COPY Folder and if not present, then copy the contents creating new folders
# folderName eg: NQ-20-02_BC702506_239
def checkIfFolderExists(folderName):
    fileNameNGS = getNameInFastQFormat(folderName)

    if os.path.isdir(PATH_TO_COPY + fileNameNGS):
        print(folderName + " already exists as " + PATH_TO_COPY + fileNameNGS)
    else:
        copyFastqFiles(ILLUMINA_FOLDER_NAME + folderName +  "/Files/", PATH_TO_COPY + fileNameNGS + "/" )

# copy all the fastQ files from src to dest
def copyFastqFiles(src, dest):
    sourceFiles = os.listdir("src")
    print("copying-->" + "&".join(sourceFiles))
    for fileName in sourceFiles:
        if os.path.isfile(fileName):
            shutil.copy(fileName, dest)


if __name__ == "__main__":
    sampleNames = getAllSampleNames()
    for sample in sampleNames:
        checkIfFolderExists(sample)
    pass
    



