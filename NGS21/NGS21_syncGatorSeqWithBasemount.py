import os
import yaml
import sys
import shutil
import datetime

print("Run start time: ", str(datetime.datetime.now()) + "\n")

script_path = os.path.dirname(os.path.abspath( __file__ ))
parent_path = os.path.abspath(os.path.join(script_path, '..'))

CONFIG_FILE = parent_path + "/linux_gatorseq.config.yaml"

config_dict=dict()
with open(CONFIG_FILE, 'r') as stream:
    try:
        config_dict=yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

LINUX_HPC_FASTQ_FOLDER = (config_dict['LINUX_PATHOLOGY_FASTQ_FOLDER_NGS21'])

# get list of all folders in basemount/gator_seq samples
# for each sample, checks if fastq are in the LINUX_HPC_FASTQ_FOLDER, if they are not present download them.
def downloadFiles(baseMountDir):
    baseMountBase = baseMountDir + "/Projects/Gatorseq_NGS/Samples/"
    #baseMountBase = baseMountDir + "/Projects/Gatorseq NGS/Samples/"
    allSamples = os.listdir(baseMountBase)

    # ToDo: check with prof, for NQ report
    eligibleSamples = [sample for sample in allSamples if sample[:2] == "NQ" and os.path.isdir(baseMountBase + sample)]
    
    print("eligibleSamples:", eligibleSamples)
    for sample in eligibleSamples:
        if "_" in sample:
            sampleFolderName = sample.split("_")[0] + "-000000/" + sample + "-000000/"
            print("sampleFolderName:", sampleFolderName)

            # if os.path.isdir(LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName):
            #     print("nothing to do")
            # else:
            #     #make the directory
            #     os.makedirs(LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName)
            
            #copy fastq
            if not os.path.isdir(LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName):
                os.makedirs(LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName)
            
            for fastqFile in os.listdir(baseMountBase + sample + "/Files/"):
                if not os.path.isfile( LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName + "/" + fastqFile ):
                    try:
                        print("copying:", baseMountBase + sample + "/Files/" + fastqFile, " to ", LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName + "/" + fastqFile)
                        shutil.copyfile(baseMountBase + sample + "/Files/" + fastqFile, LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName + "/" + fastqFile)
                    except:
                        print("copy file error")
                        os.remove(LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName +  "/" + fastqFile)
        else:
            print("file structure not clear.. please check")    
                    
if __name__ == "__main__":
    #baseMountDir = "BaseMount"
    baseMountDir = "/home/path-svc-mol/BaseSpace_Mount"
    #basemount --unmount ~/$baseMountDir
    #basemount  --config UFMOL_ENTERPRISE  ~/$baseMountDir

    downloadFiles(baseMountDir)

