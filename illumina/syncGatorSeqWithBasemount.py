import os
import yaml
import sys
import shutil

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

LINUX_HPC_FASTQ_FOLDER = (config_dict['LINUX_HPC_FASTQ_FOLDER'])

# get list of all folders in basemount/gator_seq samples
# for each sample, checks if fastq are in the LINUX_HPC_FASTQ_FOLDER, if they are not present download them.
def downloadFiles(baseMountDir):
    baseMountBase = baseMountDir + "/Projects/Gatorseq_NGS/Samples/"
    allSamples = os.listdir( baseMountBase )
    eligibleSamples = [sample for sample in allSamples if sample[:2] == "NQ" and os.path.isdir(baseMountBase + sample)]
    
    print("eligibleSamples:", eligibleSamples)
    for sample in eligibleSamples:
        if "_" in sample:
            sampleFolderName = sample.split("_")[0] + "-basemount/" + sample + "-basemount/"
            print("sampleFolderName:", sampleFolderName)

            if os.path.isdir(LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName):
                print("nothing to do")
            else:
                #make the directory
                os.makedirs(LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName)

                #copy fastq
                try:
                    for fastqFile in os.listdir(baseMountBase + sample + "/Files/"):
                        shutil.copyfile(baseMountBase + sample + "/Files/" + fastqFile, LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName + "/" + fastqFile)
                except:
                    print("copy file error")
                    os.removedirs(LINUX_HPC_FASTQ_FOLDER + "/" + sampleFolderName)
        else:
            print("file structure not clear.. please check")    
        break       
                    
if __name__ == "__main__":
    baseMountDir = "BaseMount"

    downloadFiles(baseMountDir)

