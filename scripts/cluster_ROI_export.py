""" 
@author: Joao Periquito 
iBEAt CLUSTER MAIN Scrpit
2022
Download XNAT dataset -> Name Standardization -> Execute MDR    -> Custom Moddeling (DCE, T2*, DCE)  -> T1 & T2 modelling with parallelization (done in the main)
    XNAT_cluster.py   ->  RENAME_Cluster.py   -> MDR_Cluster.py -> MODELLING_cluster.py
(T1 & T2 modelling are done in the main due to parallelization requirements)

TO RUN THE SCRIPT YOU USE: python main_cluster.py --num n (WHERE n is an integer with the value of the XNAT dataset)
"""

# To develop the application
# --------------------------
# py -3 -m venv .venv           # create virtual environment
# .venv/Scripts/activate        # activate virtual environment
import os
import datetime
import dbdicom as db

import pipelines.mapping as map
import pipelines.export_ROI_stats as export_ROIs
import zipfile
import shutil

if __name__ == '__main__':

    # Define the directory you want to list folders from
    directory_path = "G:\\My Drive\\CLUSTER_RESULTS\\Exeter\\Successful"
    #directory_path = "C:\\Users\\md1jdsp\\Desktop\\G_drive"
    mask_path = "C:\\Users\\md1jdsp\\Desktop\\iBEAt_KidneyMasks"

    dataset_list=[]
    items = os.listdir(directory_path)
    for element in items:        
        if len(element)<30 and element != 'desktop.ini' and element != 'Done'and element != 'Failed':
            dataset_list.append(element)

    mask_list = os.listdir(mask_path)

    extract_to_directory = "C:\\Users\\md1jdsp\\Desktop\\iBEAT_results_temp"

    for dataset in dataset_list:
        ExperimentName = dataset

        pathScan = directory_path + "\\" + ExperimentName

        with zipfile.ZipFile(pathScan, 'r') as zip_ref:
            zip_ref.extractall(extract_to_directory)

        ExperimentName = ExperimentName[:-4]
        pathScan = extract_to_directory + "\\" + ExperimentName

        # try:
        #     MaskExperimentName_LK = ExperimentName.split("-")[1] + "_" + ExperimentName.split("-")[2] + "_LK"
        #     MaskExperimentName_RK = ExperimentName.split("-")[1] + "_" + ExperimentName.split("-")[2] + "_RK"
        # except:
        #     MaskExperimentName_LK = ExperimentName.split("_")[-1][:4] + "_" + ExperimentName.split("_")[-1][4:] + "_LK"
        #     MaskExperimentName_RK = ExperimentName.split("_")[-1][:4] + "_" + ExperimentName.split("_")[-1][4:] + "_RK"

        try:
            MaskExperimentName_LK = ExperimentName.split("_")[0].split("-")[1] + "_" + ExperimentName.split("_")[0].split("-")[2] + "_LK"
            MaskExperimentName_RK = ExperimentName.split("_")[0].split("-")[1] + "_" + ExperimentName.split("_")[0].split("-")[2] + "_RK"
        except:
            MaskExperimentName_LK = ExperimentName.split("_")[0].split("-")[1] + "_" + ExperimentName.split("_")[1] + "_LK"
            MaskExperimentName_RK = ExperimentName.split("_")[0].split("-")[1] + "_" + ExperimentName.split("_")[1] + "_RK"

        shutil.copytree(mask_path+"\\"+MaskExperimentName_LK, pathScan+"\\"+"Mask_LK")
        shutil.copytree(mask_path+"\\"+MaskExperimentName_RK, pathScan+"\\"+"Mask_RK")

        filename_log = pathScan + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "ROI_stats_LogFile.txt" #TODO FIND ANOTHER WAY TO GET A PATH
        
        try: 
            UsedCores = int(len(os.sched_getaffinity(0)))
        except: 
            UsedCores = int(os.cpu_count())

        folder = db.database(path=pathScan)
        folder.scan()

        folder.set_log(filename_log)
        folder.log("Analysis of " + pathScan.split('//')[-1] + " has started!")
        folder.log("CPU cores: " + str(UsedCores))
        
        # try:
        #     map.main(folder)
        # except Exception as e:
        #     folder.log("Modelling was NOT completed; error: " + str(e))

        try:
            export_ROIs.main(folder,ExperimentName)
        except Exception as e:
            folder.log("Export ROI was NOT completed; error: " + str(e))
        
        shutil.rmtree(pathScan)

