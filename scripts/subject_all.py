import os
import datetime
import dbdicom as db

import pipelines.rename as rename
import pipelines.mdr as mdr
import pipelines.mapping as map
import pipelines.export_ROI_stats as export_ROIs
import pipelines.segment as seg

import scripts.upload as upload
from scripts import xnat

def single_subject(username, password, path, dataset):
    
    #import data from XNAT
    ExperimentName = xnat.main(username, password, path, dataset)
    #ExperimentName = "Leeds_REP_VOL_001_01"
    pathScan = path + "//" + ExperimentName
    filename_log = pathScan + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "MDRauto_LogFile.txt" #TODO FIND ANOTHER WAY TO GET A PATH
    
    #available CPU cores
    try: 
        UsedCores = int(len(os.sched_getaffinity(0)))
    except: 
        UsedCores = int(os.cpu_count())

    folder = db.database(path=pathScan)
    folder.set_log(filename_log)
    folder.log("Analysis of " + pathScan.split('//')[-1] + " has started!")
    folder.log("CPU cores: " + str(UsedCores))
    
    #name standardization 
    try:
        print("starting renaming")
        rename.main(folder)
    except Exception as e:
        folder.log("Renaming was NOT completed; error: " + str(e))

    #Apply motion correction using MDR
    try:
        print("starting mdr")
        mdr.main(folder)
    except Exception as e:
        folder.log("Renaming was NOT completed; error: " + str(e))

    #Custom modelling
    try:
        print("staring mapping")
        map.main(folder)
    except Exception as e:
        folder.log("Modelling was NOT completed; error: " + str(e))

    #Generate masks using unetr, apply alignment, extract biomarkers to a .csv
    try:
        print('starting parameter extraction')
        filename_csv = export_ROIs.main(folder,ExperimentName)
    except Exception as e:
        folder.log("Parameter extraction was NOT completed; error: " + str(e))
    
    #upload images, logfile and csv to google drive
    upload.main(pathScan, filename_log, filename_csv)
