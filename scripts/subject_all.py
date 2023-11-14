import os
import datetime
import dbdicom as db

import pipelines.rename as rename
import pipelines.mdr as mdr
import pipelines.mapping as map
import scripts.upload as upload
import scripts.segmentation as segmentation
from scripts import xnat


def single_subject(username, password, path, dataset):
    
    ExperimentName = xnat.main(username, password, path, dataset)
    pathScan = path + "//" + ExperimentName
    filename_log = pathScan + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "MDRauto_LogFile.txt" #TODO FIND ANOTHER WAY TO GET A PATH
    
    try: 
        UsedCores = int(len(os.sched_getaffinity(0)))
    except: 
        UsedCores = int(os.cpu_count())

    folder = db.database(path=pathScan)
    folder.set_log(filename_log)
    folder.log("Analysis of " + pathScan.split('//')[-1] + " has started!")
    folder.log("CPU cores: " + str(UsedCores))
    
    try:
        print("starting renaming")
        rename.main(folder)
    except Exception as e:
        folder.log("Renaming was NOT completed; error: " + str(e))

    try:
        print("starting mdr")
        mdr.main(folder)
    except Exception as e:
        folder.log("Renaming was NOT completed; error: " + str(e))

    try:
        print("staring mapping")
        map.main(folder)
    except Exception as e:
        folder.log("Modelling was NOT completed; error: " + str(e))

    try:
        print('starting segmentation')
        segmentation.main(folder,pathScan)
    except Exception as e:
        folder.log("Segmentation was NOT completed; error: " + str(e))

    upload.main(pathScan, filename_log)
