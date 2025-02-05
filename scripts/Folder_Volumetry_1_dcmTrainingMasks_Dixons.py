import os
import datetime
import dbdicom as db

from utilities import xnat, upload_folder_to_drive
from scripts import steps_core
from scripts import steps_internal


def single_subject(username, password, path, dataset):
    
    #Import data from XNAT
    if isinstance(dataset,str) and '_' in dataset:
        ExperimentName = xnat.main(username, password, path, SpecificDataset=dataset)
        pathScan = os.path.join(path, ExperimentName)
    elif len(dataset)==3:
        ExperimentName = xnat.main(username, password, path, dataset)
        pathScan = os.path.join(path, ExperimentName)
    elif dataset == 'load':
        pathScan = os.path.join(path)
    
    filename_log = pathScan +"_"+ datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "LogFile.txt"

    #Available CPU cores
    try: 
        UsedCores = int(len(os.sched_getaffinity(0)))
    except: 
        UsedCores = int(os.cpu_count())

    database = db.database(path=pathScan)
    database.set_log(filename_log)
    database.log("Analysis of " + pathScan.split('//')[-1] + " has started!")
    database.log("CPU cores: " + str(UsedCores))
    
    # HARMONIZATION

    # steps_core.rename_all_series(database)
    # steps_core.harmonize_dce(database)
    # steps_internal.harmonize_subject_name(database,dataset)
    
    # SEGMENTATION
    
    # steps_core.fetch_dl_models(database)
    # steps_internal.fetch_kidney_masks(database)
    # steps_core.segment_kidneys(database)
    steps_internal.export_segmentations_folder_volumetry_1(database)

    upload_folder_to_drive( os.path.join(database.path() + '_output',database.PatientName), '1-UQX4KIVnMPKHRtVQAixNY2EgqVecZ_I')
