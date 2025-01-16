import os
import datetime
import dbdicom as db

from utilities import xnat, upload_project_sinus_fat
from scripts import steps_core, steps_internal


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
    
    filename_log = pathScan +"_"+ datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "MDRauto_LogFile.txt" #TODO FIND ANOTHER WAY TO GET A PATH

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

    steps_core.rename_all_series(database)
    steps_core.harmonize_subject_name(database)
    
    # SEGMENTATION

    steps_core.fetch_dl_models(database)
    steps_internal.fetch_kidney_masks(database)
    steps_core.segment_kidneys(database)
    steps_core.segment_renal_sinus_fat(database)
    steps_internal.export_sinus_fat_segmentations(database)

    # MEASUREMENT
    
    steps_core.measure_kidney_volumetrics(database)
    steps_core.measure_sinus_fat_volumetrics(database)

    #upload images, logfile and csv to google drive
    filename_csv = os.path.join(database.path() + '_output',database.PatientName[0] + '_biomarkers.csv')
    upload_project_sinus_fat.main(pathScan, filename_log, filename_csv)