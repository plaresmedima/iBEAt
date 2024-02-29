import os
import datetime
import dbdicom as db

import pipelines.mdr as mdr
import pipelines.mapping as map
import pipelines.export_ROI_stats as export_ROIs
import pipelines.DCE_analysis as DCE_analysis

import scripts.upload as upload
import scripts.QC_mdr as check_mdr
import scripts.QC_mapping as check_maps
from scripts import xnat, steps


def single_subject(username, password, path, dataset):
    
    #Import data from XNAT
    if isinstance(dataset,str) and '_' in dataset:
        ExperimentName = xnat.main(username, password, path, SpecificDataset=dataset)
    elif len(dataset)==3:
        ExperimentName = xnat.main(username, password, path, dataset)
    elif dataset == 'load':
        ExperimentName = os.path.basename(path)
    
    pathScan = os.path.join(path, ExperimentName)
    filename_log = pathScan +"_"+ datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "MDRauto_LogFile.txt" #TODO FIND ANOTHER WAY TO GET A PATH
    
    # THIS NEEDS ANOTHER APPROACH!!!
    unetr = 'C:\\Users\\steve\\Dropbox\\Data\\dl_models\\UNETR_kidneys_v1.pth'

    #Available CPU cores
    try: 
        UsedCores = int(len(os.sched_getaffinity(0)))
    except: 
        UsedCores = int(os.cpu_count())

    folder = db.database(path=pathScan)
    folder.set_log(filename_log)
    folder.log("Analysis of " + pathScan.split('//')[-1] + " has started!")
    folder.log("CPU cores: " + str(UsedCores))
    

    ## Harmonize series descriptions

    steps.rename_all_series(folder)


    ## Segmentation steps

    # Calculate kidney and renal sinus fat masks
    steps.fetch_kidney_masks(folder) # Dummy placeholder for now!!
    steps.segment_kidneys(folder, unetr)
    steps.segment_renal_sinus_fat(folder)

    # Export kidney masks and canvas for manual editing offline
    steps.compute_whole_kidney_canvas(folder)
    steps.export_kidney_segmentations(folder)
    
    # Generate biomarker measurements
    steps.measure_kidney_volumetrics(folder)
    steps.measure_sinus_fat_volumetrics(folder)

    #Apply motion correction using MDR
    try:
        print("starting mdr")
        mdr.main(folder)
        check_mdr.main(folder)
    except Exception as e:
        folder.log("Renaming was NOT completed; error: " + str(e))

    #Custom modelling
    try:
        print("staring mapping")
        map.main(folder)
        check_maps.main(folder)
    except Exception as e:
        folder.log("Modelling was NOT completed; error: " + str(e))

    #Custom DCE
    try:
        print("staring DCE Analysis")
        DCE_analysis.main(folder,master_table)
    except Exception as e:
        folder.log("DCE Analysis was NOT completed; error: " + str(e))


    #Generate masks using unetr, apply alignment, extract biomarkers to a .csv
    try:
        print('starting parameter extraction')
        master_table = export_ROIs.main(master_table, folder)
    except Exception as e:
        folder.log("Parameter extraction was NOT completed; error: " + str(e))



    
    #upload images, logfile and csv to google drive
    #upload.main(pathScan, filename_log, filename_csv)
    filename_csv = os.path.join(folder.path() + '_output', 'biomarkers.csv')
    upload.main(pathScan, filename_log, filename_csv)
