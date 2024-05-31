import os
import datetime
import dbdicom as db

from utilities import xnat, upload
from scripts import steps


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
    
    # ## HARMONIZATION

    # steps.rename_all_series(database)
    # steps.harmonize_pc(database)
    # steps.harmonize_t2(database)
    # steps.harmonize_mt(database)
    # steps.harmonize_dti(database)
    # steps.harmonize_ivim(database)
    # steps.harmonize_dce(database)
    # steps.harmonize_subject_name(database)
    
    # ## SEGMENTATION

    # steps.fetch_dl_models(database)
    # steps.fetch_kidney_masks(database)
    # steps.segment_kidneys(database)
    # steps.segment_renal_sinus_fat(database)
    # steps.segment_aorta_on_dce(database)
    # steps.segment_renal_artery(database)
    # steps.compute_whole_kidney_canvas(database)
    # steps.export_segmentations(database) 

    # ## MOTION CORRECTION

    # steps.mdreg_t1(database)
    # steps.mdreg_t2(database)
    # steps.mdreg_t2star(database)
    # steps.mdreg_mt(database)
    # #steps.mdreg_ivim(database)
    # steps.mdreg_dti(database)
    # steps.mdreg_dce(database)
    # steps.export_mdreg(database)

    # # MAPPING

    # steps.map_T1(database)
    # steps.map_T2(database)
    # steps.map_T2star(database)
    # steps.map_MT(database)
    # #steps.map_IVIM(database)
    # steps.map_DTI(database)
    # steps.map_DCE(database)
    # steps.export_mapping(database)

    # # ALIGNMENT

    # steps.align_dixon(database) 
    # steps.align_T1(database)
    # steps.align_T2(database)
    # steps.align_T2star(database)
    # steps.align_MT(database)
    # #steps.align_IVIM(database)
    # steps.align_DTI(database)
    # steps.align_DCE(database)
    # #steps.align_ASL(database)
    steps.export_alignment(database)

    # # MEASUREMENT

    # steps.fill_gaps(database)
    # steps.cortex_medulla(database)
    
    # steps.measure_kidney_volumetrics(database)
    # steps.measure_sinus_fat_volumetrics(database)
    # steps.measure_t1_maps(database)
    # steps.measure_t2_maps(database)
    # steps.measure_t2star_maps(database)
    # steps.measure_mt_maps(database)
    # #steps.measure_ivim_maps(database)
    # steps.measure_dti_maps(database)
    # #steps.measure_asl_maps(database)
    # steps.measure_dce_maps(database)

    # # ROI analysis

    # steps.roi_fit_T1(database)
    # steps.roi_fit_T2(database)
    # steps.roi_fit_T2star(database)
    # steps.roi_fit_PC(database)
    # steps.roi_fit_DCE(database)
    # steps.roi_fit_DCE_cm(database)
    # # TODO: IVIM, DTI

        
    #upload images, logfile and csv to google drive
    #upload.main(pathScan, filename_log, filename_csv)
    filename_csv = os.path.join(database.path() + '_output', 'biomarkers.csv')
    upload.main(pathScan, filename_log, filename_csv)
