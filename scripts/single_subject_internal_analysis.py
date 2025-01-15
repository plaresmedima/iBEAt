import os
import datetime
import dbdicom as db

from utilities import xnat, upload
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
    steps_core.harmonize_pc(database)
    steps_core.harmonize_t2(database)
    steps_core.harmonize_mt(database)
    steps_core.harmonize_dti(database)
    steps_core.harmonize_ivim(database)
    steps_core.harmonize_dce(database)
    steps_internal.harmonize_t1_t2(database)
    steps_core.harmonize_subject_name(database)
    
    # SEGMENTATION

    steps_core.fetch_dl_models(database)
    steps_internal.fetch_kidney_masks(database)
    steps_core.segment_kidneys(database)
    steps_core.segment_renal_sinus_fat(database)
    steps_core.segment_aorta_on_dce(database)
    steps_core.segment_left_renal_artery(database)
    steps_core.segment_right_renal_artery(database)
    steps_core.compute_whole_kidney_canvas(database)
    steps_core.export_segmentations(database) 

    # MOTION CORRECTION

    steps_core.mdreg_t1(database)
    steps_core.mdreg_t2(database)
    steps_core.mdreg_t2star(database)
    steps_core.mdreg_mt(database)
    steps_core.mdreg_ivim(database)
    steps_core.mdreg_dti(database)
    steps_core.mdreg_dce(database)
    steps_internal.mdreg_t1_t2(database)
    steps_core.export_mdreg(database)
    

    # MAPPING

    steps_core.map_T1(database)
    steps_internal.map_T1_from_T1_T2_mdr(database)
    steps_internal.map_T2_from_T1_T2_mdr(database)
    steps_core.map_T2(database)
    steps_core.map_T2star(database)
    steps_core.map_MT(database)
    steps_core.map_IVIM(database)
    steps_core.map_DTI(database)
    steps_core.map_DCE(database)
    steps_core.export_mapping(database)

    # ALIGNMENT

    steps_core.align_dixon(database) 
    steps_core.align_T1(database)
    steps_core.align_T2(database)
    steps_core.align_T2star(database)
    steps_core.align_MT(database)
    steps_core.align_IVIM(database)
    steps_core.align_DTI(database)
    steps_core.align_DCE(database)
    steps_core.align_ASL(database)
    steps_core.export_alignment(database)

    # MEASUREMENT

    steps_core.fill_gaps(database)
    steps_core.cortex_medulla(database)
    
    steps_core.measure_kidney_volumetrics(database)
    steps_core.measure_sinus_fat_volumetrics(database)
    steps_core.measure_t1_maps(database)
    steps_core.measure_t2_maps(database)
    steps_core.measure_t2star_maps(database)
    steps_core.measure_mt_maps(database)
    steps_core.measure_ivim_maps(database)
    steps_core.measure_dti_maps(database)
    steps_core.measure_asl_maps(database)
    steps_core.measure_dce_maps(database)

    # ROI analysis

    steps_core.roi_fit_T1(database)
    steps_core.roi_fit_T2(database)
    steps_core.roi_fit_T2star(database)
    steps_core.roi_fit_PC(database)
    steps_core.roi_fit_DCE(database)
    steps_core.roi_fit_DCE_cm(database)
    steps_core.roi_fit_IVIM(database)
    steps_core.roi_fit_DTI(database)

        
    #upload images, logfile and csv to google drive
    filename_csv = os.path.join(database.path() + '_output',database.PatientName[0] + '_biomarkers.csv')
    upload.main(pathScan, filename_log, filename_csv)
