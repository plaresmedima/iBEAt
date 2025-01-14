import os
import datetime
import dbdicom as db

from scripts import steps_external


def single_subject(path):
    
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

    steps_external.rename_all_series(database)
    steps_external.harmonize_pc(database)
    steps_external.harmonize_t2(database)
    steps_external.harmonize_mt(database)
    steps_external.harmonize_dti(database)
    steps_external.harmonize_ivim(database)
    steps_external.harmonize_dce(database)
    steps_external.harmonize_subject_name(database)
    
    # SEGMENTATION

    steps_external.fetch_dl_models(database)
    steps_external.segment_kidneys(database)
    steps_external.segment_renal_sinus_fat(database)
    steps_external.segment_aorta_on_dce(database)
    steps_external.segment_left_renal_artery(database)
    steps_external.segment_right_renal_artery(database)
    steps_external.compute_whole_kidney_canvas(database)
    steps_external.export_segmentations(database) 

    # MOTION CORRECTION

    steps_external.mdreg_t1(database)
    steps_external.mdreg_t2(database)
    steps_external.mdreg_t2star(database)
    steps_external.mdreg_mt(database)
    steps_external.mdreg_ivim(database)
    steps_external.mdreg_dti(database)
    steps_external.mdreg_dce(database)
    steps_external.export_mdreg(database)
    

    # MAPPING

    steps_external.map_T1(database)
    steps_external.map_T2(database)
    steps_external.map_T2star(database)
    steps_external.map_MT(database)
    steps_external.map_IVIM(database)
    steps_external.map_DTI(database)
    steps_external.map_DCE(database)
    steps_external.export_mapping(database)

    # ALIGNMENT

    steps_external.align_dixon(database) 
    steps_external.align_T1(database)
    steps_external.align_T2(database)
    steps_external.align_T2star(database)
    steps_external.align_MT(database)
    steps_external.align_IVIM(database)
    steps_external.align_DTI(database)
    steps_external.align_DCE(database)
    steps_external.align_ASL(database)
    steps_external.export_alignment(database)

    # MEASUREMENT

    steps_external.fill_gaps(database)
    steps_external.cortex_medulla(database) 
    steps_external.measure_kidney_volumetrics(database)
    steps_external.measure_sinus_fat_volumetrics(database)
    steps_external.measure_t1_maps(database)
    steps_external.measure_t2_maps(database)
    steps_external.measure_t2star_maps(database)
    steps_external.measure_mt_maps(database)
    steps_external.measure_ivim_maps(database)
    steps_external.measure_dti_maps(database)
    steps_external.measure_asl_maps(database)
    steps_external.measure_dce_maps(database)

    # ROI analysis

    steps_external.roi_fit_T1(database)
    steps_external.roi_fit_T2(database)
    steps_external.roi_fit_T2star(database)
    steps_external.roi_fit_PC(database)
    steps_external.roi_fit_DCE(database)
    steps_external.roi_fit_DCE_cm(database)
    steps_external.roi_fit_IVIM(database)
    steps_external.roi_fit_DTI(database)
