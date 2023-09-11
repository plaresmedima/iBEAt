import pandas as pd
from dbdicom.wrappers import vreg, scipy
from dbdicom.pipelines import input_series
import os, sys, time
import csv
import datetime

def perfusion(database,master_table):

    export_study = "ASL"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'ASL_kidneys_pCASL_cor-oblique_fb_RBF_moco',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        return None, None
    dixon = series[0]
    rbf = series[1]
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    rbf_moved = []
    rbf_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_rigid_transformation(rbf, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(rbf, params, target=dixon, description='RBF - ' + kidney[1])

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        rbf_moved.append(moved)
        rbf_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)
    rbf_table = pd.concat(rbf_stats)
    master_table = pd.concat([master_table,rbf_table])

    return master_table

def fa(database,master_table):

    export_study = "FA"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_FA_Map',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        return None, None
    dixon = series[0]
    fa_map = series[1]
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    fa_map_moved = []
    fa_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_rigid_transformation(fa_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(fa_map, params, target=dixon, description='FA - ' + kidney[1])

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        fa_map_moved.append(moved)
        fa_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    fa_map_table = pd.concat(fa_map_stats)
    master_table = pd.concat([master_table,fa_map_table])

    return master_table

def ADC(database,master_table):

    export_study = "ADC"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_MD_Map',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        return None, None
    dixon = series[0]
    ADC_map = series[1]
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    ADC_map_moved = []
    ADC_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_rigid_transformation(ADC_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(ADC_map, params, target=dixon, description='ADC - ' + kidney[1])

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        ADC_map_moved.append(moved)
        ADC_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    ADC_map_table = pd.concat(ADC_map_stats)
    master_table = pd.concat([master_table,ADC_map_table])

    return master_table

def MTR(database,master_table):

    export_study = "MTR"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'MT_ON_kidneys_cor-oblique_bh_mdr_moco_MTR',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        return None, None
    dixon = series[0]
    MTR_map = series[1]
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    MTR_map_moved = []
    MTR_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_rigid_transformation(MTR_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(MTR_map, params, target=dixon, description='MTR - ' + kidney[1])

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        MTR_map_moved.append(moved)
        MTR_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    MTR_map_table = pd.concat(MTR_map_stats)
    master_table = pd.concat([master_table,MTR_map_table])

    return master_table

def T2s(database,master_table):

    export_study = "T2s"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        '2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2s_T2s_Map',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        return None, None
    dixon = series[0]
    T2s_map = series[1]
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    T2s_map_moved = []
    T2s_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_rigid_transformation(T2s_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(T2s_map, params, target=dixon, description='T2s - ' + kidney[1])

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        T2s_map_moved.append(moved)
        T2s_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    T2s_map_table = pd.concat(T2s_map_stats)
    master_table = pd.concat([master_table,T2s_map_table])

    return master_table

def main(folder,ExperimentName):

    start_time = time.time()
    folder.log(": ROI Stats has started!")

    master_table = pd.DataFrame()
    try:
        print('ASL export started')
        master_table = perfusion(folder,master_table)
    except Exception as e: 
        folder.log("ASL export was NOT completed; error: "+str(e))

    try:
        print('FA export started')
        master_table = fa(folder,master_table)
    except Exception as e: 
        folder.log("FA export was NOT completed; error: "+str(e))

    try:
        print('ADC export started')
        master_table = ADC(folder,master_table)
    except Exception as e: 
        folder.log("ADC export was NOT completed; error: "+str(e))

    try:
        print('MTR export started')
        master_table = MTR(folder,master_table)
    except Exception as e: 
        folder.log("MTR export was NOT completed; error: "+str(e))

    try:
        print('T2s export started')
        master_table = T2s(folder,master_table)
    except Exception as e: 
        folder.log("T2s export was NOT completed; error: "+str(e))

    filename_csv = datetime.datetime.now().strftime('%Y%m%d_%H%M_') + ExperimentName+'.csv'
    master_table.to_csv(filename_csv, index=False)

    folder.save()
    folder.log("Export was completed --- %s seconds ---" % (int(time.time() - start_time)))