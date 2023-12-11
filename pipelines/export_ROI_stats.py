import pandas as pd
from dbdicom.wrappers import vreg, scipy
from dbdicom.pipelines import input_series
import os, sys, time
import csv
import datetime
import pipelines.segment as seg

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
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'ASL_kidneys_pCASL_cor-oblique_fb_RBF_moco',
        'LK [overlay]',
        'RK [overlay]',
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
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_FA_Map',
        'LK [overlay]',
        'RK [overlay]',
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
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_MD_Map',
        'LK [overlay]',
        'RK [overlay]',
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
        #params = vreg.find_sbs_rigid_transformation(ADC_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        # Apply transformation to rbf image
       # moved = vreg.apply_sbs_passive_rigid_transformation(ADC_map, params)
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
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'MT_ON_kidneys_cor-oblique_bh_mdr_moco_MTR',
        'LK [overlay]',
        'RK [overlay]',
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
        #params = vreg.find_rigid_transformation(MTR_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        params = vreg.find_sbs_rigid_transformation(MTR_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        # Apply transformation to rbf image
        moved = vreg.apply_sbs_passive_rigid_transformation(MTR_map, params)
        #moved = vreg.apply_rigid_transformation(MTR_map, params, target=dixon, description='MTR - ' + kidney[1])

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
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        '2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2s_T2s_Map',
        'LK [overlay]',
        'RK [overlay]',
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
        #params = vreg.find_rigid_transformation(T2s_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        params = vreg.find_sbs_rigid_transformation(T2s_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        moved = vreg.apply_sbs_passive_rigid_transformation(T2s_map, params)
        #moved = vreg.apply_rigid_transformation(T2s_map, params, target=dixon, description='T2s - ' + kidney[1])

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

def T2s_water_fraction(database,master_table):

    export_study = "T2s water fraction"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2s_fw_Map',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2s_fw_Map',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return None, None
    dixon = series[0]
    T2s_water_fraction_map = series[1]
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    T2s_water_fraction_map_moved = []
    T2s_water_fraction_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        #params = vreg.find_rigid_transformation(T2s_water_fraction_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        params = vreg.find_sbs_rigid_transformation(T2s_water_fraction_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        # Apply transformation to rbf image
        moved = vreg.apply_sbs_passive_rigid_transformation(T2s_water_fraction_map, params)
        #moved = vreg.apply_rigid_transformation(T2s_water_fraction_map, params, target=dixon, description='T2s water fraction - ' + kidney[1])

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        T2s_water_fraction_map_moved.append(moved)
        T2s_water_fraction_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    T2s_water_fraction_map_table = pd.concat(T2s_water_fraction_map_stats)
    master_table = pd.concat([master_table,T2s_water_fraction_map_table])

    return master_table

def DCE_FP(database,master_table):

    export_study = "DCE FP"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DCE_kidneys_cor-oblique_fb_mdr_par_FP',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DCE_kidneys_cor-oblique_fb_mdr_par_FP',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return None, None
    dixon = series[0]
    DCE_FP_map = series[1]
    DCE_FP_map = DCE_FP_map.subseries(InPlanePhaseEncodingDirection='ROW')
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    DCE_FP_map_moved = []
    DCE_FP_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_sbs_rigid_transformation(DCE_FP_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        #params = vreg.find_rigid_transformation(DCE_FP_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        #moved = vreg.apply_rigid_transformation(DCE_FP_map, params, target=dixon, description='DCE FP - ' + kidney[1])
        moved = vreg.apply_sbs_passive_rigid_transformation(DCE_FP_map, params)
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        DCE_FP_map_moved.append(moved)
        DCE_FP_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    DCE_FP_map_table = pd.concat(DCE_FP_map_stats)
    master_table = pd.concat([master_table,DCE_FP_map_table])

    return master_table

def DCE_TP(database,master_table):

    export_study = "DCE TP"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DCE_kidneys_cor-oblique_fb_mdr_par_TP',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DCE_kidneys_cor-oblique_fb_mdr_par_TP',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return None, None
    dixon = series[0]
    DCE_TP_map = series[1]
    DCE_TP_map = DCE_TP_map.subseries(InPlanePhaseEncodingDirection='ROW')
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    DCE_TP_map_moved = []
    DCE_TP_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_sbs_rigid_transformation(DCE_TP_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        #params = vreg.find_rigid_transformation(DCE_TP_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        #moved = vreg.apply_rigid_transformation(DCE_TP_map, params, target=dixon, description='DCE TP - ' + kidney[1])
        moved = vreg.apply_sbs_passive_rigid_transformation(DCE_TP_map, params)
        

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        DCE_TP_map_moved.append(moved)
        DCE_TP_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    DCE_TP_map_table = pd.concat(DCE_TP_map_stats)
    master_table = pd.concat([master_table,DCE_TP_map_table])

    return master_table

def DCE_PS(database,master_table):

    export_study = "DCE PS"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DCE_kidneys_cor-oblique_fb_mdr_par_PS',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DCE_kidneys_cor-oblique_fb_mdr_par_PS',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return None, None
    dixon = series[0]
    DCE_PS_map = series[1]
    DCE_PS_map = DCE_PS_map.subseries(InPlanePhaseEncodingDirection='ROW')
    
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    DCE_PS_map_moved = []
    DCE_PS_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_sbs_rigid_transformation(DCE_PS_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        #params = vreg.find_rigid_transformation(DCE_PS_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        #moved = vreg.apply_rigid_transformation(DCE_PS_map, params, target=dixon, description='DCE PS - ' + kidney[1])
        moved = vreg.apply_sbs_passive_rigid_transformation(DCE_PS_map, params)

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        DCE_PS_map_moved.append(moved)
        DCE_PS_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    DCE_PS_map_table = pd.concat(DCE_PS_map_stats)
    master_table = pd.concat([master_table,DCE_PS_map_table])

    return master_table

def DCE_TE(database,master_table):

    export_study = "DCE TE"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DCE_kidneys_cor-oblique_fb_mdr_par_TE',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DCE_kidneys_cor-oblique_fb_mdr_par_TE',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return None, None
    dixon = series[0]
    DCE_TE_map = series[1]
    DCE_TE_map = DCE_TE_map.subseries(InPlanePhaseEncodingDirection='ROW')
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    DCE_TE_map_moved = []
    DCE_TE_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        #params = vreg.find_rigid_transformation(DCE_TE_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        params = vreg.find_sbs_rigid_transformation(DCE_TE_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        #moved = vreg.apply_rigid_transformation(DCE_TE_map, params, target=dixon, description='DCE TE - ' + kidney[1])
        moved = vreg.apply_sbs_passive_rigid_transformation(DCE_TE_map, params)
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        DCE_TE_map_moved.append(moved)
        DCE_TE_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    DCE_TE_map_table = pd.concat(DCE_TE_map_stats)
    master_table = pd.concat([master_table,DCE_TE_map_table])

    return master_table

def T1(database,master_table):

    export_study = "T1"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_par_T1',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_par_T1',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return None, None
    dixon = series[0]
    T1_map = series[1]
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    T1_map_moved = []
    T1_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        #params = vreg.find_rigid_transformation(T1_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        params = vreg.find_sbs_rigid_transformation(T1_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        # Apply transformation to rbf image

        #moved = vreg.apply_rigid_transformation(T1_map, params, target=dixon, description='T1 - ' + kidney[1])
        moved = vreg.apply_sbs_passive_rigid_transformation(T1_map, params)


        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        T1_map_moved.append(moved)
        T1_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    T1_map_table = pd.concat(T1_map_stats)
    master_table = pd.concat([master_table,T1_map_table])

    return master_table

def T2(database,master_table):

    export_study = "T2"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'T2map_kidneys_cor-oblique_mbh_magnitude_mdr_par_T2',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'T2map_kidneys_cor-oblique_mbh_magnitude_mdr_par_T2',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return None, None
    dixon = series[0]
    T2_map = series[1]
    lk = series[2]
    rk = series[3]

    # Perform a separate registration for each target region
    T2_map_moved = []
    T2_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_sbs_rigid_transformation(T2_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        #moved = vreg.apply_rigid_transformation(T2_map, params, target=dixon, description='T2 - ' + kidney[1])
        moved = vreg.apply_sbs_passive_rigid_transformation(T2_map, params)

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved)

        # Organise results
        moved.move_to(study)
        T2_map_moved.append(moved)
        T2_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    T2_map_table = pd.concat(T2_map_stats)
    master_table = pd.concat([master_table,T2_map_table])

    return master_table

def main(folder,ExperimentName):

    start_time = time.time()
    folder.scan()

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

    try:
        print('T2s water fraction export started')
        master_table = T2s_water_fraction(folder,master_table)
    except Exception as e: 
        folder.log("T2s water fraction export was NOT completed; error: "+str(e))

    # try:
    #     print('DCE FP water fraction export started')
    #     master_table = DCE_FP(folder,master_table)
    # except Exception as e: 
    #     folder.log("DCE FP export was NOT completed; error: "+str(e))

    # try:
    #     print('DCE TP water fraction export started')
    #     master_table = DCE_TP(folder,master_table)
    # except Exception as e: 
    #     folder.log("DCE TP export was NOT completed; error: "+str(e))

    # try:
    #     print('DCE PS water fraction export started')
    #     master_table = DCE_PS(folder,master_table)
    # except Exception as e: 
    #     folder.log("DCE PS export was NOT completed; error: "+str(e))

    # try:
    #     print('DCE TE water fraction export started')
    #     master_table = DCE_TE(folder,master_table)
    # except Exception as e: 
    #     folder.log("DCE TE export was NOT completed; error: "+str(e))

    try:
        print('T1 export started')
        master_table = T1(folder,master_table)
    except Exception as e: 
        folder.log("T1 export was NOT completed; error: "+str(e))

    try:
        print('T2 export started')
        master_table = T2(folder,master_table)
    except Exception as e: 
        folder.log("T2 export was NOT completed; error: "+str(e))

    filename_csv = datetime.datetime.now().strftime('%Y%m%d_%H%M_') + ExperimentName+'.csv'
    master_table.to_csv(filename_csv, index=False)

    folder.save()
    folder.log("Export was completed --- %s seconds ---" % (int(time.time() - start_time)))

    return filename_csv