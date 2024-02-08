import pandas as pd
from dbdicom.wrappers import vreg, scipy
from dbdicom.pipelines import input_series
import os, sys, time
import csv
import datetime
import pipelines.segment as seg
import scripts.QC_alignment as export_alignemnt


def perfusion(database,master_table):

    export_study = "ASL"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'ASL_kidneys_pCASL_cor-oblique_fb_M0_moco',
        'ASL_kidneys_pCASL_cor-oblique_fb_RBF_moco',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'ASL_kidneys_pCASL_cor-oblique_fb_M0_moco',
        'ASL_kidneys_pCASL_cor-oblique_fb_RBF_moco',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return master_table
    dixon   = series[0]
    asl_ref = series[1]
    rbf     = series[2]
    lk      = series[3]
    rk      = series[4]

    # Perform a separate registration for each target region
    rbf_moved = []
    rbf_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_rigid_transformation(asl_ref, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(rbf, params, description='RBF - ' + kidney[1])

        export_alignemnt.main(moved, kidney[0],export_study,database.path())
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

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

def DTI(database,master_table):

    export_study = "RD"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_MD_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_MD_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_RD_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_AD_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_Planarity_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_Linearity_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_Sphericity_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_FA_Map',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_MD_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_MD_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_RD_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_AD_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_Planarity_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_Linearity_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_Sphericity_Map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_DTI_FA_Map',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return master_table
    dixon   = series[0]
    dti_ref = series[1]
    MD_map  = series[2]
    RD_map  = series[3]
    AD_map  = series[4]
    planarity_map  = series[5]
    linearity_map  = series[6]
    sphericity_map = series[7]
    fa_map  = series[8]
    lk      = series[9]
    rk      = series[10]

    # Perform a separate registration for each target region
    MD_map_moved = []
    MD_map_stats = []
    RD_map_moved = []
    RD_map_stats = []
    AD_map_moved = []
    AD_map_stats = []
    planarity_map_moved = []
    planarity_map_stats = []
    linearity_map_moved = []
    linearity_map_stats = []
    sphericity_map_moved = []
    sphericity_map_stats = []
    fa_map_moved = []
    fa_map_stats = []

    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        # Perform coregistration based on m0
        params = vreg.find_rigid_transformation(MD_map, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(MD_map, params, description='MD - ' + kidney[1])
        export_alignemnt.main(moved, kidney[0],export_study,database.path())
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

        # Organise results
        moved.move_to(study)
        MD_map_moved.append(moved)
        MD_map_stats.append(df)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(RD_map, params, description='RD - ' + kidney[1])
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

        # Organise results
        moved.move_to(study)
        RD_map_moved.append(moved)
        RD_map_stats.append(df)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(AD_map, params, description='AD - ' + kidney[1])
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

        # Organise results
        moved.move_to(study)
        AD_map_moved.append(moved)
        AD_map_stats.append(df)
    
        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(planarity_map, params, description='Planarity - ' + kidney[1])
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

        # Organise results
        moved.move_to(study)
        planarity_map_moved.append(moved)
        planarity_map_stats.append(df)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(linearity_map, params, description='Linearity - ' + kidney[1])
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

        # Organise results
        moved.move_to(study)
        linearity_map_moved.append(moved)
        linearity_map_stats.append(df)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(sphericity_map, params, description='Sphericity - ' + kidney[1])
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

        # Organise results
        moved.move_to(study)
        sphericity_map_moved.append(moved)
        sphericity_map_stats.append(df)

        # Apply transformation to rbf image
        moved = vreg.apply_rigid_transformation(fa_map, params, description='FA - ' + kidney[1])
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

        # Organise results
        moved.move_to(study)
        fa_map_moved.append(moved)
        fa_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)
    MD_map_table = pd.concat(RD_map_stats)
    RD_map_table = pd.concat(RD_map_stats)
    AD_map_table = pd.concat(RD_map_stats)
    planarity_map_table = pd.concat(RD_map_stats)
    linearity_map_table = pd.concat(RD_map_stats)
    sphericity_map_table = pd.concat(RD_map_stats)
    fa_map_table = pd.concat(RD_map_stats)
    master_table = pd.concat([master_table,MD_map_table, RD_map_table, AD_map_table, planarity_map_table, linearity_map_table, sphericity_map_table, fa_map_table])

    return master_table


def MTR(database,master_table):

    export_study = "MTR"

    # Get input parameters
    series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'MT_ON_kidneys_cor-oblique_bh',
        'MT_ON_kidneys_cor-oblique_bh_mdr_moco_MTR',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'MT_ON_kidneys_cor-oblique_bh',
        'MT_ON_kidneys_cor-oblique_bh_mdr_moco_MTR',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return master_table
    dixon = series[0]
    MT_ref = series[1]
    MTR_map = series[2]
    lk = series[3]
    rk = series[4]

    # Perform a separate registration for each target region
    MTR_map_moved = []
    MTR_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        #params = vreg.find_rigid_transformation(MTR_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        params = vreg.find_sbs_rigid_transformation(MT_ref, dixon, tolerance=0.1)
        # Apply transformation to rbf image
        moved = vreg.apply_sbs_passive_rigid_transformation(MTR_map, params)
        export_alignemnt.main(moved, kidney[0],export_study,database.path())
        #moved = vreg.apply_rigid_transformation(MTR_map, params, target=dixon, description='MTR - ' + kidney[1])

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

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
        'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_par_S0',
        '2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2s_T2s_Map',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2s_fw_Map',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_par_S0',
        '2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2s_T2s_Map',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2s_fw_Map',
        'LK',
        'RK',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return master_table
    dixon = series[0]
    T2s_ref = series[1]
    T2s_map = series[2]
    T2s_fw = series[3]
    lk = series[4]
    rk = series[5]

    # Perform a separate registration for each target region
    T2s_map_moved = []
    T2s_map_stats = []
    T2s_water_fraction_map_moved = []
    T2s_water_fraction_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        #params = vreg.find_rigid_transformation(T2s_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        params = vreg.find_sbs_rigid_transformation(T2s_ref, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        moved = vreg.apply_sbs_passive_rigid_transformation(T2s_map, params)
        export_alignemnt.main(moved, kidney[0],export_study,database.path())
        #moved = vreg.apply_rigid_transformation(T2s_map, params, target=dixon, description='T2s - ' + kidney[1])

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

        # Organise results
        moved.move_to(study)
        T2s_map_moved.append(moved)
        T2s_map_stats.append(df)

        # Apply transformation to rbf image
        moved = vreg.apply_sbs_passive_rigid_transformation(T2s_fw, params)
        #moved = vreg.apply_rigid_transformation(T2s_map, params, target=dixon, description='T2s - ' + kidney[1])

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

        # Organise results
        moved.move_to(study)
        T2s_water_fraction_map_moved.append(moved)
        T2s_water_fraction_map_stats.append(df)

    # Keep a copy of the kidneys in the export folder for overlay.
    #lk.copy_to(study)
    #rk.copy_to(study)

    T2s_map_table = pd.concat(T2s_map_stats)
    T2s_fw_map_table = pd.concat(T2s_water_fraction_map_stats)
    master_table = pd.concat([master_table,T2s_map_table,T2s_fw_map_table])

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
            return master_table
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
        export_alignemnt.main(moved, kidney[0],export_study,database.path())
        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

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
        'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_par_a',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T1_T1_Map_v2',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_par_a',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T1_T1_Map_v2',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return master_table
    dixon = series[0]
    T1_ref = series[1]
    T1_map = series[2]
    lk = series[3]
    rk = series[4]

    # Perform a separate registration for each target region
    T1_map_moved = []
    T1_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        #params = vreg.find_rigid_transformation(T1_map, dixon, tolerance=0.1, region=kidney[0], margin=0)
        params = vreg.find_sbs_rigid_transformation(T1_ref, dixon, tolerance=0.1, region=kidney[0], margin=0)
        # Apply transformation to rbf image

        #moved = vreg.apply_rigid_transformation(T1_map, params, target=dixon, description='T1 - ' + kidney[1])
        moved = vreg.apply_sbs_passive_rigid_transformation(T1_map, params)
        export_alignemnt.main(moved, kidney[0],export_study,database.path())

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

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
        'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast',
        'T2map_kidneys_cor-oblique_mbh_magnitude_mdr_par_S0',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2_T2_Map_v2',
        'LK',
        'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        series_desc = [   
        'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast',
        'T2map_kidneys_cor-oblique_mbh_magnitude_mdr_par_S0',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2_T2_Map_v2',
        'LK [overlay]',
        'RK [overlay]',
    ]
        series, study = input_series(database, series_desc, export_study)
        if series is None:
            return master_table
    dixon = series[0]
    T2_ref = series[1]
    T2_map = series[2]
    lk = series[3]
    rk = series[4]

    # Perform a separate registration for each target region
    T2_map_moved = []
    T2_map_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_sbs_rigid_transformation(T2_ref, dixon, tolerance=0.1, region=kidney[0], margin=0)

        # Apply transformation to rbf image
        #moved = vreg.apply_rigid_transformation(T2_map, params, target=dixon, description='T2 - ' + kidney[1])
        moved = vreg.apply_sbs_passive_rigid_transformation(T2_map, params)
        export_alignemnt.main(moved, kidney[0],export_study,database.path())

        # Get ROI statistics
        df = scipy.mask_statistics(kidney[0], moved, export_study)

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

def main(master_table, folder):

    start_time = time.time()
    folder.scan()

    try:
        folder.log("ASL export has started")
        master_table = perfusion(folder,master_table)
        folder.log("ASL export was completed")
    except Exception as e: 
        folder.log("ASL export was NOT completed; error: "+str(e))

    try:
        folder.log("FA export has started")
        master_table = DTI(folder,master_table)
        folder.log("FA export was completed")
    except Exception as e: 
        folder.log("FA export was NOT completed; error: "+str(e))

    try:
        folder.log("MTR export has started")
        master_table = MTR(folder,master_table)
        folder.log("MTR export was completed")
    except Exception as e: 
        folder.log("MTR export was NOT completed; error: "+str(e))

    try:
        folder.log("T2* export has started")
        master_table = T2s(folder,master_table)
        folder.log("T2* export was completed")
    except Exception as e: 
        folder.log("T2* export was NOT completed; error: "+str(e))

    try:
        folder.log("DCE TE export has started")
        master_table = DCE_TE(folder,master_table)
        folder.log("DCE TE export was completed")
    except Exception as e: 
        folder.log("DCE TE export was NOT completed; error: "+str(e))

    try:
        folder.log("T1 export has started")
        master_table = T1(folder,master_table)
        folder.log("T1 export was completed")
    except Exception as e: 
        folder.log("T1 export was NOT completed; error: "+str(e))

    try:
        folder.log("T2 export has started")
        master_table = T2(folder,master_table)
        folder.log("T2 export was completed")
    except Exception as e: 
        folder.log("T2 export was NOT completed; error: "+str(e))

    folder.save()
    folder.log("Export was completed --- %s seconds ---" % (int(time.time() - start_time)))

    return master_table