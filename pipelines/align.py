from dbdicom.extensions import vreg
from dbdicom.pipelines import input_series

export_study = "Alignment"


def t1(database):

    # Get input parameters
    desc = [ 
        'LK',
        'RK',  
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_S0_map',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T1_map',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_FA_map',
   ]
    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform T1 alignment: not all required data exist.')

    lk = series[0]
    rk = series[1]    
    dixon = series[2]
    ref = series[3]

    maps_moved = []
    for kidney in [(lk,'LK'), (rk,'RK')]:
        dixon.message('Coregistering to kidney ' + kidney[1])
        params = vreg.find_rigid_transformation(ref, dixon, tolerance=0.1, region=kidney[0])
        for i, map_series in enumerate(series[3:]):
            dixon.progress(i+1, len(series[3:]), 'Aligning maps to kidney ' + kidney[1])
            moved = vreg.apply_rigid_transformation(map_series, params, description=desc[i+3] + '_align_' + kidney[1])
            moved.move_to(study)
            maps_moved.append(moved)
    return maps_moved


def t2(database):

    # Get input parameters
    desc = [ 
        'LK',
        'RK',  
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'T2map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_S0_map',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2_map',
    ]
    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform T2 alignment: not all required data exist.')

    lk = series[0]
    rk = series[1]    
    dixon = series[2]
    ref = series[3]

    maps_moved = []
    for kidney in [(lk,'LK'), (rk,'RK')]:
        dixon.message('Coregistering to kidney ' + kidney[1])
        params = vreg.find_rigid_transformation(ref, dixon, tolerance=0.1, region=kidney[0])
        for i, map_series in enumerate(series[3:]):
            dixon.progress(i+1, len(series[3:]), 'Aligning maps to kidney ' + kidney[1])
            moved = vreg.apply_rigid_transformation(map_series, params, description=desc[i+3] + '_align_' + kidney[1])
            moved.move_to(study)
            maps_moved.append(moved)
    return maps_moved


def t2star(database):

    # Get input parameters
    desc = [ 
        'LK',
        'RK',  
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_S0_map',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2s_map',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_fw_map',
   ]
    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform T2* alignment: not all required data exist.')

    lk = series[0]
    rk = series[1]    
    dixon = series[2]
    ref = series[3]

    maps_moved = []
    for kidney in [(lk,'LK'), (rk,'RK')]:
        dixon.message('Coregistering to kidney ' + kidney[1])
        params = vreg.find_rigid_transformation(ref, dixon, tolerance=0.1, region=kidney[0])
        for i, map_series in enumerate(series[3:]):
            dixon.progress(i+1, len(series[3:]), 'Aligning maps to kidney ' + kidney[1])
            moved = vreg.apply_rigid_transformation(map_series, params, description=desc[i+3] + '_align_' + kidney[1])
            moved.move_to(study)
            maps_moved.append(moved)
    return maps_moved

def mt(database):

    # Get input parameters
    desc = [ 
        'LK',
        'RK',  
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'MT_kidneys_cor-oblique_bh_mdr_moco_AVR',
        'MT_kidneys_cor-oblique_bh_mdr_moco_MTR',
    ]
    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MT alignment: not all required data exist.')

    lk = series[0]
    rk = series[1]    
    dixon = series[2]
    ref = series[3]

    maps_moved = []
    for kidney in [(lk,'LK'), (rk,'RK')]:
        dixon.message('Coregistering to kidney ' + kidney[1])
        params = vreg.find_rigid_transformation(ref, dixon, tolerance=0.1, region=kidney[0])
        for i, map_series in enumerate(series[3:]):
            dixon.progress(i+1, len(series[3:]), 'Aligning maps to kidney ' + kidney[1])
            moved = vreg.apply_rigid_transformation(map_series, params, description=desc[i+3] + '_align_' + kidney[1])
            moved.move_to(study)
            maps_moved.append(moved)
    return maps_moved


def ivim(database):

    # Get input parameters
    desc = [ 
        'LK',
        'RK',  
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'IVIM_kidneys_cor-oblique_fb_mdr_moco_S0_map',
        'IVIM_kidneys_cor-oblique_fb_mdr_moco_MD_map',
        'IVIM_kidneys_cor-oblique_fb_mdr_moco_Df_map',
        'IVIM_kidneys_cor-oblique_fb_mdr_moco_ff_map',
    ]
    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform IVIM alignment: not all required data exist.')

    lk = series[0]
    rk = series[1]    
    dixon = series[2]
    ref = series[3]

    maps_moved = []
    for kidney in [(lk,'LK'), (rk,'RK')]:
        dixon.message('Coregistering to kidney ' + kidney[1])
        params = vreg.find_rigid_transformation(ref, dixon, tolerance=0.1, region=kidney[0])
        for i, map_series in enumerate(series[3:]):
            dixon.progress(i+1, len(series[3:]), 'Aligning maps to kidney ' + kidney[1])
            moved = vreg.apply_rigid_transformation(map_series, params, description=desc[i+3] + '_align_' + kidney[1])
            moved.move_to(study)
            maps_moved.append(moved)
    return maps_moved


def dti(database):

    # Get input parameters
    desc = [ 
        'LK',
        'RK',  
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_MD_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_RD_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_AD_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_Planarity_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_Linearity_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_Sphericity_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_FA_map',
    ]
    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform DTI alignment: not all required data exist.')

    lk = series[0]
    rk = series[1]    
    dixon = series[2]
    ref = series[3]

    maps_moved = []
    for kidney in [(lk,'LK'), (rk,'RK')]:
        dixon.message('Coregistering to kidney ' + kidney[1])
        params = vreg.find_rigid_transformation(ref, dixon, tolerance=0.1, region=kidney[0])
        for i, map_series in enumerate(series[3:]):
            dixon.progress(i+1, len(series[3:]), 'Aligning maps to kidney ' + kidney[1])
            moved = vreg.apply_rigid_transformation(map_series, params, description=desc[i+3] + '_align_' + kidney[1])
            moved.move_to(study)
            maps_moved.append(moved)
    return maps_moved

def dce(database):

    # Get input parameters
    desc = [ 
        'LK',
        'RK',  
        'DCE_kidneys_cor-oblique_fb_mdr_moco_AUC_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco_AVD_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco_RPF_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco_MTT_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco',
    ]
    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform DCE alignment: not all required data exist.')
  
    moving = series[2]
    maps_moved = []
    for kidney in [0,1]:
        moving.message('Coregistering to kidney ' + desc[kidney])
        params = vreg.find_sbs_rigid_transformation(moving, series[kidney], tolerance=0.1, region=series[kidney])
        for i, map_series in enumerate(series[2:]):
            moving.progress(i+1, len(series[2:]), 'Aligning maps to kidney ' + desc[kidney])
            moved = vreg.apply_sbs_passive_rigid_transformation(
                map_series, params, 
                description=desc[i+2] + '_' + desc[kidney] + '_align')
            moved.move_to(study)
            maps_moved.append(moved)
    return maps_moved


def asl(database):

    # Get input parameters
    desc = [   
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'ASL_kidneys_pCASL_cor-oblique_fb_M0_moco',
        'ASL_kidneys_pCASL_cor-oblique_fb_RBF_moco',
        'LK',
        'RK',
    ]
    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform ASL alignment: not all required data exist.')
        
    dixon   = series[0]
    asl_ref = series[1]
    rbf     = series[2]
    lk      = series[3]
    rk      = series[4]

    # Perform a separate registration for each target region
    rbf_moved = []
    for kidney in [(lk,'LK'), (rk,'RK')]:
        params = vreg.find_rigid_transformation(asl_ref, dixon, tolerance=0.1, region=kidney[0], margin=0)
        moved = vreg.apply_rigid_transformation(rbf, params, description='RBF - ' + kidney[1])
        moved.move_to(study)
        rbf_moved.append(moved)

    return rbf_moved








