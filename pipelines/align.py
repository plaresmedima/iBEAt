import numpy as np
from dbdicom.extensions import vreg, elastix
from dbdicom.pipelines import input_series
from dbdicom.utils import vreg as vr
from dipy.align.imaffine import MutualInformationMetric, AffineRegistration
from dipy.align.transforms import TranslationTransform2D


export_study = "Alignment"


def _inslice_active_translation(fixed, moving, progress, applyto=[]):
    metric = MutualInformationMetric(nbins=32, sampling_proportion=None)
    affreg = AffineRegistration(
        metric = metric,
        level_iters = [10000, 1000, 100],
        sigmas = [3.0, 1.0, 0.0],
        factors = [4, 2, 1])
    transform = TranslationTransform2D()
    for z in range(moving.shape[2]):
        progress(z+1, moving.shape[2], 'Performing in-slice coregistration..')
        params0 = None
        mapping = affreg.optimize(fixed[:,:,z,0], moving[:,:,z,0], transform, params0)
        moving[:,:,z,0] = mapping.transform(moving[:,:,z,0], 'linear')
        for arr in applyto:
            for t in range(arr.shape[3]):
                arr[:,:,z,t] = mapping.transform(arr[:,:,z,t], 'linear')
    return moving


def _sbs_rigid(array_static, affine_static, array_moving, affine_moving, progress, slice_thickness):
    _, _, static_pixel_spacing = vr.affine_components(affine_static)
    rot_gradient_step, translation_gradient_step, _ = vr.affine_resolution(array_static.shape, static_pixel_spacing)
    gradient_step = np.concatenate((1.0*rot_gradient_step, 0.5*translation_gradient_step))
    optimization = {
        'method': 'GD', 
        'options': {'gradient step': gradient_step, 'tolerance': 0.1}, 
    }
    try:
        parameters = vr.align_slice_by_slice(
            moving = array_moving, 
            moving_affine = affine_moving, 
            static = array_static, 
            static_affine = affine_static, 
            parameters = np.zeros(6, dtype=np.float32), 
            resolutions = [4,2,1],
            transformation = vr.rigid,
            metric = vr.mutual_information,
            optimization = optimization,
            slice_thickness = slice_thickness,
            progress = lambda z, nz: progress(z+1, nz, 'Coregistering slice-by-slice using rigid transformations'),
        )
    except:
        print('Failed to align volumes..')
        parameters = np.zeros(6, dtype=np.float32)

    return parameters


def _align(database, desc):

    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform '+desc[2]+' alignment: not all required data exist.')

    moving = series[2]
    dims = ('SliceLocation','InstanceNumber')
    affine_moving = moving.affine()
    slice_thickness = moving.values('SliceThickness')[0]
    
    moved = []
    for kidney in [0,1]:
        moving.message('Coregistering to kidney ' + desc[kidney])

        # Inslice translation
        applyto = [s.pixel_values(dims) for s in series[3:]]
        array_moving = moving.pixel_values(dims)
        array_static = vreg.pixel_values(series[kidney], dims, on=moving)
        array_moving = _inslice_active_translation(array_static, array_moving, moving.progress, applyto)

        # Through-slice rigid - find new affines for each slice
        array_static = series[kidney].pixel_values('SliceLocation')
        affine_static = series[kidney].affine()
        array_static, affine_static = vr.mask_volume(array_static, affine_static, array_static, affine_static, 0)
        parameters = _sbs_rigid(array_static, affine_static, array_moving[...,0], affine_moving, moving.progress, slice_thickness)
        affines_moving = vr.passive_rigid_transform_slice_by_slice(affine_moving, parameters)

        # Save as DICOM
        map_arrays = [array_moving] + applyto
        for m, map_series in enumerate(series[2:]):
            moved_map = map_series.new_sibling(SeriesDescription = desc[2+m] + '_' + desc[kidney] + '_align')
            frames = map_series.frames(dims)
            cnt=0
            for z in range(frames.shape[0]):
                for t in range(frames.shape[1]):
                    cnt+=1
                    map_series.progress(cnt, frames.size, 'Applying transformation.. ')
                    affine_zt = vr.multislice_to_singleslice_affine(affines_moving[z], frames[z,t].SliceThickness)
                    frames_zt = frames[z,t].copy_to(moved_map)
                    frames_zt.set_affine(affine_zt)
                    frames_zt.set_pixel_array(map_arrays[m][:,:,z,t])
            moved_map.move_to(study)
            moved.append(moved_map)
    return moved[:2]


def t1(database):

    desc = [ 
        'LK',
        'RK',  
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T1_map',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_S0_map',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_FA_map',
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco',
   ]
    return _align(database, desc)
    
    
def t2(database):

    desc = [ 
        'LK',
        'RK',  
        'T2map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2_map',
        'T2map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_S0_map',
        'T2map_kidneys_cor-oblique_mbh_magnitude_mdr_moco',
    ]
    return _align(database, desc)


def t2star(database):

    desc = [ 
        'LK',
        'RK',  
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_T2star_map',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco_S0_map',
        'T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco'
    ]
    return _align(database, desc)


def ivim(database):

    desc = [ 
        'LK',
        'RK',  
        'IVIM_kidneys_cor-oblique_fb_mdr_moco_S0_map',
        'IVIM_kidneys_cor-oblique_fb_mdr_moco_MD_map',
        'IVIM_kidneys_cor-oblique_fb_mdr_moco_Df_map',
        'IVIM_kidneys_cor-oblique_fb_mdr_moco_ff_map',
        'IVIM_kidneys_cor-oblique_fb_mdr_moco',
    ]
    return _align(database, desc)


def dti(database):

    desc = [ 
        'LK',
        'RK',  
        'DTI_kidneys_cor-oblique_fb_mdr_moco_MD_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_RD_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_AD_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_Planarity_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_Linearity_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_Sphericity_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco_FA_map',
        'DTI_kidneys_cor-oblique_fb_mdr_moco',
    ]
    return _align(database, desc)


def dce(database):

    desc = [ 
        'LK',
        'RK',  
        'DCE_kidneys_cor-oblique_fb_mdr_moco_AUC_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco_AVD_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco_RPF_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco_MTT_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco',
    ]
    return _align(database, desc)


# TODO: MT is 3D! - needs the same approach as ASL
def mt(database):

    desc = [ 
        'LK',
        'RK',  
        'MT_kidneys_cor-oblique_bh_mdr_moco_AVR',
        'MT_kidneys_cor-oblique_bh_mdr_moco_MTR',
    ]
    return _align(database, desc)


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


def dixon(database):

    # Get input parameters
    desc = [   
        'T1w_abdomen_dixon_cor_bh_out_phase',
        'T1w_abdomen_dixon_cor_bh_in_phase',
        'T1w_abdomen_dixon_cor_bh_water',
        'T1w_abdomen_dixon_cor_bh_fat',
        'T1w_abdomen_dixon_cor_bh_fat_post_contrast',
    ]
    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform DIXON alignment: not all required data exist.')

    coregistered, followers = elastix.coregister_3d_to_3d(
        series[3], series[4],
        transformation = "AdvancedMeanSquares",
        metric = "AdvancedMattesMutualInformation",
        final_grid_spacing = 25.0,
        apply_to = series[:3],
    )
    coregistered.move_to(study)
    for series in followers:
        series.move_to(study)
    return coregistered, followers


def fill_gaps(database):

    to_fill = [
        'DCE_kidneys_cor-oblique_fb_mdr_moco_AVD_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco_RPF_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco_MTT_map',        
    ]
    ref = ['T1w_abdomen_dixon_cor_bh_fat_post_contrast']

    output = []
    for kidney in ['LK','RK']:
        desc = ref + [kidney] + [d + '_' + kidney + '_align' for d in to_fill]
        series, study = input_series(database, desc, export_study)
        if series is None:
            raise RuntimeError('Cannot fill gaps between slices on LK: not all required data exist.')
        for ms_series in series[2:]:
            output_series = vreg.fill_slice_gaps(ms_series, series[0], slice_thickness=1.0, mask=series[1])
            output_series.move_to(study)
            output.append(output_series)

    return output











