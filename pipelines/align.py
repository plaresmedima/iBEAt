import numpy as np
from dbdicom.extensions import vreg as vreg_dicom
from dbdicom.extensions import elastix
from dbdicom.pipelines import input_series
import vreg
from dipy.align.imaffine import MutualInformationMetric, AffineRegistration
from dipy.align.transforms import TranslationTransform2D


export_study = "3: Alignment"



def _rigid(array_static, affine_static, array_moving, affine_moving, progress):
    _, _, static_pixel_spacing = vreg.affine_components(affine_static)
    rot_gradient_step, translation_gradient_step, _ = vreg.affine_resolution(array_static.shape, static_pixel_spacing)
    gradient_step = np.concatenate((1.0*rot_gradient_step, 0.5*translation_gradient_step))
    optimization = {
        'method': 'GD', 
        'options': {'gradient step': gradient_step, 'tolerance': 0.1}, 
    }
    try:
        msg = 'Coregistering using rigid transformation'
        parameters = vreg.align(
            moving = array_moving, 
            moving_affine = affine_moving, 
            static = array_static, 
            static_affine = affine_static, 
            parameters = np.zeros(6, dtype=np.float32), 
            resolutions = [4,2,1],
            transformation = vreg.rigid,
            metric = vreg.mutual_information,
            optimization = optimization,
            progress = lambda z, nz: progress(z+1, nz, msg),
        )
    except:
        print('Failed to align volumes..')
        parameters = np.zeros(6, dtype=np.float32)

    return parameters


def _sbs_translation_inslice(
        array_static, affine_static, 
        array_moving, affine_moving, 
        progress, 
        slice_thickness):
    
    _, _, static_pixel_spacing = vreg.affine_components(affine_static)
    _, translation_gradient_step, _ = vreg.affine_resolution(array_static.shape, static_pixel_spacing)
    gradient_step = 0.5*translation_gradient_step[:2]
    optimization = {
        'method': 'GD', 
        'options': {'gradient step': gradient_step, 'tolerance': 0.1}, 
    }
    # optimization = {
    #     'method': 'LS', 
    #     'options': {'bounds':(-20, 20)}
    # }
    # optimization = {
    #     'method': 'brute', 
    #     'options': {'grid':[[-50, 50, 50], [-50, 50, 50]]}
    # }
    # optimization = {
    #     'method': 'ibrute', 
    #     'options': {'bounds':[[-20, 20], [-20, 20]], 'num': 3, 'nit':4}
    # }
    try:
        msg = 'Coregistering slice-by-slice using inslice translations'
        parameters = vreg.align_slice_by_slice(
            moving = array_moving, moving_affine = affine_moving, 
            static = array_static, static_affine = affine_static, 
            parameters = np.zeros(2, dtype=np.float32), 
            resolutions = [4,2,1],
            transformation = vreg.translate_inslice,
            metric = vreg.mutual_information,
            optimization = optimization,
            slice_thickness = slice_thickness,
            progress = lambda z, nz: progress(z+1, nz, msg),
        )
    except:
        print('Failed to align volumes..')
        parameters = [np.zeros(2, dtype=np.float32)] * array_moving.shape[2]

    return parameters


def _align(database, desc):

    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform '+desc[2]+' alignment: not all required data exist.')

    moving = series[3]
    dims = ('SliceLocation','InstanceNumber')
    affine_moving = moving.affine()
    array_moving = moving.pixel_values(dims)
    slice_thickness = moving.values('SliceThickness')[0]
    
    moved = []
    #for kidney in [1,2]:
    for kidney in [1]: # for testing only LK
        moving.message('Coregistering to kidney ' + desc[kidney])

        array_static = series[0].pixel_values('SliceLocation')
        affine_static = series[0].affine()

        # array_static = series[kidney].pixel_values('SliceLocation')
        # affine_static = series[kidney].affine()

        array_mask = series[kidney].pixel_values('SliceLocation')
        affine_mask = series[kidney].affine()

        # Get bounding box around kidney mask
        # TODO This may cause the slice running away if its moved outside the masked volume
        # Include bounds on the translation to prevent this
        array_static, affine_static = vreg.mask_volume(
            array_static, affine_static, 
            array_mask, affine_mask,
            10,
        )

        # Perform inslice translation per slice
        translations = _sbs_translation_inslice(
            array_static, affine_static, 
            array_moving[..., 0], affine_moving, 
            moving.progress, slice_thickness,
        )
        affines_moving_translate = vreg.passive_inslice_translation_slice_by_slice(
            affine_moving, 
            translations,
            slice_thickness,
        ) 

        affines_moving = affines_moving_translate
        # # Perform rigid transformation per slice
        # affines_moving = []
        # for z, affine_moving_z in enumerate(affines_moving_translate):
        #     pz = _rigid(
        #         array_static, affine_static, 
        #         array_moving[..., z, 0], affine_moving_z,
        #         moving.progress,
        #     )
        #     affines_moving.append(
        #         vreg.passive_rigid_transform(affine_moving_z, pz)
        #     )

        # Apply and save as DICOM
        applyto = [s.pixel_values(dims) for s in series[4:]]
        map_arrays = [array_moving] + applyto
        for m, map_series in enumerate(series[3:]):
            desc_series = desc[3+m] + '_' + desc[kidney] + '_align'
            moved_map = map_series.new_sibling(SeriesDescription=desc_series)
            frames = map_series.frames(dims)
            cnt=0
            for z in range(frames.shape[0]):
                for t in range(frames.shape[1]):
                    cnt+=1
                    map_series.progress(cnt, frames.size, 'Saving results.. ')
                    #affine_zt = vreg.multislice_to_singleslice_affine(affines_moving[z], frames[z,t].SliceThickness)
                    affine_zt = affines_moving[z]
                    frames_zt = frames[z,t].copy_to(moved_map)
                    frames_zt.set_affine(affine_zt)
                    frames_zt.set_pixel_array(map_arrays[m][:,:,z,t])
            moved_map.move_to(study)
            moved.append(moved_map)
    return moved[:2]





def t1(database):

    desc = [ 
        'Dixon_water', # 'Dixon_in_phase', 'Dixon_out_phase', 'Dixon_fat', 'Dixon_water'
        'LK',
        'RK',  
        'T1m_magnitude_mdr_moco_T1_map',
        # 'T1m_magnitude_mdr_moco_S0_map',
        # 'T1m_magnitude_mdr_moco_T1FAcorr_map',
        # 'T1m_magnitude_mdr_moco',
   ]
    return _align(database, desc)
    
    
def t2(database):

    desc = [ 
        'LK',
        'RK',  
        'T2m_magnitude_mdr_moco_T2_map',
        'T2m_magnitude_mdr_moco_S0_map',
        'T2m_magnitude_mdr_moco',
    ]
    return _align(database, desc)


def t2star(database):

    desc = [ 
        'LK',
        'RK',  
        'T2starm_magnitude_mdr_moco_T2star_map',
        'T2starm_magnitude_mdr_moco_S0_map',
        'T2starm_magnitude_mdr_moco_f_fat_map',
        'T2starm_magnitude_mdr_moco'
    ]
    return _align(database, desc)


def ivim(database):

    desc = [ 
        'LK',
        'RK',  
        'IVIM_mdr_moco_S0_map',
        'IVIM_mdr_moco_D_map',
        'IVIM_mdr_moco_D_star_map',
        'IVIM_mdr_moco_Perfusion_fraction_map',
        'IVIM_mdr_moco',
    ]
    return _align(database, desc)


def dti(database):

    desc = [ 
        'LK',
        'RK',  
        'DTI_mdr_moco_MD_map',
        'DTI_mdr_moco_RD_map',
        'DTI_mdr_moco_AD_map',
        'DTI_mdr_moco_Planarity_map',
        'DTI_mdr_moco_Linearity_map',
        'DTI_mdr_moco_Sphericity_map',
        'DTI_mdr_moco_FA_map',
        'DTI_mdr_moco_fit',
        'DTI_mdr_moco',
    ]
    return _align(database, desc)


def dce(database):

    desc = [ 
        'LK',
        'RK',  
        'DCE_mdr_moco_AUC_map',
        'DCE_mdr_moco_AVD_map',
        'DCE_mdr_moco_RPF_map',
        'DCE_mdr_moco_MTT_map',
        'DCE_mdr_moco',
    ]
    return _align(database, desc)


# TODO: MT is 3D! - needs the same approach as ASL
def mt(database):

    desc = [ 
        'LK',
        'RK',  
        'MT_mdr_moco_AVR_map',
        'MT_mdr_moco_MTR_map',
    ]
    return _align(database, desc)


def asl(database):

    # Get input parameters
    desc = [   
        'Dixon_post_contrast_water',
        'ASL_M0_moco',
        'ASL_RBF_moco',
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
        params = vreg_dicom.find_rigid_transformation(asl_ref, dixon, tolerance=0.1, region=kidney[0], margin=0)
        moved = vreg_dicom.apply_rigid_transformation(rbf, params, description='RBF - ' + kidney[1])
        moved.move_to(study)
        rbf_moved.append(moved)

    return rbf_moved


def dixon(database):

    # Get input parameters
    desc = [   
        'Dixon_out_phase',
        'Dixon_in_phase',
        'Dixon_water',
        'Dixon_fat',
        'Dixon_post_contrast_fat',
    ]
    series, study = input_series(database, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform DIXON alignment: not all required data exist.')

    coregistered, followers = elastix.coregister_3d_to_3d(
        series[3], series[4],
        transformation = "AdvancedMeanSquares", # That does not look right - AMS is a metric not a transformation
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
        'DCE_mdr_moco_AVD_map',
        'DCE_mdr_moco_RPF_map',
        'DCE_mdr_moco_MTT_map',        
    ]
    ref = ['Dixon_post_contrast_fat']

    output = []
    for kidney in ['LK','RK']:
        desc = ref + [kidney] + [d + '_' + kidney + '_align' for d in to_fill]
        series, study = input_series(database, desc, export_study)
        if series is None:
            raise RuntimeError('Cannot fill gaps between slices on LK: not all required data exist.')
        for ms_series in series[2:]:
            output_series = vreg_dicom.fill_slice_gaps(ms_series, series[0], slice_thickness=1.0, mask=series[1])
            output_series.move_to(study)
            output.append(output_series)

    return output



# OBSOLETE

# Issue caused by this in LK
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
    _, _, static_pixel_spacing = vreg.affine_components(affine_static)
    rot_gradient_step, translation_gradient_step, _ = vreg.affine_resolution(array_static.shape, static_pixel_spacing)
    gradient_step = np.concatenate((1.0*rot_gradient_step, 0.5*translation_gradient_step))
    optimization = {
        'method': 'GD', 
        'options': {'gradient step': gradient_step, 'tolerance': 0.1}, 
    }
    try:
        parameters = vreg.align_slice_by_slice(
            moving = array_moving, 
            moving_affine = affine_moving, 
            static = array_static, 
            static_affine = affine_static, 
            parameters = np.zeros(6, dtype=np.float32), 
            resolutions = [4,2,1],
            transformation = vreg.rigid,
            metric = vreg.mutual_information,
            optimization = optimization,
            slice_thickness = slice_thickness,
            progress = lambda z, nz: progress(z+1, nz, 'Coregistering slice-by-slice using rigid transformations'),
        )
    except:
        print('Failed to align volumes..')
        parameters = np.zeros(6, dtype=np.float32)

    return parameters

def __align(database, desc):

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
        array_static = vreg_dicom.pixel_values(series[kidney], dims, on=moving)
        array_moving = _inslice_active_translation(array_static, array_moving, moving.progress, applyto)

        # Through-slice rigid - find new affines for each slice
        array_static = series[kidney].pixel_values('SliceLocation')
        affine_static = series[kidney].affine()
        array_static, affine_static = vreg.mask_volume(array_static, affine_static, array_static, affine_static, 0)
        parameters = _sbs_rigid(array_static, affine_static, array_moving[...,0], affine_moving, moving.progress, slice_thickness)
        affines_moving = vreg.passive_rigid_transform_slice_by_slice(affine_moving, parameters)

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
                    affine_zt = vreg.multislice_to_singleslice_affine(affines_moving[z], frames[z,t].SliceThickness)
                    frames_zt = frames[z,t].copy_to(moved_map)
                    frames_zt.set_affine(affine_zt)
                    frames_zt.set_pixel_array(map_arrays[m][:,:,z,t])
            moved_map.move_to(study)
            moved.append(moved_map)
    return moved[:2]











