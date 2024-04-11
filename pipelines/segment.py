import numpy as np
from dbdicom.extensions import skimage, scipy, dipy, sklearn
from dbdicom.pipelines import input_series
from models import DCE_aorta, PC, UNETR_kidneys_v1



export_study = 'Segmentations'


def kidneys(database, weights):

    # Get weights file and check if valid 
    # if not os.path.isfile(weights):
    #     msg = 'The weights file ' + weights + ' has not been found. \n'
    #     msg += 'Please check that the file with model weights is in the folder, and is named ' + UNETR_kidneys_v1.filename
    #     database.dialog.information(msg)
    #     return

    database.message('Segmenting kidneys. This could take a few minutes. Please be patient..')

    # Get appropriate series and check if valid
    #series = database.series(SeriesDescription=UNETR_kidneys_v1.trained_on)
    sery, study = input_series(database, UNETR_kidneys_v1.trained_on, export_study)
    if sery is None:
        return    

    # Loop over the series and create the mask
    #desc = sery.instance().SeriesDescription
    array, header = sery.array(['SliceLocation'], pixels_first=True, first_volume=True)

    # Calculate predictions 
    masks = UNETR_kidneys_v1.apply(array, weights)
    left_kidney, right_kidney = UNETR_kidneys_v1.kidney_masks(masks)

    # Save UNETR output
    result = study.new_child(SeriesDescription = 'BK')
    result.set_array(masks, header, pixels_first=True)
    # result[['WindowCenter','WindowWidth']] = [1.0, 2.0]

    # Save and display left kidney data
    left = study.new_child(SeriesDescription = 'LK')
    left.set_array(left_kidney, header, pixels_first=True)
    # left[['WindowCenter','WindowWidth']] = [1.0, 2.0]
    
    # Save and display right kidney data
    right = study.new_child(SeriesDescription = 'RK')
    right.set_array(right_kidney, header, pixels_first=True)
    # right[['WindowCenter','WindowWidth']] = [1.0, 2.0]

    return left, right


def renal_sinus_fat(folder):

    fat = folder.series(SeriesDescription='T1w_abdomen_dixon_cor_bh_fat_post_contrast')
    lk  = folder.series(SeriesDescription='LK')
    rk  = folder.series(SeriesDescription='RK')

    kidneys = lk+rk
    sf_series = []

    if len(kidneys)==[]:
        msg = 'Cannot perform renal sinus fat segmentation: no kidney masks are available.'
        raise RuntimeError(msg)

    fat_image_masked, fat_mask = dipy.median_otsu(fat[0], median_radius=1, numpass=1)

    for kidney in kidneys:
        kidney_hull = skimage.convex_hull_image_3d(kidney)
        sinus_fat = scipy.image_calculator(fat_mask, kidney_hull, 'series 1 * series 2', integer=True)
        #sinus_fat_open = skimage.opening_3d(sinus_fat)
        sinus_fat_largest = scipy.extract_largest_cluster_3d(sinus_fat)
        sinus_fat_largest.SeriesDescription = kidney.instance().SeriesDescription + 'SF'
        sf_series.append(sinus_fat_largest)
        # Cleanup
        kidney_hull.remove()
        #sinus.remove()
        sinus_fat.remove()
    
    fat_image_masked.remove()
    fat_mask.remove()   

    return sf_series


def compute_whole_kidney_canvas(database):
    series_desc = [
        'T1w_abdomen_dixon_cor_bh_fat_post_contrast',
        'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
    ] 
    features, study = input_series(database, series_desc, export_study)
    if features is None:
        return
    clusters = sklearn.sequential_kmeans(features, n_clusters=2, multiple_series=True)
    for c in clusters:
        c.move_to(study)
    return clusters


def aorta_on_dce(folder):

    desc = "DCE_aorta_axial_fb"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot create DCE-AIF mask: series ' + desc + ' does not exist. ')

    axial, header = series.array(['AcquisitionTime'], pixels_first=True, first_volume=True)

    aif_mask = DCE_aorta.segment(axial)

    aif_mask_series = study.new_series(SeriesDescription='DCE-AIF')
    aif_mask_series.set_array(aif_mask, header[0], pixels_first=True)

    return aif_mask_series


def renal_artery(folder):

    desc = [
        'PC_RenalArtery_Right_EcgTrig_fb_120_magnitude',
        'PC_RenalArtery_Right_EcgTrig_fb_120_velocity',
        'PC_RenalArtery_Left_EcgTrig_fb_120_magnitude',
        'PC_RenalArtery_Left_EcgTrig_fb_120_velocity',
    ]
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot create PC-RA mask: some series are missing.')
    
    dims = ['InstanceNumber']
    dx = series[0].values('PixelSpacing')[0][0]

    # Calculate right artery mask:
    mag, mag_hdr_right = series[0].array(dims, pixels_first=True, first_volume=True)
    vel, vel_hdr_right = series[1].array(dims, pixels_first=True, first_volume=True) 
    right_mask = PC.renal_artery_mask(mag, -vel, pixel_spacing=dx)

    # Calculate left artery mask:
    mag, mag_hdr_left = series[2].array(dims, pixels_first=True, first_volume=True)
    vel, vel_hdr_left = series[3].array(dims, pixels_first=True, first_volume=True)
    left_mask = PC.renal_artery_mask(mag, vel, pixel_spacing=dx)

    # Save as DICOM
    right_mask_series = study.new_series(SeriesDescription='PC-RRA')
    right_mask_series.set_array(right_mask, mag_hdr_right[0], pixels_first=True)

    left_mask_series = study.new_series(SeriesDescription='PC-LRA')
    left_mask_series.set_array(left_mask, mag_hdr_left[0], pixels_first=True)

    return right_mask_series, left_mask_series


def cortex_medulla(database):

    features = [
        'DCE_kidneys_cor-oblique_fb_mdr_moco_AVD_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco_RPF_map',
        'DCE_kidneys_cor-oblique_fb_mdr_moco_MTT_map',  
    ]

    output = []
    for kidney in ['LK','RK']:
        desc = [kidney] + [f + '_' + kidney + '_align_fill' for f in features]
        series, study = input_series(database, desc, export_study)
        if series is None:
            raise ValueError('Cannot separate cortex and medulla for kidney '+kidney+': required sequences are missing.')
        clusters, cluster_features = sklearn.kmeans(series[1:], series[0], n_clusters=3, multiple_series=True, return_features=True)
        # Background = cluster with smallest AVD
        background = np.argmin([c[0] for c in cluster_features])
        # Cortex = cluster with largest RPF
        cortex = np.argmax([c[1] for c in cluster_features]) 
        # Medulla = cluster with largest MTT
        medulla = np.argmax([c[2] for c in cluster_features])
        # Check
        remainder = {0,1,2} - {background, cortex, medulla}
        if len(remainder) > 0:
            raise ValueError('Problem separating cortex and medulla: identified clusters do not have the expected values.')
        clusters[cortex].SeriesDescription = kidney + 'C'
        clusters[cortex].move_to(study)
        clusters[medulla].SeriesDescription = kidney + 'M'
        clusters[medulla].move_to(study)
        clusters[background].SeriesDescription = kidney + 'B'
        clusters[background].move_to(study)
        cm = scipy.image_calculator(clusters[cortex], clusters[medulla], '+')
        cm.SeriesDescription = kidney + 'CM'
        output += clusters + [cm]
    return output