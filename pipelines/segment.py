import os
import os.path

from dbdicom.extensions import skimage, scipy, dipy, sklearn
from dbdicom.pipelines import input_series
import mapping.UNETR_kidneys_v1 as unetr


export_study = 'Segmentations'


def kidneys(database, weights):

    # Get weights file and check if valid 
    # if not os.path.isfile(weights):
    #     msg = 'The weights file ' + weights + ' has not been found. \n'
    #     msg += 'Please check that the file with model weights is in the folder, and is named ' + unetr.filename
    #     database.dialog.information(msg)
    #     return

    database.message('Segmenting kidneys. This could take a few minutes. Please be patient..')

    # Get appropriate series and check if valid
    #series = database.series(SeriesDescription=unetr.trained_on)
    sery, study = input_series(database, unetr.trained_on, export_study)
    if sery is None:
        return    

    # Loop over the series and create the mask
    #desc = sery.instance().SeriesDescription
    array, header = sery.array(['SliceLocation'], pixels_first=True, first_volume=True)

    # Calculate predictions 
    masks = unetr.apply(array, weights)
    left_kidney, right_kidney = unetr.kidney_masks(masks)

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

