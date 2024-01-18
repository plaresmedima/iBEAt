import os
import time

from dbdicom.wrappers import skimage, scipy, dipy, sklearn
from dbdicom.pipelines import input_series
import mapping.UNETR_kidneys_v1 as unetr
import matplotlib.pyplot as plt


export_study = 'Segmentations'

def segment_kidneys(database, weights):

    start_time = time.time()
    database.log('Kidney segmentation has started.')

    # Get weights file and check if valid 
    # if not os.path.isfile(weights):
    #     msg = 'The weights file ' + weights + ' has not been found. \n'
    #     msg += 'Please check that the file with model weights is in the folder, and is named ' + unetr.filename
    #     database.dialog.information(msg)
    #     return

    # Get appropriate series and check if valid
    series = database.series(SeriesDescription=unetr.trained_on)
    sery, study = input_series(database, unetr.trained_on, export_study)
    if sery is None:
        return    

    # Loop over the series and create the mask
    desc = sery.instance().SeriesDescription
    array, header = sery.array(['SliceLocation'], pixels_first=True, first_volume=True)

    # Calculate predictions 
    masks = unetr.apply(array, weights)
    left_kidney, right_kidney = unetr.kidney_masks(masks)

    # Save UNETR output
    result = study.new_child(SeriesDescription = 'UNETR kidneys v1')
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

    # kidneys = [left, right]
    # features = skimage.volume_features(kidneys)

    database.save()

    database.log("Kidney segmentation was completed --- %s seconds ---" % (int(time.time() - start_time)))

def compute_whole_kidney_canvas(database):
    start_time = time.time()
    database.log('Sequential K-means has started.')
    series_desc = [
        'T1w_abdomen_dixon_cor_bh_fat_post_contrast',
        'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
    ] 
    # series_desc = [
    #     'T1w_abdomen_dixon_cor_bh_fat',
    #     'T1w_abdomen_dixon_cor_bh_out_phase'
    # ] 

    features, study = input_series(database, series_desc, export_study)
    if features is None:
        return None
    clusters = sklearn.sequential_kmeans(features, n_clusters=2, multiple_series=True)
    for c in clusters:
        c.move_to(study)
    database.log("Sequential kmeans was completed --- %s seconds ---" % (int(time.time() - start_time)))
    return clusters


def export_whole_kidney_canvas(folder):
    start_time = time.time()
    folder.log("Export segmentation results has started")

    resultsFolder = 'segmentation_canvas'

    path = folder.path()
    resultsPath = os.path.join(path,resultsFolder)

    os.mkdir(resultsPath)

    fat_desc        = 'T1w_abdomen_dixon_cor_bh_fat_post_contrast' 
    out_desc        = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
    in_desc         = 'T1w_abdomen_dixon_cor_bh_in_phase_post_contrast'
    water_desc      = 'T1w_abdomen_dixon_cor_bh_water_post_contrast'
    k_means1_desc   = 'KMeans cluster 1'
    k_means2_desc   = 'KMeans cluster 2'

    fat             = folder.series(SeriesDescription=fat_desc)
    out_ph          = folder.series(SeriesDescription=out_desc)
    in_ph           = folder.series(SeriesDescription=in_desc)
    water           = folder.series(SeriesDescription=water_desc)
    k_means1        = folder.series(SeriesDescription=k_means1_desc)
    k_means2        = folder.series(SeriesDescription=k_means2_desc)

    exportToFolder = fat + out_ph + in_ph + water + k_means1 + k_means2
    
    for series in exportToFolder:
        print(series.SeriesDescription)    
        series.export_as_dicom(resultsPath)

    folder.log("Export segmentation results was completed --- %s seconds ---" % (int(time.time() - start_time)))

def export_masks(folder):
    start_time = time.time()
    folder.log("Export segmentation results has started")

    path = folder.path()

    resultsFolder = 'masks'
    resultsPath = os.path.join(path,resultsFolder)

    os.mkdir(resultsPath)

    lk_mask        = 'LK' 
    rk_mask        = 'RK'
    fat_desc        = 'T1w_abdomen_dixon_cor_bh_fat_post_contrast' 
    out_desc        = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
    in_desc         = 'T1w_abdomen_dixon_cor_bh_in_phase_post_contrast'
    water_desc      = 'T1w_abdomen_dixon_cor_bh_water_post_contrast'


    LK             = folder.series(SeriesDescription=lk_mask)
    RK             = folder.series(SeriesDescription=rk_mask)
    fat            = folder.series(SeriesDescription=fat_desc)
    out_ph         = folder.series(SeriesDescription=out_desc)
    in_ph          = folder.series(SeriesDescription=in_desc)
    water          = folder.series(SeriesDescription=water_desc)

    exportToFolder = LK + RK + fat + out_ph + in_ph + water
    
    for series in exportToFolder:
        print(series.SeriesDescription)    
        series.export_as_dicom(resultsPath)

    folder.log("Export maks results was completed --- %s seconds ---" % (int(time.time() - start_time)))


def compute_renal_sinus_fat(database):
    start_time = time.time()
    database.log("Compute renal sinus fat has started")
    series_desc = [
        'T1w_abdomen_dixon_cor_bh_fat_post_contrast',
        'LK', 'RK',
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        return None, None
    fat_image = series[0]
    kidneys = series[1:]

    sf_series = []
    fat_image_masked, fat_mask = dipy.median_otsu(fat_image, median_radius=1, numpass=1)

    for kidney in kidneys:

        # Pipeline calculation
        kidney_hull = skimage.convex_hull_image_3d(kidney)
        sinus_fat = scipy.image_calculator(fat_mask, kidney_hull, 'series 1 * series 2', integer=True)
        sinus_fat_open = skimage.opening_3d(sinus_fat)
        sinus_fat_largest = scipy.extract_largest_cluster_3d(sinus_fat_open)
        sinus_fat_largest.SeriesDescription = kidney.instance().SeriesDescription + 'SF'

        # Remove intermediate results and append
        kidney_hull.remove()
        sinus_fat.remove()
        sinus_fat_open.remove()
        sinus_fat_largest.move_to(study)
        sf_series.append(sinus_fat_largest)

    fat_image_masked.remove()
    fat_mask.remove()

    # Collect features & display
    df = skimage.volume_features(sf_series)

    database.log("Renal sinus fat computation was completed --- %s seconds ---" % (int(time.time() - start_time)))

    return sf_series, df