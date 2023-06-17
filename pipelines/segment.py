import os

from dbdicom.wrappers import skimage, scipy, dipy, sklearn
from dbdicom.pipelines import input_series

export_study = 'SegmentationResults'


def compute_whole_kidney_canvas(database):
    series_desc = [
        'T1w_abdomen_dixon_cor_bh_fat_post_contrast',
        'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
    ] 
    features, study = input_series(database, series_desc, export_study)
    if features is None:
        return None
    clusters = sklearn.sequential_kmeans(features, n_clusters=2, multiple_series=True)
    for c in clusters:
        c.move_to(study)
    return clusters


def export_whole_kidney_canvas(folder, path):

    resultsFolder = 'segmentation_results'
    resultsPath = path + '_' + resultsFolder

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


def compute_renal_sinus_fat(database):
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

    return sf_series, df