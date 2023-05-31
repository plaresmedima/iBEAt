"""
Batch script to overlay right- and left kidney masks on an out of phase DIXON image.

The result is saved in a new directory along with a copy of the out of phase images.
"""

import os
import dbdicom as db
from dbdicom.wrappers import skimage, scipy, dipy

# Define import and export directories
data_dir = 'C:\\Users\\md1jdsp\\Desktop\\'
source_dir = os.path.join(data_dir, 'sinus_sample')
results_dir = os.path.join(data_dir, 'output_sinus_fat')

# Series descriptions of all data
out_desc        = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'

left_mask_desc  = 'LK'
right_mask_desc = 'RK'

folder = db.database(path=source_dir)
df = []

for patient in folder.patients():

    out = patient.series(SeriesDescription=out_desc)
    lk = patient.series(SeriesDescription='LK')
    rk = patient.series(SeriesDescription='RK')
    kidneys = lk+rk

    cleanup = True

    sf_series = []
    fat_image_masked, fat_mask = dipy.median_otsu(out[0], median_radius=1, numpass=1)
    for kidney in kidneys:
        # Pipeline calculation
        kidney_hull = skimage.convex_hull_image_3d(kidney)
        # sinus = scipy.image_calculator(kidney_hull, kidney, 'series 1 - series 2', integer=True)
        # sinus_fat = scipy.image_calculator(fat_mask, sinus, 'series 1 * series 2', integer=True)
        sinus_fat = scipy.image_calculator(fat_mask, kidney_hull, 'series 1 * series 2', integer=True)
        sinus_fat_largest = scipy.extract_largest_cluster_3d(sinus_fat)
        sinus_fat_largest.SeriesDescription = kidney.instance().SeriesDescription + 'SF'
        # Append and display
        sf_series.append(sinus_fat_largest)

        # Remove intermediate results
        if cleanup:
            kidney_hull.remove()
            #sinus.remove()
            sinus_fat.remove()
    fat_image_masked.remove()
    if cleanup:
        fat_mask.remove()
    #   Collect features & display
    df.append(skimage.volume_features(sf_series))

    for series in sf_series:
        series.export_as_dicom(results_dir)
