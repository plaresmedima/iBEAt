"""
Batch script to overlay right- and left kidney masks on an out of phase DIXON image.

The result is saved in a new directory along with a copy of the out of phase images.
"""

import os
import dbdicom as db
from dbdicom.wrappers import scipy

#source_dir = 'C://Users//UOS//Desktop//13_04_2023_masks_original'
#source_dir = 'C://Users//md1jdsp//Desktop//test_source'
#results_dir = 'C://Users//md1jdsp//Desktop//test_results'

# Define import and export directories
data_dir = 'C:\\Users\\steve\\Dropbox\\Data'
source_dir = os.path.join(data_dir, 'test_source')
results_dir = os.path.join(data_dir, 'test_results')

# Series descriptions of all data
out_desc        = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
left_mask_desc  = 'LK'
right_mask_desc = 'RK'


# Loop over all subfolders of the source directory
for subfolder in os.listdir(source_dir):
    
    # Define paths for the subfolders
    import_path = os.path.join(source_dir, subfolder)

    # Open the source folder
    folder = db.database(path=import_path)

    # Find the right series
    out_ph = folder.series(SeriesDescription=out_desc)
    mask_l = folder.series(SeriesDescription=left_mask_desc)
    mask_r = folder.series(SeriesDescription=right_mask_desc)

    # Export the out of phase series
    out_ph[0].export_as_dicom(results_dir)

    # Overlay each mask on out-phase and export the result
    for mask in mask_l+mask_r:
        overlay = scipy.map_to(mask, out_ph[0])
        overlay.export_as_dicom(results_dir)

    # Remove temporary files
    folder.restore()