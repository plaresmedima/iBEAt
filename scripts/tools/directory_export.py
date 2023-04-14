import os
import dbdicom as db
from dbdicom.wrappers import scipy
import matplotlib.pyplot as plt
import numpy as np

#source_dir = 'C://Users//UOS//Desktop//13_04_2023_masks_original'
source_dir = 'C://Users//md1jdsp//Desktop//test_source'
resultsPath = 'C://Users//md1jdsp//Desktop//test_results'

file_list = os.listdir(source_dir)

out_desc        = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
left_mask_desc  = 'LK'
right_mask_desc = 'RK'


for file_name in file_list:
    
    studie = source_dir + "//" + file_name

    folder = db.database(path=studie)

    out_ph = folder.series(SeriesDescription=out_desc)
    mask_l = folder.series(SeriesDescription=left_mask_desc)
    mask_r = folder.series(SeriesDescription=right_mask_desc)

    print(out_ph[0].SeriesDescription)
    print(mask_l[0].SeriesDescription)
    print(mask_r[0].SeriesDescription)

    mask_overlaid = []

    for mask in mask_l+mask_r:
        overlaid = scipy.map_to(mask, out_ph[0])
        mask_overlaid.append(overlaid)

    # array_mask_l_overlaid, header_mask_l_overlaid = mask_overlaid[0].array(['SliceLocation'], pixels_first=True)
    # plt.imshow(np.squeeze(array_mask_l_overlaid[:,:,90,0]))
    # plt.show()

    exportToFolder = out_ph + mask_overlaid

    for series in exportToFolder:
        print(series.SeriesDescription)    
        series.export_as_dicom(resultsPath+"//"+ file_name)


