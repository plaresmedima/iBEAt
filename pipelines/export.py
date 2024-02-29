import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from dbdicom.extensions import vreg


def kidney_masks_as_dicom(folder):

    folder.message('Exporting whole kidney masks as dicom..')
    results_path = folder.path() + '_masks'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

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
        series.export_as_dicom(results_path)


def kidney_masks_as_png(database,backgroud_series = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast',RK_mask = 'RK', LK_mask = 'LK' ):

    database.message('Exporting masks as png..')
    results_path = database.path() + '_QC'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    series_img = database.series(SeriesDescription=backgroud_series)
    series_mask_LK = database.series(SeriesDescription=LK_mask)
    series_mask_RK = database.series(SeriesDescription=RK_mask)

    overlay_mask_LK  = vreg.map_to(series_mask_LK[0],series_img[0])
    overlay_mask_RK  = vreg.map_to(series_mask_RK[0],series_img[0])

    array_series_img, _ = series_img[0].array(['SliceLocation','AcquisitionTime'], pixels_first=True)
    array_overlay_mask_LK , _  = overlay_mask_LK.array(['SliceLocation','AcquisitionTime'], pixels_first=True)
    array_overlay_mask_RK , _  = overlay_mask_RK.array(['SliceLocation','AcquisitionTime'], pixels_first=True)

    overlay_mask_LK.remove()
    overlay_mask_RK.remove()

    array_series_img = np.squeeze(array_series_img)
    array_overlay_mask_LK = np.squeeze(array_overlay_mask_LK)
    array_overlay_mask_RK = np.squeeze(array_overlay_mask_RK)

    array_series_img = array_series_img.transpose((1,0,2))
    array_overlay_mask_LK = array_overlay_mask_LK.transpose((1,0,2))
    array_overlay_mask_RK = array_overlay_mask_RK.transpose((1,0,2))

    array_overlay_mask = array_overlay_mask_LK + array_overlay_mask_RK

    num_row_cols = int(np.ceil (np.sqrt(array_overlay_mask.shape[2])))

    fig, ax = plt.subplots(nrows=num_row_cols, ncols=num_row_cols,gridspec_kw = {'wspace':0, 'hspace':0},figsize=(num_row_cols,num_row_cols))
    i=0
    for row in ax:
        for col in row:
            if i>=array_overlay_mask.shape[2]:
                col.set_xticklabels([])
                col.set_yticklabels([])
                col.set_aspect('equal')
                col.axis("off")
            else:  
            
                col.imshow(array_series_img[:,:,i], 'gray', interpolation='none',vmin=0,vmax=np.mean(array_series_img)+np.std(array_series_img))
                col.imshow(array_overlay_mask[:,:,i], 'jet' , interpolation='none', alpha=0.5)
                col.set_xticklabels([])
                col.set_yticklabels([])
                col.set_aspect('equal')
                col.axis("off")
            i = i +1 
    fig.suptitle('Kidney Masks', fontsize=14)
    fig.savefig(os.path.join(results_path, 'Masks.png'), dpi=600)


def whole_kidney_canvas(folder):

    folder.message('Exporting whole kidney canvas as dicom..')
    results_path = folder.path() + '_segmentation_canvas'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

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
        series.export_as_dicom(results_path)

    




