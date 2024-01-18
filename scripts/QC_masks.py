import numpy as np
import matplotlib.pyplot as plt
from dbdicom.wrappers import scipy
import os

def main(database,backgroud_series = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast',RK_mask = 'RK', LK_mask = 'LK' ):

    series_img = database.series(SeriesDescription=backgroud_series)
    series_mask_LK = database.series(SeriesDescription=LK_mask)
    series_mask_RK = database.series(SeriesDescription=RK_mask)

    overlay_mask_LK  = scipy.map_to(series_mask_LK[0],series_img[0])
    overlay_mask_RK  = scipy.map_to(series_mask_RK[0],series_img[0])

    array_series_img, header_series_img       = series_img[0].array       (['SliceLocation','AcquisitionTime'], pixels_first=True)
    array_overlay_mask_LK , header_overlay_mask_LK  = overlay_mask_LK.array(['SliceLocation','AcquisitionTime'], pixels_first=True)
    array_overlay_mask_RK , header_overlay_mask_RK  = overlay_mask_RK.array(['SliceLocation','AcquisitionTime'], pixels_first=True)

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
                col.imshow(array_overlay_mask[:,:,i], 'jet' , interpolation='none', alpha=0.2)
                col.set_xticklabels([])
                col.set_yticklabels([])
                col.set_aspect('equal')
                col.axis("off")
            i = i +1 
    fig.suptitle('Kidney Masks', fontsize=14)
    fig.savefig(os.path.join(database.path(),'QC', 'Masks' + '.png'), dpi=600)

