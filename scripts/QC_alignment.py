import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dbdicom.wrappers import scipy
import os

def main(background_series,mask_series,label,path_to_save):

    overlay_mask  = scipy.map_to(mask_series,background_series)

    array_background, _       = background_series.array (['SliceLocation','AcquisitionTime'], pixels_first=True)
    array_overlay_mask, _     = overlay_mask.array      (['SliceLocation','AcquisitionTime'], pixels_first=True)

    array_background  = np.squeeze(array_background)
    array_overlay_mask = np.squeeze(array_overlay_mask)

    print(array_background)
    print(overlay_mask)

    if len(array_background.shape) > 3:
        array_background  = np.squeeze(array_background[...,0])

    # print(array_background.shape)
    # print(array_overlay_mask.shape)

    array_background   = np.transpose(array_background,(1,0,2))
    array_overlay_mask = np.transpose(array_overlay_mask,(1,0,2))

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
            
                col.imshow(array_background[:,:,i], 'gray', interpolation='none',vmin=0,vmax=np.median(array_background)+np.std(array_background))
                col.imshow(array_overlay_mask[:,:,i], 'jet' , interpolation='none', alpha=0.2)
                col.set_xticklabels([])
                col.set_yticklabels([])
                col.set_aspect('equal')
                col.axis("off")
            i = i +1 
    fig.suptitle(label + "alignment", fontsize=14)
    fig.savefig(os.path.join(path_to_save, label + '_alignment_' + mask_series.SeriesDescription + '.png'), dpi=600)
    print('Done')
    
