import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import imageio

from dbdicom.extensions import vreg
from pipelines.roi_fit import load_aif


def kidney_masks_as_dicom(folder):

    folder.message('Exporting whole kidney masks as dicom..')
    results_path = os.path.join(folder.path() + '_output', 'masks')
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    fat_desc = 'T1w_abdomen_dixon_cor_bh_fat_post_contrast' 
    out_desc = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
    in_desc = 'T1w_abdomen_dixon_cor_bh_in_phase_post_contrast'
    water_desc = 'T1w_abdomen_dixon_cor_bh_water_post_contrast'
    k_means1_desc = 'KMeans cluster 1'
    k_means2_desc = 'KMeans cluster 2'
    lk_mask = 'LK' 
    rk_mask = 'RK'

    fat = folder.series(SeriesDescription=fat_desc)
    out_ph = folder.series(SeriesDescription=out_desc)
    in_ph = folder.series(SeriesDescription=in_desc)
    water = folder.series(SeriesDescription=water_desc)
    k_means1 = folder.series(SeriesDescription=k_means1_desc)
    k_means2 = folder.series(SeriesDescription=k_means2_desc)
    LK = folder.series(SeriesDescription=lk_mask)
    RK = folder.series(SeriesDescription=rk_mask)

    export_to_folder = LK + RK + fat + out_ph + in_ph + water + k_means1 + k_means2
    
    for series in export_to_folder:
        series.export_as_dicom(results_path)



def kidney_masks_as_png(database,backgroud_series = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast',RK_mask = 'RK', LK_mask = 'LK' ):

    database.message('Exporting masks as png..')
    results_path = database.path() + '_output'
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



def aif_as_png(folder):

    folder.message('Exporting AIF as png..')
    results_path = folder.path() + '_output'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    time, aif = load_aif(folder)
    time = time/60

    fig, ax = plt.subplots(1,1,figsize=(5,5))
    ax.plot(time, aif, 'b-', label='Aorta', linewidth=3.0)
    ax.plot(time, 0*time, color='gray')
    ax.set(xlabel='Time (min)', ylabel='DCE signal (a.u.)')
    ax.legend()

    plt.savefig(os.path.join(results_path, 'aif.png'), dpi=600)
    plt.close()
    

def mapping(database):

    results_path = database.path() + '_output'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    for series in database.series():
        desc = series['SeriesDescription']
        if desc[-4:] == '_map':

            array, _ = series.array(['SliceLocation'], pixels_first=True, first_volume=True)
            
            # Create frames for gif movie
            frames = []
            for z in range(array.shape[2]):
                series.progress(z+1, array.shape[2], 'Exporting '+desc+' as gif..')
                smin = np.amin(array[:,:,z])
                smax = np.amax(array[:,:,z])
                frame = 255*(array[:,:,z]-smin)/(smax-smin)
                frame = frame.astype(np.uint8)
                frames.append(frame.T)

            # Save the frames as a GIF
            imageio.mimsave(os.path.join(results_path, 'map_'+ desc + '.gif'), frames, duration=100)



def mdreg(database):

    results_path = database.path() + '_output'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    for series in database.series():
        desc = series['SeriesDescription']
        if desc[-9:] == '_mdr_moco':

            # Get data
            array, _ = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True, first_volume=True)
            array = np.reshape(array, (array.shape[0], array.shape[1], -1))

            # Create frames for gif movie
            frames = []
            for z in range(array.shape[2]):
                series.progress(z+1, array.shape[2], 'Exporting '+desc+' as gif..')
                smin = np.amin(array[:,:,z])
                smax = np.amax(array[:,:,z])
                frame = 255*(array[:,:,z]-smin)/(smax-smin)
                frame = frame.astype(np.uint8)
                frames.append(frame.T)

            # Save frames as .gif
            # Faster speed for the big series
            if desc == "IVIM_kidneys_cor-oblique_fb_mdr_moco":
                duration=25
            elif desc == "DTI_kidneys_cor-oblique_fb_mdr_moco":
                duration=25
            elif desc == "DCE_kidneys_cor-oblique_fb_mdr_moco":
                duration=20
            else:
                duration=100

            imageio.mimsave(os.path.join(results_path, 'mdr_'+desc + '.gif'), frames, duration=duration)
    

def alignment(database):

    database.message('Exporting alignments as png..')
    results_path = database.path() + '_output'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    lk = database.series(SeriesDescription='LK')
    rk = database.series(SeriesDescription='RK')

    for s, background_series in enumerate(database.series(StudyDescription='Alignment')):

        desc = background_series.instance().SeriesDescription
        if desc[-2:] == 'LK':
            overlay_mask  = vreg.map_to(lk, background_series)
        else:
            overlay_mask  = vreg.map_to(rk, background_series)

        array_background, _ = background_series.array(['SliceLocation'], pixels_first=True, first_volume=True)
        array_overlay_mask, _ = overlay_mask.array(['SliceLocation'], pixels_first=True, first_volume=True)

        array_background  = np.transpose(array_background,(1,0,2))
        array_overlay_mask = np.transpose(array_overlay_mask,(1,0,2))

        num_row_cols = int(np.ceil(np.sqrt(array_overlay_mask.shape[2])))
        fig, ax = plt.subplots(nrows=num_row_cols, ncols=num_row_cols, gridspec_kw = {'wspace':0, 'hspace':0}, figsize=(num_row_cols,num_row_cols))
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
        fig.suptitle(desc, fontsize=14)
        fig.savefig(os.path.join(results_path, desc + '.png'), dpi=600)




