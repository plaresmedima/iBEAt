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

    fat_desc = 'Dixon_post_contrast_fat_' 
    out_desc = 'Dixon_post_contrast_out_phase'
    in_desc = 'Dixon_post_contrast_in_phase'
    water_desc = 'Dixon_post_contrast_water'
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



def kidney_masks_as_png(database,backgroud_series = 'Dixon_post_contrast_out_phase',RK_mask = 'RK', LK_mask = 'LK' ):

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

    #overlay_mask_LK.remove()
    #overlay_mask_RK.remove()

    array_series_img = np.squeeze(array_series_img)
    array_overlay_mask_LK = np.squeeze(array_overlay_mask_LK)
    array_overlay_mask_RK = np.squeeze(array_overlay_mask_RK)

    array_series_img = array_series_img.transpose((1,0,2))
    array_overlay_mask_LK = array_overlay_mask_LK.transpose((1,0,2))
    array_overlay_mask_RK = array_overlay_mask_RK.transpose((1,0,2))

    array_overlay_mask = array_overlay_mask_LK + array_overlay_mask_RK
    array_overlay_mask[array_overlay_mask >0.5] = 1
    array_overlay_mask[array_overlay_mask <0.5] = np.nan
    


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
            
                # Display the background image
                col.imshow(array_series_img[:,:,i], cmap='gray', interpolation='none', vmin=0, vmax=np.mean(array_series_img) + 2 * np.std(array_series_img))
            
                # Overlay the mask with transparency
                col.imshow(array_overlay_mask[:,:,i], cmap='autumn', interpolation='none', alpha=0.5)

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


    scale_dict = {

    'T1m_magnitude_mdr_moco_err_map'       : (0,50), 
    'T1m_magnitude_mdr_moco_T1_map'        : (1000,2000), 
    'T1m_magnitude_mdr_moco_T1FAcorr_map'  : (0,20),
    'T2m_magnitude_mdr_moco_err_map'       : (0,100), 
    'T2m_magnitude_mdr_moco_T2_map'        : (40,100), 
    'T2starm_magnitude_mdr_moco_err_map'   : (0,50),
    'T2starm_magnitude_mdr_moco_T2star_map': (20,80),
    'T2starm_magnitude_mdr_moco_f_fat_map' : (0,0.4),
    'MT_mdr_moco_MTR_map'                  : (0,50), 
    'MT_mdr_moco_AVR_map'                  : (0,120), 
    'DTI_mdr_moco_FA_map'                  : (0,0.6), 
    'DTI_mdr_moco_MD_map'                  : (0.001,0.003),
    'DTI_mdr_moco_Sphericity_map'          : (0.5,1),
    'DTI_mdr_moco_Linearity_map'           : (0,0.5),
    'DTI_mdr_moco_Planarity_map'           : (0,0.5),
    'DTI_mdr_moco_AD_map'                  : (0.001,0.003),
    'DTI_mdr_moco_RD_map'                  : (0.001,0.003), 
    'DCE_mdr_moco_ATT_map'                 : (0,15),
    'DCE_mdr_moco_RPF_map'                 : (100,400),
    'DCE_mdr_moco_AVD_map'                 : (0,300),
    'DCE_mdr_moco_MTT_map'                 : (0,200), 
    'DCE_mdr_moco_FP_map'                  : (50,300),
    'DCE_mdr_moco_TP_map'                  : (0,100), 
    'DCE_mdr_moco_VP_map'                  : (50,150),
    'DCE_mdr_moco_FT_map'                  : (10,150),
    'DCE_mdr_moco_TT_map'                  : (0,150)
    }

    results_path = database.path() + '_output'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    for series in database.series():
        desc = series['SeriesDescription']

        if desc[-4:] == '_map':
            if '_AD_' in desc or '_RD_' in desc or '_Sphericity_' in desc or '_Linearity_' in desc or '_Planarity_' in desc:
                continue
            else:    
                print(desc)
                array, _ = series.array(['SliceLocation'], pixels_first=True, first_volume=True)
                array  = np.transpose(array,(1,0,2))
                
                vmin, vmax = scale_dict.get(desc, (np.median(array) - np.std(array), np.median(array) + np.std(array)))

                num_row_cols = int(np.ceil (np.sqrt(array.shape[2])))

                fig, ax = plt.subplots(nrows=num_row_cols, ncols=num_row_cols,gridspec_kw = {'wspace':0, 'hspace':0},figsize=(num_row_cols,num_row_cols))
                i=0
                for row in ax:
                    for col in row:
                        if i>=array.shape[2]:
                            col.set_xticklabels([])
                            col.set_yticklabels([])
                            col.set_aspect('equal')
                            col.axis("off")
                        else:  
                            
                            if "_S0_" in desc:
                                im = col.imshow(array[:,:,i], cmap='gray', interpolation='none', vmin=vmin, vmax=vmax)
                            else:
                                im = col.imshow(array[:,:,i], cmap='jet', interpolation='none', vmin=vmin, vmax=vmax)

                            col.set_xticklabels([])
                            col.set_yticklabels([])
                            col.set_aspect('equal')
                            col.axis("off")
                        i = i +1 
                
                #cbar = fig.colorbar(im, ax=ax[-1, -1], orientation='vertical', fraction=0.046, pad=0.04)
                fig.suptitle(desc, fontsize=10)
                fig.savefig(os.path.join(results_path, 'map_'+ desc +'.png'), dpi=600)




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
            if desc == "IVIM_mdr_moco":
                duration=25
            elif desc == "DTI_mdr_moco":
                duration=25
            elif desc == "DCE_mdr_moco":
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

    desc_list = [
    
    'MT_mdr_moco_MTR_map_LK_align',
    'MT_mdr_moco_MTR_map_RK_align',
    'T1m_magnitude_mdr_moco_T1_map_LK_align',
    'T1m_magnitude_mdr_moco_T1_map_RK_align',
    'T2m_magnitude_mdr_moco_T2_map_LK_align',
    'T2m_magnitude_mdr_moco_T2_map_RK_align',
    'T2starm_magnitude_mdr_moco_T2star_map_LK_align',
    'T2starm_magnitude_mdr_moco_T2star_map_RK_align',
    'DTI_mdr_moco_MD_map_LK_align',
    'DTI_mdr_moco_MD_map_RK_align',
    'DCE_mdr_moco_AUC_map_LK_align',
    'DCE_mdr_moco_AUC_map_RK_align'
    ]

    scale_dict = {

    'MT_mdr_moco_MTR_map_LK_align'                  : (0,80), 
    'MT_mdr_moco_MTR_map_RK_align'                  : (0,80),
    'T1m_magnitude_mdr_moco_T1_map_LK_align'        : (1000,2000),
    'T1m_magnitude_mdr_moco_T1_map_RK_align'        : (1000,2000),
    'T2m_magnitude_mdr_moco_T2_map_LK_align'        : (10,100),
    'T2m_magnitude_mdr_moco_T2_map_RK_align'        : (10,100),
    'T2starm_magnitude_mdr_moco_T2star_map_LK_align': (0,80),
    'T2starm_magnitude_mdr_moco_T2star_map_RK_align': (0,80),
    'DTI_mdr_moco_MD_map_LK_align'                  : (0.001,0.003),
    'DTI_mdr_moco_MD_map_RK_align'                  : (0.001,0.003),
    }



    for s, background_series in enumerate(database.series(StudyDescription='3: Alignment')):
        
        desc = background_series.instance().SeriesDescription

        if desc in desc_list:

            if 'LK' in desc:
                overlay_mask  = vreg.map_to(lk[0], background_series)
            else:
                overlay_mask  = vreg.map_to(rk[0], background_series)

            array_background, _ = background_series.array(['SliceLocation'], pixels_first=True, first_volume=True)
            array_overlay_mask, _ = overlay_mask.array(['SliceLocation'], pixels_first=True, first_volume=True)

            array_background  = np.transpose(array_background,(1,0,2))
            array_overlay_mask = np.transpose(array_overlay_mask,(1,0,2))

            array_overlay_mask[array_overlay_mask >0.5] = 1
            array_overlay_mask[array_overlay_mask <0.5] = np.nan

            vmin, vmax = scale_dict.get(desc, (0, np.median(array_background) + 2 * np.std(array_background)))

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
                    
                        col.imshow(array_background[:,:,i], 'gray', interpolation='none',vmin=vmin,vmax=vmax)
                        col.imshow(array_overlay_mask[:,:,i], 'autumn' , interpolation='none', alpha=0.5)
                        col.set_xticklabels([])
                        col.set_yticklabels([])
                        col.set_aspect('equal')
                        col.axis("off")
                    i = i +1 
            fig.suptitle(desc, fontsize=10)
            fig.savefig(os.path.join(results_path, desc + '.png'), dpi=600)




