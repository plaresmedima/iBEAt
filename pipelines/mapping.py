import os
import time
import numpy as np
from tqdm import tqdm
import multiprocessing

from dipy.core.gradients import gradient_table
import dipy.reconst.dti as dti
from dipy.reconst.dti import fractional_anisotropy
from scipy.integrate import trapz
from mapping import t1_t2_alone
from mapping import t1_alone

import mapping.t1_exp
import mapping.t1_philips
import mapping.t1
import mapping.t2
import mapping.T2s

import pipelines.segment as seg

from utilities import map_to_dixon_mask
from utilities import fill_kidney_holes_interp_v2


def T1_Philips(series, study, mask=None):

    start_time = time.time()
    series.log("T1 mapping (Philips) has started")

    array, header = series.array(['SliceLocation',(0x2005, 0x1572)], pixels_first=True)

    if slice is not None:
        magnitude_array_T2s_slice = magnitude_array_T2s[:,:,int(slice-1),:]
        magnitude_array_T2s = magnitude_array_T2s_slice
    
    if mask is not None:
        mask = np.transpose(mask)
        for i_slice in range (np.shape(magnitude_array_T2s)[2]):
            for i_w in range (np.shape(magnitude_array_T2s)[3]):
                magnitude_array_T2s[:,:,i_slice,i_w]=magnitude_array_T2s[:,:,i_slice,i_w]*mask

    fittedMaps = mapping.t1_philips.main(array, header)

    M0map, T1_app_map, T1_map, rsquaremap = fittedMaps

    M0_map_series = series.SeriesDescription + "_T1_" + "M0_Map"
    M0_map_series = study.new_series(SeriesDescription=M0_map_series)
    M0_map_series.set_array(M0map,np.squeeze(header[:,0]),pixels_first=True)

    T1_app_map_series = series.SeriesDescription + "_T1_" + "T1_app_Map"
    T1_app_map_series = study.new_series(SeriesDescription=T1_app_map_series)
    T1_app_map_series.set_array(T1_app_map,np.squeeze(header[:,0]),pixels_first=True)

    T1_map_series = series.SeriesDescription + "_T1_" + "T1_Map"
    T1_map_series = study.new_series(SeriesDescription=T1_map_series)
    T1_map_series.set_array(T1_map,np.squeeze(header[:,0]),pixels_first=True)

    rsquare_map_series = series.SeriesDescription + "_T1_" + "rsquare_Map"
    rsquare_map_series = study.new_series(SeriesDescription=rsquare_map_series)
    rsquare_map_series.set_array(rsquaremap,np.squeeze(header[:,0]),pixels_first=True)

    series.log("T1 mapping (Philips) was completed --- %s seconds ---" % (int(time.time() - start_time)))

    return M0_map_series, T1_app_map_series, T1_map_series


def T1_MOLLI(series, study):

    start_time = time.time()
    series.log("T1 mapping (Siemens) has started")
    array, header = series.array(['SliceLocation','InversionTime'], pixels_first=True, first_volume=True)

    # PARAMETER VARIABLES INITIALIZATION
    model_fit = np.empty(array.shape)
    pars = np.empty(array.shape[:3] + (5,) )

    # LOOP THROUGH SLICE
    for z in range(array.shape[2]):
        msg = 'Fitting T1 model: slice ' + str(1+z) + ' out of ' + str(array.shape[2])
        TI = np.array([hdr['InversionTime'] for hdr in header[z,:]])
        for x in range(array.shape[0]):
            series.progress(1+x, array.shape[0], msg)
            for y in range(array.shape[1]):
                fit, S0, S0_App, T1_app, T1, r_square = mapping.t1_exp.main(array[x,y,z,:], TI)
                model_fit[x,y,z,:] = fit
                pars[x,y,z,0] = S0
                pars[x,y,z,1] = S0_App
                pars[x,y,z,2] = T1_app
                pars[x,y,z,3] = T1
                pars[x,y,z,4] = r_square

    M0_map_series = study.new_series(SeriesDescription="M0")
    M0_map_series.set_array(pars[...,0], header[:,0], pixels_first=True)

    T1_app_map_series = study.new_series(SeriesDescription="T1app")
    T1_app_map_series.set_array(pars[...,2], header[:,0], pixels_first=True)

    T1_map_series = study.new_series(SeriesDescription="T1")
    T1_map_series.set_array(pars[...,3], header[:,0], pixels_first=True)

    rsquare_map_series = study.new_series(SeriesDescription="Rsquare")
    rsquare_map_series.set_array(pars[...,4], header[:,0], pixels_first=True)

    series.log("T1 mapping (Siemens) was completed --- %s seconds ---" % (int(time.time() - start_time)))

    return M0_map_series, T1_app_map_series, T1_map_series

def T2s(series=None, mask=None,export_ROI=False,slice=None,Fat_export=False,study = None):

    start_time = time.time()
    series.log("T2* mapping has started")

    series_T2s = series

    array, header = series.array(['SliceLocation', 'EchoTime'], pixels_first=True)

    TE_list = [hdr["EchoTime"] for hdr in (header[0,:,0])] 

    #Check if the data corresponds to the Siemens protocol (12 TEs)        
    if len(TE_list) == 12 and np.max(TE_list) < 50:

        #app.dialog.information("T2* Mapping has started")

        magnitude_array_T2s = array

        if slice is not None:
            magnitude_array_T2s_slice = magnitude_array_T2s[:,:,int(slice-1),:]
            magnitude_array_T2s = magnitude_array_T2s_slice

        if mask is not None:
            mask = np.transpose(mask)
            for i_slice in range (np.shape(magnitude_array_T2s)[2]):
                for i_w in range (np.shape(magnitude_array_T2s)[3]):
                    magnitude_array_T2s[:,:,i_slice,i_w]=magnitude_array_T2s[:,:,i_slice,i_w]*mask

        #T2* mapping input: T2*-weighted images (x,y,z,TE), echo times, wezel as optional argument to create progress bars in to wezel interface
        M0map, fwmap, T2smap, rsquaremap = mapping.T2s.main(magnitude_array_T2s, TE_list)

        #wezel vizualitation of T2* mapping parameters: M0 map, Water Fraction map, T2* map,T2* r square (goodness of fit)
        M0_map_series = series_T2s.SeriesDescription + "_T2s_" + "M0_Map"
        M0_map_series = study.new_series(SeriesDescription=M0_map_series)
        M0_map_series.set_array(M0map,np.squeeze(header[:,0]),pixels_first=True)

        fw_map_series = series_T2s.SeriesDescription + "_T2s_" + "fw_Map"
        fw_map_series = study.new_series(SeriesDescription=fw_map_series)
        fw_map_series.set_array(fwmap,np.squeeze(header[:,0]),pixels_first=True)

        T2s_map_series = series_T2s.SeriesDescription + "_T2s_" + "T2s_Map"
        T2s_map_series = study.new_series(SeriesDescription=T2s_map_series)
        T2s_map_series.set_array(T2smap,np.squeeze(header[:,0]),pixels_first=True)

        rsquare_map_series = series_T2s.SeriesDescription + "_T2s_" + "rsquare_Map"
        rsquare_map_series = study.new_series(SeriesDescription=rsquare_map_series)
        rsquare_map_series.set_array(rsquaremap,np.squeeze(header[:,0]),pixels_first=True)

        series.log("T2* mapping was completed --- %s seconds ---" % (int(time.time() - start_time))) 

        return M0_map_series, fw_map_series, T2s_map_series


def T1T2_Modelling(series_T1_T2, study=None):

    series_T1 = series_T1_T2[0]
    series_T2 = series_T1_T2[1]

    array_T1, header_T1 = series_T1.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    array_T2, header_T2 = series_T2.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)

    array_T1 = np.squeeze(array_T1[:,:,:,:,0])
    array_T2 = np.squeeze(array_T2[:,:,:,:,0])

    header_T1 = np.squeeze(header_T1[:,...])
    header_T2 = np.squeeze(header_T2[:,...])

    TR = 4.6                            #in ms
    FA = 12    #in degrees
    FA_rad = FA/360*(2*np.pi)           #convert to rads
    N_T1 = 66                           #number of k-space lines
    FA_Cat  = [(-FA/5)/360*(2*np.pi), (2*FA/5)/360*(2*np.pi), (-3*FA/5)/360*(2*np.pi), (4*FA/5)/360*(2*np.pi), (-5*FA/5)/360*(2*np.pi)] #cat module

    FA_Cat  = [(-FA/5)/360*(2*np.pi), (2*FA/5)/360*(2*np.pi), (-3*FA/5)/360*(2*np.pi), (4*FA/5)/360*(2*np.pi), (-5*FA/5)/360*(2*np.pi)] #cat module

    TE = [0,30,40,50,60,70,80,90,100,110,120]
    Tspoil = 1
    N_T2 = 72
    Trec = 463*2
    FA_eff = 0.6

    T1_S0_map = np.zeros(np.shape(array_T1)[0:3])
    T1_map = np.zeros(np.shape(array_T1)[0:3])
    FA_Eff_map = np.zeros(np.shape(array_T1)[0:3])
    Ref_Eff_map = np.zeros(np.shape(array_T1)[0:3])
    T2_S0_map = np.zeros(np.shape(array_T1)[0:3])
    T2_map = np.zeros(np.shape(array_T1)[0:3])
    T1_rsquare_map = np.zeros(np.shape(array_T1)[0:3])
    T2_rsquare_map = np.zeros(np.shape(array_T1)[0:3])

    for i in range(np.shape(array_T1)[2]):
        Kidney_pixel_T1 = np.squeeze(array_T1[...,i,:])
        Kidney_pixel_T2 = np.squeeze(array_T2[...,i,:])

        if np.size(np.shape(np.squeeze(header_T1)))==2:
            TI_temp =  [float(hdr['InversionTime']) for hdr in header_T1[i,:]]
        elif np.size(np.shape(np.squeeze(header_T1)))==3:
            TI_temp =  [float(hdr['InversionTime']) for hdr in header_T1[i,:,0]]

        pool = multiprocessing.Pool(processes=os.cpu_count())

        arguments =[]
        pool = multiprocessing.Pool(initializer=multiprocessing.freeze_support,processes=os.cpu_count())
        for (x, y), _ in np.ndenumerate(Kidney_pixel_T1[..., 0]):
            t1_value = Kidney_pixel_T1[x, y, :]
            t2_value = Kidney_pixel_T2[x, y, :]

            arguments.append((x,y,t1_value,t2_value,TI_temp,TE,FA_rad,TR,N_T1,N_T2,FA_Cat,Trec,FA_eff,Tspoil))

        results = list(tqdm(pool.imap(t1_t2_alone.main, arguments), total=len(arguments), desc='Processing pixels of slice ' + str(i)))

        for result in results:
            xi = result[0]
            yi = result[1]
            T1 = result[2]
            T2 = result[3]
            S0_T1 = result[4]
            S0_T2 = result[5]
            FA_eff = result[6]
            r_squared_T1 = result[7]
            r_squared_T2 = result[8]

            r_squared_T1 = result[7]
            r_squared_T2 = result[8]
            T1_map[xi,yi,i] = T1
            T2_map[xi,yi,i] = T2
            T1_S0_map[xi,yi,i] = S0_T1
            T2_S0_map[xi,yi,i] = S0_T2
            FA_Eff_map[xi,yi,i] = FA_eff
            T1_rsquare_map[xi,yi,i] = r_squared_T1
            T2_rsquare_map[xi,yi,i] = r_squared_T2

    T1_S0_map_series = series_T1.SeriesDescription + "_T1_" + "S0_Map_v2"
    T1_S0_map_series = series_T1.new_series(SeriesDescription=T1_S0_map_series)
    T1_S0_map_series.set_array(np.squeeze(T1_S0_map),np.squeeze(header_T1[:,0]),pixels_first=True)

    T1_map_series = series_T1.SeriesDescription + "_T1_" + "T1_Map_v2"
    T1_map_series = series_T1.new_series(SeriesDescription=T1_map_series)
    T1_map_series.set_array(np.squeeze(T1_map),np.squeeze(header_T1[:,0]),pixels_first=True)

    FA_Eff_map_series = series_T1.SeriesDescription + "_T1_" + "FA_Eff_Map_v2"
    FA_Eff_map_series = series_T1.new_series(SeriesDescription=FA_Eff_map_series)
    FA_Eff_map_series.set_array(np.squeeze(FA_Eff_map),np.squeeze(header_T1[:,0]),pixels_first=True)

    T2_S0_map_series = series_T1.SeriesDescription + "_T2_" + "S0_Map_v2"
    T2_S0_map_series = series_T1.new_series(SeriesDescription=T2_S0_map_series)
    T2_S0_map_series.set_array(np.squeeze(T2_S0_map),np.squeeze(header_T2[:,0]),pixels_first=True)

    T2_map_series = series_T1.SeriesDescription + "_T2_" + "T2_Map_v2"
    T2_map_series = series_T1.new_series(SeriesDescription=T2_map_series)
    T2_map_series.set_array(np.squeeze(T2_map),np.squeeze(header_T2[:,0]),pixels_first=True)

    T1_rsquare_map_series = series_T1.SeriesDescription + "_T1_" + "rsquare_Map_v2"
    T1_rsquare_map_series = series_T1.new_series(SeriesDescription=T1_rsquare_map_series)
    T1_rsquare_map_series.set_array(np.squeeze(T1_rsquare_map),np.squeeze(header_T1[:,0]),pixels_first=True)

    T2_rsquare_map_series = series_T1.SeriesDescription + "_T2_" + "rsquare_Map_v2"
    T2_rsquare_map_series = series_T1.new_series(SeriesDescription=T2_rsquare_map_series)
    T2_rsquare_map_series.set_array(np.squeeze(T2_rsquare_map),np.squeeze(header_T2[:,0]),pixels_first=True)

def T1_then_T2(series_T1_T2, series_mask, study=None):

    #Prepare T1w/T2w arrays
    series_T1 = series_T1_T2[0]
    series_T2 = series_T1_T2[1]

    array_T1, header_T1 = series_T1.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    array_T2, header_T2 = series_T2.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)

    array_T1 = np.squeeze(array_T1[:,:,:,:,0])
    array_T2 = np.squeeze(array_T2[:,:,:,:,0])

    header_T1 = np.squeeze(header_T1[:,...])
    header_T2 = np.squeeze(header_T2[:,...])
    #######################

    #Setup parameters
    TR = 4.6                            #in ms
    FA = 12    #in degrees
    FA_rad = FA/360*(2*np.pi)           #convert to rads
    N_T1 = 66                           #number of k-space lines
    FA_Cat  = [(-FA/5)/360*(2*np.pi), (2*FA/5)/360*(2*np.pi), (-3*FA/5)/360*(2*np.pi), (4*FA/5)/360*(2*np.pi), (-5*FA/5)/360*(2*np.pi)] #cat module

    TE = [0,30,40,50,60,70,80,90,100,110,120]
    Tspoil = 1
    N_T2 = 72
    Trec = 463*2
    FA_eff = 0.6

    #initialize result arrays
    T1_S0_map = np.zeros(np.shape(array_T1)[0:3])
    T1_map = np.zeros(np.shape(array_T1)[0:3])
    FA_Eff_map = np.zeros(np.shape(array_T1)[0:3])
    Ref_Eff_map = np.zeros(np.shape(array_T1)[0:3])
    T2_S0_map = np.zeros(np.shape(array_T1)[0:3])
    T2_map = np.zeros(np.shape(array_T1)[0:3])
    T1_rsquare_map = np.zeros(np.shape(array_T1)[0:3])
    T2_rsquare_map = np.zeros(np.shape(array_T1)[0:3])

    #perform T1 mapping slice by slice
    for i in range(np.shape(array_T1)[2]):
        Kidney_pixel_T1 = np.squeeze(array_T1[...,i,:])

        if np.size(np.shape(np.squeeze(header_T1)))==2:
            TI_temp =  [float(hdr['InversionTime']) for hdr in header_T1[i,:]]
        elif np.size(np.shape(np.squeeze(header_T1)))==3:
            TI_temp =  [float(hdr['InversionTime']) for hdr in header_T1[i,:,0]]

        pool = multiprocessing.Pool(processes=os.cpu_count())

        arguments =[]
        pool = multiprocessing.Pool(initializer=multiprocessing.freeze_support,processes=os.cpu_count())
        for (x, y), _ in np.ndenumerate(Kidney_pixel_T1[..., 0]):
            t1_value = Kidney_pixel_T1[x, y, :]

            arguments.append((x,y,t1_value,TI_temp,FA_rad,TR,N_T1,FA_Cat,FA_eff))
            
        results = list(tqdm(pool.imap(t1_alone.main, arguments), total=len(arguments), desc='Processing pixels of slice ' + str(i)))

        for result in results:
            xi = result[0]
            yi = result[1]
            T1 = result[2]
            S0_T1 = result[3]
            FA_eff = result[4]
            r_squared_T1 = result[5]

            T1_map[xi,yi,i] = T1
            T1_S0_map[xi,yi,i] = S0_T1
            FA_Eff_map[xi,yi,i] = FA_eff
            T1_rsquare_map[xi,yi,i] = r_squared_T1

    #Coregiser T1 map to dixon mask using active alignment 
    T1_map_coreg_series,params_T1map = map_to_dixon_mask.main(T1_map,series_mask)
    T1_map_coreg_array, T1_map_coreg_header = T1_map_coreg_series.array(['SliceLocation'], pixels_first=True)
    T1_map_coreg_array = np.squeeze(T1_map_coreg_array)

    #fill coreg T1 map array usint scipy interpret
    mask_array, mask_header = series_mask.array(['SliceLocation'], pixels_first=True)
    mask_array = np.squeeze(mask_array)
    T1_map_coreg_array_filled = fill_kidney_holes_interp_v2.main(T1_map_coreg_array, mask_array)

    #Split T2w images by Echo Time
    T2w_splitted_echoes = series_T2.split_by("Echo Time")
    
    #### JUST TO CONFIRM THIS IS THE FIRST ECHO ####
    #array_T2_TE0, header_T2_TE0 = T2w_splitted_echoes[0].array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    #print(header_T2_TE0['Echo Time'])
    ################################################

    #Coregister T2w images to dixon mask using active alignment
    T2w,params_T2w = map_to_dixon_mask.main(T2w_splitted_echoes[0],series_mask) 
    T2w_imgs_coreg_series = []
    for echo_series in T2w_splitted_echoes:
        T2w_imgs_coreg_series.append(map_to_dixon_mask.main(echo_series,series_mask,params_T2w))

    #fill coreg T2w images array usint scipy interpret (by T2_prep)
    T2w_imgs_coreg_array_filled = np.zeros(T2w_imgs_coreg_array.shape)
    for series in T2w_imgs_coreg_series:
        T2w_imgs_coreg_array, T2w_imgs_coreg_header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
        T2w_imgs_coreg_array = np.squeeze(T2w_imgs_coreg_array)
        for T2_prep in range(np.shape(T2w_imgs_coreg_array)[3]):
            T2w_imgs_coreg_array_temp = np.squeeze(T2w_imgs_coreg_array[:,:,:,T2_prep])
            T2w_imgs_coreg_array_filled_temp = fill_kidney_holes_interp_v2.main(T2w_imgs_coreg_array_temp, mask_array)
            T2w_imgs_coreg_array_filled[:,:,:,T2_prep] = T2w_imgs_coreg_array_filled_temp

    #perform T2 mapping slice by slice
    for i in range(np.shape(T2w_imgs_coreg_array_filled)[2]):

        Kidney_pixel_T2 = np.squeeze(T2w_imgs_coreg_array_filled[...,i,:])

        arguments =[]
        pool = multiprocessing.Pool(initializer=multiprocessing.freeze_support,processes=os.cpu_count())
        for (x, y), _ in np.ndenumerate(Kidney_pixel_T1[..., 0]):
            t1_map =   T1_map_coreg_array_filled[x,y,i]
            t2_value = Kidney_pixel_T2[x, y, :]

            arguments.append((x,y,t1_map,t2_value,TE,FA_rad,TR,N_T2,Trec,FA_eff,Tspoil))

        results = list(tqdm(pool.imap(t1_t2_alone.main, arguments), total=len(arguments), desc='Processing pixels of slice ' + str(i)))

        for result in results:
            xi = result[0]
            yi = result[1]
            T2 = result[2]
            S0_T2 = result[3]
            FA_eff = result[4]
            r_squared_T2 = result[5]
            T2_map[xi,yi,i] = T2
            T2_S0_map[xi,yi,i] = S0_T2
            FA_Eff_map[xi,yi,i] = FA_eff
            T2_rsquare_map[xi,yi,i] = r_squared_T2

    #save calculated series
    T1_S0_map_series = series_T1.SeriesDescription + "_T1_" + "S0_Map_v2"
    T1_S0_map_series = series_T1.new_series(SeriesDescription=T1_S0_map_series)
    T1_S0_map_series.set_array(np.squeeze(T1_S0_map),np.squeeze(header_T1[:,0]),pixels_first=True)

    T1_map_series = series_T1.SeriesDescription + "_T1_" + "T1_Map_v2"
    T1_map_series = series_T1.new_series(SeriesDescription=T1_map_series)
    T1_map_series.set_array(np.squeeze(T1_map),np.squeeze(header_T1[:,0]),pixels_first=True)

    FA_Eff_map_series = series_T1.SeriesDescription + "_T1_" + "FA_Eff_Map_v2"
    FA_Eff_map_series = series_T1.new_series(SeriesDescription=FA_Eff_map_series)
    FA_Eff_map_series.set_array(np.squeeze(FA_Eff_map),np.squeeze(header_T1[:,0]),pixels_first=True)

    T2_S0_map_series = series_T1.SeriesDescription + "_T2_" + "S0_Map_v2"
    T2_S0_map_series = series_T1.new_series(SeriesDescription=T2_S0_map_series)
    T2_S0_map_series.set_array(np.squeeze(T2_S0_map),np.squeeze(header_T2[:,0]),pixels_first=True)

    T2_map_series = series_T1.SeriesDescription + "_T2_" + "T2_Map_v2"
    T2_map_series = series_T1.new_series(SeriesDescription=T2_map_series)
    T2_map_series.set_array(np.squeeze(T2_map),np.squeeze(header_T2[:,0]),pixels_first=True)

    T1_rsquare_map_series = series_T1.SeriesDescription + "_T1_" + "rsquare_Map_v2"
    T1_rsquare_map_series = series_T1.new_series(SeriesDescription=T1_rsquare_map_series)
    T1_rsquare_map_series.set_array(np.squeeze(T1_rsquare_map),np.squeeze(header_T1[:,0]),pixels_first=True)

    T2_rsquare_map_series = series_T1.SeriesDescription + "_T2_" + "rsquare_Map_v2"
    T2_rsquare_map_series = series_T1.new_series(SeriesDescription=T2_rsquare_map_series)
    T2_rsquare_map_series.set_array(np.squeeze(T2_rsquare_map),np.squeeze(header_T2[:,0]),pixels_first=True)

def IVIM(series=None, mask=None,export_ROI=False, study = None):

        series_IVIM = series

        array, header = series_IVIM.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
        b_vals = [0,10.000086, 19.99908294, 30.00085926, 50.00168544, 80.007135, 100.0008375, 199.9998135, 300.0027313, 600.0]

        pixel_array_IVIM = array

        if mask is not None:
                mask=np.transpose(mask)
                for i_slice in range (np.shape(pixel_array_IVIM)[2]):
                    for i_w in range (np.shape(pixel_array_IVIM)[3]):
                        pixel_array_IVIM[:,:,i_slice,i_w]=pixel_array_IVIM[:,:,i_slice,i_w]*mask

        S0map, Dmap,rsquaremap = mapping.IVIM.main(pixel_array_IVIM,b_vals)

        S0_map_series = series_IVIM.SeriesDescription + "_IVIM_" + "S0_Map"
        S0_map_series = study.new_series(SeriesDescription=S0_map_series)
        S0_map_series.set_array(np.squeeze(S0map),np.squeeze(header[:,0]),pixels_first=True)
        
        D_map_series = series_IVIM.SeriesDescription + "_IVIM_" + "D_Map"
        D_map_series = study.new_series(SeriesDescription=D_map_series)
        D_map_series.set_array(np.squeeze(Dmap),np.squeeze(header[:,0]),pixels_first=True)

        rsquare_map_series = series_IVIM.SeriesDescription + "_IVIM_" + "rsquare_Map"
        rsquare_map_series = study.new_series(SeriesDescription=rsquare_map_series)
        rsquare_map_series.set_array(np.squeeze(rsquaremap),np.squeeze(header[:,0]),pixels_first=True)

def DTI(series=None, mask=None,export_ROI=False, study = None):

    # # Added this as the original version (below) gave an error
    # fit, par = fit_DTI(series)
    # return par[0], par[1]

    # Previous code

    start_time = time.time()
    series.log("DTI-FA & ADC mapping has started")

    series_DTI = series

    array, header = series_DTI.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    pixel_array_DTI = array
    header = np.squeeze(header)
    
    b_vals_check = [float(hdr[(0x19, 0x100c)]) for hdr in header[0,:]]
    b_vecs_check = [hdr[(0x19, 0x100e)] for hdr in header[0,:]]

    #Check if the data corresponds to the Siemens protocol (more than 1 unique b-value)     
    return_vals = None, None   
    if len(b_vals_check) >= 1 and np.shape(b_vecs_check)[0] >=6:

######FROM DIPY
        gtab = gradient_table(np.squeeze(b_vals_check), np.squeeze(b_vecs_check))
        tenmodel = dti.TensorModel(gtab)
        tenfit = tenmodel.fit(np.squeeze(pixel_array_DTI))

        FAmap = fractional_anisotropy(tenfit.evals)
        MDmap = dti.mean_diffusivity(tenfit.evals)
######FROM DIPY          

        MDmap[MDmap>0.005]=0.005
        MDmap[MDmap<0]=0

        FA_map_series = series_DTI.SeriesDescription + "_DTI_" + "FA_Map"
        FA_map_series = study.new_series(SeriesDescription=FA_map_series)
        FA_map_series.set_array(np.squeeze(FAmap),np.squeeze(header[:,0]),pixels_first=True)

        MD_map_series = series_DTI.SeriesDescription + "_DTI_" + "MD_Map"
        MD_map_series = study.new_series(SeriesDescription=MD_map_series)
        MD_map_series.set_array(np.squeeze(MDmap),np.squeeze(header[:,0]),pixels_first=True)

        return_vals = FA_map_series, MD_map_series

    series.log("DTI-FA & ADC mapping was completed --- %s seconds ---" % (int(time.time() - start_time))) 
    
    return return_vals

def MTR(series=None, mask=None, export_ROI=False, study=None):

    start_time = time.time()
    series.log("MTR mapping has started")
        
    array_mt_moco, header_mt_moco = series.array(['SliceLocation', 'AcquisitionTime'],pixels_first=True)
    header_on = header_mt_moco[:,0,0]

    array_mt_off = np.squeeze(array_mt_moco[:,:,:,0])
    array_mt_on = np.squeeze(array_mt_moco[:,:,:,1])
    array_mtr = np.zeros((np.shape(array_mt_off)[0:3]))
    
    for s in range (np.shape(array_mt_off)[2]):
        temp_off_moco = np.squeeze(array_mt_off[:,:,s])
        temp_on_moco = np.squeeze(array_mt_on[:,:,s])

        array_mtr[:,:,s] = np.divide((temp_off_moco - temp_on_moco),temp_off_moco, out=np.zeros_like(temp_off_moco - temp_on_moco), where=temp_off_moco!=0) * 100
    
    study = series.parent()
    mtr = series.SeriesDescription + '_MTR'
    mtr = study.new_series(SeriesDescription=mtr)
    mtr.set_array(array_mtr, np.squeeze(header_on), pixels_first=True)

    series.log("MTR mapping was completed --- %s seconds ---" % (int(time.time() - start_time))) 
    return mtr


def DCE_MAX(series=None, mask=None,export_ROI=False, study=None):

    start_time = time.time()
    series.log("DCE-MAX mapping has started")

    series_DCE = series

    array, header = series_DCE.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    
    header = np.squeeze(header)
    if np.shape(array)[-1] == 2:
        array = array[...,0]

    if np.shape(header)[-1] == 2:
        header = header[...,0:1]
        header = np.squeeze(header)

    pixel_array_DCE = array
    number_slices = np.shape(pixel_array_DCE)[2]

    DCE_Max_map = np.empty(np.shape(pixel_array_DCE)[0:3])
    DCE_Area_map = np.empty(np.shape(pixel_array_DCE)[0:3])

    timeDCE = np.zeros(header.shape[1])

    for slice in range(number_slices):   
        for i_2 in range(header.shape[1]):
            if len (header.shape) == 2:
                tempTime = str(header[slice,i_2]['AcquisitionTime'])
            elif len (header.shape==3):
                tempTime = str(header[slice,i_2,0]['AcquisitionTime'])

            beforepoint = tempTime.split(".")[0]
            afterpoint = tempTime.split(".")[1]
            tempH = int(beforepoint[0:2])
            tempM = int(beforepoint[2:4])
            tempS = int(beforepoint[4:])
            tempRest = float("0." + afterpoint)
            timeDCE[i_2] = tempH*3600+tempM*60+tempS+tempRest
        timeDCE -=timeDCE[0]

        array_DCE_temp = np.squeeze(pixel_array_DCE[:,:,slice,:])

        for xi in range((np.size(array_DCE_temp,0))):
            for yi in range((np.size(array_DCE_temp,1))):

                
                Kidney_pixel_DCE = np.squeeze(np.array(array_DCE_temp[xi,yi,:]))
                DCE_Max_map[xi,yi,slice] = np.max(Kidney_pixel_DCE-np.mean(Kidney_pixel_DCE[0:11]))
                DCE_Area_map[xi,yi,slice] = trapz(Kidney_pixel_DCE,timeDCE)

    DCEMax_map_series = series_DCE.SeriesDescription + "_DCE_" + "Max_Map"
    DCEMax_map_series = study.new_series(SeriesDescription=DCEMax_map_series)
    DCEMax_map_series.set_array(np.squeeze(DCE_Max_map),np.squeeze(header[:,0]),pixels_first=True)

    DCEArea_map_series = series_DCE.SeriesDescription + "_DCE_" + "AUC_Map"
    DCEArea_map_series = study.new_series(SeriesDescription=DCEArea_map_series)
    DCEArea_map_series.set_array(np.squeeze(DCE_Area_map),np.squeeze(header[:,0]),pixels_first=True)

    series.log("DCE-MAX mapping was completed --- %s seconds ---" % (int(time.time() - start_time)))

    return DCEMax_map_series, DCEArea_map_series


def main(folder):

    start_time = time.time()
    folder.log("AI segmentation has started!")
    weights = 'UNETR_kidneys_v1.pth'
    seg.segment_kidneys(folder, weights)
    folder.scan()

    folder.log("Modelling has started!")

    list_of_series = folder.series()

    current_study = list_of_series[0].parent()
    study = list_of_series[0].new_pibling(StudyDescription=current_study.StudyDescription + '_ModellingResults')

    for i,series in enumerate(list_of_series):
        print(series['SeriesDescription'])
        if series["SequenceName"] is not None:

            if series['SeriesDescription'] == "T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco":
                try:
                    print('Starting T2s')
                    series.log("T2s mapping has started")
                    #T2s(series, study=study)
                except Exception as e: 
                    series.log("T2* mapping was NOT completed; error: "+str(e))

            elif series['SeriesDescription'] == "DTI_kidneys_cor-oblique_fb_mdr_moco":
                try:
                    print('Starting DTI')
                    series.log("DTI mapping has started")
                    #DTI(series, study=study)
                except Exception as e: 
                    series.log("DTI-FA & ADC mapping was NOT completed; error: "+str(e))

            elif series['SeriesDescription'] == "DCE_kidneys_cor-oblique_fb_mdr_moco":
                try:
                    print('Starting DCE_MAX')
                    series.log("DCE mapping has started")
                    #DCE_MAX(series, study=study)
                except Exception as e: 
                    series.log("DCE-MAX mapping was NOT completed; error: "+str(e))

            elif series.SeriesDescription == 'MT_ON_kidneys_cor-oblique_bh_mdr_moco':
                try:
                    print('Starting MTR')
                    series.log("MTR mapping has started")
                    #MTR(series, study=study)
                except Exception as e: 
                    series.log("MTR mapping was NOT completed; error: "+str(e))

            elif series.SeriesDescription == 'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco':
                print(series.Manufacturer)
                if series.Manufacturer != 'SIEMENS': # TODO: Check this, something not right
                    try:
                        T1_Philips(series, study=study)  
                    except Exception as e: 
                        series.log("T1 mapping was NOT completed; error: "+str(e))
                else:
                    T1 = series
                    for i_2,series in enumerate(list_of_series):
                        print(series['SeriesDescription'])
                        if series['SeriesDescription'] == "T2map_kidneys_cor-oblique_mbh_magnitude_mdr_moco":
                            T2 = series
                            for i_2,series in enumerate(list_of_series):
                                print(series['SeriesDescription'])
                                if series['SeriesDescription'] == "LK":
                                    Kidney_mask = series
                                    series.log("T1, T2 mapping has started")
                                    T1_then_T2([T1,T2],Kidney_mask, study=study)
                            break

    folder.save()
    folder.log("Modelling was completed --- %s seconds ---" % (int(time.time() - start_time)))
