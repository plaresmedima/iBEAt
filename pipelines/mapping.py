import numpy as np


from dbdicom.pipelines import input_series

from pipelines.roi_fit import load_aif

import models
from models import (
    T1_look_locker_spoiled,
    T2_mono_exp,
    T2star_mono_exp,
    t1_t2_combined, 
    IVIM_nonlinear, 
    IVIM_semilinear,
    DTI_dipy, 
    DWI_linear,
    DWI_nonlinear,
    DCE_descriptive, 
    DCE_deconv, 
    DCE_2CM,
)



# from utilities import map_to_dixon_mask
# from utilities import fill_kidney_holes_interp_v2

export_study = 'Parameter maps'


def T1map(folder):

    # Find source DICOM data
    desc = "T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        desc = "T1map_kidneys_cor-oblique_mbh_magnitude"
        series, study = input_series(folder, desc, export_study)
        if series is None:
            raise RuntimeError('Cannot perform T1 mapping: series ' + desc + 'does not exist. ')

    # Load data - note each slice can have different TIs
    array, header = series.array(['SliceLocation', 'InversionTime'], pixels_first=True, first_volume=True)
    TI = [np.array([hdr['InversionTime'] for hdr in header[z,:]]) for z in range(array.shape[2])] 

    # Calculate fit slice by slice
    fit = np.empty(array.shape)
    par = np.empty(array.shape[:3] + (3,) )
    for z in range(array.shape[2]):
        series.progress(z+1, array.shape[2], 'Fitting T1 model')
        fit[:,:,z,:], par[:,:,z,:] = T1_look_locker_spoiled.fit(array[:,:,z,:], TI[z], xtol=1e-3, bounds=True)

    # Derive error map
    ref = np.linalg.norm(array, axis=-1)
    err = 100*np.linalg.norm(fit-array, axis=-1)/ref
    err[ref==0] = 0

    # Save results to DICOM database
    series = study.new_series(SeriesDescription=desc + "_" + 'fit')
    series.set_array(fit, header, pixels_first=True)
    series = study.new_series(SeriesDescription=desc+'_err_map')
    series.set_array(err, header[:,0], pixels_first=True)
    for i, p in enumerate(T1_look_locker_spoiled.pars()):
        series = study.new_series(SeriesDescription=desc+'_' + p + '_map')
        series.set_array(par[...,i], header[:,0], pixels_first=True)

    return (series,)


def T2map(folder):

    # Find source DICOM data
    desc = "T2map_kidneys_cor-oblique_mbh_magnitude_mdr_moco"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        desc = "T2map_kidneys_cor-oblique_mbh_magnitude"
        series, study = input_series(folder, desc, export_study)
        if series is None:
            raise RuntimeError('Cannot perform T2 mapping: series ' + desc + 'does not exist. ')

    # Load data 
    array, header = series.array(['SliceLocation', 'InversionTime'], pixels_first=True, first_volume=True)
    TI = series.values('InversionTime', dims=('SliceLocation', 'InversionTime'))

    # Calculate fit slice by slice
    fit = np.empty(array.shape)
    par = np.empty(array.shape[:3] + (3,) )
    for z in range(array.shape[2]):
        series.progress(z+1, array.shape[2], 'Fitting T2 model')
        fit[:,:,z,:], par[:,:,z,:] = T2_mono_exp.fit(array[:,:,z,:], TI[z,:], xtol=1e-3, bounds=True)

    # Derive error map
    ref = np.linalg.norm(array, axis=-1)
    err = 100*np.linalg.norm(fit-array, axis=-1)/ref
    err[ref==0] = 0

    # Save results to DICOM database
    series = study.new_series(SeriesDescription=desc + "_" + 'fit')
    series.set_array(fit, header, pixels_first=True)
    series = study.new_series(SeriesDescription=desc+'_err_map')
    series.set_array(err, header[:,0], pixels_first=True)
    for i, p in enumerate(T2_mono_exp.pars()):
        series = study.new_series(SeriesDescription=desc+'_' + p +'_map')
        series.set_array(par[...,i], header[:,0], pixels_first=True)
    return series


def T1T2(folder):

    if folder.Manufacturer != 'SIEMENS': 
        raise ValueError('Only Siemens implementation available at the moment')

    desc = [
        'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco', 
        'T2map_kidneys_cor-oblique_mbh_magnitude_mdr_moco']
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform T1 and T2 mapping: series ' + str(desc) + 'do not exist. ')

    array_T1, header_T1 = series[0].array(['SliceLocation', 'AcquisitionTime'], pixels_first=True, first_volume=True)
    array_T2, header_T2 = series[1].array(['SliceLocation', 'AcquisitionTime'], pixels_first=True, first_volume=True)
    TI_slices = [ [hdr['InversionTime'] for hdr in header_T1[z,:]] for z in range(header_T1.shape[0])]

    par = t1_t2_combined.fit((array_T1, array_T2), TI_slices)
    pnames = t1_t2_combined.pars()
    maps = []
    for p in range(len(pnames)):
        series = study.new_series(SeriesDescription='T1T2_' + pnames[p] + '_map')
        series.set_array(par[...,p], header_T1[:,0], pixels_first=True)
        maps.append(series)
    return maps[0], maps[1]


def T2starmap(folder):
    
    # Find source DICOM data
    desc = "T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        desc = "T2star_map_kidneys_cor-oblique_mbh_magnitude"
        series, study = input_series(folder, desc, export_study)
        if series is None:
            raise RuntimeError('Cannot perform T2* mapping: series ' + desc + 'does not exist. ')

    array, header = series.array(['SliceLocation', 'EchoTime'], pixels_first=True, first_volume=True)
    echo_times = [hdr["EchoTime"] for hdr in (header[0,:])] 
    
    # Calculate fit slice by slice
    fit = np.empty(array.shape)
    par = np.empty(array.shape[:3] + (3,) )
    for z in range(array.shape[2]):
        series.progress(z+1, array.shape[2], 'Fitting T2* model')
        fit[:,:,z,:], par[:,:,z,:] = T2star_mono_exp.fit(array[:,:,z,:], echo_times, xtol=1e-3, bounds=True)

    # Derive error map
    ref = np.linalg.norm(array, axis=-1)
    err = 100*np.linalg.norm(fit-array, axis=-1)/ref
    err[ref==0] = 0

    # Save results to DICOM database
    series = study.new_series(SeriesDescription=desc + "_" + 'fit')
    series.set_array(fit, header, pixels_first=True)
    series = study.new_series(SeriesDescription=desc+'_err_map')
    series.set_array(err, header[:,0], pixels_first=True)
    for i, p in enumerate(T2star_mono_exp.pars()):
        series = study.new_series(SeriesDescription=desc+'_' + p + '_map')
        series.set_array(par[...,i], header[:,0], pixels_first=True)
    return series


def MT(folder):
    
    # Find source DICOM data
    desc = "MT_kidneys_cor-oblique_bh_mdr_moco"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        desc = "MT_kidneys_cor-oblique_bh"
        series, study = input_series(folder, desc, export_study)
        if series is None:
            raise RuntimeError('Cannot perform MT mapping: series ' + desc + 'does not exist. ')
        
    array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True, first_volume=True)

    # Calculate
    array_mtr = 100*(array[...,0]-array[...,1])/array[...,0]
    array_mtr[array[...,0]==0] = 0
    array_mtr[array_mtr>100] = 100
    array_mtr[array_mtr<-100] = -100
    array_avr = np.mean(array, axis=-1)
    
    # Save as DICOM
    mtr = study.new_series(SeriesDescription=desc + '_MTR_map')
    mtr.set_array(array_mtr, header[:,0], pixels_first=True)
    avr = study.new_series(SeriesDescription=desc + '_AVR_map')
    avr.set_array(array_avr, header[:,0], pixels_first=True)
    
    return mtr


def DTI(folder):

    # Find appropriate series
    desc = "DTI_kidneys_cor-oblique_fb_mdr_moco"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        desc = "DTI_kidneys_cor-oblique_fb"
        series, study = input_series(folder, desc, export_study)
        if series is None:
            raise RuntimeError('Cannot perform mapping on DTI: series ' + desc + 'does not exist. ')

    # Read data
    array, header = series.array(['SliceLocation', 'InstanceNumber'], pixels_first=True, first_volume=True)
    bvals, bvecs = series.values('DiffusionBValue', 'DiffusionGradientOrientation', dims=('SliceLocation', 'InstanceNumber'))

    # Compute
    series.message('Fitting DTI model..')
    fit, pars = DTI_dipy.fit(array, bvals[0,:], np.stack(bvecs[0,:]), fit_method='NLLS')
    #fit, pars = DTI_dipy.fit(array, bvals[0,:], np.stack(bvecs[0,:]), fit_method='WLS') # for debugging
    ref = np.linalg.norm(array, axis=-1)
    err = 100*np.linalg.norm(fit-array, axis=-1)/ref
    err[ref==0] = 0

    # Save as DICOM
    series = study.new_series(SeriesDescription=desc + '_fit')
    series.set_array(fit, header, pixels_first=True)
    series = study.new_series(SeriesDescription=desc + '_fiterr_map')
    series.set_array(err, header[:,0], pixels_first=True)
    maps = []
    for i, p in enumerate(DTI_dipy.pars()):
        series = study.new_series(SeriesDescription=desc + '_' + p +'_map')
        series.set_array(pars[...,i], header[:,0], pixels_first=True)
        maps.append(series)

    return maps[0], maps[1]


def ivim(folder):

    # Find appropriate series
    desc = "IVIM_kidneys_cor-oblique_fb_mdr_moco"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        desc = 'IVIM_kidneys_cor-oblique_fb'
        series, study = input_series(folder, desc, export_study)
        if series is None:
            raise RuntimeError('Cannot perform mapping on IVIM: series ' + desc + 'does not exist. ')

    array, header = series.array(['SliceLocation', 'InstanceNumber'], pixels_first=True, first_volume=True)
    bvals, bvecs = series.values('DiffusionBValue', 'DiffusionGradientOrientation', dims=('SliceLocation', 'InstanceNumber'))

    # Calculate fit in one line
    # fit, par = fit_image(IVIM_nonlinear, array, bvals[0,:], xtol=1e-3, bounds=True)

    # Calculate fit slice by slice so progress bar can be shown
    fit = np.empty(array.shape)
    par = np.empty(array.shape[:3] + (len(IVIM_nonlinear.pars()),) )
    for z in range(array.shape[2]):
        series.progress(z+1, array.shape[2], 'Fitting IVIM model')

        fit[:,:,z,:], par[:,:,z,:] = IVIM_nonlinear.fit(array[:,:,z,:], 
            bvals[0,:10].astype(np.float32), xtol=1e-3, bounds=True, parallel=True)
        
    S0, Df, MD, ff = IVIM_nonlinear.derived(par)

    # Save as DICOM
    fit_series = study.new_series(SeriesDescription=desc + "_fit")
    fit_series.set_array(fit, header, pixels_first=True)
    Df_series = study.new_series(SeriesDescription=desc + "_Df_map")
    Df_series.set_array(Df, header[:,0], pixels_first=True)
    MD_series = study.new_series(SeriesDescription=desc + "_MD_map")
    MD_series.set_array(MD, header[:,0], pixels_first=True)
    S0_series = study.new_series(SeriesDescription=desc + "_S0_map")
    S0_series.set_array(S0, header[:,0], pixels_first=True)
    ff_series = study.new_series(SeriesDescription=desc + "_ff_map")
    ff_series.set_array(ff, header[:,0], pixels_first=True)

    return Df_series, ff_series


def DCE(folder):

    # Find appropriate series
    desc = "DCE_kidneys_cor-oblique_fb_mdr_moco"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        desc = "DCE_kidneys_cor-oblique_fb"
        series, study = input_series(folder, desc, export_study)
        if series is None:
            raise RuntimeError('Cannot perform mapping on DCE: series ' + desc + ' does not exist. ')

    # Extract data
    time, aif = load_aif(folder)
    array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True, first_volume=True)
    
    # Calculate maps
    series.message('Calculating descriptive parameters..')
    MAX, AUC, ATT = DCE_descriptive.fit(array, aif, time[1]-time[0], baseline=15)
    series.message('Performing linear fit..')
    fit, par = DCE_2CM.fit(array, aif, time, baseline=15)
    series.message('Deconvolving..')
    RPF, AVD, MTT = DCE_deconv.fit(array, aif, time[1]-time[0], baseline=15, regpar=0.1)
    
    # Save maps as DICOM
    fit_series = study.new_series(SeriesDescription=desc + "_fit")
    fit_series.set_array(fit, header, pixels_first=True)
    
    MAX_series = study.new_series(SeriesDescription=desc + "_MAX_map")
    AUC_series = study.new_series(SeriesDescription=desc + "_AUC_map")
    ATT_series = study.new_series(SeriesDescription=desc + "_ATT_map")
    RPF_series = study.new_series(SeriesDescription=desc + "_RPF_map")
    AVD_series = study.new_series(SeriesDescription=desc + "_AVD_map")
    MTT_series = study.new_series(SeriesDescription=desc + "_MTT_map")
    FP_series = study.new_series(SeriesDescription=desc + "_FP_map")
    TP_series = study.new_series(SeriesDescription=desc + "_TP_map")
    VP_series = study.new_series(SeriesDescription=desc + "_VP_map")
    FT_series = study.new_series(SeriesDescription=desc + "_FT_map")
    TT_series = study.new_series(SeriesDescription=desc + "_TT_map")

    MAX_series.set_array(MAX, header[:,0], pixels_first=True)
    AUC_series.set_array(AUC, header[:,0], pixels_first=True)
    ATT_series.set_array(ATT, header[:,0], pixels_first=True)
    RPF_series.set_array(RPF, header[:,0], pixels_first=True)
    AVD_series.set_array(AVD, header[:,0], pixels_first=True)
    MTT_series.set_array(MTT, header[:,0], pixels_first=True)
    FP_series.set_array(par[...,0], header[:,0], pixels_first=True)
    TP_series.set_array(par[...,1], header[:,0], pixels_first=True)
    VP_series.set_array(par[...,0]*par[...,1]/60, header[:,0], pixels_first=True)
    FT_series.set_array(par[...,2], header[:,0], pixels_first=True)
    TT_series.set_array(par[...,3], header[:,0], pixels_first=True)

    return RPF_series, AVD_series, MTT_series



