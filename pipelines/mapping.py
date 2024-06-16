import numpy as np

from dbdicom.pipelines import input_series
import dcmri

from pipelines.roi_fit import load_aif

import models.t1
import models.t2
import models.T2star
import models.DTI
import models.IVIM

export_study = '2: Parameter maps'


def T1(folder):

    # Find source DICOM data
    series, study, desc = _map_input(folder, "T1m_magnitude")

    # Model parameters (Siemens)
    TR = 4.6 # Echo Spacing is msec but not in header -> Set as TR in harmonize
    FA_cat = [-1, 2, -3, 4, -5] # Catalization module confirmed by Siemens (Peter Schmitt): Magn Reson Med 2003 Jan;49(1):151-7. doi: 10.1002/mrm.10337
    N_T1 = 66 # Number of k-space lines (hardcoded from Siemens protocol)
    FA_nom = 12 # Flip angle in degrees (hardcoded from Siemens protocol)

    model = models.t1.Bloch()
    kwargs = {'TR':TR, 'FA_cat':FA_cat, 'N_T1':N_T1, 'FA':FA_nom}
    dims = ['SliceLocation', 'AcquisitionTime']

    # Load data - note each slice can have different TIs
    array, header = series.array(dims, pixels_first=True, first_volume=True)
    TI = series.values('InversionTime', dims=tuple(dims)).astype(np.float32)

    _map(series, array, header, model, TI, study, desc, **kwargs)



def T2(folder):

    # Find source DICOM data
    series, study, desc = _map_input(folder, "T2m_magnitude")

    # Load data
    model = models.t2.MonoExp()
    dims = ['SliceLocation', 'AcquisitionTime']
    array, header = series.array(dims, pixels_first=True, first_volume=True)
    TI = series.values('InversionTime', dims=tuple(dims)).astype(np.float32)

    _map(series, array, header, model, TI, study, desc)



def T2star(folder):
    
    # Find source DICOM data
    series, study, desc = _map_input(folder, "T2starm_magnitude")

    model = models.T2star.BiExp()
    dims = ['SliceLocation', 'EchoTime']
    array, header = series.array(dims, pixels_first=True, first_volume=True)
    echo_times = series.values('EchoTime', dims=tuple(dims)).astype(np.float32)

    _map(series, array, header, model, echo_times, study, desc)


def MT(folder):
    
    # Find source DICOM data
    series, study, desc = _map_input(folder, "MT")
        
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
    
    return (mtr,)


def DTI(folder):

    # Find appropriate series
    series, study, desc = _map_input(folder, "DTI")

    # Read data
    model = models.DTI.DiPy()
    array, header = series.array(['SliceLocation', 'InstanceNumber'], pixels_first=True, first_volume=True)
    bvals, bvecs = series.values('DiffusionBValue', 'DiffusionGradientOrientation', dims=('SliceLocation', 'InstanceNumber'))

    # Compute
    series.message('Fitting DTI model..')
    fit, pars = model.fit(array, bvals[0,:], np.stack(bvecs[0,:]), fit_method='NLLS')
    #err = model.error(array, fit)

    # Save as DICOM
    series = study.new_series(SeriesDescription=desc + '_fit')
    series.set_array(fit, header, pixels_first=True)
    #series = study.new_series(SeriesDescription=desc + '_fiterr_map')
    #series.set_array(err, header[:,0], pixels_first=True)
    maps = []
    for i, p in enumerate(model.pars()):
        series = study.new_series(SeriesDescription=desc + '_' + p +'_map')
        series.set_array(pars[...,i], header[:,0], pixels_first=True)
        maps.append(series)

    return maps[0], maps[1]


def IVIM(folder):

    # Find appropriate series
    series, study, desc = _map_input(folder, "IVIM")

    model = models.IVIM.DiPy()
    array, header = series.array(['SliceLocation', 'InstanceNumber'], pixels_first=True, first_volume=True)
    bvals, bvecs = series.values('DiffusionBValue', 'DiffusionGradientOrientation', dims=('SliceLocation', 'InstanceNumber'))

    # Compute
    series.message('Fitting IVIM model..')
    fit, pars = model.fit(array, bvals[0,:], np.stack(bvecs[0,:]), fit_method='TRR')
    #err = model.error(array, fit)

    # Save as DICOM
    series = study.new_series(SeriesDescription=desc + '_fit')
    series.set_array(fit, header, pixels_first=True)
    #series = study.new_series(SeriesDescription=desc + '_fiterr_map')
    #series.set_array(err, header[:,0], pixels_first=True)



    # Calculate fit slice by slice so progress bar can be shown
    # fit = np.empty(array.shape)
    # par = np.empty(array.shape[:3] + (len(model.pars()),) )
    # for z in range(array.shape[2]):
    #     series.progress(z+1, array.shape[2], 'Fitting IVIM model')
    #     fit[:,:,z,:], par[:,:,z,:] = model.fit(array[:,:,z,:], bvals[0,:10].astype(np.float32), xtol=1e-3, bounds=True, parallel=True) 
    # S0, Df, D, ff = model.derived(par)

    # Save as DICOM
    # fit_series = study.new_series(SeriesDescription=desc + "_fit")
    # fit_series.set_array(fit, header, pixels_first=True)
    # Df_series = study.new_series(SeriesDescription=desc + "_Df_map")
    # Df_series.set_array(Df, header[:,0], pixels_first=True)
    # MD_series = study.new_series(SeriesDescription=desc + "_D_map")
    # MD_series.set_array(MD, header[:,0], pixels_first=True)
    # S0_series = study.new_series(SeriesDescription=desc + "_S0_map")
    # S0_series.set_array(S0, header[:,0], pixels_first=True)
    # ff_series = study.new_series(SeriesDescription=desc + "_ff_map")
    # ff_series.set_array(ff, header[:,0], pixels_first=True)

    # return Df_series, ff_series

    maps = []
    for i, p in enumerate(model.pars()):
        series = study.new_series(SeriesDescription=desc + '_' + p +'_map')
        series.set_array(pars[...,i], header[:,0], pixels_first=True)
        maps.append(series)

    return maps[0], maps[1]


def DCE(folder):

    # Find appropriate series       
    series, study, desc = _map_input(folder, "DCE")

    # Extract data
    time, aif = load_aif(folder)
    array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True, first_volume=True)
    
    # Calculate maps
    series.message('Calculating descriptive parameters..')
    MAX, AUC, ATT, SO = dcmri.pixel_descriptives(array, aif, time[1]-time[0], baseline=15)
    series.message('Performing linear fit..')
    fit, par = dcmri.pixel_2cfm_linfit(array, aif, time, baseline=15)
    series.message('Deconvolving..')
    RPF, AVD, MTT = dcmri.pixel_deconvolve(array, aif, time[1]-time[0], baseline=15, regpar=0.1)
    
    # Save maps as DICOM
    fit_series = study.new_series(SeriesDescription=desc + "_fit")
    fit_series.set_array(fit, header, pixels_first=True)

    MAX[MAX<0]=0
    MAX[MAX>10000]=10000
    AUC[AUC<0]=0
    AUC[AUC>10000]=10000
    ATT[ATT<0]=0
    ATT[ATT>10000]=10000
    RPF[RPF<0]=0
    RPF[RPF>10000]=10000
    AVD[AVD<0]=0
    AVD[AVD>10000]=10000
    MTT[MTT<0]=0
    MTT[MTT>10000]=10000

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



def _map_input(folder, desc):
    series, study = input_series(folder, desc + '_mdr_moco', export_study)
    if series is not None:
        return series, study, desc + '_mdr_moco'
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform mapping: series ' + desc + 'does not exist. ')
    return series, study, desc


def _map(series, array, header, model, pars, study, desc, **kwargs):

    # Calculate fit slice by slice
    fit = np.empty(array.shape)
    par = np.empty(array.shape[:3] + (len(model.pars()),) )
    for z in range(array.shape[2]):
        series.progress(z+1, array.shape[2], 'Fitting model')
        fit[:,:,z,:], par[:,:,z,:] = model.fit(array[:,:,z,:], pars[z,:], xtol=1e-3, bounds=True, **kwargs)
    err = model.error(array, fit)

    # Save results to DICOM database
    series = study.new_series(SeriesDescription=desc + "_" + 'fit')
    series.set_array(fit, header, pixels_first=True)
    series = study.new_series(SeriesDescription=desc+'_err_map')
    series.set_array(err, header[:,0], pixels_first=True)
    for i, p in enumerate(model.pars()):
        series = study.new_series(SeriesDescription=desc+'_' + p + '_map')
        series.set_array(par[...,i], header[:,0], pixels_first=True)

    return (series,)




