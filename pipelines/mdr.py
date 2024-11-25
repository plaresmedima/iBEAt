import os
import numpy as np

from dbdicom.pipelines import input_series

import itk
import mdreg
import dcmri
from mdreg.models import constant

import models.DWI
import models.t1
import models.t1_t2
import models.t2
import models.T2star
import models.DTI
import models.IVIM

from pipelines.roi_fit import load_aif

export_study = '1: MDR results'


def T1(folder):

    desc = "T1m_magnitude"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on T1: series ' + desc + 'does not exist. ')

    TR = 4.6 # Echo Spacing is msec but not in header -> Set as TR in harmonize
    array, header = series.array(['SliceLocation', 'InversionTime'], pixels_first=True, first_volume=True)
    signal_pars = [{'xdata': np.array([hdr['InversionTime'] for hdr in header[z,:]]), 'TR':TR} for z in range(array.shape[2])]

    series.message('Setting up MDR..')
    signal_model = models.t1.MonoExp().fit
    
    return _mdr(series, array, header, signal_model, signal_pars, study)


def T2(folder):

    desc = "T2m_magnitude"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on T2: series ' + desc + 'does not exist. ')

    array, header = series.array(['SliceLocation', 'InversionTime'], pixels_first=True, first_volume=True)
    TE = series.values('InversionTime',  dims=('SliceLocation', 'InversionTime')).astype(np.float16)

    series.message('Setting up MDR..')
    signal_pars = [{'xdata':TE[z,:]} for z in range(TE.shape[0])]
    signal_model = models.t2.MonoExpOffset().fit
    
    return _mdr(series, array, header, signal_model, signal_pars, study)

def T1_T2(folder):

    desc_T1_T2 = "T1m_T2m_magnitude"
    series_T1_T2, study_T1_T2 = input_series(folder, desc_T1_T2, export_study)
    if series_T1_T2 is None:
        raise RuntimeError('Cannot perform MDR on T1: series ' + desc_T1 + 'does not exist. ')

    desc_T1 = "T1m_magnitude"
    series_T1, study_T1 = input_series(folder, desc_T1, export_study)
    if series_T1 is None:
        raise RuntimeError('Cannot perform MDR on T1: series ' + desc_T1 + 'does not exist. ')
    
    desc_T2 = "T2m_magnitude"
    series_T2, study_T2 = input_series(folder, desc_T2, export_study)
    if series_T2 is None:
        raise RuntimeError('Cannot perform MDR on T1: series ' + desc_T2 + 'does not exist. ')

    array_T1, header_T1 = series_T1.array(['SliceLocation', 'InversionTime'], pixels_first=True, first_volume=True)
    array_T2, header_T2 = series_T2.array(['SliceLocation', 'InversionTime'], pixels_first=True, first_volume=True)
    
    TE = series_T2.values('InversionTime',  dims=('SliceLocation', 'InversionTime')).astype(np.float16)
    TI = series_T1.values('InversionTime',  dims=('SliceLocation', 'InversionTime')).astype(np.float16)

    array   = np.concatenate((array_T1, array_T2), axis=3)
    header  = np.concatenate((header_T1,header_T2),axis=1)
    signal_pars_T1_T2 = np.concatenate((TI,TE),axis=1)
    signal_pars = [{'xdata':signal_pars_T1_T2[z,:]} for z in range(signal_pars_T1_T2.shape[0])]

    series_T1.message('Setting up MDR..')
    signal_model = models.t1_t2.MonoExp_T1_T2().fit
    
    return _mdr(series_T1_T2, array, header, signal_model, signal_pars, study_T1_T2)

def T2star(folder):

    desc = "T2starm_magnitude"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on T2*: series ' + desc + 'does not exist. ')

    array, header = series.array(['SliceLocation', 'EchoTime'], pixels_first=True, first_volume=True)
    signal_pars = [{'xdata':np.array([hdr.EchoTime for hdr in header[z,:]])} for z in range(array.shape[2])]
    
    series.message('Setting up MDR..')
    signal_model = models.T2star.MonoExp().fit

    return _mdr(series, array, header, signal_model, signal_pars, study)


def MT(folder): # TODO: Note this is a 3D sequence - do not coreg slice by slice - needs 3D registration

    desc = 'MT'
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on MT: series ' + desc + ' does not exist. ')

    array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True, first_volume=True)
    signal_pars = [{} for _ in range(array.shape[2])]

    series.message('Setting up MDR..')
    signal_model = constant.main

    return _mdr(series, array, header, signal_model, signal_pars, study)


def DTI(folder):

    desc = "DTI"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on DTI: series ' + desc + 'does not exist. ')

    dims = ['SliceLocation', 'InstanceNumber']
    bvals, bvecs = series.values('DiffusionBValue', 'DiffusionGradientOrientation', dims=dims)
    array, header = series.array(dims, pixels_first=True, first_volume=True)
    
    series.message('Setting up MDR..')
    signal_pars = [{'bvals':bvals[z,:], 'bvecs':np.stack(bvecs[z,:]), 'fit_method':'WLS'} for z in range(array.shape[2])]
    signal_model = models.DTI.DiPy().fit

    return _mdr(series, array, header, signal_model, signal_pars, study)


def IVIM(folder, series=None,study=None):

    desc = 'IVIM'
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on IVIM: series ' + desc + 'does not exist. ')

    dims = ['SliceLocation', 'InstanceNumber']
    bvals, bvecs = series.values('DiffusionBValue', 'DiffusionGradientOrientation', dims=dims)
    array, header = series.array(dims, pixels_first=True, first_volume=True)

    series.message('Setting up MDR..')
    signal_pars = [{'bvals':bvals[z,:], 'bvecs':np.stack(bvecs[z,:])} for z in range(array.shape[2])]
    signal_model = models.DWI.Lin().fit
    return _mdr(series, array, header, signal_model, signal_pars, study)


def DCE(folder):

    desc = "DCE"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on DCE: not all series are there.')
    
    time, aif = load_aif(folder)
    array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True, first_volume=True)
    
    series.message('Setting up MDR..')
    signal_pars = [{'aif':aif, 'time':time, 'baseline':15} for _ in range(array.shape[2])]
    signal_model = dcmri.pixel_2cfm_linfit

    return _mdr(series, array, header, signal_model, signal_pars, study)


def default_elastix_parameters():
    # See here for default bspline settings and explanation of parameters
    # https://github.com/SuperElastix/ElastixModelZoo/tree/master/models%2Fdefault
    param_obj = itk.ParameterObject.New()
    parameter_map_bspline = param_obj.GetDefaultParameterMap('bspline')
    param_obj.AddParameterMap(parameter_map_bspline) 
    param_obj.SetParameter("FixedImagePyramid", "FixedRecursiveImagePyramid") # "FixedSmoothingImagePyramid"
    param_obj.SetParameter("MovingImagePyramid", "MovingRecursiveImagePyramid") # "MovingSmoothingImagePyramid"
    param_obj.SetParameter("Metric", "AdvancedMeanSquares")
    param_obj.SetParameter("FinalGridSpacingInPhysicalUnits", "50.0")
    param_obj.SetParameter("ErodeMask", "false")
    param_obj.SetParameter("ErodeFixedMask", "false")
    #param_obj.SetParameter("NumberOfResolutions", "4") 
    #param_obj.SetParameter("MaximumNumberOfIterations", "500") # down from 500
    param_obj.SetParameter("MaximumStepLength", "0.1") 
    #param_obj.SetParameter("NumberOfSpatialSamples", "2048")
    #param_obj.SetParameter("BSplineInterpolationOrder", "1")
    #param_obj.SetParameter("FinalBSplineInterpolationOrder", "3")
    #param_obj.SetParameter("DefaultPixelValue", "0")
    param_obj.SetParameter("WriteResultImage", "false")
    return param_obj


def _mdr(series, array, header, signal_model, signal_pars, study, downsample=2, elastix_parameters=default_elastix_parameters()):

    # Define 3D output arrays
    model_fit = np.zeros(array.shape)
    coreg = np.zeros(array.shape)
    def_max = np.zeros(array.shape[:3])
    
    # Perform MDR slice by slice
    number_slices = array.shape[2]
    for z in range(number_slices):
        series.progress(z+1, number_slices, 'Performing model-driven registration')
        mdr = mdreg.MDReg()
        mdr.log = False # True for debugging
        mdr.parallel = False # Parallellization is too slow
        mdr.max_iterations = 5
        mdr.downsample = downsample
        mdr.signal_model = signal_model
        mdr.signal_parameters = signal_pars[z]
        mdr.elastix = elastix_parameters
        mdr.set_array(array[:,:,z,:])
        mdr.pixel_spacing = header[z,0].PixelSpacing
        mdr.fit() 
        model_fit[:,:,z,:] = mdr.model_fit
        coreg[:,:,z,:] = mdr.coreg
        def_max[:,:,z] = np.amax(np.linalg.norm(mdr.deformation, axis=2),axis=2)

    # Save results as DICOM
    desc = series.SeriesDescription
    fit = study.new_series(SeriesDescription=desc + '_mdr_fit')
    fit.set_array(model_fit, header, pixels_first=True)
    moco = study.new_series(SeriesDescription = desc + '_mdr_moco')
    moco.set_array(coreg, header, pixels_first=True)
    defo = study.new_series(SeriesDescription = desc + '_mdr_defo')
    defo.set_array(def_max, header[:,0], pixels_first=True)

    return fit, moco





