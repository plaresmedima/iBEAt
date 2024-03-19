import os
import numpy as np
import imageio

from dbdicom.pipelines import input_series

import itk
import mdreg
from mdreg.models import constant
from models import (
    T1_look_locker_spoiled,
    T2_mono_exp,
    T2star_mono_exp, 
    DCE_2CM,
    DTI_dipy,
    DWI_linear,
)
from pipelines.roi_fit import load_aif


export_study = 'MDR results'


def T1(folder):

    desc = "T1map_kidneys_cor-oblique_mbh_magnitude"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on T1: series ' + desc + 'does not exist. ')

    array, header = series.array(['SliceLocation', 'InversionTime'], pixels_first=True, first_volume=True)
    signal_pars = [{'xdata': np.array([hdr['InversionTime'] for hdr in header[z,:]])} for z in range(array.shape[2])]

    series.message('Loading elastix parameters..')
    signal_model = T1_look_locker_spoiled.fit
    elastix_parameters = default_elastix_parameters()
    downsample = 2
    
    return _mdr(series, array, header, signal_model, elastix_parameters, signal_pars, study, downsample)


def T2(folder):

    desc = "T2map_kidneys_cor-oblique_mbh_magnitude"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on T2: series ' + desc + 'does not exist. ')

    array, header = series.array(['SliceLocation', 'InversionTime'], pixels_first=True, first_volume=True)
    TE = series.values('InversionTime',  dims=('SliceLocation', 'InversionTime'))

    series.message('Loading elastix parameters..')
    signal_pars = [{'xdata':TE[z,:]} for z in range(TE.shape[1])]
    signal_model = T2_mono_exp.fit
    elastix_parameters = default_elastix_parameters()
    downsample = 2
    
    return _mdr(series, array, header, signal_model, elastix_parameters, signal_pars, study, downsample)


def T2star(folder):

    desc = "T2star_map_kidneys_cor-oblique_mbh_magnitude"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on T2*: series ' + desc + 'does not exist. ')

    array, header = series.array(['SliceLocation', 'EchoTime'], pixels_first=True, first_volume=True)
    signal_pars = [{'xdata':np.array([hdr.EchoTime for hdr in header[z,:]])} for z in range(array.shape[2])]
    
    series.message('Loading elastix parameters..')
    signal_model = T2star_mono_exp.fit
    elastix_parameters = default_elastix_parameters()
    downsample = 2

    return _mdr(series, array, header, signal_model, elastix_parameters, signal_pars, study, downsample)


def MT(folder): # Note this is a 3D sequence - do not coreg slice by slice - needs 3D registration

    desc = 'MT_kidneys_cor-oblique_bh'
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on MT: series ' + desc + ' does not exist. ')

    array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True, first_volume=True)
    signal_pars = [{} for _ in range(array.shape[2])]

    series.message('Loading elastix parameters..')
    signal_model = constant.main
    elastix_parameters = default_elastix_parameters()
    downsample = 2

    return _mdr(series, array, header, signal_model, elastix_parameters, signal_pars, study, downsample)


def DTI(folder):

    desc = "DTI_kidneys_cor-oblique_fb"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on DTI: series ' + desc + 'does not exist. ')

    dims = ['SliceLocation', 'InstanceNumber']
    bvals, bvecs = series.values('DiffusionBValue', 'DiffusionGradientOrientation', dims=dims)
    array, header = series.array(dims, pixels_first=True, first_volume=True)
    
    signal_pars = [{'bvals':bvals[z,:], 'bvecs':np.stack(bvecs[z,:]), 'fit_method':'WLS'} for z in range(array.shape[2])]
    signal_model = DTI_dipy.fit

    series.message('Loading elastix parameters..')
    elastix_parameters = default_elastix_parameters()
    downsample = 2

    return _mdr(series, array, header, signal_model, elastix_parameters, signal_pars, study, downsample)


def IVIM(folder, series=None,study=None):

    desc = 'IVIM_kidneys_cor-oblique_fb'
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on IVIM: series ' + desc + 'does not exist. ')

    array, header = series.array(['SliceLocation', 'InstanceNumber'], pixels_first=True, first_volume=True)
    bvals, bvecs = series.values('DiffusionBValue', 'DiffusionGradientOrientation', dims=('SliceLocation', 'InstanceNumber'))

    signal_pars = [{'bvals':bvals[z,:], 'bvecs':bvecs[z,:]} for z in range(array.shape[2])]
    signal_model = DWI_linear.fit
    elastix_parameters = default_elastix_parameters()
    downsample = 2

    return _mdr(series, array, header, signal_model, elastix_parameters, signal_pars, study, downsample)


def DCE(folder):

    desc = "DCE_kidneys_cor-oblique_fb"
    series, study = input_series(folder, desc, export_study)
    if series is None:
        raise RuntimeError('Cannot perform MDR on DCE: not all series are there.')
    
    time, aif = load_aif(folder)
    array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True, first_volume=True)
    
    signal_pars = [{'aif':aif, 'time':time, 'baseline':15} for _ in range(array.shape[2])]
    signal_model = DCE_2CM.fit
    elastix_parameters = default_elastix_parameters()
    downsample = 2

    return _mdr(series, array, header, signal_model, elastix_parameters, signal_pars, study, downsample)



def _mdr(series, array, header, signal_model, elastix_parameters, signal_pars, study, downsample):

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


