import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dbdicom.extensions import vreg

from pipelines import measure
from models import (
    T1_look_locker_spoiled,
    DCE_2CM_nonlinear,
)


def T1(folder):

    results_path = folder.path() + '_QC'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    dyn_desc = 'T1map_kidneys_cor-oblique_mbh_magnitude_mdr_moco'
    dims = ('SliceLocation', 'InversionTime')
    model = T1_look_locker_spoiled

    vals_kidneys = []
    for kidney in ['LK','RK']:

        # Check if the required series are there and raise an error if not
        dyn_kidney = dyn_desc + '_' + kidney + '_align' 
        folder.message('Finding ' + dyn_kidney)
        dyn_kidney = folder.series(SeriesDescription=dyn_kidney)
        if dyn_kidney == []:
            raise ValueError('Cannot perform T1 ROI analysis: missing dynamic series aligned to kidney ' + kidney)
        folder.message('Finding ' + kidney)
        kidney_mask = folder.series(SeriesDescription=kidney)
        if kidney_mask == []:
            raise ValueError('Cannot perform T1 ROI analysis: missing mask for kidney' + kidney)
        map_kidney = dyn_desc + '_T1_map_' + kidney + '_align' 
        folder.message('Finding ' + map_kidney)
        map_kidney = folder.series(SeriesDescription=map_kidney)
        if map_kidney == []:
            raise ValueError('Cannot perform T1 ROI analysis: missing T1 aligned to kidney ' + kidney)
        
        # Load curve and T1 values for the kidney
        TI = dyn_kidney[0].values('InversionTime', dims=dims)
        TI = TI[0,:].astype(np.float32)
        dyn_kidney, vals_kidney = load_roi_curve(dyn_kidney[0], kidney_mask[0], map_kidney[0], dims=dims)
        vals_kidneys.append(vals_kidney)

        # Calculate the fit
        pars = model.fit_signal((TI, dyn_kidney, 1e-6, True))
        fit = model.signal(TI, *pars)

        # Export the results
        fig, ax = plt.subplots(1,1,figsize=(5,5))
        ax.plot(TI, dyn_kidney, 'ro', label='Signal for kidney ' + kidney, linewidth=3.0)
        ax.plot(TI, fit, 'b-', label='T1 model fit for kidney ' + kidney, linewidth=3.0)
        ax.plot(TI, 0*TI, color='gray')
        ax.set(xlabel='Inversion time (msec)', ylabel='Signal (a.u.)')
        ax.legend()

        plt.savefig(os.path.join(results_path, 'model fit T1 ('+ kidney +').png'), dpi=600)
        plt.close()

        # Update master table
        p = model.pars()
        measure.add_rows(folder, [
            [folder.PatientID, kidney, 'Kidney', 'ROI', pars[1], 'msec', kidney + '-' + p[1] + '-ROI', 'ROI fit'],
            [folder.PatientID, kidney, 'Kidney', 'ROI', pars[2], 'deg', kidney + '-' + p[2] + '-ROI', 'ROI fit'],
        ])

    # Build dynamic for BK
    kidney = 'BK'
    nt = len(vals_kidneys[0])
    dyn_kidney = np.zeros(nt)
    for t in range(nt):
        vals_t = list(vals_kidneys[0][t]) + list(vals_kidneys[1][t])
        dyn_kidney[t] = np.mean(vals_t)

    # Calculate the fit
    pars = model.fit_signal((TI, dyn_kidney, 1e-6, True))
    fit = model.signal(TI, *pars)

    # Export the results
    fig, ax = plt.subplots(1,1,figsize=(5,5))
    ax.plot(TI, dyn_kidney, 'ro', label='Signal for kidney ' + kidney, linewidth=3.0)
    ax.plot(TI, fit, 'b-', label='T1 model fit for kidney ' + kidney, linewidth=3.0)
    ax.plot(TI, 0*TI, color='gray')
    ax.set(xlabel='Inversion time (msec)', ylabel='Signal (a.u.)')
    ax.legend()

    plt.savefig(os.path.join(results_path, 'model fit T1 ('+ kidney +').png'), dpi=600)
    plt.close()

    # Update master table
    p = model.pars()
    measure.add_rows(folder, [
        [folder.PatientID, kidney, 'Kidney', 'ROI', pars[1], 'msec', kidney + '-' + p[1] + '-ROI', 'ROI fit'],
        [folder.PatientID, kidney, 'Kidney', 'ROI', pars[2], 'deg', kidney + '-' + p[2] + '-ROI', 'ROI fit'],
    ])



def dce(folder):

    dyn_desc = "DCE_kidneys_cor-oblique_fb_mdr_moco"
    time, aif = load_aif(folder)

    dyn_kidneys = []
    t1_kidneys = []
    for kidney in ['LK','RK']:

        # Check if the required series are there and raise an error if not
        dyn_kidney = dyn_desc + '_' + kidney + '_align' 
        dyn_kidney = folder.series(SeriesDescription=dyn_kidney)
        if dyn_kidney == []:
            raise ValueError('Cannot perform DCE ROI analysis: missing dynamic series aligned to kidney ' + kidney)
        kidney_mask = folder.series(SeriesDescription=kidney)
        if kidney_mask == []:
            raise ValueError('Cannot perform DCE ROI analysis: missing mask for kidney' + kidney)
        map_kidney = dyn_desc + '_AUC_map_' + kidney + '_align' 
        map_kidney = folder.series(SeriesDescription=map_kidney)
        if map_kidney == []:
            raise ValueError('Cannot perform DCE ROI analysis: missing AUC aligned to kidney ' + kidney)
        
        # Load curve and T1 values for the kidney
        dyn_kidney = load_roi_curve(dyn_kidney[0], kidney_mask[0], map_kidney[0])
        dyn_kidneys.append(dyn_kidney)
        t1_kidney = measure.read_master_table(folder, kidney + '-T1-Median')
        t1_kidneys.append(t1_kidney/1000)


    # Get sequence parameters
    #FA     = float(header[0,0]['FlipAngle'])
    #Weight = float(header[0,0]['PatientWeight'])
    #TR = 2.2/1000

    fit, pars = DCE_2CM_nonlinear.fit_roi(time, aif, dyn_kidneys, t1_kidneys)

    # Save pars to master table
    # Export fits on top of curves to QC folder
    # Do this first for T1, then come back



def load_roi_curve(dynamic, mask, map, dims=('SliceLocation', 'AcquisitionTime')):
    
    # Because each slice has a different orientation,
    # the mask is list of 2D images
    mask_array, _ = vreg.mask_array(mask, on=map, dim=dims[1])
    array = dynamic.pixel_values(dims)
    # old API
    # array, _ = dynamic.array(list(dims), pixels_first=True, first_volume=True)

    curve = np.zeros(array.shape[3])
    vals = []
    for t in range(array.shape[3]):
        vals_t = []
        for z in range(array.shape[2]):
            array_zt = array[:,:,z,t]
            mask_z = mask_array[z][:,:,0,0] > 0.5
            vals_t += list(array_zt[mask_z])
        curve[t] = np.mean(vals_t)
        vals.append(vals_t)
    return curve, vals


def load_aif(folder):

    dce = folder.series(SeriesDescription="DCE_aorta_axial_fb")
    aif = folder.series(SeriesDescription='DCE-AIF')
    if dce==[]:
        raise RuntimeError('No DCE data found: cannot export AIF')
    if aif==[]:
        raise RuntimeError('No AIF data found: cannot export AIF')
    array, header = dce[0].array(['AcquisitionTime'], pixels_first=True, first_volume=True)
    aif_mask, _ = aif[0].array(pixels_first=True, first_volume=True)
    dce[0].message('Loading AIF curve..')
    time = np.array([header[k]['AcquisitionTime'] for k in range(header.shape[0])])
    time -= time[0]
    loc = aif_mask != 0
    aif = np.array([np.mean(array[:,:,k][loc]) for k in range(array.shape[2])])

    return time, aif




    



    
