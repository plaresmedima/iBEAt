import os
import time
import numpy as np

import mdreg
import mdreg.models.T2star_parallel
import mdreg.models.T2_parallel
import mdreg.models.T1_parallel
import mdreg.models.DWI_monoexponential_parallel
import mdreg.models.DTI
import mdreg.models.DCE_2CFM

import utilities.autoaif

elastix_pars = os.path.join(os.path.join(os.path.dirname(os.path.dirname(__file__)),'utilities'), 'elastix')

def MDRegT1(series, study):
    start_time = time.time()
    series.log("T1 motion correction has started")
    if series.Manufacturer == 'SIEMENS':
        sort_by = 'InversionTime'

    else:
        sort_by = (0x2005, 0x1572)

    array, header = series.array(['SliceLocation', sort_by], pixels_first=True)
    
    signal_pars =0
    signal_model = mdreg.models.T1_parallel
    elastix_file = 'BSplines_T1.txt'
    number_slices = array.shape[2]
    
    vals = _mdr(series, number_slices, array, header, signal_model, elastix_file, signal_pars, sort_by=sort_by, study=study)
    series.log("T1 motion correction was completed --- %s seconds ---" % (int(time.time() - start_time)))
    return vals

def MDRegT2star(series=None,study=None):
    start_time = time.time()
    series.log("T2* motion correction has started")

    array, header = series.array(['SliceLocation', 'EchoTime'], pixels_first=True)

    signal_pars = 0
    signal_model = mdreg.models.T2star_parallel
    elastix_file = 'BSplines_T2star.txt'
    number_slices = array.shape[2]

    vals = _mdr(series, number_slices, array, header, signal_model, elastix_file, signal_pars, sort_by='EchoTime',study=study)
    series.log("T2* motion correction was completed --- %s seconds ---" % (int(time.time() - start_time)))
    return vals


def MDRegT2(series=None, study=None):
    """Perform MDR on all slices using a T2 mono-exp model"""
    start_time = time.time()
    series.log("T2 motion correction has started")
    
    #series = zoom(series, 0.5)

    if series.Manufacturer == 'SIEMENS':
        array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
        signal_pars = [0,30,40,50,60,70,80,90,100,110,120]
    else:
        array, header = series.array(['SliceLocation', 'EchoTime'], pixels_first=True)
        signal_pars = series.EchoTime

    signal_model = mdreg.models.T2_parallel
    elastix_file = 'BSplines_T2.txt'
    number_slices = array.shape[2]
    
    vals = _mdr(series, number_slices, array, header, signal_model, elastix_file, signal_pars, sort_by='None', study=study)
    series.log("T2 motion correction was completed --- %s seconds ---" % (int(time.time() - start_time)))
    return vals


def MDRegIVIM(series=None,study=None):
    """Perform MDR on all slices using a DWI mono-exp model"""
    start_time = time.time()
    series.log("IVIM motion correction has started")

    #series = zoom(series, 0.5)
    array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)

    signal_pars = [0,10.000086, 19.99908294, 30.00085926, 50.00168544, 80.007135, 100.0008375, 199.9998135, 300.0027313, 600.0]
    signal_model = mdreg.models.DWI_monoexponential_parallel
    elastix_file = 'BSplines_IVIM.txt'

    number_slices = array.shape[2]
    vals = _mdr(series, number_slices, array, header, signal_model, elastix_file, signal_pars, sort_by='None',study=study)
    series.log("IVIM motion correction was completed --- %s seconds ---" % (int(time.time() - start_time)))
    return vals

def MDRegDTI(series=None,study=None):
    """Perform MDR on all slices using a DTI model"""
    start_time = time.time()
    series.log("DTI motion correction has started")

    if series.Manufacturer == 'SIEMENS':
        array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    else:
        array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
        array = np.transpose(array, (0, 1, 2, 4, 3))
        header = np.transpose(header, (0, 2, 1))

    signal_pars = 0
    signal_model = mdreg.models.DTI
    elastix_file = 'BSplines_DTI.txt'
    number_slices = array.shape[2]

    vals = _mdr(series, number_slices, array, header, signal_model, elastix_file, signal_pars, sort_by='DTI',study=study)
    series.log("DTI motion correction was completed --- %s seconds ---" % (int(time.time() - start_time)))
    return vals

def MDRegMT(series=None,study=None):
    """Perform MDR on all slices using a MT model"""
    start_time = time.time()
    series[0].log("MT motion correction has started")

    mt_off =series[0]
    mt_on =series[1]

    array_off, header_off = mt_off.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    array_on, header_on = mt_on.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)

    #array_off = np.reshape(array_off,np.shape(array_off)+(1,))
    #array_on = np.reshape(array_on,np.shape(array_on)+(1,))
    array = np.concatenate((array_off, array_on), axis=3)
    header = np.concatenate((header_off, header_on), axis=1)
    
    #array = np.reshape(array,np.shape(array)+(1,))

    signal_pars = []
    signal_model = mdreg.models.constant
    elastix_file = 'BSplines_MT.txt'

    number_slices = array.shape[2]
    vals = _mdr(mt_on, number_slices, array, header, signal_model, elastix_file, signal_pars, sort_by='None',study=study)
    series[0].log("MT motion correction was completed --- %s seconds ---" % (int(time.time() - start_time)))
    return vals

def MDRegDCE(series=None, study=None):
    """Perform MDR on all slices using a DCE linear model"""
    start_time = time.time()
    series.log("DCE motion correction has started")

    #series = zoom(series, 0.5)
    array, header = series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)

    signal_pars = 0
    signal_model = mdreg.models.DCE_2CFM
    elastix_file = 'BSplines_DCE.txt'

    number_slices = array.shape[2]
    vals = _mdr(series, number_slices, array, header, signal_model, elastix_file, signal_pars, sort_by='DCE', study=study)
    series.log("DCE motion correction was completed --- %s seconds ---" % (int(time.time() - start_time)))
    return vals


def _mdr(series, number_slices, array, header, signal_model, elastix_file, signal_pars,sort_by, study=None):
    """ MDR fit function.  

    Args:
    ----
    app             (wezel)
    series          (wezel series): contains the user selected series in wezel GUI
    number_slices   (numpy.1darray): number of total slices (third dimention of the parameter "array")
    array           (numpy.ndarray): DICOM data with shape [x-dim, y-dim, z-dim (slice), t-dim (time series)]
    header          list containing all DICOM headers  
    signal_model    python script with the model fit for the differents mapping techniques
    elastix_file    (.txt file): elaxtix text file: BSplines_*mapping technique*.txt 
    signal_pars     (numpy.ndarray): either contains an hardcoded parameter not visible in DICOM headers (e.g.: T2 - EchoTime or DWI - b-values), or is 0 triggers the extraction of the relevant mapping parameter from DICOM headers
    sort by         (string array): either specifies the needed variable in DICOM headers (e.g.: EchoTime, IversionTime) or triggers a specific pre-process step (e.g.: DTI, DCE) 


    Returns
    -------
    model_fit       (numpy.ndarray): signal model fit at all time-series with shape [x-dim, y-dim, z-dim (slice), t-dim (time series), parameter (e.g.: S0, T2)].
    par             (numpy.ndarray): signal model fit at all time-series with shape [x-dim, y-dim, z-dim (slice), t-dim (time series), parameter (e.g.: S0, T2)].
    moco            (numpy.ndarray): signal model fit at all time-series with shape [x-dim, y-dim, z-dim (slice), t-dim (time series), parameter (e.g.: S0, T2)].
    """

    # PERFORM MDR
    parameters = signal_model.pars()

    # PARAMETER VARIABLES INITIALIZATION
    model_fit = np.empty(array.shape)
    coreg = np.empty(array.shape)
    pars = np.empty(array.shape[:3] + (len(parameters),) )

    # LOOP THROUGH SLICES
    for slice in range(number_slices):
        #print("Part -8  before mdreg: " + str(psutil.virtual_memory()[3]/1000000000))
        mdr = mdreg.MDReg()
        #print("Part -7  after mdreg: " + str(psutil.virtual_memory()[3]/1000000000))

        if signal_pars!=0:                                                                          #if condition for hardcoded parameters e.g.: T2 "EchoTime"
            mdr.signal_parameters = signal_pars

        else:
            if sort_by == "DTI":                                                                    #extracting DTI relevant parameters from DICOM headers                                              
                        
                        if series.Manufacturer == 'SIEMENS':
                            b_values = [float(hdr[(0x19, 0x100c)]) for hdr in header[slice,:,0]]
                            b_vectors = [hdr[(0x19, 0x100e)] for hdr in header[slice,:,0]]
                            orientation = [hdr.ImageOrientationPatient for hdr in header[slice,:,0]]
                        else:
                            b_values = [float(hdr[(0x18, 0x9087)]) for hdr in header[slice,:,0]]
                            b_vectors = [hdr[(0x18, 0x9089)] for hdr in header[slice,:,0]]
                            orientation = [hdr.ImageOrientationPatient for hdr in header[slice,:,0]]
                        mdr.signal_parameters = [b_values, b_vectors, orientation]

            elif sort_by =="DCE" and signal_pars==0:                                                
                        
                        # GET AIF
                        cutRatio=0.25             #create a window around the center of the image where the aorta is
                        filter_kernel=(15,15)     #gaussian kernel for smoothing the image to destroy any noisy single high intensity filter
                        regGrow_threshold = 2     #threshold for the region growing algorithm

                        if series.Manufacturer == 'SIEMENS':
                            for i in range(header.shape[0]):
                                if (header[i,0,0]["ImageOrientationPatient"]== [1, 0, 0, 0, 1, 0]):
                                    aortaslice = int(i + 1)
                                    break
                                else:
                                    aortaslice = int(9)
                        else:
                            aortaslice = int(1)

                        aif = utilities.autoaif.DCEautoAIF(array, header, series, aortaslice, cutRatio, filter_kernel, regGrow_threshold)

                        time = np.zeros(header.shape[1])

                        for i_2 in range(header.shape[1]):
                            tempTime = str(header[slice,i_2,0]['AcquisitionTime'])
                            beforepoint = tempTime.split(".")[0]
                            afterpoint = tempTime.split(".")[1]
                            tempH = int(beforepoint[0:2])
                            tempM = int(beforepoint[2:4])
                            tempS = int(beforepoint[4:])
                            tempRest = float("0." + afterpoint)
                            time[i_2] = tempH*3600+tempM*60+tempS+tempRest
                        time -=time[0]
               
                        baseline = 15
                        hematocrit = 0.45
                        signal_pars = [aif, time, baseline, hematocrit]
                        mdr.signal_parameters = signal_pars
                        print("DCE parameters were calculated ")

            else:
                mdr.signal_parameters = [hdr[sort_by] for hdr in (header[slice,:,0])]                #extracting relevant parameters from DICOM headers using the string sort_by
                print(mdr.signal_parameters)
                print("sort by")
        #STORING RESULTS
        
        #print("Part -6  before set array: " + str(psutil.virtual_memory()[3]/1000000000))
        mdr.set_array(array[:,:,slice,:,0])    
        mdr.pixel_spacing = header[slice,0,0].PixelSpacing
        mdr.signal_model = signal_model
        mdr.read_elastix(os.path.join(elastix_pars, elastix_file))
        #mdr.pinned_message = 'MDR for slice ' + str(slice+1)
        
        #print("Part -5  before mdr.fit: " + str(psutil.virtual_memory()[3]/1000000000))
        mdr.fit()

        model_fit[:,:,slice,:,0] = mdr.model_fit
        coreg[:,:,slice,:,0] = mdr.coreg
        pars[:,:,slice,:] = mdr.pars

   # array = [x,y,z,TE]
   # pars = [x,y,z,2]
   # model_fit = [x,y,z,TE]
   # coreg = [x,y,z,TE]
    if study is None:
        study = series.parent()
    #EXPORT RESULTS TO wezel GUI USING DICOM
    for p in range(len(parameters)):
        par = series.SeriesDescription + '_mdr_par_' + parameters[p]
        par = study.new_series(SeriesDescription=par)
        #par = series.new_sibling(SeriesDescription=par)
        par.set_array(np.squeeze(pars[...,p]), np.squeeze(header[:,0]), pixels_first=True)
    fit = series.SeriesDescription + '_mdr_fit'
    fit = study.new_series(SeriesDescription=fit)
    fit.set_array(model_fit, np.squeeze(header[:,:]), pixels_first=True)
    moco = series.SeriesDescription + '_mdr_moco'
    moco = study.new_series(SeriesDescription = moco)
    moco.set_array(coreg, np.squeeze(header[:,:]), pixels_first=True)

    return fit, moco

def main(folder):
    start_time = time.time()
    folder.log("MDR has started!")

    current_study = folder.series()[0].parent()
    study = folder.series()[0].new_pibling(StudyDescription=current_study.StudyDescription + '_MDRresults')

    for series in folder.series():

        SeqName = series["SeriesDescription"]

        if SeqName is not None:
            print(SeqName)

            if SeqName == "T2star_map_kidneys_cor-oblique_mbh_magnitude":
                try:
                    print('starting t2s')
                    MDRegT2star(series, study=study)
                except Exception as e: 
                    folder.log("T2* motion correction was NOT completed; error: "+str(e))

            elif SeqName == "T1map_kidneys_cor-oblique_mbh_magnitude":
                try:
                    print("starting T1 mapping")
                    MDRegT1(series, study)
                except Exception as e: 
                    folder.log("T1 motion correction was NOT completed; error: "+str(e))

            elif SeqName == "T2map_kidneys_cor-oblique_mbh_magnitude":
                try:
                    print("starting T2 mapping")
                    MDRegT2(series, study=study)
                except Exception as e: 
                    folder.log("T2 motion correction was NOT completed; error: "+str(e))   

            elif SeqName == "DTI_kidneys_cor-oblique_fb":
                try:
                    print('starting dti')
                    MDRegDTI(series, study=study) 
                except Exception as e: 
                    folder.log("DTI motion correction was NOT completed; error: "+str(e))
            
            elif SeqName == "MT_OFF_kidneys_cor-oblique_bh":
                try:
                    MT_OFF = series
                    for series in folder.series():
                        if series['SeriesDescription'] == "MT_ON_kidneys_cor-oblique_bh":
                            MT_ON = series
                            break
                    MDRegMT([MT_OFF, MT_ON], study=study) 
                except Exception as e: 
                    folder.log("MT motion correction was NOT completed; error: "+str(e))

            elif SeqName == "DCE_kidneys_cor-oblique_fb":
                try:
                    print('starting DCE')
                    MDRegDCE(series, study=study)   
                except Exception as e: 
                    folder.log("DCE motion correction was NOT completed; error: "+str(e))

    folder.save()
    folder.log("MDR was completed --- %s seconds ---" % (int(time.time() - start_time)))

def mdr_slice_by_slice(series, study, sort_by, signal_model, elastix_file):

    array, header = series.array(['SliceLocation',sort_by], pixels_first=True)
    
    # PARAMETER VARIABLES INITIALIZATION
    parameters = signal_model.pars()
    model_fit = np.empty(array.shape)
    coreg = np.empty(array.shape)
    pars = np.empty(array.shape[:3] + (len(parameters),) )

    # LOOP THROUGH SLICES
    s=0
    number_slices = array.shape[2]
    for slice in range(number_slices):
        s+=1
        series.progress(s, number_slices, 'Performing model-driven registration')
        mdr = mdreg.MDReg()
        mdr.signal_parameters = [hdr[sort_by] for hdr in (header[slice,:,0])]
        mdr.set_array(array[:,:,slice,:,0])    
        mdr.pixel_spacing = header[slice,0,0].PixelSpacing
        mdr.signal_model = signal_model
        mdr.read_elastix(os.path.join(elastix_pars, elastix_file))
        mdr.fit()
        model_fit[:,:,slice,:,0] = mdr.model_fit
        coreg[:,:,slice,:,0] = mdr.coreg
        pars[:,:,slice,:] = mdr.pars

    for p in range(len(parameters)):
        par = study.new_series(SeriesDescription='mdr_par_' + parameters[p])
        par.set_array(np.squeeze(pars[...,p]), np.squeeze(header[:,0]), pixels_first=True)
    fit = study.new_series(SeriesDescription='mdr_fit')
    fit.set_array(model_fit, np.squeeze(header[:,:]), pixels_first=True)
    moco = study.new_series(SeriesDescription = 'mdr_moco')
    moco.set_array(coreg, np.squeeze(header[:,:]), pixels_first=True)

    return fit, moco
