import time
import numpy as np
from dbdicom.wrappers import scipy
from dbdicom.wrappers import vreg, scipy
import scripts.QC_alignment as export_alignemnt
from pipelines.kidney_dce import fit as fit_DCE
import pandas as pd
import os

def main(folder,master_table):

    DCE_moco_imgs_description = "DCE_kidneys_cor-oblique_fb_mdr_moco"
    DCE_aorta_mask_description = "DCE_kidneys_cor-oblique_fb_DCE_ART"
    T1_map_moco_description = "T1map_kidneys_cor-oblique_mbh_magnitude_mdr_par_T1"
    RK_description = 'RK'
    LK_description = 'LK'

    folder.log("DCE analysis has started!")

    DCE_moco_series = folder.series(SeriesDescription=DCE_moco_imgs_description)
    T1_moco_series = folder.series(SeriesDescription=T1_map_moco_description)
    Aorta_mask_series = folder.series(SeriesDescription=DCE_aorta_mask_description)
    LK_mask_series = folder.series(SeriesDescription=LK_description)
    RK_mask_series = folder.series(SeriesDescription=RK_description)

    DCE_kidneys = DCE_moco_series[0].split_by('ImageOrientationPatient')

    time_aorta, signal_aorta = get_signal_aorta(DCE_moco_series,Aorta_mask_series)
    #plt.plot(time_aorta,signal_aorta,'.')
    #plt.show()
    # print(np.mean(signal_aorta))

    signal_lk = get_signal_kidney(DCE_kidneys,LK_mask_series,folder)
    #plt.close()
    #plt.plot(time_aorta,signal_lk,'.')
    #plt.show()
    # print(np.mean(signal_lk))

    signal_rk = get_signal_kidney(DCE_kidneys,RK_mask_series,folder)
    #plt.close()
    #plt.plot(time_aorta,signal_rk,'.')
    #plt.show()
    # print(np.mean(signal_rk))



    T1_lk = get_T1_kidney(T1_moco_series,LK_mask_series,folder)
    T1_rk = get_T1_kidney(T1_moco_series,RK_mask_series,folder)

    T1_lk = T1_lk/1000
    T1_rk = T1_rk/1000

    R1_lk = 1/T1_lk
    R1_rk = 1/T1_rk

    _, header = DCE_moco_series[0].array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    header = np.squeeze(header)
    #TR = float(header[0,0]['Repetition Time'])
    #TR     = float(header[0,0]['RepetitionTime'])/1000
    FA     = float(header[0,0]['FlipAngle'])
    Weight = float(header[0,0]['PatientWeight'])

    TR = 2.2/1000
    # FA = 10
    # Weight = 80

    # # Internal time resolution & acquisition time
    # dt = 1.0                # sec
    # tmax = 40*60.0          # Total acquisition time (sec)

    # # Default values for experimental parameters
    # tacq = 1.61             # Time to acquire a single datapoint (sec)
    # field_strength = 3.0    # Field strength (T)
    # weight = 70.0           # Patient weight in kg
    # conc = 0.5             # mmol/mL (https://www.bayer.com/sites/default/files/2020-11/primovist-pm-en.pdf) #DOTOREM = 0.5mmol/ml, Gadovist = 1.0mmol/ml 
    # dose = 0.05            # mL per kg bodyweight (quarter dose)
    # rate = 2                # Injection rate (mL/sec)
    # TR = 2.2/1000.0        # Repetition time (sec)
    # FA = 10.0               # Nominal flip angle (degrees)
    # TI = 85/1000
    # TSAT = 25.5/1000

    # # Physiological parameters
    # Hct = 0.45
    
    DCE_data_csv = pd.DataFrame({"time": np.array(time_aorta), "aorta": np.array(signal_aorta), "lk": np.array(signal_lk), "rk": np.array(signal_rk)})
    DCE_data_csv.to_csv(os.path.join(folder.path(), 'DCE_Data' +'.csv'),index=False)

    lines = [

    "RK_T1 (s): " + str(T1_rk),
    "LK_T1 (s): " + str(T1_lk),
    "Weight: " + str(header[0,0]['PatientWeight']),

    ]
    txt_path = os.path.join(folder.path(), 'DCE_Data' +'.txt')

    with open(txt_path, "w") as file:
        for line in lines:
            file.write(line + "\n")

    fit_DCE(
        2.2/1000,  #sec
        FA, #degrees
        3.0, #T
        0.05, #mL/kg  
        0.5, # mmol/ml 
        2, # ml/sec
        Weight, #kg
        0.45, 
        time_aorta, #sec
        signal_aorta,  #au
        1, #1/sec
        signal_rk, #au
        R1_rk, #1/sec
        150, #mL
        path = folder.path(), # full path for exporting diagnostics
        #name = header[0,0]['PatientID'],
        name = 'subject',
        ROI = 'rk',
        export_aorta = True,
        )
    
    fit_DCE(
        2.2/1000,  #sec
        FA, #degrees
        3.0, #T
        0.05, #mL/kg  
        0.5, # mmol/ml 
        2, # ml/sec
        Weight, #kg
        0.45, 
        time_aorta, #sec
        signal_aorta,  #au
        1, #1/sec
        signal_lk, #au
        R1_lk, #1/sec
        150, #mL
        path = folder.path(), # full path for exporting diagnostics
        #name = header[0,0]['PatientID'],
        name = 'subject',
        ROI = 'lk',
        export_aorta = True,
        )
    
    return master_table

def get_T1_kidney(T1_moco_series,kidney_mask_series,folder):
    
    T1_moco_series = T1_moco_series[0]
    
    overlay_mask  = scipy.map_to(kidney_mask_series[0], T1_moco_series)

    export_alignemnt.main(T1_moco_series, overlay_mask,'DCE_T1_alignment_TEST',folder.path())

    array_T1_moco, _ = T1_moco_series.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    array_T1_moco    = np.squeeze(array_T1_moco)

    kidney_Imgs = array_T1_moco
    print(kidney_Imgs.shape)

    array_kidney_mask, _ = overlay_mask.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    array_kidney_mask = np.squeeze(array_kidney_mask)

    T1_mask = np.squeeze(kidney_Imgs) * np.squeeze(array_kidney_mask)
    T1 =np.median(T1_mask[T1_mask>=1])
    return T1

def get_signal_aorta(DCE_moco_series,Aorta_mask_series):

    aortaslice = int(9)

    array_DCE_moco, header_DCE_moco = DCE_moco_series[0].array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    array_DCE_moco    = np.squeeze(array_DCE_moco)

    aorta_Imgs = array_DCE_moco[:,:,aortaslice-1,...]

    array_aorta_mask, _ = Aorta_mask_series[0].array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    array_aorta_mask = np.squeeze(array_aorta_mask)

    aif =[]
    for z in range(aorta_Imgs.shape[2]):
        tmask = np.squeeze(aorta_Imgs[:,:,z]) * np.squeeze(array_aorta_mask)
        aif.append(np.mean(tmask[tmask!=0]))

    time = np.zeros(header_DCE_moco.shape[1])

    for i_2 in range(header_DCE_moco.shape[1]):
        tempTime = str(header_DCE_moco[aortaslice-1,i_2,0]['AcquisitionTime'])
        time[i_2] = float(tempTime)

        # beforepoint = tempTime.split(".")[0]
        # afterpoint = tempTime.split(".")[1]
        # tempH = int(beforepoint[0:2])
        # tempM = int(beforepoint[2:4])
        # tempS = int(beforepoint[4:])
        # tempRest = float("0." + afterpoint)
        # time[i_2] = tempH*3600+tempM*60+tempS+tempRest
    time -=time[0]

    time_aorta = np.array(time)
    signal_aorta= np.array(aif)

    return time_aorta, signal_aorta

def get_signal_kidney(DCE_moco_series,kidney_mask_series,folder):
    
    DCE_kidneys = DCE_moco_series[0]
    
    overlay_mask  = scipy.map_to(kidney_mask_series[0], DCE_kidneys)

    export_alignemnt.main(DCE_kidneys, overlay_mask,'DCE_alignment_TEST',folder.path())

    array_DCE_moco, _ = DCE_kidneys.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    array_DCE_moco    = np.squeeze(array_DCE_moco)

    kidney_Imgs = array_DCE_moco

    array_kidney_mask, _ = overlay_mask.array(['SliceLocation', 'AcquisitionTime'], pixels_first=True)
    array_kidney_mask = np.squeeze(array_kidney_mask)

    aif =[]
    print(array_kidney_mask.shape)
    print(kidney_Imgs.shape)

    for z in range(kidney_Imgs.shape[3]):
        tmask = np.squeeze(kidney_Imgs[:,:,:,z]) * np.squeeze(array_kidney_mask)
        aif.append(np.mean(tmask[tmask!=0]))

    signal_kidney= np.array(aif)

    return signal_kidney



    
