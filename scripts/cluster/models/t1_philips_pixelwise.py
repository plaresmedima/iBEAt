"""
@author: Joao Periquito
iBEAt study T2* model-fit 
Siemens 3T PRISMA - Leeds (T2* sequence)
2021
"""

import numpy as np
from tqdm import tqdm

from .model_library.single_pixel_forward_models import t1_exp_fm


#from iBEAt_Model_Library.single_pixel_forward_models import iBEAT_T2s_FM

def main(T1_images_to_be_fitted, sequenceParam):
    """ main function that performs the T2* model-fit with shared parameters at single pixel level. 

    Args
    ----
    T1_images_to_be_fitted (numpy.ndarray): pixel value for time-series (i.e. at each T2 prep time) with shape [x,:]
    
    sequenceParam (list): [TE_list]


    Returns
    -------
    fitted_parameters: Map with signal model fitted parameters: 'S0', 'T1','T2','Flip Efficency','180 Efficency'.  
    """
    

    T1_map = np.zeros((np.size(T1_images_to_be_fitted,0), np.size(T1_images_to_be_fitted,1), np.size(T1_images_to_be_fitted,2)))
    T1_app_map = np.zeros((np.size(T1_images_to_be_fitted,0), np.size(T1_images_to_be_fitted,1), np.size(T1_images_to_be_fitted,2)))
    M0map  = np.zeros((np.size(T1_images_to_be_fitted,0), np.size(T1_images_to_be_fitted,1), np.size(T1_images_to_be_fitted,2)))
    rsquaremap  = np.zeros((np.size(T1_images_to_be_fitted,0), np.size(T1_images_to_be_fitted,1), np.size(T1_images_to_be_fitted,2)))

    for i in tqdm (range(np.shape(T1_images_to_be_fitted)[2]),desc="Slice Completed..."):

        tempImpageSlice_T1 = np.squeeze(T1_images_to_be_fitted[:,:,i,:])
        TI = np.array([hdr[(0x2005, 0x1572)] for hdr in sequenceParam[i,:,0]])

        for xi in tqdm(range((np.size(tempImpageSlice_T1,0))),desc="Rows Completed..."):
            for yi in range((np.size(tempImpageSlice_T1,1))):
                
                Kidney_pixel_T1 = np.squeeze(np.array(tempImpageSlice_T1[xi,yi,:]))

                if Kidney_pixel_T1[0] == 0:
                    continue

                fit, S0, S0_App, T1_app, T1, r_square = t1_exp_fm.main(Kidney_pixel_T1, TI)

                residuals = Kidney_pixel_T1 - fit
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((Kidney_pixel_T1-np.mean(Kidney_pixel_T1))**2)
                r_squared = 1 - (ss_res / ss_tot)

                if (np.isnan(r_squared)): r_squared = 0

                M0map[xi, yi,i] = S0
                T1_app_map[xi, yi,i] = T1_app
                T1_map[xi,yi,i] = T1
                rsquaremap[xi,yi,i] = r_squared


    fittedMaps = M0map, T1_app_map, T1_map, rsquaremap
    #print('iupi')
    return fittedMaps
