"""
@author: Joao Periquito
iBEAt study T2* model-fit 
Siemens 3T PRISMA - Leeds (T2* sequence)
2021
"""

import numpy as np
from tqdm import tqdm
import mapping

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
    

    T2_map = np.zeros((np.size(T1_images_to_be_fitted,0), np.size(T1_images_to_be_fitted,1), np.size(T1_images_to_be_fitted,2)))
    M0map  = np.zeros((np.size(T1_images_to_be_fitted,0), np.size(T1_images_to_be_fitted,1), np.size(T1_images_to_be_fitted,2)))
    rsquaremap  = np.zeros((np.size(T1_images_to_be_fitted,0), np.size(T1_images_to_be_fitted,1), np.size(T1_images_to_be_fitted,2)))
    Cmap  = np.zeros((np.size(T1_images_to_be_fitted,0), np.size(T1_images_to_be_fitted,1), np.size(T1_images_to_be_fitted,2)))

    for i in tqdm (range(np.shape(T1_images_to_be_fitted)[2]),desc="Slice Completed..."):

        tempImpageSlice_T2 = np.squeeze(T1_images_to_be_fitted[:,:,i,:])
        TE = sequenceParam


        for xi in tqdm(range((np.size(tempImpageSlice_T2,0))),desc="Rows Completed..."):
            for yi in range((np.size(tempImpageSlice_T2,1))):
                
                Kidney_pixel_T2 = np.squeeze(np.array(tempImpageSlice_T2[xi,yi,:]))

                if Kidney_pixel_T2[0] == 0:
                    continue

                fit, S0, T2,C, r_square = mapping.t2_exp_w_offset.main(Kidney_pixel_T2, TE)

                residuals = Kidney_pixel_T2 - fit
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((Kidney_pixel_T2-np.mean(Kidney_pixel_T2))**2)
                r_squared = 1 - (ss_res / ss_tot)

                if (np.isnan(r_squared)): r_squared = 0

                M0map[xi, yi,i] = S0
                T2_map[xi,yi,i] = T2
                rsquaremap[xi,yi,i] = r_squared
                Cmap[xi,yi,i] = C


    fittedMaps = M0map, T2_map, C, rsquaremap
    #print('iupi')
    return fittedMaps
