"""
@author: Joao Periquito
iBEAt study T2* model-fit 
Siemens 3T PRISMA - Leeds (T2* sequence)
2021
"""

import numpy as np
from tqdm import tqdm

from scipy.optimize import curve_fit

def pars():
    return [
        "T2s_M0_Map",
        "T2s_fw_Map",
        "T2s_T2s_Map",
        "T2s_rsquare_Map",
    ]

def Mono_Exp_T2s_with_Water_Fat(x,M_eq, fw, T2sw):
    """ MonoExponential Fit of a T2* sequence takin in account in/out Fat.

    M_eq, T2sw, T2sf and fw (fraction of water): tissue parameters, T2sf = 9.3ms (Le Ster, C. et al. doi:10.1002/jmri.25205)
    x : List of echo times used between 3.7 and 44.3ms)
    """
    S_T2s = np.zeros(np.size(x))
    for m in range(np.size(x)):
        if m % 2 == 0:
            n=-1
            #print(n)
        else:
            n=1
            #print(n)

        S_T2s[m] = M_eq*(fw*np.exp(-x[m]/T2sw) + n*(1-fw)*np.exp(-x[m]/9.3))
    return S_T2s


def T2s_fitting(images_to_be_fitted, TE_list, ModelParam):
    """ curve_fit function for T2*-mapping.

    Args
    ----
    images_to_be_fitted (numpy.ndarray): pixel value for time-series (i.e. at each TE time) with shape [x,:]
    TE_list (list): list of echo times
    ModelParam: T2sf: T2* of Fat at 3.0T = 9.3 (Le Ster, C. et al. doi:10.1002/jmri.25205)

    Returns
    -------
    fit (list): signal model fit per pixel
    S0 (numpy.float64): fitted parameter 'S0' per pixel 
    T2s (numpy.float64): fitted parameter 'T2*' (ms) per pixel.
    fw (numpy.float64): fitted parameter 'water fraction per pixel [0-100%]
    """
    T2sf = ModelParam

    lb = [0,     0,     0]
    ub = [10000,     1,   200]
    initial_guess = [np.max(images_to_be_fitted),1,60] 


    try:
        fittedParameters, pcov = curve_fit(
            Mono_Exp_T2s_with_Water_Fat, TE_list, images_to_be_fitted, 
            initial_guess,
            bounds=(lb,ub),
            method='trf',
            maxfev=1000)    
    
        fit = []

        fit.append(Mono_Exp_T2s_with_Water_Fat(TE_list,fittedParameters[0],fittedParameters[1],fittedParameters[2]))

        S0 = fittedParameters[0]
        fw = fittedParameters[1]
        T2s = fittedParameters[2]

    except:
        fit = 0
        S0  = 0
        fw  = 0
        T2s = 0

    return fit, S0,fw, T2s


def main(T2s_images_to_be_fitted, sequenceParam):
    """ main function that performs the T2* model-fit with shared parameters at single pixel level. 

    Args
    ----
    T2s_images_to_be_fitted (numpy.ndarray): pixel value for time-series (i.e. at each T2 prep time) with shape [x,:]
    
    sequenceParam (list): [TE_list]


    Returns
    -------
    fitted_parameters: Map with signal model fitted parameters: 'S0', 'T1','T2','Flip Efficency','180 Efficency'.  
    """
    TE_list = np.array(sequenceParam)

    T2smap = np.zeros((np.size(T2s_images_to_be_fitted,0), np.size(T2s_images_to_be_fitted,1), np.size(T2s_images_to_be_fitted,2)))
    M0map  = np.zeros((np.size(T2s_images_to_be_fitted,0), np.size(T2s_images_to_be_fitted,1), np.size(T2s_images_to_be_fitted,2)))
    fwmap  = np.zeros((np.size(T2s_images_to_be_fitted,0), np.size(T2s_images_to_be_fitted,1), np.size(T2s_images_to_be_fitted,2)))
    rsquaremap  = np.zeros((np.size(T2s_images_to_be_fitted,0), np.size(T2s_images_to_be_fitted,1), np.size(T2s_images_to_be_fitted,2)))

    for i in tqdm (range(np.shape(T2s_images_to_be_fitted)[2]),desc="Slice Completed..."):

        tempImpageSlice_T2s = np.squeeze(T2s_images_to_be_fitted[:,:,i,:])

        for xi in tqdm(range((np.size(tempImpageSlice_T2s,0))),desc="Rows Completed..."):
            
            for yi in range((np.size(tempImpageSlice_T2s,1))):
                
                Kidney_pixel_T2s = np.squeeze(np.array(tempImpageSlice_T2s[xi,yi,:]))

                if Kidney_pixel_T2s[0] == 0:
                    continue

                T2s_fat = 9.3     #T2* of Fat at 3.0T = 9.3 (Le Ster, C. et al. doi:10.1002/jmri.25205)
                results = T2s_fitting(Kidney_pixel_T2s, TE_list, T2s_fat)
                Fit = results[0]
                Fitted_Parameters = [results[1],results[2],results[3]]

                residuals = Kidney_pixel_T2s - Fit
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((Kidney_pixel_T2s-np.mean(Kidney_pixel_T2s))**2)
                r_squared = 1 - (ss_res / ss_tot)

                if (np.isnan(r_squared)): r_squared = 0

                M0map[xi, yi,i] = Fitted_Parameters[0]
                fwmap[xi, yi,i] = Fitted_Parameters[1]
                T2smap[xi,yi,i] = Fitted_Parameters[2]
                rsquaremap[xi,yi,i] = r_squared

    return M0map, fwmap, T2smap, rsquaremap