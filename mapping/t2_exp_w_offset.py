"""
@author: Joao Periquito
iBEAt study T2-mapping Philips TSE
2023
"""

import numpy as np
from scipy.optimize import curve_fit

def T2_exp_Eq(TE,S0,T2,C):
    """ T1 Calculation using MOLLI-FISP
    """
    S = np.sqrt((S0*np.exp(-TE/T2))**2 + (C**2))

    return S


def T2_fitting(images_to_be_fitted, TE):
    """ curve_fit function for T1-mapping.

    Args
    ----
    images_to_be_fitted (numpy.ndarray): pixel value for time-series (i.e. at each TI time) with shape [x,:]
    TI (list): list of inversion times
    sequenceParam (list): [FA, TR, N]

    Returns
    -------
    fit (list): signal model fit per pixel
    S0 (numpy.float64): fitted parameter 'S0' per pixel 
    T1 (numpy.float64): fitted parameter 'T1' (ms) per pixel.
    """   


    lb =            [0     , 0  ,0]
    ub =            [np.inf, 1500,np.inf]

    S0 = np.max(images_to_be_fitted)
    T2 = 80
    C = np.min(images_to_be_fitted)

    initial_guess = [S0, T2,C] 
    try:
        popt, pcov = curve_fit(T2_exp_Eq, TE, images_to_be_fitted, initial_guess, bounds=(lb,ub), method='trf',maxfev=500)
        S0    = popt[0]
        T2    = popt[1]
        C     = popt[2]

    except:
        S0    = 0
        T2    = 0
        C     = 0

    residuals = images_to_be_fitted - T2_exp_Eq(TE,S0,T2,C)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((images_to_be_fitted-np.mean(images_to_be_fitted))**2)
    r_squared = 1 - (ss_res / ss_tot)

    
    fit=T2_exp_Eq(TE,S0,T2,C)

    return fit, S0,T2,C,r_squared



def main(images_to_be_fitted, TE):
    """ main function that performs the T2 model-fit at single pixel level. 

    Args
    ----
    images_to_be_fitted (numpy.ndarray): pixel value for time-series (i.e. at each echo time) with shape [x,:]
    TE_list (list): list of echo times

    Returns
    -------
    fit (list): signal model fit per pixel
    fitted_parameters (list): list with signal model fitted parameters 'S0','fw' (fraction of water) and 'T2sw'.  
    """

    results = T2_fitting(images_to_be_fitted, TE)

    fit     = results[0]
    S0      = results[1]
    T2      = results[2]
    C      = results[3]
    r_square= results[4]

    return fit, S0, T2,C, r_square