"""
@author: Joao Periquito
iBEAt study T1-mapping forward model-fit
Siemens 3T PRISMA - Leeds (T1-mapping sequence)
2021
"""

import numpy as np
from scipy.optimize import curve_fit

def T1_FLASH_MOLLI_Eq(TI,M_eq, M_eq_App,T1_App):
    """ T1 Calculation using MOLLI.

    TI : list of inversion times (between 100 and 7700ms)
    M_eq, M_eq_App, T1_App: tissue parameters
    Inv_Eff: Efficency of the 180 inversion pulse
    """
    S_t = M_eq_App -(M_eq+M_eq_App)*np.exp(-TI/T1_App)

    return np.abs(S_t)

def T1_corrected(M_eq, M_eq_App,T1_App):

    A = M_eq_App
    B = M_eq_App + M_eq
    if (A ==0):
        T1 = 0
    else:
        T1 = (B/A - 1)*T1_App

    return T1

def T1_fitting(images_to_be_fitted, TI):
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
    lb =            [0     , 0     , 0     ]
    ub =            [np.inf, np.inf, np.inf]

    MeqA = np.max(images_to_be_fitted)
    Meq = np.max(images_to_be_fitted)*1.9-np.max(images_to_be_fitted)
    T1app = 1500/(1.9-1)

    initial_guess = [MeqA, Meq, T1app] 
    try:
        popt, pcov = curve_fit(T1_FLASH_MOLLI_Eq, TI, images_to_be_fitted, initial_guess, bounds=(lb,ub), method='trf',maxfev=500)
        S0      = popt[0]
        S0_App  = popt[1]
        T1_app  = popt[2]
    except:
        S0      = 0
        S0_App  = 0
        T1_app  = 0

    fit = T1_FLASH_MOLLI_Eq(TI,S0, S0_App,T1_app)
    residuals = images_to_be_fitted - fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((images_to_be_fitted-np.mean(images_to_be_fitted))**2)
    r_squared = 1 - (ss_res / ss_tot)

    T1 = T1_corrected(S0, S0_App,T1_app)

    return fit, S0, S0_App,T1_app,T1,r_squared



def main(images_to_be_fitted, TI):

    results = T1_fitting(images_to_be_fitted, TI)

    fit     = results[0]
    S0      = results[1]
    S0_App  = results[2]
    T1_app  = results[3]
    T1      = results[4]
    r_square= results[5]

    return fit, S0, S0_App, T1_app, T1, r_square