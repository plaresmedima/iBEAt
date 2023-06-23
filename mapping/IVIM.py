"""
@author: Joao Periquito
iBEAt study IVIM model-fit
Siemens 3T PRISMA - Leeds (DW-EPI sequence)
2021
"""

import numpy as np
from scipy.optimize import curve_fit
from tqdm import tqdm



def Mono_Exp_IVIM(x,S0, D):
    """ MonoExponential Fit for IVIM data.

    S0, D, Ds and f (perfsion fraction): tissue parameters
    x : List of b-values used between 0 and 600 s/mm2)
    """
    S_IVIM = S0*(np.exp(-x*D))

    return S_IVIM

def Bi_Exp_IVIM(x,S0, D,Ds,f):
    """ BiExponential Fit for IVIM data.

    S0, D, Ds and f (perfsion fraction): tissue parameters
    x : List of b-values used between 0 and 600 s/mm2)
    """
    S_IVIM = S0*(f*np.exp(-x*Ds) + (1-f)*np.exp(-x*D))

    return S_IVIM


def Bi_Exp_IVIM_fitting(images_to_be_fitted, Bval_list):
    """ curve_fit function for IVIM-mapping.

    Args
    ----
    images_to_be_fitted (numpy.ndarray): pixel value for time-series (i.e. at each TE time) with shape [x,:]
    Bval_list (list): list of b-values


    Returns
    -------
    fit (list): signal model fit per pixel
    S0 (numpy.float64): fitted parameter 'S0' per pixel 
    D (numpy.float64): fitted parameter 'tissue diffusion' (mm2/s) per pixel.
    Ds (numpy.float64): fitted parameter 'pseudo diffusion' (mm2/s) per pixel.
    f (numpy.float64): fitted parameter 'pseudo diffusion fraction per pixel [0-100%]
    """

    lb_mono = [0,      0]
    ub_mono = [1,      1]
    initial_guess_mono = [1,0.002] 

    lb_bi = [0,      0]
    ub_bi = [1,      1]
    initial_guess_bi = [1,0.02] 

    try:
        
        fittedParameters_mono, pcov = curve_fit(Mono_Exp_IVIM, Bval_list[len(images_to_be_fitted)-3:len(images_to_be_fitted)], images_to_be_fitted[len(images_to_be_fitted)-3:len(images_to_be_fitted)], initial_guess_mono,bounds=(lb_mono,ub_mono),method='trf',maxfev=5000)

        fittedParameters, pcov = curve_fit(lambda x, S0,Ds: Bi_Exp_IVIM(x,S0, fittedParameters_mono[1],Ds,1-fittedParameters_mono[0]), Bval_list, images_to_be_fitted,initial_guess_bi,bounds=(lb_bi,ub_bi),method='trf',maxfev=5000)    
    
        fit = []

        fit.append(Bi_Exp_IVIM(Bval_list,fittedParameters[0],fittedParameters_mono[1],fittedParameters[1],1-fittedParameters_mono[0]))

        S0 = fittedParameters[0]
        D = fittedParameters_mono[1]
        Ds = fittedParameters[1]
        f = 1-fittedParameters_mono[0]



    except:
        fit = np.zeros(len(images_to_be_fitted))
        S0  = 0
        D  = 0
        Ds = 0
        f = 0

    return fit, S0,D, Ds,f


def Mono_Exp_IVIM_fitting(images_to_be_fitted, Bval_list):
    """ curve_fit function for mono-exp IVIM-mapping (ADC only)

    Args
    ----
    images_to_be_fitted (numpy.ndarray): pixel value for time-series (i.e. at each TE time) with shape [x,:]
    Bval_list (list): list of b-values


    Returns
    -------
    fit (list): signal model fit per pixel
    S0 (numpy.float64): fitted parameter 'S0' per pixel 
    D (numpy.float64): fitted parameter 'tissue diffusion' (mm2/s) per pixel.
    Ds (numpy.float64): fitted parameter 'pseudo diffusion' (mm2/s) per pixel.
    f (numpy.float64): fitted parameter 'pseudo diffusion fraction per pixel [0-100%]
    """

    lb_mono = [0,      0]
    ub_mono = [1,      1]
    initial_guess_mono = [1,0.002] 

    try:
        
        fittedParameters_mono, pcov = curve_fit(Mono_Exp_IVIM, Bval_list, images_to_be_fitted, initial_guess_mono,bounds=(lb_mono,ub_mono),method='trf',maxfev=5000)
    
        fit = []

        fit.append(Mono_Exp_IVIM(Bval_list,fittedParameters_mono[0],fittedParameters_mono[1]))

        S0 = fittedParameters_mono[0]
        D = fittedParameters_mono[1]


    except:
        fit = np.zeros(len(images_to_be_fitted))
        S0  = 0
        D  = 0

    return fit, S0, D



def main(IVIM_images_to_be_fitted, sequenceParam):
    
    bvals_list = np.array(sequenceParam)
    #print(bvals_list)

    bvals_unique = bvals_list[0:10]
    inds = bvals_list.argsort()
    bvals_list = bvals_list[inds]
    #print(bvals_list)

    S0map = np.zeros((np.size(IVIM_images_to_be_fitted,0), np.size(IVIM_images_to_be_fitted,1), np.size(IVIM_images_to_be_fitted,2)))
    Dmap = np.zeros((np.size(IVIM_images_to_be_fitted,0), np.size(IVIM_images_to_be_fitted,1), np.size(IVIM_images_to_be_fitted,2)))
    Dsmap  = np.zeros((np.size(IVIM_images_to_be_fitted,0), np.size(IVIM_images_to_be_fitted,1), np.size(IVIM_images_to_be_fitted,2)))
    fmap  = np.zeros((np.size(IVIM_images_to_be_fitted,0), np.size(IVIM_images_to_be_fitted,1), np.size(IVIM_images_to_be_fitted,2)))
    rsquaremap  = np.zeros((np.size(IVIM_images_to_be_fitted,0), np.size(IVIM_images_to_be_fitted,1), np.size(IVIM_images_to_be_fitted,2)))
    IVIM_images_avg = np.zeros((np.size(IVIM_images_to_be_fitted,0), np.size(IVIM_images_to_be_fitted,1), len(bvals_unique)))

    #IVIM_images_to_be_fitted = IVIM_images_to_be_fitted[:,:,2:4,:]
    #IVIM_images_to_be_fitted = IVIM_images_to_be_fitted[172:292,311:397,:,:]
    

    for i in tqdm (range(np.shape(IVIM_images_to_be_fitted)[2]),desc="Slice Completed..."):

        tempImageSlice_IVIM = np.squeeze(IVIM_images_to_be_fitted[:,:,i,:])
        tempImageSlice_IVIM_sorted = np.squeeze(tempImageSlice_IVIM[:,:,[inds]])

        bval_counter=0
        for k in range(len(bvals_unique)):
            b_temp = bvals_unique[k]
            b_vals_temp = bvals_list[bvals_list == b_temp]
            
            #print(bval_counter)
            #print(bval_counter+len(b_vals_temp))
            #print(bvals_list[bval_counter:bval_counter+len(b_vals_temp)])

            IVIM_images_avg[:,:,k] = np.average(tempImageSlice_IVIM_sorted[:,:,bval_counter:bval_counter+len(b_vals_temp)],axis =2)
            
            bval_counter = bval_counter + len(b_vals_temp)
        
        for xi in tqdm(range((np.size(IVIM_images_avg,0))),desc="Rows Completed..."):
            
            for yi in range((np.size(IVIM_images_avg,1))):
                
                Kidney_pixel_IVIM = np.squeeze(np.array(IVIM_images_avg[xi,yi,:]))

                if (Kidney_pixel_IVIM[0]==0):
                    continue
                 
                #print(Kidney_pixel_IVIM)
                Kidney_pixel_IVIM = Kidney_pixel_IVIM/Kidney_pixel_IVIM[0]

                results = Mono_Exp_IVIM_fitting(Kidney_pixel_IVIM, bvals_unique)
                Fit = results[0]
                Fitted_Parameters = [results[1], results[2]]

                residuals =  Kidney_pixel_IVIM - Fit
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((Kidney_pixel_IVIM-np.mean(Kidney_pixel_IVIM))**2)
                r_squared = 1 - (ss_res / ss_tot)
                if (np.isnan(r_squared)): r_squared = 0

                #print(r_squared)
                #print(Fitted_Parameters[0])
                #print(Fitted_Parameters[1])
                #print(Fitted_Parameters[2])
                #print(Fitted_Parameters[3])

                #plt.plot(bvals_unique,Kidney_pixel_IVIM,'.',label ="data")
                #plt.plot(bvals_unique,np.squeeze(Fit),'--', label="fitted")
                #plt.show()

                S0map[xi, yi,i] = Fitted_Parameters[0]
                Dmap[xi, yi,i] = Fitted_Parameters[1]*1000
                #Dsmap[xi,yi,i] = Fitted_Parameters[2]*1000
                #fmap[xi,yi,i] = Fitted_Parameters[3]
                rsquaremap[xi,yi,i] = r_squared
                #else:
                #S0map[xi, yi,i] = 0
                #Dmap[xi, yi,i] = 0
                #Dsmap[xi,yi,i] = 0
                #fmap[xi,yi,i] = 0
                #rsquaremap[xi,yi,i] = 0

    #fittedMaps = np.squeeze(S0map), np.squeeze(Dmap), np.squeeze(Dsmap),np.squeeze(fmap), np.squeeze(rsquaremap)
    fittedMaps = np.squeeze(S0map), np.squeeze(Dmap),np.squeeze(rsquaremap)
    #print('iupi')
    return fittedMaps