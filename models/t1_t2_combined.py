import os
import multiprocessing
import tqdm
import numpy as np

import models.t1
import models.t2


# Hardwired sequence parameters
TR = 4.6                            #in ms
FA = 12    #in degrees
FA_rad = FA/360*(2*np.pi)           #convert to rads
N_T1 = 66                           #number of k-space lines
FA_Cat  = [
    (-FA/5)/360*(2*np.pi), 
    (2*FA/5)/360*(2*np.pi), 
    (-3*FA/5)/360*(2*np.pi), 
    (4*FA/5)/360*(2*np.pi), 
    (-5*FA/5)/360*(2*np.pi)] #cat module
TE = [0,30,40,50,60,70,80,90,100,110,120]
Tspoil = 1
N_T2 = 72
Trec = 463*2
FA_eff = 0.6



def pars():
    return [
        "T1_T1_Map_v2",
        "T2_T2_Map_v2",
        "T1_S0_Map_v2",
        "T2_S0_Map_v2",
        "T1_FA_Eff_Map_v2",
        "T1_rsquare_Map_v2",
        "T2_rsquare_Map_v2",
     ]


def t1_t2_alone(arguments):

    try:
        sig_T1, sig_T2, TI_temp = arguments

        fit_T1, pars_T1 = models.t1.main(sig_T1, TI_temp, [FA_rad, TR, N_T1, FA_Cat])
                                                                                                            
        S0_T1,T1,FA_eff = pars_T1

        fit_T2, pars_T2 = models.t2.main(sig_T2, TE, [T1,Tspoil,FA_rad,TR, N_T2,Trec,FA_eff])

        S0_T2, T2, FA_eff_2 = pars_T2

        #r squared calculation 

        ss_res_T1 = np.sum(np.nan_to_num((sig_T1-fit_T1)**2))
        ss_res_T2 = np.sum(np.nan_to_num((sig_T2-fit_T2)**2))
        ss_tot_T1 = np.sum(np.nan_to_num((sig_T1-np.nanmean(sig_T1))**2))
        ss_tot_T2 = np.sum(np.nan_to_num((sig_T2-np.nanmean(sig_T2))**2))

        if ss_tot_T1  == 0:
            r_squared_T1 = 0
        else:
            r_squared_T1 = np.nan_to_num(1 - (ss_res_T1 / ss_tot_T1))
        
        if ss_tot_T2  == 0:
            r_squared_T2 = 0
        else:
            r_squared_T2 = np.nan_to_num(1 - (ss_res_T2 / ss_tot_T2))

        #replace possible nan (from division by 0: ss_res_T1 / ss_tot_T1) to 0
        if (np.isnan(r_squared_T1)): r_squared_T1 = 0
        if (np.isnan(r_squared_T2)): r_squared_T2 = 0
        
    except:
        T1 = T2 = S0_T1 = S0_T2 = FA_eff = r_squared_T1 = r_squared_T2 = 0
        
    return T1, T2, S0_T1, S0_T2, FA_eff, r_squared_T1, r_squared_T2


def fit(arrays, TI_slice):

    shape = arrays[0].shape
    array_T1 = arrays[0].reshape((shape[0]*shape[1], shape[2], shape[3]))
    array_T2 = arrays[1].reshape((shape[0]*shape[1], shape[2], shape[3]))
    par = np.empty( (shape[0]*shape[1], shape[2], shape[3], len(pars())) )

    for z in range(array_T1.shape[1]):

        pool = multiprocessing.Pool(initializer=multiprocessing.freeze_support,processes=os.cpu_count())

        arguments = [(array_T1[i,z,:], array_T2[i,z,:], TI_slice[z]) for i in range(shape[0]*shape[1])]
        params = list(tqdm(pool.imap(t1_t2_alone, arguments), total=len(arguments), desc='Processing pixels of slice ' + str(i)))

        for i, p in enumerate(params):
            par[i,z,:] = p

    return par.reshape((shape[0], shape[1], shape[2], shape[3], len(pars())))



