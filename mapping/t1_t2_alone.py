""" 
@author: Joao Periquito 
iBEAt T1 & T2 Custom modelling
2022
Apply forward modelling to T1 and T2 data
"""

import numpy as np
import mapping.t2


def main(arguments):

    try:
        x,y,t1_map,t2_value,TE,FA_rad,TR,N_T2,Trec,FA_eff,Tspoil = arguments

        fit_T2, fitted_parameters_T2 = mapping.t2.main (t2_value, TE,[t1_map,Tspoil,FA_rad,TR, N_T2,Trec,FA_eff])

        S0_T2, T2, FA_eff_2 =  fitted_parameters_T2
        
        residuals_T2 = t2_value-np.squeeze(fit_T2) 

        #r squared calculation 
 
        ss_res_T2 = np.sum(np.nan_to_num(residuals_T2**2))
        ss_tot_T2 = np.sum(np.nan_to_num((t2_value-np.nanmean(t2_value))**2))
        
        if ss_tot_T2  == 0:
            r_squared_T2 = 0
        else:
            r_squared_T2 = np.nan_to_num(1 - (ss_res_T2 / ss_tot_T2))


        #replace possible nan (from division by 0: ss_res_T1 / ss_tot_T1) to 0
        if (np.isnan(r_squared_T2)): r_squared_T2 = 0
        
    except:
        T2 = S0_T2 = FA_eff = r_squared_T2 = 0
        
    return x,y, T2, S0_T2, FA_eff, r_squared_T2