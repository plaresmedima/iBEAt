""" 
@author: Joao Periquito 
iBEAt T1 & T2 Custom modelling
2022
Apply forward modelling to T1 and T2 data
"""

import numpy as np
import mapping.t1

def main(arguments):

    try:
        x,y,t1_value,TI_temp,FA_rad,TR,N_T1,FA_Cat,FA_eff  = arguments

        fit_T1, fitted_parameters_T1 = mapping.t1.main (t1_value, TI_temp, [FA_rad, TR, N_T1,FA_Cat])
                                                                                                            
        S0_T1,T1,FA_eff = fitted_parameters_T1
        
        residuals_T1 = t1_value-np.squeeze(fit_T1) 

        #r squared calculation 

        ss_res_T1 = np.sum(np.nan_to_num(residuals_T1**2))
        ss_tot_T1 = np.sum(np.nan_to_num((t1_value-np.nanmean(t1_value))**2))

        if ss_tot_T1  == 0:
            r_squared_T1 = 0
        else:
            r_squared_T1 = np.nan_to_num(1 - (ss_res_T1 / ss_tot_T1))


        #replace possible nan (from division by 0: ss_res_T1 / ss_tot_T1) to 0
        if (np.isnan(r_squared_T1)): r_squared_T1 = 0
        
    except:
        T1 = S0_T1 = FA_eff = r_squared_T1 =  0
        
    return x,y,T1, S0_T1, FA_eff, r_squared_T1