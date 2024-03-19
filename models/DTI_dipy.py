import numpy as np

from dipy.core.gradients import gradient_table
import dipy.reconst.dti as dti


def pars():
    return [
        'S0',
        'FA',
        'MD',
        'Sphericity',
        'Linearity',
        'Planarity',
        'AD',
        'RD',
    ]

def fit(array, bvals:np.ndarray=None, bvecs:np.ndarray=None, fit_method="WLS"):
    # https://github.com/dipy/dipy/blob/master/dipy/reconst/dti.py

    # Fit DTI model
    # Alternative fit methods are:
    #  'WLS' weighted least squares;
    #  'LS' ordinary least squares;
    #  'NLLS' non-linear least quares; 
    #  'RT' for RESTORE.
    gtab = gradient_table(bvals, bvecs)
    tenmodel = dti.TensorModel(gtab, fit_method=fit_method, return_S0_hat=True)
    tenfit = tenmodel.fit(array)
    
    # MD[MD>0.005]=0.005
    # MD[MD<0]=0

    fit = tenfit.predict(gtab)
    pars = ( 
        tenfit.S0_hat, 
        tenfit.fa, 
        tenfit.md, 
        tenfit.sphericity, 
        tenfit.linearity, 
        tenfit.planarity, 
        tenfit.ad, 
        tenfit.rd,
    )
    pars = np.stack(pars, axis=-1)

    return fit, pars