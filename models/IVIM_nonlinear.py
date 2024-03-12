import numpy as np

from scipy.optimize import curve_fit

import models.DWI_nonlinear as dwi
from utilities.fit import fit_image


def init(sig):
    smax = np.amax(sig)
    if smax<=0:
        return [0.0025, 0]
    else:
        return [0.0025, smax]

def bnds(sig):
    smax = np.amax(sig)
    if smax<=0:
        upper = [1, np.inf]
    else:
        upper = [1.0, 2*smax]
    lower = [0,0]
    return lower, upper

def signal(bvals, Df, Sf, Dxx, Dyy, Dzz, Sd):
    Sflow = Sf*np.exp(-bvals*Df)
    Sflow = np.concatenate([Sflow,Sflow,Sflow])
    Sdif = dwi.signal(bvals, Dxx, Dyy, Dzz, Sd)
    return Sflow + Sdif

def fit_signal(args):
    bvals, ydata, xtol, bounds = args
    if bounds==True:
        bounds_dwi = dwi.bnds(ydata)
        bounds_ivim = bnds(ydata)
    else:
        bounds_dwi = (-np.inf, np.inf)
        bounds_ivim = (-np.inf, np.inf)
    pars_dwi = dwi.init(ydata)
    pars_ivim = init(ydata)
    try:
        # # Fit the DWI model to the last 3 b-values
        n = int(len(ydata)/3)
        ydif = np.concatenate([ydata[n-3:n], ydata[2*n-3:2*n], ydata[3*n-3:3*n]])
        pars_dwi, _ = curve_fit(
            dwi.signal, 
            bvals[-3:], ydif, p0=pars_dwi, xtol=xtol, bounds=bounds_dwi)
        # Fit the IVIM signal on all data points, but fixing the DWI part
        pars_ivim, _ = curve_fit(
            lambda b, Df, Sf: signal(b, Df, Sf, *tuple(pars_dwi)),
            bvals, ydata, p0=pars_ivim, xtol=xtol, bounds=bounds_ivim)
    except RuntimeError:
        # If the model does not fit, return initial values
        pass
    return tuple(pars_ivim) + tuple(pars_dwi) # Df, Sf, Dxx, Dyy, Dzz, Sd


def pars():
    return ['Df', 'Sf', 'Dxx', 'Dyy', 'Dzz', 'Sd']


def fit(*args, **kwargs):
    return fit_image(signal, fit_signal, *args, **kwargs)


def derived(par):
    Df = par[...,0]
    S0f = par[...,1] 
    MD = np.mean(par[...,2:5], axis=-1)
    S0d = par[...,5]
    S0 = S0f+S0d
    ff = S0f/S0
    ff[S0==0]=0
    ff[ff>1]=1
    ff[ff<0]=0
    return S0, Df, MD, ff
