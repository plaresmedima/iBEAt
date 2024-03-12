import os 
import numpy as np
from scipy.optimize import curve_fit
from utilities.fit import fit_image


def pars():
    return ['Dxx', 'Dyy', 'Dzz', 'S0']

def init(sig):
    smax=np.amax(sig)
    if smax<=0:
        return [0.0025, 0.0025, 0.0025, 0]
    else:
        return [0.0025, 0.0025, 0.0025, smax]

def bnds(sig):
    smax = np.amax(sig)
    if smax<=0:
        upper = [0.1, 0.1, 0.1, np.inf]
    else:
        upper = [0.1, 0.1, 0.1, 2*smax]
    lower = [0,0,0,0]
    return lower, upper

def signal(bvals, Dxx, Dyy, Dzz, s0):
    return np.concatenate([
        s0*np.exp(-bvals*Dxx),
        s0*np.exp(-bvals*Dyy),
        s0*np.exp(-bvals*Dzz),
    ])


def fit_signal(args):
    xdata, ydata, xtol, bounds = args
    if bounds==True:
        bounds = bnds(ydata)
    else:
        bounds = (-np.inf, np.inf)
    pars = init(ydata)
    try:
        pars, _ = curve_fit(signal, xdata, ydata, p0=pars, xtol=xtol, bounds=bounds)
    except:
        pass
    return pars


def fit(*args, **kwargs):
    return fit_image(signal, fit_signal, *args, **kwargs)

