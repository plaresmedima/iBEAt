import numpy as np
from scipy.optimize import curve_fit
from utilities.fit import fit_image


def pars():
    return ['S0', 'T2star']

def init(sig):
    return [np.amax(sig), 50]

def bnds(sig):
    return [0,5], [5*np.amax(sig), 200]

def signal(TE, S0, T2star):
    if T2star==0:
        return TE*0
    return S0*np.exp(-TE/T2star)

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

