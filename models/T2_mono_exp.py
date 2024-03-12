import numpy as np
from scipy.optimize import curve_fit
from utilities.fit import fit_image


def pars():
    return ['S0', 'T2']

def init(sig):
    smax = np.amax(sig)
    if smax<=0:
        return [0, 80]
    else:
        return [smax, 80]

def bnds(sig):
    smax = np.amax(sig)
    if smax<=0:
        return [0,10,0], [np.inf, 500]
    else:
        return [0,10,0], [5*smax, 500]

def signal(TE,S0,T2):
    if T2==0:
        return TE*0
    return S0*np.exp(-TE/T2)

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