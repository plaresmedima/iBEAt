import numpy as np
from scipy.optimize import curve_fit
from utilities.fit import fit_image



# M = Mss + (Minit-Mss) exp(-TI*R1app)
# Minit = -M0
# M = Mss - (M0+Mss) exp(-TI*R1app)
# B/A = (M0+Mss)/Mss
# B/A = (1-cosFA*E) / (1-E)  +  1
# (1 - (B/A-1)*(1-E)) / E = cosFA

# This works but less eleagnt - using the direct fitting
# def der(par, TR, FAapp):
#     M0, Mss, T1app = par[...,0], par[...,1], par[...,2]
#     zeros = (M0==0)
#     T1 = (Mss/M0) * T1app
#     T1[zeros] = 0
#     T1[T1>3000] = 3000
#     T1[T1<0] = 0
#     E = np.exp(-TR/T1)
#     E[zeros] = 0
#     cosFA = (1 - (Mss/M0)*(1-E)) / E
#     cosFA[cosFA>1] = 1
#     cosFA[cosFA<-1] = -1
#     cosFA[zeros] = np.cos(FAapp*np.pi/180)
#     FA = np.arccos(cosFA)*180/np.pi
#     return T1, FA

# def pars():
#     return ['S0', 'Sss', 'T1_apparent']

# def init(sig):
#     smax = np.amax(sig)
#     return [smax, smax, 1000]

# def bnds(sig):
#     smax = np.amax(sig)
#     return [0, 0, 0], [2*smax, 2*smax, 3000]

# def signal(TI, M0, Mss, T1app):
#     if T1app==0:
#         return np.abs(M0)
#     return np.abs(M0 - (M0+Mss) * np.exp(-TI/T1app))
    
    
def pars():
    return ['S0', 'T1', 'FAcorr']

def init(sig):
    return [np.amax(sig), 1500, 12]

def bnds(sig):
    smax = np.amax(sig)
    return ([0, 50, 0], [smax, 3000, 25])

def signal(TI, M0, T1, FA):
    TR = 4.6 # Needs a keyword argument
    if T1==0:
        return np.abs(M0)
    cFA = np.cos(FA*np.pi/180)
    E = np.exp(-TR/T1)
    E = (1-E) / (1-cFA*E)
    Mss = M0 * E
    T1app = T1 * E
    return np.abs(Mss - (M0+Mss) * np.exp(-TI/T1app))
    
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