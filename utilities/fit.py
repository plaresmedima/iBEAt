import os 
import multiprocessing
import numpy as np

try: 
    num_workers = int(len(os.sched_getaffinity(0)))
except: 
    num_workers = int(os.cpu_count())


def fit_image(signal, fit_signal, imgs:np.ndarray, xdata=None, xtol=1e-3, bounds=False, parallel=True):
    """Fit a single-pixel model pixel-by-pixel to a 2D or 3D image"""
    
    # Reshape to (x,t)
    shape = imgs.shape
    imgs = imgs.reshape((-1,shape[-1]))
    
    # Perform the fit pixelwise
    if parallel:
        args = [(xdata, imgs[x,:], xtol, bounds) for x in range(imgs.shape[0])]
        pool = multiprocessing.Pool(processes=num_workers)
        fit_pars = pool.map(fit_signal, args)
        pool.close()
        pool.join()
    else: # for debugging
        fit_pars = [fit_signal((xdata, imgs[x,:], xtol, bounds)) for x in range(imgs.shape[0])]

    # Create output arrays
    npars = len(fit_pars[0])
    fit = np.empty(imgs.shape)
    par = np.empty((imgs.shape[0], npars))
    for x, p in enumerate(fit_pars):
        fit[x,:] = signal(xdata, *p)
        par[x,:] = p

    # Return in original shape
    fit = fit.reshape(shape)
    par = par.reshape(shape[:-1] + (npars,))
       
    return fit, par
