import os 
import multiprocessing
import numpy as np
from scipy.optimize import curve_fit

try: 
    num_workers = int(len(os.sched_getaffinity(0)))
except: 
    num_workers = int(os.cpu_count())




def fit_image(signal, fit_signal, imgs:np.ndarray, xdata=None, xtol=1e-3, bounds=False, parallel=True, const=None, **kwargs):
    """Fit a single-pixel model pixel-by-pixel to a 2D or 3D image"""
    
    # Reshape to (x,t)
    shape = imgs.shape
    imgs = imgs.reshape((-1,shape[-1]))
    if const is not None:
        const = np.ravel(const)
    
    # Perform the fit pixelwise
    if parallel:
        #args = [(xdata, imgs[x,:], xtol, bounds, kwargs) for x in range(imgs.shape[0])]
        args = []
        for x in range(imgs.shape[0]):
            if const is None:
                xdata_x = np.array(xdata)
            else:
                const_pixel = np.array([const[x]])
                xdata_x = np.concatenate((const_pixel, xdata))
            args_x = (xdata_x, imgs[x,:], xtol, bounds, kwargs)
            args.append(args_x)
        pool = multiprocessing.Pool(processes=num_workers)
        fit_pars = pool.map(fit_signal, args)
        pool.close()
        pool.join()

    else: # for debugging
        fit_pars = []
        for x in range(imgs.shape[0]):
            if const is None:
                xdata_x = np.array(xdata)
            else:
                const_pixel = np.array([const[x]])
                xdata_x = np.concatenate((const_pixel, xdata))
            fit_pars_x = fit_signal((xdata_x, imgs[x,:], xtol, bounds, kwargs))
            fit_pars.append(fit_pars_x)

        #fit_pars = [fit_signal((xdata_x, imgs[x,:], xtol, bounds, kwargs)) for x in range(imgs.shape[0])]

    # Create output arrays
    npars = len(fit_pars[0])
    fit = np.empty(imgs.shape)
    par = np.empty((imgs.shape[0], npars))
    for x, p in enumerate(fit_pars):
        fit[x,:] = signal(xdata, *p, **kwargs)
        par[x,:] = p

    # Return in original shape
    fit = fit.reshape(shape)
    par = par.reshape(shape[:-1] + (npars,))
       
    return fit, par



class Model:

    def pars(self):
        raise NotImplementedError('No pars method implemented - this is required.')
    
    def bnds(self, *args):
        return (-np.inf, np.inf)
    
    def init(self, *args):
        raise NotImplementedError('No init method implemented - this is required.')
    
    def signal(self, *args, **kwargs):
        raise NotImplementedError('No signal method implemented - this is required.')
    
    def fit_signal(self, args):
        xdata, ydata, xtol, bounds, kwargs = args
        if bounds==True:
            bounds = self.bnds(ydata)
        else:
            bounds = (-np.inf, np.inf)
        pars = self.init(ydata)
        def fit_func(x,*p):
            return self.signal(x, *p, **kwargs)   
        try:
            pars, _ = curve_fit(fit_func, xdata, ydata, p0=pars, xtol=xtol, bounds=bounds)
        except:
            pass
        return pars
    
    def fit(self, *args, **kwargs):
        return fit_image(self.signal, self.fit_signal, *args, **kwargs)
    
    def error(self, array, fit):
        ref = np.linalg.norm(array, axis=-1)
        err = 100*np.linalg.norm(fit-array, axis=-1)/ref
        err[ref==0] = 0
        return err