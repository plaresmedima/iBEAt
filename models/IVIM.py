import numpy as np
from utilities.fit import Model
from dipy.core.gradients import gradient_table
from dipy.reconst.ivim import IvimModel
from scipy.optimize import curve_fit
from models.DWI import NonLin as DWI
import multiprocessing
import os

try:
    num_workers = int(len(os.sched_getaffinity(0)))
except:
    num_workers = int(os.cpu_count())


class DiPy():

    def pars(self):
        return [
            'S0',
            'Perfusion_fraction',
            'D_star',
            'D',
            
        ]

    def fit(self, array, bvals:np.ndarray=None, bvecs:np.ndarray=None, fit_method="TRR"):
        # array (x,y,z,t)
        # bvals (t)
        # bvecs (t) 

        # 2-stage Trust-Region Reflective
        # based NLLS fitting method (also containing the linear fit): `trr` and the
        # Variable Projections based fitting method: `varpro`.
        gtab = gradient_table(bvals, bvecs, b0_threshold=0)
        ivimmodel_vp = IvimModel(gtab, fit_method=fit_method, bounds=([0., 0., 0.,0.], [2*np.max(array), 1., 1., 1.]),split_b_S0=100,split_b_D=100)
        ivim_fit = ivimmodel_vp.fit(array)
        
        fit = ivim_fit.predict(gtab)

        ivim_fit.S0_predicted[ivim_fit.S0_predicted<0] = 0
        ivim_fit.S0_predicted[ivim_fit.S0_predicted>100000] = 100000
        ivim_fit.perfusion_fraction[ivim_fit.perfusion_fraction<0] = 0
        ivim_fit.perfusion_fraction[ivim_fit.perfusion_fraction>1] = 1 
        ivim_fit.D_star[ivim_fit.D_star<0] = 0
        ivim_fit.D_star[ivim_fit.D_star>1] = 1
        ivim_fit.D[ivim_fit.D<0] = 0
        ivim_fit.D[ivim_fit.D>1] = 1 

        pars = ( 
            ivim_fit.S0_predicted,
            ivim_fit.perfusion_fraction, 
            ivim_fit.D_star, 
            ivim_fit.D,   
        )

        pars = np.stack(pars, axis=-1)

        return fit, pars


class SemiLin(Model):

    def pars(self):
        return ['Df', 'Sf', 'Dxx', 'Dyy', 'Dzz', 'Sd']

    def init(self, sig):
        smax = np.amax(sig)
        if smax<=0:
            return (0.0025, 0)
        else:
            return (0.0025, smax)

    def bnds(self, sig):
        smax = np.amax(sig)
        if smax<=0:
            upper = [1, np.inf]
        else:
            upper = [1.0, 2*smax]
        lower = [0,0]
        return lower, upper

    def signal(self, bvals, Df, Sf, Dxx, Dyy, Dzz, Sd):
        Sflow = Sf*np.exp(-bvals*Df)
        return np.concatenate([
            Sflow + Sd*np.exp(-bvals*Dxx),
            Sflow + Sd*np.exp(-bvals*Dyy),
            Sflow + Sd*np.exp(-bvals*Dzz),
        ])

    def fit_signal(self, args):
        bvals, ydata, pars_dwi, xtol, bounds, kwargs = args
        if bounds==True:
            bounds = self.bnds(ydata)
        else:
            bounds = (-np.inf, np.inf)
        pars = self.init(ydata)
        try:
            # Fit the IVIM signal on all data points, but fixing the DWI part
            pars, _ = curve_fit(
                lambda b, Df, Sf: self.signal(b, Df, Sf, *tuple(pars_dwi)),
                bvals, ydata, p0=pars, xtol=xtol, bounds=bounds)
        except RuntimeError:
            # If the model does not fit, return initial values
            pass
        return tuple(pars) + tuple(pars_dwi) # Df, Sf, Dxx, Dyy, Dzz, Sd


    def fit(self, imgs:np.ndarray, bvals=None, bvecs=None, xtol=1e-3, bounds=False, parallel=True):

        # Reshape to (x,t)
        shape = imgs.shape
        imgs = imgs.reshape((-1,shape[-1]))

        # Fit the DWI parameters linearly to the last 3 b-values
        k, n = 3, len(bvals)
        imgs_dwi = np.concatenate([imgs[:,n-k:n], imgs[:,2*n-k:2*n], imgs[:,3*n-k:3*n]], axis=-1)
        bvals_dwi = np.concatenate([bvals[n-k:], bvals[n-k:], bvals[n-k:]])
        bvecs_dwi = np.concatenate([bvecs[n-k:], bvecs[n-k:], bvecs[n-k:]])
        _, pars_dwi = DWI().fit(imgs_dwi, bvals_dwi, bvecs_dwi)

        # Perform the fit pixelwise
        if parallel:
            args = [(bvals, imgs[x,:], pars_dwi[x,:], xtol, bounds) for x in range(imgs.shape[0])]
            pool = multiprocessing.Pool(processes=num_workers)
            fit_pars = pool.map(self.fit_signal, args)
            pool.close()
            pool.join()
        else: # for debugging
            fit_pars = [self.fit_signal((bvals, imgs[x,:], pars_dwi[x,:], xtol, bounds)) for x in range(imgs.shape[0])]

        # Create output arrays
        npars = len(self.pars())
        fit = np.empty(imgs.shape)
        par = np.empty((imgs.shape[0], npars))
        for x, p in enumerate(fit_pars):
            fit[x,:] = self.signal(bvals, *p)
            par[x,:] = p

        # Return in original shape
        fit = fit.reshape(shape)
        par = par.reshape(shape[:-1] + (npars,))
        
        return fit, par

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



class NonLin(Model):

    def pars(self):
        return ['Df', 'Sf', 'Dxx', 'Dyy', 'Dzz', 'Sd']

    def init(self, sig):
        smax = np.amax(sig)
        if smax<=0:
            return [0.0025, 0]
        else:
            return [0.0025, smax]

    def bnds(self, sig):
        smax = np.amax(sig)
        if smax<=0:
            upper = [1, np.inf]
        else:
            upper = [1.0, 2*smax]
        lower = [0,0]
        return lower, upper

    def signal(self, bvals, Df, Sf, Dxx, Dyy, Dzz, Sd):
        Sflow = Sf*np.exp(-bvals*Df)
        Sflow = np.concatenate([Sflow,Sflow,Sflow])
        Sdif = DWI().signal(bvals, Dxx, Dyy, Dzz, Sd)
        return Sflow + Sdif

    def fit_signal(self, args):
        bvals, ydata, xtol, bounds, kwargs = args
        if bounds==True:
            bounds_dwi = DWI().bnds(ydata)
            bounds_ivim = self.bnds(ydata)
        else:
            bounds_dwi = (-np.inf, np.inf)
            bounds_ivim = (-np.inf, np.inf)
        pars_dwi = DWI().init(ydata)
        pars_ivim = self.init(ydata)
        try:
            # # Fit the DWI model to the last 3 b-values
            n = int(len(ydata)/3)
            ydif = np.concatenate([ydata[n-3:n], ydata[2*n-3:2*n], ydata[3*n-3:3*n]])
            pars_dwi, _ = curve_fit(
                DWI().signal, 
                bvals[-3:], ydif, p0=pars_dwi, xtol=xtol, bounds=bounds_dwi)
            # Fit the IVIM signal on all data points, but fixing the DWI part
            pars_ivim, _ = curve_fit(
                lambda b, Df, Sf: self.signal(b, Df, Sf, *tuple(pars_dwi)),
                bvals, ydata, p0=pars_ivim, xtol=xtol, bounds=bounds_ivim)
        except RuntimeError:
            # If the model does not fit, return initial values
            pass
        return tuple(pars_ivim) + tuple(pars_dwi) # Df, Sf, Dxx, Dyy, Dzz, Sd

    def derived(self, par):
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
