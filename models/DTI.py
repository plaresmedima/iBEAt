import numpy as np

from dipy.core.gradients import gradient_table
import dipy.reconst.dti as dti


class DiPy():

    def pars(self):
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

    def fit(self, array, bvals:np.ndarray=None, bvecs:np.ndarray=None, fit_method="LS"):
        # array (x,y,z,t)
        # bvals (t)
        # bvecs (t) 
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

        tenfit.S0_hat[tenfit.S0_hat<0] = 0
        tenfit.S0_hat[tenfit.S0_hat>100000] = 100000
        tenfit.fa[tenfit.fa<0] = 0
        tenfit.fa[tenfit.fa>1] = 1
        tenfit.md[tenfit.md<0] = 0
        tenfit.md[tenfit.md>0.05] = 0.05
        tenfit.sphericity[tenfit.sphericity<0] = 0
        tenfit.sphericity[tenfit.sphericity>1] = 1
        tenfit.linearity[tenfit.linearity<0] = 0
        tenfit.linearity[tenfit.linearity>1] = 1
        tenfit.planarity[tenfit.planarity<0] = 0
        tenfit.planarity[tenfit.planarity>1] = 1
        tenfit.ad[tenfit.ad<0] = 0
        tenfit.ad[tenfit.ad>0.05] = 0.05
        tenfit.rd[tenfit.rd<0] = 0
        tenfit.rd[tenfit.rd>0.05] = 0.05


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
    

class Lin(): #Obsolete - replace by DiPy

    def pars(self):
        return ['Dxx', 'Dxy', 'Dxz', 'Dyy', 'Dyz', 'Dzz', 'S0']


    def fit(self, im:np.ndarray, bvals:np.ndarray=None, bvecs:np.ndarray=None):
        # im = (x,t) or (x,y,t) or (x,y,z,t)

        # Reshape to 2D
        shape = im.shape
        im = im.reshape((-1,shape[-1]))

        # b matrix
        b = np.zeros((3, 3, len(bvals)))
        for i in range(len(bvals)):
            # bvec at bval=0 is undefined so better not use - could be None
            if bvals[i] > 0: 
                bvec_i = np.array(bvecs[i])
                b[:, :, i] = np.outer(bvals[i]*bvec_i, bvec_i)
                
        # Sort all b matrices in to a vector b=[bxx,2*bxy,2*bxz,byy,2*byz,bzz];
        b = np.vstack((b[0,0,:],
                        2*b[0,1,:],
                        2*b[0,2,:],
                        b[1,1,:],
                        2*b[1,2,:],
                        b[2,2,:])).T
            
        # take the log of the image to linearise the equation
        imlog = np.log(im)

        # Where the image is zero, the log is infinite. 
        zeros = np.isinf(imlog)
        imlog[zeros] = 0
        
        # Add another column to handle the constant term:
        # S = S0 exp(-bD)
        # log(S) = log(S0) - b*D
        # becomes:
        # Slog = [-b, 1] * [D; log(S0)] = [-b, 1] * X
        mat = np.c_[-b, np.ones(b.shape[0])]
        X = np.linalg.lstsq(mat, imlog.T, rcond=None)[0]

        fit = np.exp(np.dot(mat, X).T)
        fit[zeros] = 0
        
        par = np.empty((im.shape[0], 7))
        par[:,:6] = X[:6,:].T
        par[:,6] = np.exp(X[6,:]).T

        # Return in original shape
        fit = fit.reshape(shape)
        par = par.reshape(shape[:-1] + (7,))
        
        return fit, par