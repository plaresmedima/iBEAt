import numpy as np
from utilities.fit import Model


class Lin():
    # Can be replaced by dipy implementation

    def pars(self):
        return ['Dxx', 'Dyy', 'Dzz', 'S0']

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

        # Values corresponding to off diagonal elements are all zero
        # Take out and fit only Dxx, Dyy, Dzz
        # Sort all b matrices in to a vector b=[bxx,byy,bzz];
        b = np.vstack((b[0,0,:],
                        b[1,1,:],
                        b[2,2,:])).T
        
        # take the log of the image to linearise the equation
        # Use absolute value in case some essentially zero signals are negative due to rounding errors
        imlog = np.log(np.abs(im))

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
        
        par = np.empty((im.shape[0], 4))
        par[:,:3] = X[:3,:].T
        par[:,3] = np.exp(X[3,:]).T

        # Return in original shape
        fit = fit.reshape(shape)
        par = par.reshape(shape[:-1] + (4,))
        
        return fit, par
    

class NonLin(Model):

    # Can be replaced by dipy implementation

    def pars(self):
        return ['Dxx', 'Dyy', 'Dzz', 'S0']

    def init(self, sig):
        smax=np.amax(sig)
        if smax<=0:
            return [0.0025, 0.0025, 0.0025, 0]
        else:
            return [0.0025, 0.0025, 0.0025, smax]

    def bnds(self, sig):
        smax = np.amax(sig)
        if smax<=0:
            upper = [0.1, 0.1, 0.1, np.inf]
        else:
            upper = [0.1, 0.1, 0.1, 2*smax]
        lower = [0,0,0,0]
        return lower, upper

    def signal(self, bvals, Dxx, Dyy, Dzz, s0):
        return np.concatenate([
            s0*np.exp(-bvals*Dxx),
            s0*np.exp(-bvals*Dyy),
            s0*np.exp(-bvals*Dzz),
        ])


