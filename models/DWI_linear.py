# Alternative to DTI_dipy
# Produces the same results so this has become superfluous.
# Keeping it for its transparency for now

import numpy as np

def pars():
    return ['Dxx', 'Dyy', 'Dzz', 'S0']


def fit(im:np.ndarray, bvals:np.ndarray=None, bvecs:np.ndarray=None):
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

