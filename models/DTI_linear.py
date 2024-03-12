# Alternative to DTI_dipy
# Produces the same results so this has become superfluous.
# Keeping it for its transparency for now

import numpy as np

def pars():
    return ['Dxx', 'Dxy', 'Dxz', 'Dyy', 'Dyz', 'Dzz', 'S0']


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

