import numpy as np

def fit(imgs:np.ndarray, aif, dt, baseline=1, regpar=0.1, Hct=0.45):

    # Reshape to 2D (x,t)
    shape = np.shape(imgs)
    imgs = imgs.reshape((-1,shape[-1]))

    # Calculate baseline
    S0 = np.mean(imgs[:,:baseline], axis=1)
    Sa0 = np.mean(aif[:baseline])
    
    # Subtract baseline
    ca = (aif-Sa0)/(1-Hct)
    for t in range(imgs.shape[1]):
        imgs[:,t] = imgs[:,t]-S0

    # Create AIF matrix
    nt = len(ca)
    A = np.zeros((nt,nt))
    for k in range(nt):
        A[k:,k] = dt*ca[:nt-k]

    # Invert A with SVD
    U, S, Vh = np.linalg.svd(A, full_matrices=False) # A = np.dot(U*S, Vh)
    Sinv = np.zeros(len(S))
    Smin = regpar*np.amax(S)
    Spos = S[S>Smin]
    Sinv[:len(Spos)] = 1/S[:len(Spos)]
    Ainv = np.dot(Vh.T * Sinv, U.T)

    # Deconvolve
    IRF = np.dot(Ainv, imgs.T).T

    # Derive parameters
    # Return in original shape
    RPF = np.amax(IRF, axis=1).reshape(shape[:-1])
    AVD = dt*np.sum(IRF, axis=1).reshape(shape[:-1])
    MTT = AVD/RPF
    MTT[RPF==0] = 0

    # Convert to conventional units and apply bounds
    RPF*=6000 # mL/min/100mL
    RPF[RPF<0]=0
    RPF[RPF>2000]=2000
    AVD*=100 # mL/100mL
    AVD[AVD<0]=0
    AVD[AVD>500]=500
    MTT[MTT<0] = 0 #sec
    MTT[MTT>600] = 600

    return RPF, AVD, MTT
