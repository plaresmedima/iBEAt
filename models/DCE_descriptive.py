import numpy as np

def fit(array, aif, dt, baseline=1):

    # Calculate baseline
    Sa0 = np.mean(aif[:baseline])
    S0 = np.mean(array[...,:baseline], axis=-1)

    # Subtract baseline
    ca = aif-Sa0
    # for t in range(array.shape[-1]):
    #     array[...,t] = array[...,t]-S0

    # Calculate parameters
    max_enh = np.amax(array, axis=-1) - S0
    sum_enh = np.sum(array, axis=-1) - S0*array.shape[-1]
    MAX = max_enh/np.amax(ca)
    AUC = sum_enh/np.sum(ca)
    MTT = dt*sum_enh/max_enh
    MTT[max_enh==0] = 0

    # Apply bounds
    MAX*=100
    MAX[MAX<0]=0
    MAX[MAX>200]=200
    AUC*=100 # mL/100mL
    AUC[AUC<0]=0
    AUC[AUC>500]=500
    MTT[MTT<0] = 0
    MTT[MTT>600] = 600

    return MAX, AUC, MTT