import numpy as np

def r_square(signal,fit):

    residuals = signal - fit
    # ss_res = np.sum(residuals**2)
    # ss_tot = np.sum((signal-np.mean(signal))**2)
    # r_squared = 1 - (ss_res / ss_tot)

    return 100*np.sum(residuals**2)/np.sum(signal**2)

