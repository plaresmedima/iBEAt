import numpy as np
from utilities.fit import Model
from models import relaxometry


class MonoExp(Model):

    def pars(self):
        return ['S0', 'T1', 'T1FAcorr']

    def init(self, sig):
        return [np.amax(sig), 1500, 12]

    def bnds(self, sig):
        smax = np.amax(sig)
        return ([0, 50, 0], [2*smax, 3000, 25])

    def signal(self, TI, M0, T1, FA, TR=5.0):
        if T1==0:
            return np.abs(M0)
        cFA = np.cos(FA*np.pi/180)
        E = np.exp(-TR/T1)
        E = (1-E) / (1-cFA*E)
        Mss = M0 * E
        T1app = T1 * E
        return np.abs(Mss - (M0+Mss) * np.exp(-TI/T1app))



class Bloch(Model):

    def pars(self):
        return ['S0', 'T1', 'T1FAcorr']

    def init(self, sig):
        return [np.amax(sig), 1200, 12]

    def bnds(self, sig):
        smax = np.amax(sig)
        return ([0, 50, 0], [2*smax, 3000, 25])

    def signal(self, TI, S0, T1, FAcorr, 
            # Defaults from Siemens iBEAt protocol
            TR = 4.6, # msec
            FA_cat = [-1, 2, -3, 4, -5], # Catalization module confirmed by Siemens (Peter Schmitt): Magn Reson Med 2003 Jan;49(1):151-7. doi: 10.1002/mrm.10337
            N_T1 = 66, # Number of k-space lines (hardcoded from Siemens protocol)
            FA = 12, # Flip angle in degrees (hardcoded from Siemens protocol)
            ):
        FA_cat = [(f*FA/len(FA_cat))/180*np.pi for f in FA_cat]
        FA_eff = FAcorr/FA
        FA = FA/180*np.pi
        return relaxometry.signalSequenceT1_FLASH(S0, T1, TI, FA, FA_eff, TR, N_T1, FA_cat)
