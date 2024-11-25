import numpy as np
from utilities.fit import Model



class MonoExp_T1_T2(Model):

    def pars(self):
        return ['S0', 'T1', 'T2','FA_corr']

    def init(self, sig):
        return [np.amax(sig), 1500, 80, 12]

    def bnds(self, sig):
        smax = np.amax(sig)
        return ([0, 50, 0, 0 ], [2*smax, 3000, 200, 25])

    def signal(self, TI_TE, S0, T1, T2, FA, TR=4.6):

        TI = TI_TE[0:28]
        TE = TI_TE[28:]
        
        cFA = np.cos(FA*np.pi/180)
        E = np.exp(-TR/T1)
        E = (1-E) / (1-cFA*E)
        Mss = S0 * E
        T1app = T1 * E

        y_T1 = np.abs(Mss - (S0+Mss) * np.exp(-TI/T1app))

        y_T2 = S0*np.exp(-TE/T2)

        return  np.concatenate((y_T1,y_T2),axis=0)