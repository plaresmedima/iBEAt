
import numpy as np
from utilities.fit import Model
from models import relaxometry


class MonoExpOffset(Model):

    def pars(self):
        return ['S0', 'T2', 'T2offset']

    def init(self, sig):
        smax = max([np.amax(sig),1])
        return [smax, 80, 0.01]

    def bnds(self, sig):
        smax = np.amax(sig)
        if smax<=0:
            return [0,10,0], [np.inf, 500, 1]
        else:
            return [0,10,0], [5*smax, 500, 1]

    def signal(self,TE,S0,T2,C):
        if T2==0:
            return TE*0
        return np.sqrt((S0*np.exp(-TE/T2))**2 + ((S0*C)**2)) 
       

class MonoExp(Model):
    
    def pars(self):
        return ['S0', 'T2']

    def init(self, sig):
        smax = max([np.amax(sig),1])
        return [smax, 80.0]

    def bnds(self, sig):
        smax = np.amax(sig)
        if smax<=0:
            return [0.0,10.0], [np.inf, 500.0]
        else:
            return [0.0,10.0], [5*smax, 500.0]

    def signal(self, TE,S0,T2):
        if T2==0:
            return TE*0
        return S0*np.exp(-TE/T2)   


# class Bloch(Model):

#     def pars(self):
#         return ['S0', 'T2', 'T2FAcorr']
    
#     def init(self, sig):
#         return [np.amax(sig), 80, 12]
    
#     def bnds(self, sig):
#         smax = np.amax(sig)
#         return ([0, 0, 0], [2*smax, 200, 25])
    
#     def signal(self, Tprep, S0, T2, FAcorr, 
#             T1 = 1400, # Get from data
#             # Defaults from Siemens iBEAt protocol
#             Tspoil = 1,# Spoil time in ms
#             N_T2 = 72,# Number of k-space lines (hardcoded from Siemens protocol)
#             Trec = 463*2,# Recovery time in ms (hardcoded from Siemens protocol)
#             TR = 4.6,# TR in ms (hardcoded from Siemens protocol)
#             FA = 12, # Flip angle in degrees (hardcoded from Siemens protocol) converted to radians
#         ):
#         FA_eff = FAcorr/FA
#         FA = FA/180*np.pi
#         return relaxometry.signalSequenceT2prep(Tprep, S0, T2, T1, Tspoil, FA, FA_eff, TR, N_T2, Trec)
    
class Bloch(Model):

    def pars(self):
        return ['S0', 'T2', 'T2FAcorr']
    
    def init(self, sig):
        return [np.amax(sig), 80, 12]
    
    def bnds(self, sig):
        smax = np.amax(sig)
        return ([0, 0, 0], [2*smax, 200, 25])
    
    def signal(self, xdata, S0, T2, FAcorr,
            # Defaults from Siemens iBEAt protocol
            Tspoil = 1,# Spoil time in ms
            N_T2 = 72,# Number of k-space lines (hardcoded from Siemens protocol)
            Trec = 463*2,# Recovery time in ms (hardcoded from Siemens protocol)
            TR = 4.6,# TR in ms (hardcoded from Siemens protocol)
            FA = 12, # Flip angle in degrees (hardcoded from Siemens protocol) converted to radians
        ):
        T1 = xdata[0]
        Tprep = xdata[1:]
        FA_eff = FAcorr/FA
        FA = FA/180*np.pi
        return relaxometry.signalSequenceT2prep(Tprep, S0, T2, T1, Tspoil, FA, FA_eff, TR, N_T2, Trec)