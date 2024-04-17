import numpy as np
from utilities.fit import Model


class BiExp2D(Model):

    def pars(self):
        return ["S0", "T2star", "f_fat", "phase", "freq_shift"]

    def init(self, sig):
        return [np.amax(sig), 60, 0, 0, -np.inf]

    def bnds(self, sig):
        return [0,0,0,0,0], [5*np.amax(sig), 200, 1, 2*np.pi, np.inf]

    def signal(self, TE:np.ndarray, S0, T2star, f_fat, phase, freq_shift):
        T2star_fat = 9.3    # Le Ster, C. et al. doi:10.1002/jmri.25205
        chem_shift = 3.4e-6
        larmor_freq = 128e+6
        wat_phase = phase + TE*larmor_freq*(1+freq_shift)
        fat_phase = phase + TE*larmor_freq*(1+freq_shift)*(1+chem_shift)
        wat_magn = S0*(1-f_fat)*np.exp(-TE/T2star)
        fat_magn = S0*f_fat*np.exp(-TE/T2star_fat)
        wat_x = wat_magn * np.cos(wat_phase)
        wat_y = wat_magn * np.sin(wat_phase)
        fat_x = fat_magn * np.cos(fat_phase)
        fat_y = fat_magn * np.sin(fat_phase)
        return np.concatenate((wat_x + fat_x, wat_y + fat_y))


class BiExp(Model):

    def pars(self):
        return ["S0", "T2star", "f_fat"]

    def init(self, sig):
        return [np.amax(sig), 60, 0]

    def bnds(self, sig):
        return [0,0,0], [5*np.amax(sig), 200, 1]

    def signal(self, TE:np.ndarray, S0, T2star, f_fat):
        T2star_fat = 9.3    # Le Ster, C. et al. doi:10.1002/jmri.25205
        phase = np.zeros(TE.size)
        for m in range(TE.size):
            if m % 2 == 0:
                phase[m]=-1
            else:
                phase[m]=1
        wat = np.exp(-TE/T2star)
        fat = np.exp(-TE/T2star_fat)
        return S0*((1-f_fat)*wat + phase*f_fat*fat)


class MonoExp(Model):
    
    def pars(self):
        return ['S0', 'T2star']

    def init(self, sig):
        return [np.amax(sig), 50]

    def bnds(self, sig):
        return [0,5], [5*np.amax(sig), 200]

    def signal(self, TE, S0, T2star):
        if T2star==0:
            return TE*0
        return S0*np.exp(-TE/T2star)