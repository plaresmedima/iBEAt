import numpy as np

def freeRecoveryMagnetization(M_init, t, M_eq, T1):
    """ Free Longitudinal Recovery

    M_init: initial magnetization state
    M_eq, T1: tissue parameters
    t: recovery time
    """
    E = np.exp(-t/T1)
    Mt = M_init * E  + M_eq * (1-E)
    return Mt

def pulse(M_init, FA):
    """ Alpha Pulse

    M_init: initial magnetization state
    FA: used flip-angle in radians
    """
    return np.cos(FA) * M_init


def FLASHreadout(M_init, M_eq, T1, FA, TR, N):
    """ FLASH readout

    M_init: initial magnetization state
    M_eq, T1: tissue parameters
    FA : Flip angle in radians
    TR : time between FA pulses in ms
    N : number of readout pulses
    """
    M_current = M_init
    
    for i in range(int(N)):
        M_current = pulse(M_current, FA)
        M_current = freeRecoveryMagnetization(M_current, TR, M_eq, T1)
    return M_current

def FLASHreadout_CatModule(M_init, M_eq, T1, FA, TR, N):
    """ Catalization Module

    M_init: initial magnetization state
    M_eq, T1: tissue parameters
    FA : list of flip angle in radians
    TR : time between FA pulses in ms
    N : number of readout pulses
    """
    M_current = M_init
    
    for i in range(int(N)):
        M_current = pulse(M_current, FA[i])
        M_current = freeRecoveryMagnetization(M_current, TR, M_eq, T1)
    return M_current

def freeDecayMagnetization(M_init, t, T2):
    """ T2 prep

    M_init: initial magnetization state
    t: T2 prep duration
    T2: tissue parameters
    """
    E = np.exp(-t/T2)
    return E * M_init


def signalSequenceT2prepOneShot(M_init, M_eq, T1, T2, Tprep, Tspoil, FA, TR, N):
    """ One shot of a T2-prep sequence

    M_eq, T1, T2 : tissue parameters
    Tprep : preparation time in ms (between 0 and 120ms)
    Tspoil: time for spoiling after prep pulse (1ms)
    FA : Flip angle in radians (12 degrees)
    TR : time between FA pulses in ms (about 4.6ms)
    N : number of readout pulses (72)
    Trec : recovery time after readout (2 * 463ms)
    k-space center: 48 (72/2 + 12 due to partial Fourier) 
    """

    Mcurrent = freeDecayMagnetization(M_init, Tprep, T2) # prep pulse
    Mcurrent = freeRecoveryMagnetization(Mcurrent, Tspoil, M_eq, T1) # during spoiling
    Mcurrent = FLASHreadout(Mcurrent, M_eq, T1, FA, TR, N/2-12) # readout

    return Mcurrent 


def signalSequenceT2prep(Tprep, M_eq, T2, T1 , Tspoil, FA, FA_Eff, TR, N, Trec):
    """ All shots of a T2 prep sequence.

    M_eq, T1, T2 : tissue parameters
    Tprep : preparation time in ms (between 0 and 120ms)
    Tspoil: time for spoiling after prep pulse (1ms)
    FA : Flip angle in radians (12 degrees)
    TR : time between FA pulses in ms (about 4.6ms)
    N : number of readout pulses (72)
    Trec : recovery time after readout (2 * 463ms)
    k-space center: 48 (72/2 + 12 due to partial Fourier) 
    """
    FA = FA * FA_Eff
    M_result = np.zeros(11) 
    M_current = M_eq
    M_current = signalSequenceT2prepOneShot(M_current, M_eq, T1, T2, Tprep[0], Tspoil, FA, TR, N) # signal at first Tprep
    M_result[0]=M_current

    for t in range(1, np.size(Tprep)):
        
        M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2+12) # readout
        M_current = freeRecoveryMagnetization(M_current, 120-Tprep[t-1], M_eq, T1) # recovery after readout
        M_current = freeRecoveryMagnetization(M_current, Trec, M_eq, T1) # recovery after readout
        M_current = signalSequenceT2prepOneShot(M_current, M_eq, T1, T2, Tprep[t], Tspoil, FA, TR, N) # signal at first Tprep
        M_result[t] =M_current
    return M_result


def signalSequenceT1_FLASH(M_eq, T1, TI, FA, FA_Eff, TR, N, FA_Cat):
    """ All shots of a T2 prep sequence.

    M_eq, T1: tissue parameters
    TI : list of inversion times (between 100 and 7700ms)
    FA : Flip angle in radians (12 degrees)
    TR : time between FA pulses in ms (about 4.6ms)
    N : number of readout pulses (66)

    k-space center: 13 (66/2 - 20 due to 5/8 partial Fourier) 
    """
    FA = FA*FA_Eff
    FA_Cat = np.array(FA_Cat)*FA_Eff
    M_result = np.zeros(28) 
    ####### 1st SET: 16 TI's ########


#SLICE 1
    #M_current = M_eq*(-1)                                          # 180 in z-
    #M_current = freeRecoveryMagnetization(M_current, 9120, M_eq, T1)
    #M_current = M_current*(-1)
    #M_current = freeRecoveryMagnetization(M_current, 5060, M_eq, T1)
    #M_current = M_current*(-1)
    #M_current = freeRecoveryMagnetization(M_current, 10200, M_eq, T1)

#SLICE 2
    #M_current = M_current*(-1)
    #M_current = freeRecoveryMagnetization(M_current, 9120, M_eq, T1)
    #M_current = M_current*(-1)
    #M_current = freeRecoveryMagnetization(M_current, 5060, M_eq, T1)
    #M_current = M_current*(-1)
    #M_current = freeRecoveryMagnetization(M_current, 10200, M_eq, T1)
    #M_current = M_current*(-1)

    # Inversion                                              
    M_current = M_eq*(-1)                                          # 180 in z-
    M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1)     # half of inversion pulse delay
    
    # TIfill  
    M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1)     # TIfill delay
    
    # Catalization Module (25ms)
    M_current = FLASHreadout_CatModule(M_current, M_eq, T1, FA_Cat, TR, 5)
    #M_current = freeRecoveryMagnetization(M_current, 2, M_eq, T1)   # 2ms to make 25ms
    M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2-20)   
    #Acquisition (66 lines: 13 + 53) 
    M_result[0] =M_current                                          # save result (13 lines)
    M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2+20)       # rest of the readout (53 lines)

    for t in range(1, np.size(TI)-12):

        M_current = freeRecoveryMagnetization(M_current, 80, M_eq, T1)  # fixed 80ms delay Siemens
        M_current = freeRecoveryMagnetization(M_current, 80, M_eq, T1)  # fixed 80ms delay Siemens
        M_current = freeRecoveryMagnetization(M_current, 13, M_eq, T1)  # inversion pulse duration
        M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1) # TI fill delay
        M_current = FLASHreadout_CatModule(M_current, M_eq, T1, FA_Cat, TR, 5) 
       # M_current = freeRecoveryMagnetization(M_current, 2, M_eq, T1)   # 2ms to make 25ms
        M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2-20)   # k-space center
        M_result[t] =M_current                                      # save result
        M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2+20)   # rest of the readout
    
    ####### Recovery ########
    #Beat 1  2
    M_current = freeRecoveryMagnetization(M_current, 507*2, M_eq, T1)   # recovery time (2 beats: ms)
    #M_current = freeRecoveryMagnetization(M_current, 80+80+13, M_eq, T1)  
    #M_current = FLASHreadout_CatModule(M_current, M_eq, T1, FA_Cat, TR, 5)
    #M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N)
    #M_current = freeRecoveryMagnetization(M_current, 507-80-80-13-23, M_eq, T1)

    #M_current = freeRecoveryMagnetization(M_current, 80+80+13, M_eq, T1)  
    #M_current = FLASHreadout_CatModule(M_current, M_eq, T1, FA_Cat, TR, 5)
    #M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N)
    #M_current = freeRecoveryMagnetization(M_current, 507-80-80-13-23, M_eq, T1)


    
    ####### 2nd SET: 8 TI's ########
    M_current = freeRecoveryMagnetization(M_current, 80, M_eq, T1)
    M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1)
    M_current = M_current*(-1)
    M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1)
    M_current = freeRecoveryMagnetization(M_current, 80, M_eq, T1)
    M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1)
    M_current = FLASHreadout_CatModule(M_current, M_eq, T1, FA_Cat, TR, 5)           
    #M_current = freeRecoveryMagnetization(M_current, 2, M_eq, T1)       
    M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2-20)       
    M_result[16] = M_current                                          
    M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2+20)

    for t in range(1, np.size(TI)-20):

        M_current = freeRecoveryMagnetization(M_current, 80, M_eq, T1)
        M_current = freeRecoveryMagnetization(M_current, 13, M_eq, T1)
        M_current = freeRecoveryMagnetization(M_current, 80, M_eq, T1)
        M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1)
        M_current = FLASHreadout_CatModule(M_current, M_eq, T1, FA_Cat, TR, 5)             # 5 flash readouts
        #M_current = freeRecoveryMagnetization(M_current, 2, M_eq, T1)   # 2m spoiler gradient in z
        M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2-20)   # k-space center
        M_result[16+t] = M_current                                      # save result
        M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2+20)
     
    ####### Recovery 2 ########
    #Beat 1 & 2
    M_current = freeRecoveryMagnetization(M_current, 507*2, M_eq, T1)   # recovery time (2 beats: ms)
    #M_current = freeRecoveryMagnetization(M_current, 80+80+13, M_eq, T1)  
    #M_current = FLASHreadout_CatModule(M_current, M_eq, T1, FA_Cat, TR, 5)
    #M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N)
    #M_current = freeRecoveryMagnetization(M_current, 507-80-80-13-23, M_eq, T1)

    #M_current = freeRecoveryMagnetization(M_current, 80+80+13, M_eq, T1)  
    #M_current = FLASHreadout_CatModule(M_current, M_eq, T1, FA_Cat, TR, 5)
    #M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N)
    #M_current = freeRecoveryMagnetization(M_current, 507-80-80-13-23, M_eq, T1)
    
    ####### 3rd SET: 4 TI's ########
    M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1)
    M_current = M_current*(-1)
    M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1)
    M_current = freeRecoveryMagnetization(M_current, 80, M_eq, T1)
    M_current = freeRecoveryMagnetization(M_current, 80, M_eq, T1)
    M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1)
    M_current = FLASHreadout_CatModule(M_current, M_eq, T1, FA_Cat, TR, 5)             # 5 flash readouts
    #M_current = freeRecoveryMagnetization(M_current, 2, M_eq, T1)       # 2m spoiler gradient in z
    M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2-20)       # k-space center
    M_result[24] = M_current                                          # save result
    M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2+20)

    for t in range(1, np.size(TI)-24):

        M_current = freeRecoveryMagnetization(M_current, 13, M_eq, T1)
        M_current = freeRecoveryMagnetization(M_current, 80, M_eq, T1)
        M_current = freeRecoveryMagnetization(M_current, 80, M_eq, T1)
        M_current = freeRecoveryMagnetization(M_current, 6.5, M_eq, T1)
        M_current = FLASHreadout_CatModule(M_current, M_eq, T1, FA_Cat, TR, 5)             # 5 flash readouts
        #M_current = freeRecoveryMagnetization(M_current, 2, M_eq, T1)       # 2m spoiler gradient in z
        M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2-20)       # k-space center
        M_result[24+t] = M_current                                          # save result
        M_current = FLASHreadout(M_current, M_eq, T1, FA, TR, N/2+20)

    M_result = np.abs(M_result)

    return M_result