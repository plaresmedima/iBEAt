import pandas as pd
from models.dce.pbpk_aorta import OneScan as AortaModel
from models.dce.pbpk_kidney_nephron_short import OneScan as KidneyModel


def _fit_aorta(TR, FA, field_strength, time, signal, weight, R10):
    print('Fitting aorta...')
    aorta = AortaModel(
        weight = weight,
        dose = 0.2, #mL/kg  
        conc = 0.5, # mmol/ml 
        rate = 1, # ml/sec
        TR = TR, #sec
        FA = FA, #deg
        field_strength = field_strength, 
        tdce = time,   
        Sdce = signal,
        R1molli = R10,
        callback = False,
        ptol = 1e-3,
    )
    aorta.estimate_p()
    print('Aorta goodness of fit: ', aorta.goodness())
    aorta.fit_p()
    print('Aorta goodness of fit: ', aorta.goodness())
    return aorta


def _fit_kidney(TR, FA, field_strength, Hct, time, signal, R1, kidney_volume, aorta):
    print('Fitting kidney...')
    kidney = KidneyModel(
        TR = TR, #sec
        FA = FA, #deg
        dt = aorta.dt,
        tmax = aorta.tmax,
        field_strength = field_strength,       
        Hct = Hct,   
        J_aorta = aorta.p.value.CO*aorta.cb/1000,
        tdce = time,
        BAT = aorta.p.value.BAT,
        Sdce = signal,
        R1molli = R1,
        ptol = 1e-6,
        kidney_volume = kidney_volume, 
        CO = aorta.p.value.CO,
    )
    kidney.estimate_p()
    print('Goodness of fit (%): ', kidney.goodness())
    kidney.fit_p()
    print('Goodness of fit (%): ', kidney.goodness())
    return kidney


def fit(
        TR,  #sec
        FA, #degrees
        field_strength, #T
        dose, #mL/kg
        conc, # mmol/ml 
        rate, # ml/sec
        weight, #kg
        Hct, 
        time, #sec
        signal_aorta,  #au
        R1_aorta, #1/sec
        signal_kidney, #au
        R1_kidney, #1/sec
        volume_kidney, #mL
        path = None, # full path for exporting diagnostics
        name = 'subject',
        ROI = 'kidney',
        export_aorta = True,
        ):
    
    # Perform the fits
    aorta = _fit_aorta(TR, FA, field_strength, time, signal_aorta, weight, R1_aorta)
    kidney = _fit_kidney(TR, FA, field_strength, Hct, time, signal_kidney, R1_kidney, volume_kidney, aorta)

    # Create export parameters
    kidney_pars = kidney.export_p()
    kidney_pars['structure'] = ROI
    if export_aorta:
        aorta_pars = aorta.export_p()
        aorta_pars['structure'] = 'Aorta'
        kidney_pars = pd.concat([aorta_pars, kidney_pars])
    kidney_pars = kidney_pars[['structure','name','value','unit']]

    # Save diagnostics
    if path is not None:
        kidney.plot_fit(save=True, show=False, path=path, prefix=name+'_'+ROI)
        if export_aorta:
            aorta.plot_fit(save=True, show=False, path=path, prefix=name+'_aorta')

    return kidney_pars # dataframe



