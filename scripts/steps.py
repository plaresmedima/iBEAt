import time
from pipelines import (
    segment, 
    fetch, 
    measure, 
    export, 
    rename, 
    mdr, 
    mapping, 
    harmonize, 
    align,
)



## HARMONIZATION



def rename_all_series(database):
    start_time = time.time()
    database.log("Renaming has started!")
    try:
        rename.all_series(database)
        rename.check(database)
        database.log("Renaming was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Renaming was NOT completed; error: " + str(e))  


def harmonize_all_series(database):
    start_time = time.time()
    database.log("Harmonizing series has started!")
    try:
        harmonize.all_series(database)
        database.log("Harmonizing series was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Harmonizing series was NOT completed; error: " + str(e)) 





## SEGMENTATION
        


def fetch_kidney_masks(database):
    start_time = time.time()
    database.log("Fetching kidney masks has started")
    try:
        fetch.kidney_masks(database)
        database.log("Fetching kidney masks was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Fetching kidney masks was NOT completed; error: "+str(e))

def fetch_dl_models(database):
    start_time = time.time()
    database.log("Fetching deep-learning models has started")
    try:
        fetch.dl_models(database)
        database.log("Fetching deep-learning models was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("Fetching deep-learning models was NOT completed; error: "+str(e))


def segment_kidneys(database, weights):
    start_time = time.time()
    database.log("Kidney segmentation has started")
    try:
        lk = database.series(SeriesDescription='LK')
        if len(lk) == 0:
            segment.kidneys(database, weights)
            database.log("Kidney segmentation was completed --- %s seconds ---" % (int(time.time() - start_time)))
        else:
            database.log('Kidney masks were already present - no automated kidney segmentation was performed.')
        database.save()
    except Exception as e:
        database.log("Kidney segmentation was NOT completed; error: "+str(e))
        raise RuntimeError('Critical step failed - ending pipeline.')


def segment_renal_sinus_fat(database):
    start_time = time.time()
    database.log("Sinus fat segmentation has started")
    try:
        segment.renal_sinus_fat(database)
        database.log("Sinus segmentation was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Sinus fat segmentation was NOT completed; error: "+str(e))


def segment_aorta_on_dce(database):
    start_time = time.time()
    database.log("Aorta segmentation has started")
    try:
        segment.aorta_on_dce(database)
        database.log("Aorta segmentation was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Aorta segmentation was NOT completed; error: "+str(e))


def compute_whole_kidney_canvas(database):
    start_time = time.time()
    database.log('Starting computing canvas')
    try:
        segment.compute_whole_kidney_canvas(database)
        database.log("Computing canvas was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("Computing canvas was NOT computed; error: "+str(e))


# Set up exports as individual steps
def export_segmentations(database):
    start_time = time.time()
    database.log("Export kidney segmentations has started")
    try:
        export.kidney_masks_as_dicom(database)
        export.kidney_masks_as_png(database)
        export.whole_kidney_canvas(database) # should this be part of the dicom masks export?
        export.aif_as_png(database)
        database.log("Export kidney segmentations was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("Export kidney segmentations was NOT completed; error: "+str(e))



## MODEL-DRIVEN MOTION CORRECTION



def mdreg_t1(database):
    start_time = time.time()
    database.log("Model-driven registration for T1 has started")
    try:
        mdr.T1(database)
        database.log("Model-driven registration for T1 was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Model-driven registration for T1 was NOT completed; error: "+str(e))

def mdreg_t2(database):
    start_time = time.time()
    database.log("Model-driven registration for T2 has started")
    try:
        mdr.T2(database)
        database.log("Model-driven registration for T2 was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Model-driven registration for T2 was NOT completed; error: "+str(e))


def mdreg_t2star(database):
    start_time = time.time()
    database.log("Model-driven registration for T2* has started")
    try:
        mdr.T2star(database)
        database.log("Model-driven registration for T2* was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Model-driven registration for T2* was NOT completed; error: "+str(e))


def mdreg_mt(database):
    start_time = time.time()
    database.log("Model-driven registration for MT has started")
    try:
        mdr.MT(database)
        database.log("Model-driven registration for MT was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Model-driven registration for MT was NOT completed; error: "+str(e))


def mdreg_dti(database):
    start_time = time.time()
    database.log("Model-driven registration for DTI has started")
    try:
        mdr.DTI(database)
        database.log("Model-driven registration for DTI was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Model-driven registration for DTI was NOT completed; error: "+str(e))


def mdreg_ivim(database):
    start_time = time.time()
    database.log("Model-driven registration for IVIM has started")
    try:
        mdr.IVIM(database)
        database.log("Model-driven registration for IVIM was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Model-driven registration for IVIM was NOT completed; error: "+str(e))


def mdreg_dce(database):
    start_time = time.time()
    database.log("Model-driven registration for DCE has started")
    try:
        mdr.DCE(database)
        database.log("Model-driven registration for DCE was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Model-driven registration for DCE was NOT completed; error: "+str(e))

def export_mdreg(database):
    start_time = time.time()
    database.log("Exporting MDR results has started")
    try:
        export.mdreg(database)
        database.log("Exporting MDR results was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("Exporting MDR results was NOT completed; error: "+str(e))

        

## MAPPING


  
def map_T1(database):
    start_time = time.time()
    print('Starting T1 mapping')
    database.log("T1 mapping has started")
    try:
        mapping.T1map(database)
        database.log("T1 mapping was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("T1 mapping was NOT completed; error: "+str(e))

def map_T2(database):
    start_time = time.time()
    print('Starting T2 mapping')
    database.log("T2 mapping has started")
    try:
        mapping.T2map(database)
        database.log("T2 mapping was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("T2 mapping was NOT completed; error: "+str(e))


def map_T1T2(database):
    start_time = time.time()
    print('Starting T1 and T2 mapping')
    database.log("T1 and T2 mapping has started")
    try:
        mapping.T1T2(database)
        database.log("T1 and T2 mapping was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("T1 and T2 mapping was NOT completed; error: "+str(e))

        
def map_T2star(database):
    start_time = time.time()
    print('Starting T2*')
    database.log("T2* mapping has started")
    try:
        mapping.T2starmap(database)
        database.log("T2* mapping was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("T2* mapping was NOT completed; error: "+str(e))


def map_MT(database):
    start_time = time.time()
    print('Starting MTR')
    database.log("MTR mapping has started")
    try:
        mapping.MT(database)
        database.log("MTR mapping was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("MTR mapping was NOT completed; error: "+str(e))


def map_DTI(database):
    start_time = time.time()
    print('Starting DTI')
    database.log("DTI mapping has started")
    try:
        mapping.DTI(database)
        database.log("DTI mapping was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("DTI-FA & ADC mapping was NOT completed; error: "+str(e))


def map_IVIM(database):
    start_time = time.time()
    print('Starting IVIM')
    database.log("IVIM mapping has started")
    try:
        mapping.ivim(database)
        database.log("IVIM mapping was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("IVIM mapping was NOT completed; error: "+str(e))


def map_DCE(database):
    start_time = time.time()
    print('Starting DCE_MAX')
    database.log("DCE mapping has started")
    try:
        mapping.DCE(database)
        database.log("DCE mapping was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("DCE mapping was NOT completed; error: "+str(e))

def export_mapping(database):
    start_time = time.time()
    database.log("Exporting mapping results has started")
    try:
        export.mapping(database)
        database.log("Exporting maping results was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("Exporting mapping results was NOT completed; error: "+str(e))



## ALIGNMENT
        

def align_T1(database):
    start_time = time.time()
    print('Starting T1 alignment')
    database.log("T1 alignment has started")
    try:
        align.t1(database)
        database.log("T1 alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("T1 alignment was NOT completed; error: "+str(e))

def align_T2(database):
    start_time = time.time()
    print('Starting T2 alignment')
    database.log("T2 alignment has started")
    try:
        align.t2(database)
        database.log("T2 alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("T2 alignment was NOT completed; error: "+str(e))

def align_T2star(database):
    start_time = time.time()
    print('Starting T2* alignment')
    database.log("T2* alignment has started")
    try:
        align.t2star(database)
        database.log("T2* alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("T2* alignment was NOT completed; error: "+str(e))

def align_MT(database):
    start_time = time.time()
    print('Starting MT alignment')
    database.log("MT alignment has started")
    try:
        align.mt(database)
        database.log("MT alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("MT alignment was NOT completed; error: "+str(e))

def align_IVIM(database):
    start_time = time.time()
    print('Starting IVIM alignment')
    database.log("IVIM alignment has started")
    try:
        align.ivim(database)
        database.log("IVIM alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("IVIM alignment was NOT completed; error: "+str(e))

def align_DTI(database):
    start_time = time.time()
    print('Starting DTI alignment')
    database.log("DTI alignment has started")
    try:
        align.dti(database)
        database.log("DTI alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("DTI alignment was NOT completed; error: "+str(e))

def align_DCE(database):
    start_time = time.time()
    print('Starting DCE alignment')
    database.log("DCE alignment has started")
    try:
        align.dce(database)
        database.log("DCE alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("DCE alignment was NOT completed; error: "+str(e))

def align_ASL(database):
    start_time = time.time()
    print('Starting ASL alignment')
    database.log("ASL alignment has started")
    try:
        align.asl(database)
        database.log("ASL alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("ASL alignment was NOT completed; error: "+str(e))

def export_alignment(database):
    start_time = time.time()
    print('Exporting alignment')
    database.log("Export alignment has started")
    try:
        export.alignment(database)
        database.log("Export alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("Export alignment was NOT completed; error: "+str(e))



## MEASURE
        


def measure_kidney_volumetrics(database):
    start_time = time.time()
    database.log("Kidney volumetrics has started")
    try:
        measure.kidney_volumetrics(database)
        database.log("Kidney volumetrics was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("Kidney volumetrics was NOT completed; error: "+str(e))

def measure_sinus_fat_volumetrics(database):
    start_time = time.time()
    database.log("Sinus fat volumetrics has started")
    try:
        measure.sinus_fat_volumetrics(database)
        database.log("Sinus fat volumetrics was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("Sinus fat volumetrics was NOT completed; error: "+str(e))

def measure_t1_maps(database):
    start_time = time.time()
    database.log("T1 map measurement has started")
    try:
        measure.t1_maps(database)
        database.log("T1 map measurement was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("T1 map measurement was NOT completed; error: "+str(e))  

def measure_t2_maps(database):
    start_time = time.time()
    database.log("T2 map measurement has started")
    try:
        measure.t2_maps(database)
        database.log("T2 map measurement was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("T2 map measurement was NOT completed; error: "+str(e))     

def measure_t2star_maps(database):
    start_time = time.time()
    database.log("T2* map measurement has started")
    try:
        measure.t2star_maps(database)
        database.log("T2* map measurement was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("T2* map measurement was NOT completed; error: "+str(e))   

def measure_mt_maps(database):
    start_time = time.time()
    database.log("MT map measurement has started")
    try:
        measure.mt_maps(database)
        database.log("MT map measurement was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("MT map measurement was NOT completed; error: "+str(e)) 

def measure_ivim_maps(database):
    start_time = time.time()
    database.log("IVIM map measurement has started")
    try:
        measure.ivim_maps(database)
        database.log("IVIM map measurement was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("IVIM map measurement was NOT completed; error: "+str(e)) 

def measure_dti_maps(database):
    start_time = time.time()
    database.log("DTI map measurement has started")
    try:
        measure.dti_maps(database)
        database.log("DTI map measurement was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("DTI map measurement was NOT completed; error: "+str(e))  

def measure_dce_maps(database):
    start_time = time.time()
    database.log("DCE map measurement has started")
    try:
        measure.dce_maps(database)
        database.log("DCE map measurement was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("DCE map measurement was NOT completed; error: "+str(e))  

def measure_asl_maps(database):
    start_time = time.time()
    database.log("ASL perfusion measurement has started")
    try:
        measure.asl_maps(database)
        database.log("ASL perfusion measurement was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        database.log("ASL perfusion measurement was NOT completed; error: "+str(e))  

        

