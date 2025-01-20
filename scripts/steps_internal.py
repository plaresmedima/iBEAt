import time
from pipelines import (
    fetch_Drive_mask,
    mdr, 
    mapping, 
    harmonize,
    segment,
    export, 
)


## HARMONIZATION

def harmonize_t1_t2(database):
    start_time = time.time()
    database.log("Harmonizing T1 and T2 series (merged) has started!")
    try:
        harmonize.t1_t2_merge(database)
        database.log("Harmonizing T1 and T2 series (merged) was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Harmonizing T1 and T2 series (merged) was NOT completed; error: " + str(e)) 
        database.restore()


## SEGMENTATION

def fetch_kidney_masks(database):
    start_time = time.time()
    database.log("Fetching kidney masks has started")
    try:
        fetch_Drive_mask.kidney_masks(database)
        database.log("Fetching kidney masks was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Fetching kidney masks was NOT completed; error: "+str(e))
        database.restore()

def compute_whole_kidney_canvas(database):
    start_time = time.time()
    database.log('Starting computing canvas')
    try:
        segment.compute_whole_kidney_canvas(database)
        database.log("Computing canvas was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("Computing canvas was NOT computed; error: "+str(e))
        database.restore()

def export_sinus_fat_segmentations(database):
    start_time = time.time()
    database.log("Export kidney segmentations has started")
    try:
        export.kidney_masks_as_png(database)
    except Exception as e:
        database.log("Export kidney segmentations was NOT completed; error: "+str(e))
    
    try:
        export.kidney_masks_sinus_fat_dixon_as_dicom(database)
    except Exception as e:
        database.log("Export kidney masks sinus fat masks and dixon was NOT completed; error: "+str(e))

    database.log("Export kidney segmentations was completed --- %s seconds ---" % (int(time.time() - start_time)))


## MODEL-DRIVEN MOTION CORRECTION

def mdreg_t1_t2(database):
    start_time = time.time()
    database.log("Model-driven registration for T1 T2 has started")
    try:
        mdr.T1_T2(database)
        database.log("Model-driven registration for T1 T2 was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e:
        database.log("Model-driven registration for T1 T2 was NOT completed; error: "+str(e))
        database.restore()
    
## MAPPING


def map_T1_from_T1_T2_mdr(database):
    start_time = time.time()
    print('Starting T1 mapping from T1_T2 MDR')
    database.log("T1 mapping from T1_T2 MDR has started")
    try:
        mapping.T1_from_T1_T2_mdr(database)
        database.log("T1 mapping from T1_T2 MDR was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("T1 mapping from T1_T2 MDR was NOT completed; error: "+str(e))
        database.restore()

def map_T2_from_T1_T2_mdr(database):
    start_time = time.time()
    print('Starting T2 mapping from T1_T2 MDR')
    database.log("T2 mapping from T1_T2 MDR has started")
    try:
        mapping.T2_from_T1_T2_mdr(database)
        database.log("T2 mapping from T1_T2 MDR was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("T2 mapping from T1_T2 MDR was NOT completed; error: "+str(e))
        database.restore()


## ALIGNMENT

def export_alignment_t2s_project(database):
    start_time = time.time()
    print('Exporting alignment')
    database.log("Export alignment has started")
    try:
        export.alignment_t2s_project(database)
        database.log("Export alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("Export alignment was NOT completed; error: "+str(e))
        database.restore()

## MEASURE
    
## ROI ANALYSIS

## DATA EXPORT

def export_dicom_t2s_project(database):
    start_time = time.time()
    print('Exporting alignment')
    database.log("Export alignment has started")
    try:
        export.project_t2s_prepare_dicom(database)
        database.log("Export alignment was completed --- %s seconds ---" % (int(time.time() - start_time)))
        database.save()
    except Exception as e: 
        database.log("Export alignment was NOT completed; error: "+str(e))
        database.restore()

