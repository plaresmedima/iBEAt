import time
from pipelines import (
    fetch_Drive_mask,
    mdr, 
    mapping, 
    harmonize, 
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
        
## MEASURE
    
## ROI ANALYSIS
