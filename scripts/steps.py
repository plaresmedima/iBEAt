import time
from pipelines import segment, fetch, measure, export, rename


def rename_all_series(folder):
    start_time = time.time()
    folder.log("Renaming has started!")
    try:
        rename.all_series(folder)
        rename.check(folder)
        folder.log("Renaming was completed --- %s seconds ---" % (int(time.time() - start_time)))
        folder.save()
    except Exception as e:
        folder.log("Renaming was NOT completed; error: " + str(e))  


def fetch_kidney_masks(folder):
    start_time = time.time()
    folder.log("Fetching kidney masks has started")
    try:
        fetch.kidney_masks(folder)
        folder.log("Fetching kidney masks was completed --- %s seconds ---" % (int(time.time() - start_time)))
        folder.save()
    except Exception as e:
        folder.log("Fetching kidney masks was NOT completed; error: "+str(e))


def segment_kidneys(folder, weights):
    start_time = time.time()
    folder.log("Kidney segmentation has started")
    try:
        lk = folder.series(SeriesDescription='LK')
        if len(lk) == 0:
            segment.kidneys(folder, weights)
            folder.log("Kidney segmentation was completed --- %s seconds ---" % (int(time.time() - start_time)))
        else:
            folder.log('Kidney masks were already present - no automated kidney segmentation was performed.')
        folder.save()
    except Exception as e:
        # Is this the right way to handle? other steps depend on this one.
        folder.log("Kidney segmentation was NOT completed; error: "+str(e))


def segment_renal_sinus_fat(folder):
    start_time = time.time()
    folder.log("Sinus fat segmentation has started")
    try:
        segment.renal_sinus_fat(folder)
        folder.log("Sinus segmentation was completed --- %s seconds ---" % (int(time.time() - start_time)))
        folder.save()
    except Exception as e:
        folder.log("Sinus fat segmentation was NOT completed; error: "+str(e))


def compute_whole_kidney_canvas(folder):
    start_time = time.time()
    folder.log('Starting computing canvas')
    try:
        segment.compute_whole_kidney_canvas(folder)
        folder.log("Computing canvas was completed --- %s seconds ---" % (int(time.time() - start_time)))
        folder.save()
    except Exception as e: 
        folder.log("Computing canvas was NOT computed; error: "+str(e))


def measure_kidney_volumetrics(folder):
    start_time = time.time()
    folder.log("Kidney volumetrics has started")
    try:
        measure.kidney_volumetrics(folder)
        folder.log("Kidney volumetrics was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        folder.log("Kidney volumetrics was NOT completed; error: "+str(e))


def measure_sinus_fat_volumetrics(folder):
    start_time = time.time()
    folder.log("Sinus fat volumetrics has started")
    try:
        measure.sinus_fat_volumetrics(folder)
        folder.log("Sinus fat volumetrics was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        folder.log("Sinus fat volumetrics was NOT completed; error: "+str(e))


def export_kidney_segmentations(folder):
    start_time = time.time()
    folder.log("Export kidney segmentations has started")
    try:
        export.kidney_masks_as_dicom(folder)
        export.kidney_masks_as_png(folder)
        export.whole_kidney_canvas(folder) # should this be part of the dicom masks export?
        folder.log("Export kidney segmentations was completed --- %s seconds ---" % (int(time.time() - start_time)))
    except Exception as e:
        folder.log("Export kidney segmentations was NOT completed; error: "+str(e))



