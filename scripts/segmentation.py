""" 
@author: Joao Periquito 
iBEAt SEGMENTATION Scrpit
2023
Find iBEAt fat and opposed phase scan and execute kmeans segmentation
"""
import time
import pipelines.segment as seg


def main(folder, path):

    start_time = time.time()
    folder.log("Segmentation has started!")

    try:
        seg.compute_whole_kidney_canvas(folder)
        folder.save()
    except Exception as e: 
        folder.log("Sequential kmeans was NOT completed; error: "+str(e))

    try:
        seg.export_whole_kidney_canvas(folder, path)
    except Exception as e: 
        folder.log("Export segmentation results was NOT completed; error: "+str(e))

    folder.log("Segmentation was completed --- %s seconds ---" % (int(time.time() - start_time)))


    
