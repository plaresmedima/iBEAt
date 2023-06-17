""" 
@author: Joao Periquito 
iBEAt SEGMENTATION Scrpit
2023
Find iBEAt fat and opposed phase scan and execute kmeans segmentation
"""
import datetime
import time

import pipelines.segment as seg


def main(folder, path, filename_log):

    try:
        start_time = time.time()
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Sequential kmeans has started")
        file.close()
        
        seg.compute_whole_kidney_canvas(folder)
        folder.save()

        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Sequential kmeans was completed --- %s seconds ---" % (int(time.time() - start_time))) 
        file.close() 
    
    except Exception as e: 

        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Sequential kmeans was NOT completed; error: "+str(e)) 
        file.close()

    try:
        start_time = time.time()
        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Export segmentation results has started")
        file.close()
        
        seg.export_whole_kidney_canvas(folder, path)

        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Export segmentation results was completed --- %s seconds ---" % (int(time.time() - start_time))) 
        file.close() 
    
    except Exception as e: 

        file = open(filename_log, 'a')
        file.write("\n"+str(datetime.datetime.now())[0:19] + ": Export segmentation results was NOT completed; error: "+str(e)) 
        file.close()

    
