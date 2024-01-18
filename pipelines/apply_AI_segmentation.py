import pipelines.segment as seg
import pandas as pd

def main(folder):
    folder.log("AI segmentation has started!")
    weights = 'UNETR_kidneys_v1.pth'

    seg.segment_kidneys(folder, weights)
    #seg.compute_renal_sinus_fat(folder)
    seg.export_masks(folder)
    seg.compute_whole_kidney_canvas(folder)
    seg.export_whole_kidney_canvas(folder)





