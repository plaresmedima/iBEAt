import pipelines.segment as seg

def main(folder):
    folder.log("AI segmentation has started!")
    weights = 'UNETR_kidneys_v1.pth'
    seg.segment_kidneys(folder, weights)