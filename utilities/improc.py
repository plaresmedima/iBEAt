"""General image processing utilities"""

import numpy as np


def region_grow_range(img:np.ndarray, seed:list, min:float, max:float):
    # Taken from wezel.widgets.canvas.utils
    # Grows a region from a list of seed points
    # Selecting all neighbours with a value in the range [min, max]
    selected = np.zeros(img.shape, dtype=np.bool8)
    checked = np.zeros(img.shape, dtype=np.bool8)
    width, height = img.shape
    neighbours = [ 
        [0, -1], [1, 0], 
        [0, 1], [-1, 0],
    ]
    while seed != []:
        p = seed.pop()
        if min <= img[p[0], p[1]] <= max:
            selected[p[0], p[1]] = True
        for next in neighbours:
            x = p[0] + next[0]
            y = p[1] + next[1]
            if x < 0 or y < 0 or x >= width or y >= height:
                continue
            if not checked[x,y]:
                checked[x,y] = True
                if min <= img[x,y] <= max:
                    seed.append([x,y])
    return selected


def region_grow_thresh(img:np.ndarray, seed:list, threshold:float):
    # Adapted from wezel.widgets.canvas.utils
    # Grows a region from a list of seed points
    # Selecting all neighbours with a value with a threshold of the current seed
    selected = np.zeros(img.shape, dtype=np.bool8)
    #checked = np.zeros(img.shape, dtype=np.bool8)
    width, height = img.shape
    neighbours = [ 
        [1,1], [1,-1], [-1,-1], [-1,1],
        [0, -1], [1, 0], [0, 1], [-1, 0],
    ]
    while seed != []:
        p = seed.pop()
        selected[p[0], p[1]] = True
        for next in neighbours:
            x = p[0] + next[0]
            y = p[1] + next[1]
            if x < 0 or y < 0 or x >= width or y >= height:
                continue
            # if not checked[x,y]:
            #     checked[x,y] = True
            if not selected[x,y]:
                if np.abs(img[x,y]-img[p[0],p[1]]) < threshold:
                    selected[x,y] = True
                    seed.append([x,y])
    return selected