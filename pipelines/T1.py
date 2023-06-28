from dbdicom.pipelines import input_series
from pipelines.mdr import MDRegT1
from pipelines.mapping import T1
from dbdicom.wrappers import vreg

export_study = "T1"
magn_series = "T1map_kidneys_cor-oblique_mbh_magnitude"


def mdr(database):
    series, study = input_series(database, magn_series, export_study)
    if series is None:
        return
    return MDRegT1(series, study)

def map(database):
    series, study = input_series(database, magn_series, export_study)
    if series is None:
        return
    return T1(series, study)

def mdr_and_map(database):
    series, study = input_series(database, magn_series, export_study)
    if series is None:
        return
    _, moco = MDRegT1(series, study)
    return T1(moco, study)

def coreg_with_qa_mask(database):
    series_desc = [
        'mdr_par_T1',
        'T1w_abdomen_dixon_cor_bh_water_post_contrast',
        'mdr_par_a',
        'RK',
        'LK',
        'T1-RK',
        'T1-LK', 
    ]
    series, study = input_series(database, series_desc, export_study)
    if series is None:
        return
    series_return = []
    for r, kidney in enumerate(series[3:5]):
        params = vreg.find_sbs_rigid_transformation(series[0], series[1], tolerance=0.1, metric='mutual information', region=kidney, margin=0)
        T1_moved = vreg.apply_sbs_rigid_transformation(series[0], params, target=series[1])
        a_moved = vreg.apply_sbs_rigid_transformation(series[2], params, target=series[1])
        kidney_moved = vreg.apply_sbs_rigid_transformation(series[5+r], params, target=series[1])
        series_return.append((kidney, kidney_moved))
    return series_return



