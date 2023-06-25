from dbdicom.pipelines import input_series
from pipelines.mdr import MDRegT1

export_study = "T1"

def mdr(database):
    series, study = input_series(database, "T1map_kidneys_cor-oblique_mbh_magnitude", export_study)
    return MDRegT1(series, study)
