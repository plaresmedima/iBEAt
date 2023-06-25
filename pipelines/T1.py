from dbdicom.pipelines import input_series
from pipelines.mdr import MDRegT1
from pipelines.mapping import T1

export_study = "T1"
magn_series = "T1map_kidneys_cor-oblique_mbh_magnitude"


def mdr(database):
    series, study = input_series(database, magn_series, export_study)
    return MDRegT1(series, study)

def map(database):
    series, study = input_series(database, magn_series, export_study)
    return T1(series, study)

def mdr_and_map(database):
    series, study = input_series(database, magn_series, export_study)
    _, moco = MDRegT1(series, study)
    return T1(moco, study)
