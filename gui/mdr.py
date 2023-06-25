from wezel.gui import Action
from dbdicom.pipelines import input_series
import pipelines.mdr as mdr

study_desc = 'MDRresults'


def _if_a_database_is_open(app): 
    return app.database() is not None


def _T2star(app):
    series, study = input_series(app.database(), "T2star_map_kidneys_cor-oblique_mbh_magnitude", study_desc)
    if series is not None:
        fit, _ = mdr.MDRegT2star(series, study=study)
        app.display(fit)
    app.refresh()

def _T1(app):
    series, study = input_series(app.database(), "T1map_kidneys_cor-oblique_mbh_magnitude", study_desc)
    if series is not None:
        fit, _ = mdr.MDRegT1(series, study=study)
        app.display(fit)
    app.refresh()

def _T2(app):
    series, study = input_series(app.database(), "T2map_kidneys_cor-oblique_mbh_magnitude", study_desc)
    if series is not None:
        fit, _ = mdr.MDRegT2(series, study=study)
        app.display(fit)
    app.refresh()

def _DTI(app):
    series, study = input_series(app.database(), "DTI_kidneys_cor-oblique_fb", study_desc)
    if series is not None:
        fit, _ = mdr.MDRegDTI(series, study=study)
        app.display(fit)
    app.refresh()

def _MT(app):
    series, study = input_series(app.database(), ["MT_OFF_kidneys_cor-oblique_bh", "MT_ON_kidneys_cor-oblique_bh"], study_desc)
    if series is not None:
        fit, _ = mdr.MDRegMT(series, study=study)
        app.display(fit)
    app.refresh()

def _DCE(app):
    series, study = input_series(app.database(), "DCE_kidneys_cor-oblique_fb", study_desc)
    if series is not None:
        fit, _ = mdr.MDRegDCE(series, study=study)
        app.display(fit)
    app.refresh()


action_T2star = Action("T2* - motion correction", on_clicked=_T2star, is_clickable=_if_a_database_is_open)
action_T1 = Action("T1 - motion correction", on_clicked=_T1, is_clickable=_if_a_database_is_open)
action_T2 = Action("T2 - motion correction", on_clicked=_T2, is_clickable=_if_a_database_is_open)
action_DTI = Action("DTI - motion correction", on_clicked=_DTI, is_clickable=_if_a_database_is_open)
action_MT = Action("MT - motion correction", on_clicked=_MT, is_clickable=_if_a_database_is_open)
action_DCE = Action("DCE - motion correction", on_clicked=_DCE, is_clickable=_if_a_database_is_open)