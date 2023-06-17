from wezel.gui import Action
from dbdicom.pipelines import input_series

import pipelines.mapping as map

study_desc = 'ModellingResults'


def _if_a_database_is_open(app): 
    return app.database() is not None


def _T2star(app):
    series, study = input_series(app.database(), "T2star_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco", study_desc)
    if series is not None:
        M0, fw, T2s = map.T2s_Modelling(series, study=study)
        app.display(M0)
        app.display(fw)
        app.display(T2s)
    app.refresh()


def _DTI(app):
    series, study = input_series(app.database(), "DTI_kidneys_cor-oblique_fb_mdr_moco", study_desc)
    if series is not None:
        FA, MD = map.DTI_Modelling(series, study=study)
        app.display(FA)
        app.display(MD)
    app.refresh()


def _DCE(app):
    series, study = input_series(app.database(), "DCE_kidneys_cor-oblique_fb_mdr_moco", study_desc)
    if series is not None:
        MAX, AUC = map.DCE_MAX_Modelling(series, study=study)
        app.display(MAX)
        app.display(AUC)
    app.refresh()


def _MTR(app):
    series, study = input_series(app.database(), 'MT_ON_kidneys_cor-oblique_bh_mdr_moco', study_desc)
    if series is not None:
        MTR = map.MTR_Modelling(series, study=study)
        app.display(MTR)
    app.refresh()


def _T1(app):
    series, study = input_series(app.database(), 'T1_map_kidneys_cor-oblique_mbh_magnitude_mdr_moco', study_desc)
    if series is not None:
        M0, T1app, T1 = map.T1_Modelling_Philips(series, study=study)
        app.display(M0)
        app.display(T1app)
        app.display(T1)
    app.refresh()


action_T2star = Action("Parameter mapping (T2*)", on_clicked=_T2star, is_clickable=_if_a_database_is_open)
action_T1 = Action("Parameter mapping (T1)", on_clicked=_T1, is_clickable=_if_a_database_is_open)
action_DTI = Action("Parameter mapping (DTI)", on_clicked=_DTI, is_clickable=_if_a_database_is_open)
action_MT = Action("Parameter mapping (MT)", on_clicked=_MTR, is_clickable=_if_a_database_is_open)
action_DCE = Action("Parameter mapping (DCE)", on_clicked=_DCE, is_clickable=_if_a_database_is_open)

