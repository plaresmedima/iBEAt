from wezel.gui import Action
import scripts.cluster.cluster_mdr as mdr


def _if_a_database_is_open(app): 
    return app.database() is not None


def _input(database, series_desc):
    series = database.series(SeriesDescription=series_desc)
    if series == []:
        database.dialog.information("Cannot find series " + series_desc)
        raise FileNotFoundError
    series = series[-1]
    study_desc = 'MDRresults'
    studies = database.studies(StudyDescription=study_desc)
    if studies == []:
        study = series.new_pibling(StudyDescription=study_desc)
    else:
        study = studies[-1]
    return series, study


def _T2star(app):
    try:
        series, study = _input(app.database(), "T2star_map_kidneys_cor-oblique_mbh_magnitude")
    except:
        return
    fit, _ = mdr.MDRegT2star(series, study=study)
    app.display(fit)
    app.refresh()

def _T1(app):
    try:
        series, study = _input(app.database(), "T1map_kidneys_cor-oblique_mbh_magnitude")
    except:
        return
    fit, _ = mdr.MDRegT1(series, study=study)
    app.display(fit)
    app.refresh()

def _T2(app):
    try:
        series, study = _input(app.database(), "T2map_kidneys_cor-oblique_mbh_magnitude")
    except:
        return
    fit, _ = mdr.MDRegT2(series, study=study)
    app.display(fit)
    app.refresh()

def _DTI(app):
    try:
        series, study = _input(app.database(),  "DTI_kidneys_cor-oblique_fb")
    except:
        return
    fit, _ = mdr.MDRegDTI(series, study=study)
    app.display(fit)
    app.refresh()

def _MT(app):
    try:
        mt_off, study = _input(app.database(),  "MT_OFF_kidneys_cor-oblique_bh")
        mt_on, study = _input(app.database(),  "MT_ON_kidneys_cor-oblique_bh")
    except:
        return
    fit, _ = mdr.MDRegMT([mt_off, mt_on], study=study)
    app.display(fit)
    app.refresh()

def _DCE(app):
    try:
        series, study = _input(app.database(),  "DCE_kidneys_cor-oblique_fb")
    except:
        return
    fit, _ = mdr.MDRegDCE(series, study=study)
    app.display(fit)
    app.refresh()


action_T2star = Action("Remove breathing (T2*)", on_clicked=_T2star, is_clickable=_if_a_database_is_open)
action_T1 = Action("Remove breathing (T1)", on_clicked=_T1, is_clickable=_if_a_database_is_open)
action_T2 = Action("Remove breathing (T2)", on_clicked=_T2, is_clickable=_if_a_database_is_open)
action_DTI = Action("Remove breathing (DTI)", on_clicked=_DTI, is_clickable=_if_a_database_is_open)
action_MT = Action("Remove breathing (MT)", on_clicked=_MT, is_clickable=_if_a_database_is_open)
action_DCE = Action("Remove breathing (DCE)", on_clicked=_DCE, is_clickable=_if_a_database_is_open)