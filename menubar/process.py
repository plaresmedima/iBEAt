from os.path import dirname
import datetime
from wezel.gui import Action
import scripts.cluster.cluster_rename as rename
import scripts.cluster.cluster_mdr as mdr
import scripts.cluster.cluster_modelling as modelling
import scripts.cluster.cluster_segmentation as segmentation
import scripts.cluster.cluster_upload as upload


def _if_a_database_is_open(app): 
    return app.database() is not None


def _rename(app):
    folder = app.database()
    rename.main(folder)
    app.refresh()


def _mdr(app):
    folder = app.database()
    filename_log = folder.path() + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "log_MDR.txt"
    mdr.main(folder, filename_log)
    app.refresh()


def _mapping(app):
    folder = app.database()
    filename_log = folder.path() + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "log_mapping.txt"
    modelling.main(folder, filename_log)
    app.refresh()


def _segmentation(app):
    folder = app.database()
    filename_log = folder.path() + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "log_segmentation.txt"
    segmentation.main(folder, filename_log)
    app.refresh()


def _upload(app):
    folder = app.database()
    path = folder.path()
    filename_log = path + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "log_segmentation.txt"
    upload.main(path, filename_log)
    app.refresh()



action_rename = Action("Rename series..", on_clicked=_rename, is_clickable=_if_a_database_is_open)
action_mdr = Action("Model-driven registration..", on_clicked=_mdr, is_clickable=_if_a_database_is_open)
action_mapping = Action("Parameter mapping..", on_clicked=_mapping, is_clickable=_if_a_database_is_open)
action_segmentation = Action("Segmentation..", on_clicked=_segmentation, is_clickable=_if_a_database_is_open)
action_upload = Action("Upload to drive..", on_clicked=_upload, is_clickable=_if_a_database_is_open)

