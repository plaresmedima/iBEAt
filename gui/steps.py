import datetime
from wezel.gui import Action, Menu
from wezel.displays import TableDisplay
from wezel.plugins.pyvista import SurfaceDisplay

from pipelines import rename, fetch, segment, measure, mapping, mdr, export
import scripts.upload as upload

# Needs to be in a central location
weights = 'C:\\Users\\steve\\Dropbox\\Data\\dl_models\\UNETR_kidneys_v1.pth'

def _never(app):
    return False

def _if_a_database_is_open(app): 
    return app.database() is not None


def rename_all_series(app):
    folder = app.database()
    rename.all_series(folder)
    rename.check(folder)
    app.refresh()


def fetch_kidney_masks(app):
    folder = app.database()
    fetch.kidney_masks(folder)
    app.refresh()


def segment_kidneys(app):
    database = app.database()
    kidneys = segment.kidneys(database, weights)
    for kidney in kidneys:
        app.addWidget(SurfaceDisplay(kidney), title=kidney.label())
    app.refresh()


def segment_renal_sinus_fat(app):
    database = app.database()
    sinus_fats = segment.renal_sinus_fat(database)
    if sinus_fats is not None:
        for sinus_fat in sinus_fats:
            app.addWidget(SurfaceDisplay(sinus_fat), title=sinus_fat.label())
    app.refresh()


def whole_kidney_canvas(app):
    folder = app.database()
    clusters = segment.compute_whole_kidney_canvas(folder)
    if clusters is not None:
        app.display(clusters)
    app.refresh()


def measure_kidney_volumetrics(app):
    database = app.database()
    features = measure.kidney_volumetrics(database)
    app.addWidget(TableDisplay(features), 'Kidney - Volume features')
    app.refresh()


def measure_sinus_fat_volumetrics(app):
    database = app.database()
    features = measure.sinus_fat_volumetrics(database)
    app.addWidget(TableDisplay(features), 'Renal sinus fat - Volume features')
    app.refresh()


def export_kidney_segmentations(app):
    database = app.database()
    export.kidney_masks_as_dicom(database)
    export.kidney_masks_as_png(database)
    export.whole_kidney_canvas(database)
    app.refresh()


def all_segmentation_steps(app):
    answer = app.dialog.question('This is going to take a while. Do you want to continue?')
    if 'Yes' != answer:
        return
    fetch_kidney_masks(app)
    segment_kidneys(app)
    segment_renal_sinus_fat(app)
    whole_kidney_canvas(app)
    measure_kidney_volumetrics(app)
    measure_sinus_fat_volumetrics(app)
    export_kidney_segmentations(app)


def _mdr(app):
    folder = app.database()
    filename_log = folder.path() + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "log_MDR.txt"
    mdr.main(folder, filename_log)
    app.refresh()


def _mapping(app):
    folder = app.database()
    filename_log = folder.path() + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "log_mapping.txt"
    mapping.main(folder, filename_log)
    app.refresh()


def _upload(app):
    folder = app.database()
    path = folder.path()
    filename_log = path + datetime.datetime.now().strftime('%Y%m%d_%H%M_') + "log_segmentation.txt"
    upload.main(path, filename_log)
    app.refresh()



action_rename = Action("Rename series..", on_clicked=rename_all_series, is_clickable=_if_a_database_is_open)

action_all_segmentation_steps = Action("All segmentation steps..", on_clicked=all_segmentation_steps, is_clickable=_if_a_database_is_open)
action_fetch_kidney_masks = Action("Fetch kidney masks..", on_clicked=fetch_kidney_masks, is_clickable=_never)
action_segment_kidneys = Action("Auto-segment kidneys...", on_clicked=segment_kidneys, is_clickable=_if_a_database_is_open)
action_renal_sinus_fat = Action("Segment renal sinus fat..", on_clicked=segment_renal_sinus_fat, is_clickable=_if_a_database_is_open)
action_whole_kidney_canvas = Action("Calculate segmentation canvas..", on_clicked=whole_kidney_canvas, is_clickable=_if_a_database_is_open)
action_export_kidney_segmentations = Action("Export kidney segmentations..", on_clicked=export_kidney_segmentations, is_clickable=_if_a_database_is_open)
action_kidney_volumetrics = Action("Measure kidney volumetrics..", on_clicked=measure_kidney_volumetrics, is_clickable=_if_a_database_is_open)
action_sinus_fat_volumetrics = Action("Measure sinus fat volumetrics..", on_clicked=measure_sinus_fat_volumetrics, is_clickable=_if_a_database_is_open)


action_mdr = Action("Model-driven registration..", on_clicked=_mdr, is_clickable=_if_a_database_is_open)
action_mapping = Action("Parameter mapping..", on_clicked=_mapping, is_clickable=_if_a_database_is_open)
action_upload = Action("Upload to drive..", on_clicked=_upload, is_clickable=_if_a_database_is_open)


menu_segment = Menu('Segment..')
menu_segment.add(action_all_segmentation_steps)
menu_segment.add_separator()
menu_segment.add(action_fetch_kidney_masks)
menu_segment.add(action_segment_kidneys)
menu_segment.add(action_renal_sinus_fat)
menu_segment.add_separator()
menu_segment.add(action_whole_kidney_canvas)
menu_segment.add(action_export_kidney_segmentations)
menu_segment.add_separator()
menu_segment.add(action_kidney_volumetrics)
menu_segment.add(action_sinus_fat_volumetrics)