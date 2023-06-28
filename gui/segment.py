
from wezel.displays import TableDisplay
from wezel.gui import Action
from wezel.plugins.pyvista import SurfaceDisplay
import pipelines.segment as seg 


weights = 'C:\\Users\\steve\\Dropbox\\DLmodels\\UNETR_kidneys_v1.pth'


def _if_a_database_is_open(app): 
    return app.database() is not None


def _whole_kidney_canvas(app):
    clusters = seg.compute_whole_kidney_canvas(app.database())
    if clusters is not None:
        app.display(clusters)
    app.refresh()


def _renal_sinus_fat(app):
    sf_series, sf_biomarkers = seg.compute_renal_sinus_fat(app.database())
    if sf_series is not None:
        for sinus_fat in sf_series:
            app.addWidget(SurfaceDisplay(sinus_fat), title=sinus_fat.label())
        app.addWidget(TableDisplay(sf_biomarkers), 'RSF - Volume features')
    app.refresh()


def _UNETR_kidneys_v1(app):
    database = app.database()
    kidneys, features = seg.segment_kidneys(database, weights)
    if kidneys is None:
        return
    for kidney in kidneys:
        app.addWidget(SurfaceDisplay(kidney), title=kidney.label())
    app.addWidget(TableDisplay(features), 'Kidney - Volume features')
    app.refresh()


whole_kidney_canvas = Action("Segmentation canvas", on_clicked=_whole_kidney_canvas, is_clickable=_if_a_database_is_open)
renal_sinus_fat = Action("Renal Sinus Fat", on_clicked=_renal_sinus_fat, is_clickable=_if_a_database_is_open)
action_unetr = Action("Auto-segment kidneys", on_clicked=_UNETR_kidneys_v1, is_clickable=_if_a_database_is_open)



