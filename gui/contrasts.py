from wezel.gui import Action
from wezel.displays import TableDisplay
import pipelines.asl as asl
import pipelines.T1 as T1


def _if_a_database_is_open(app): 
    return app.database() is not None
    
def _asl_perfusion(app):
    database = app.database()
    rbf_moved, rbf_stats = asl.perfusion(database)
    if rbf_moved is not None:
        for im in rbf_moved:
            app.display(im)
        app.addWidget(TableDisplay(rbf_stats), 'RBF values')
    app.refresh()

def _T1_mdr(app):
    database = app.database()
    fit, moco = T1.mdr(database)
    if moco is not None:
        app.display(moco)
    app.refresh()

def _T1_map(app):
    database = app.database()
    pars = T1.map(database)
    if pars is not None:
        for p in pars:
            app.display(p)
    app.refresh()


asl_perfusion = Action("Renal Blood Flow (ASL)", on_clicked=_asl_perfusion, is_clickable=_if_a_database_is_open)
T1_mdr = Action("T1 - motion correction", on_clicked=_T1_mdr, is_clickable=_if_a_database_is_open)
T1_map = Action("T1 - mapping", on_clicked=_T1_map, is_clickable=_if_a_database_is_open)

