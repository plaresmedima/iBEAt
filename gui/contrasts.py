from wezel.gui import Action, Menu
from wezel.displays import TableDisplay
from wezel.plugins.pyvista import SurfaceDisplay
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

def _T1_mdr_and_map(app):
    database = app.database()
    pars = T1.mdr_and_map(database)
    if pars is not None:
        for p in pars:
            app.display(p)
    app.refresh()

def _T1_coreg_with_qa_mask(app):
    database = app.database()
    pars = T1.coreg_with_qa_mask(database)
    if pars is not None:
        for p in pars:
            viewer = SurfaceDisplay(p[1], reference=p[0], triangulate=False)
            app.addWidget(viewer, title=p[1].label())
    app.refresh()


asl_perfusion = Action("Renal Blood Flow (ASL)", on_clicked=_asl_perfusion, is_clickable=_if_a_database_is_open)
T1_mdr = Action("T1 - motion correction", on_clicked=_T1_mdr, is_clickable=_if_a_database_is_open)
T1_map = Action("T1 - mapping", on_clicked=_T1_map, is_clickable=_if_a_database_is_open)
T1_mdr_and_map = Action("T1 - motion correction and mapping", on_clicked=_T1_mdr_and_map, is_clickable=_if_a_database_is_open)
T1_coreg_qa = Action("T1 - coregistration (with QA mask)", on_clicked=_T1_coreg_with_qa_mask, is_clickable=_if_a_database_is_open)


menu_T1 = Menu('T1')
menu_T1.add(T1_mdr)
menu_T1.add(T1_map)
menu_T1.add(T1_mdr_and_map)
menu_T1.add(T1_coreg_qa)

