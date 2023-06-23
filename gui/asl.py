from wezel.gui import Action
from wezel.displays import TableDisplay
import pipelines.asl as asl


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


perfusion = Action("Renal Blood Flow (ASL)", on_clicked=_asl_perfusion, is_clickable=_if_a_database_is_open)
