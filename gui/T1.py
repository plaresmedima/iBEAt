from wezel.gui import Action
import pipelines.T1 as T1


def _if_a_database_is_open(app): 
    return app.database() is not None
    

def _mdr(app):
    database = app.database()
    fit, moco = T1.mdr(database)
    if moco is not None:
        app.display(moco)
    app.refresh()


action_mdr = Action("T1 - motion correction", on_clicked=_mdr, is_clickable=_if_a_database_is_open)