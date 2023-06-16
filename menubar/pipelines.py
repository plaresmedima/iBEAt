import pandas as pd
from wezel.gui import Action
from dbdicom.wrappers import vreg, scipy
from wezel.displays import TableDisplay


def _if_a_database_is_open(app): 
    return app.database() is not None


def _get_series(db, desc):
    # Get series and check that not empty. 
    # If there are multiple, use the last one.
    db.message('Looking for '+ desc + ' in the list')
    series = db.series(SeriesDescription=desc)
    if series != []:
        return series[-1]
    else:
        msg = 'Series ' + desc + ' is missing. \n'
        msg += 'Cannot auto-calculate ASL perfusion.'
        db.dialog.information(msg)
        raise ValueError(msg)
    

def pipeline_asl_perfusion(database):

    # Get input parameters
    m0 = _get_series(database, 'ASL_kidneys_pCASL_cor-oblique_fb_M0_moco')
    dixon = _get_series(database, 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast')
    rbf = _get_series(database, 'ASL_kidneys_pCASL_cor-oblique_fb_RBF_moco')
    lk = _get_series(database, 'LK')
    rk = _get_series(database, 'RK')

    # Perform a separate registration for each target region
    rbf_moved = []
    rbf_stats = []
    for kidney in [(lk,'LK'), (rk,'RK')]:

        # Perform coregistration based on m0
        params = vreg.find_rigid_transformation(m0, dixon, tolerance=0.1, region=kidney[0])

        # Apply transformation to rbf image
        try:
            moved = vreg.apply_rigid_transformation(rbf, params, target=dixon, description='RBF [' + kidney[1] + ']')
            df = scipy.mask_statistics(kidney[0], moved)
            rbf_moved.append(moved)
            rbf_stats.append(df)
        except ValueError:
            msg = 'RBF series cannot be aligned.'
            database.dialog.information(msg)

    return rbf_moved, pd.concat(rbf_stats)


def _asl_perfusion(app):
    database = app.database()
    try:
        rbf_moved, rbf_stats = pipeline_asl_perfusion(database)
    except:
        return
    for im in rbf_moved:
        app.display(im)
    app.addWidget(TableDisplay(rbf_stats), 'RBF values')
    app.refresh()


asl_perfusion = Action("Renal Blood Flow (ASL)", on_clicked=_asl_perfusion, is_clickable=_if_a_database_is_open)
