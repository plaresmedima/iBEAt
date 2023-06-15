from dbdicom.wrappers import sklearn, skimage, scipy, dipy
from wezel.displays import TableDisplay
from wezel.gui import Action
from wezel.plugins.pyvista import SurfaceDisplay


def _if_a_database_is_open(app): 
    return app.database() is not None


def _whole_kidney_mask(app):

    fat_desc = 'T1w_abdomen_dixon_cor_bh_fat_post_contrast' 
    out_desc = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'

    app.status.message('Finding correct sequences in the list..')

    fat = app.database().series(SeriesDescription=fat_desc)
    out = app.database().series(SeriesDescription=out_desc)
    msg = 'OK'
    if fat==[] or out==[]:
        msg = 'Not all sequences have been found - please check the list'
    elif len(fat)>1 or len(out)>1:
        msg = 'Some sequences have been repeated - please select the appropriate for analysis'
    if msg == 'OK':
        features = fat + out
    else:
        all = app.database().series()
        cancel, features = app.dialog.input(
            {'label': fat_desc, 'type':'select record', 'options': all, 'default': fat},
            {'label': out_desc, 'type':'select record', 'options': all, 'default': out},
            title = msg)
        if cancel:
            return
        
    app.status.message('Preparing for kidney preprocessing..')

    clusters = sklearn.sequential_kmeans(features, n_clusters=2, multiple_series=True)
    app.display(clusters)
    app.refresh()



def _renal_sinus_fat(app):

    app.status.message('Finding correct sequences in the list..')

    fat_image = app.database().series(SeriesDescription='T1w_abdomen_dixon_cor_bh_fat_post_contrast')
    lk = app.database().series(SeriesDescription='LK')
    rk = app.database().series(SeriesDescription='RK')

    all = app.database().series()
    cancel, selected = app.dialog.input(
        {"label":"Fat image (have you checked for fat-water swap?)", "type":"select record", "options": all, 'default': fat_image},
        {"label":"Whole-kidney regions", "type":"select records", "options": all, 'default': lk+rk},
        title = "Please select input for renal sinus fat processing")
    if cancel:
        return
    fat_image = selected[0]
    kidneys = selected[1]

    cleanup = True

    sf_series = []
    fat_image_masked, fat_mask = dipy.median_otsu(fat_image, median_radius=1, numpass=1)
    for kidney in kidneys:

        # Pipeline calculation
        kidney_hull = skimage.convex_hull_image_3d(kidney)
        sinus_fat = scipy.image_calculator(fat_mask, kidney_hull, 'series 1 * series 2', integer=True)
        sinus_fat_open = skimage.opening_3d(sinus_fat)
        sinus_fat_largest = scipy.extract_largest_cluster_3d(sinus_fat_open)
        sinus_fat_largest.SeriesDescription = kidney.instance().SeriesDescription + 'SF'

        # Remove intermediate results
        if cleanup:
            kidney_hull.remove()
            sinus_fat.remove()
            sinus_fat_open.remove()
    
        # Append and display
        sf_series.append(sinus_fat_largest)
        viewer = SurfaceDisplay(sinus_fat_largest)
        app.addWidget(viewer, title=sinus_fat_largest.label())

    if cleanup:
        fat_image_masked.remove()
        fat_mask.remove()

    # Collect features & display
    df = skimage.volume_features(sf_series)
    app.addWidget(TableDisplay(df), 'Volume features')
    app.refresh()



whole_kidney_mask = Action("Presegment whole kidneys", on_clicked=_whole_kidney_mask, is_clickable=_if_a_database_is_open)
renal_sinus_fat = Action("Renal sinus fat", on_clicked=_renal_sinus_fat, is_clickable=_if_a_database_is_open)



