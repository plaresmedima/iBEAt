import numpy as np
import scipy.ndimage as ndi
from wezel.gui import Action
from wezel.plugins.pyvista import SurfaceDisplay
from wezel.displays import TableDisplay
from dl_models import UNETR_kidneys_v1


def _if_a_database_is_open(app): 
    return app.database() is not None


def largest_cluster(array):
    # Label all features in the array
    label_img, cnt = ndi.label(array)
    # Find the label of the largest feature
    labels = range(1,cnt+1)
    size = [np.count_nonzero(label_img==l) for l in labels]
    max_label = labels[size.index(np.amax(size))]
    # Return a mask corresponding to the largest feature
    return label_img==max_label



def _UNETR_kidneys_v1(app):

    # Get file with weights and check if valid
    path = app.dialog.directory("Please select the folder with your kidney model")
    if path == '':
        app.status.hide() 
        return
    try:
        file = UNETR_kidneys_v1.weights_file(path)
    except FileNotFoundError as e:
        app.dialog.information(str(e))
        return

    # Find all series of the appropriate type
    app.status.message('Looking for ' + UNETR_kidneys_v1.trained_on)
    series = app.database().series(SeriesDescription=UNETR_kidneys_v1.trained_on)

    # If the correct sequences are not there, allow a manual selection (with warning)
    if series==[]:
        msg = 'The correct input data for this model have not been found. \n'
        msg += 'The results are not likely to be correct when applied to other data.'
        app.dialog.information(msg)
        series = app.database().series()
        sel = app.selected('Series')
        cancel, f = app.dialog.input( 
            {"label":"Apply model to:", "type":"select records", "options": series, 'default':sel},
            title = "Please select source data for segmentation manually")
        if cancel:
            return
        series = f[0]
    
    # Autosegment all series and save results in the database.
    for sery in series:
        desc = sery.instance().SeriesDescription
        array, header = sery.array(['SliceLocation','EchoTime'], pixels_first=True)
        sery.message('Calculating kidney masks for series ' + desc)

        # Calculate predictions 
        masks = UNETR_kidneys_v1.apply(array, file)

        # Clean masks
        left_kidney = largest_cluster(masks == 2)
        right_kidney = largest_cluster(masks == 1)

        # Save UNETR output
        result = sery.new_sibling(SeriesDescription = 'UNETR kidneys v1')
        result.set_array(masks, header, pixels_first=True)
        result[['WindowCenter','WindowWidth']] = [1.0,2.0]

        # Save and display left kidney data
        left = sery.new_sibling(SeriesDescription = 'LK')
        left.set_array(left_kidney, header, pixels_first=True)
        left[['WindowCenter','WindowWidth']] = [1.0,2.0]
        app.addWidget(SurfaceDisplay(left), title=left.label())

        # Save and display right kidney data
        right = sery.new_sibling(SeriesDescription = 'RK')
        right.set_array(right_kidney, header, pixels_first=True)
        right[['WindowCenter','WindowWidth']] = [1.0,2.0]
        app.addWidget(SurfaceDisplay(right), title=right.label())

        app.refresh()


segment = Action("Auto-segment kidneys", on_clicked=_UNETR_kidneys_v1, is_clickable=_if_a_database_is_open)


