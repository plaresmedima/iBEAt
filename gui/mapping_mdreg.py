from wezel.gui import Menu
import pipelines.mapping_mdreg as mapping_mdreg


def _if_a_series_is_selected(app):
    return app.nr_selected('Series') != 0


def _dti_model_fit(app):
    # Get input parameters including target for coregistration
    for series in app.selected('Series'):
        fit, par = mapping_mdreg.fit_DTI(series)
        app.display(fit)
        for p in par:
            app.display(p)
    app.refresh()


def _t1_model_fit(app):
    # Get input parameters including target for coregistration
    for series in app.selected('Series'):
        fit, par = mapping_mdreg.fit_T1(series)
        app.display(fit)
        for p in par:
            app.display(p)
    app.refresh()


def _t2_model_fit(app):
    # Get input parameters including target for coregistration
    for series in app.selected('Series'):
        fit, par = mapping_mdreg.fit_T2(series)
        app.display(fit)
        for p in par:
            app.display(p)
    app.refresh()


def _t2_star_model_fit(app):
    # Get input parameters including target for coregistration
    for series in app.selected('Series'):
        fit, par = mapping_mdreg.fit_T2star(series)
        app.display(fit)
        for p in par:
            app.display(p)
    app.refresh()


def _ivim_model_fit(app):
    # Get input parameters including target for coregistration
    for series in app.selected('Series'):
        fit, par = mapping_mdreg.fit_IVIM(series)
        app.display(fit)
        for p in par:
            app.display(p)
    app.refresh()


def _mt_model_fit(app):
    selected = app.selected('Series')
    if selected == []:
        return
    seriesList = selected[0].parent().parent().series()#selected[0].parent().children()
    value1 = seriesList.index(selected[0])
    try:
        value2 = seriesList.index(selected[1])
    except:
        value2 = value1
    
    seriesLabels = [s.label() for s in seriesList]
    cancel, f = app.dialog.input(
        {"label":"MT_OFF Series", "type":"dropdownlist", "list": seriesLabels, 'value':value1},
        {"label":"MT_ON Series", "type":"dropdownlist", "list": seriesLabels, 'value':value2},
        title = "Please select MT-OFF and MT-ON Series")
    if cancel:
        return
    series1 = seriesList[f[0]["value"]]
    series2 = seriesList[f[1]["value"]]
    fit, par = mapping_mdreg.fit_MT(series1, series2)
    app.display(fit)
    for p in par:
        app.display(p)
    app.refresh()


def _dce_model_fit(app):
    # Get input parameters including target for coregistration
    for series in app.selected('Series'):
        fit, par = mapping_mdreg.fit_DCE(series)
        app.display(fit)
        for p in par:
            app.display(p)
    app.refresh()


menu = Menu('iBEAt-reg')
menu.add_action('DTI model fit', on_clicked=_dti_model_fit, is_clickable=_if_a_series_is_selected)
menu.add_action('T1 model fit', on_clicked=_t1_model_fit, is_clickable=_if_a_series_is_selected)
menu.add_action('T2 model fit', on_clicked=_t2_model_fit, is_clickable=_if_a_series_is_selected)
menu.add_action('T2* model fit', on_clicked=_t2_star_model_fit, is_clickable=_if_a_series_is_selected)
menu.add_action('IVIM model fit', on_clicked=_ivim_model_fit, is_clickable=_if_a_series_is_selected)
menu.add_action('MT model fit', on_clicked=_mt_model_fit, is_clickable=_if_a_series_is_selected)
menu.add_action('DCE model fit', on_clicked=_dce_model_fit, is_clickable=_if_a_series_is_selected)