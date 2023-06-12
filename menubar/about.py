import webbrowser
from wezel.gui import Action


def _about_ibeat(app):
    webbrowser.open("https://bmcnephrol.biomedcentral.com/articles/10.1186/s12882-020-01901-x")


def _about_beat_dkd(app):
    webbrowser.open("https://www.beat-dkd.eu/")


ibeat = Action('iBEAt', on_clicked=_about_ibeat)
beat_dkd = Action('BEAt-DKD', on_clicked=_about_beat_dkd)



