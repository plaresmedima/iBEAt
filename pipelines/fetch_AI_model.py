import os
import utilities.zenodo_link as UNETR_zenodo
import requests

def dl_models(database):

    unetr, unetr_link= UNETR_zenodo.main()
    record_id = unetr_link.split('.')[-1]

    unetr_path = os.path.join(os.path.dirname(database.path()),unetr)

    if os.path.exists(unetr_path):
        database.log("UNETR was found in the local folder")
        return
    else:   

        zenodo_url = f"https://zenodo.org/records/{record_id}/files/{unetr}?download=1"

        with requests.get(zenodo_url) as req:
                    with open(os.path.join(os.path.dirname(database.path()),unetr), 'wb') as f:
                        for chunk in req.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)





