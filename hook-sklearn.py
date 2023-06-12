from PyInstaller.utils.hooks import collect_data_files
datas = collect_data_files('sklearn')

hiddenimports = [
    "sklearn.metrics._pairwise_distances_reduction._middle_term_computer", 
    "sklearn.metrics._pairwise_distances_reduction._datasets_pair"
]