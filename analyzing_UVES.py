# Analyzing UVES Spectroscopy with Astropy

# from : https://learn.astropy.org/rst-tutorials/UVES.html?highlight=filtertutorials

import matplotlib.pyplot as plt
%matplotlib inline

# download tar files and extract the files' data
import tarfile
from astropy.utils.data import download_file
url = 'http://data.astropy.org/tutorials/UVES/data_UVES.tar.gz'
f = tarfile.open(download_file(url, cache = True), mode = 'r|*')
working_dir_path = '.' # change the path
f.extractall(path = working_dir_path)

