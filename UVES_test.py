import matplotlib.pyplot as plt
%matplotlib inline

import tarfile
from astropy.utils.data import download_file
url = 'http://data.astropy.org/tutorials/UVES/data_UVES.tar.gz'
f = tarfile.open(download_file(url, cache = True), mode = 'r|*')
working_dir_path = 'C:/Users/Administrator/Desktop/donghun/AA'
f.extractall(path = working_dir_path)




