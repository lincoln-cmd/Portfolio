import matplotlib.pyplot as plt
%matplotlib inline

import tarfile
from astropy.utils.data import download_file
url = 'http://data.astropy.org/tutorials/UVES/data_UVES.tar.gz'
f = tarfile.open(download_file(url, cache = True), mode = 'r|*')
#working_dir_path = 'C:/Users/Administrator/Desktop/donghun/AA'
#f.extractall(path = working_dir_path)

from glob import glob
import os
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits

globpath = os.path.join(working_dir_path, 'UVES/*.fits')
print(globpath)
filelist = glob(globpath)
filelist.sort()

sp = fits.open(filelist[-1])
sp.info()
