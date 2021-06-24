import matplotlib.pyplot as plt
%matplotlib inline

import tarfile
from astropy.utils.data import download_file

from glob import glob
import os
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits

globpath = os.path.join(working_dir_path, 'UVES/*.fits')
#print('globpath : ', globpath)
filelist = glob(globpath)
filelist.sort()

sp = fits.open(filelist[-1])
