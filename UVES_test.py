import matplotlib.pyplot as plt
%matplotlib inline

import tarfile
from astropy.utils.data import download_file
url = 'http://data.astropy.org/tutorials/UVES/data_UVES.tar.gz'
f = tarfile.open(download_file(url, cache = True), mode = 'r|*')
working_dir_path = 'C:/Users/Administrator/Desktop/donghun/AA'
#f.extractall(path = working_dir_path)

from glob import glob
import os
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits

globpath = os.path.join(working_dir_path, 'UVES/*.fits')
#print('globpath : ', globpath)
filelist = glob(globpath)
filelist.sort()
#print(filelist)

sp = fits.open(filelist[2])
#sp.info()

header = sp[0].header
#print('header : ', header)

wcs = WCS(header)
index = np.arange(header['NAXIS1'])
#print('index : ', index)

wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
#print('wavelength : ', wavelength)
#print('shape : ',wavelength.shape)
wavelength = wavelength.flatten()
#print('new wavelength : ', wavelength)
#print('new shape2 : ', wavelength.shape)

flux = sp[0].data

date_obs = header['Date-OBS']
#print(date_obs)
#print('header : ', header)



