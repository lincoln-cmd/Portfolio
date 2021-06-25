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

def read_spec(filename):

    
    sp = fits.open(filename)
    
    header = sp[0].header
    
    wcs = WCS(header)
    #print('header : ', header)
    #print('wcs : ', wcs)
    
    index = np.arange(header['NAXIS1'])
    
    wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
    wavelength = wavelength.flatten()
    
    flux = sp[0]
    
    date_obs = header['Date-OBS']
    return wavelength, flux, date_obs

#print(read_spec(filelist[0]))

def read_setup(filename):
    sp = fits.open(filelist[0])
    header = sp[0].header
    
    return header['EXPTIME'], header['CRVAL1'], header['HIERARCH ESO INS PATH']

#for i in filelist:
 #   print(read_setup(i))
'''
excersice
'''
def test(filename):
    sp = fits.open(filename)
    header = sp[0].header
    
    return header['CDELT1'], header['DATAMIN'], header['DATAMAX']

print(test(filelist[11]))





'''
excersice
'''