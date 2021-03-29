# Analyzing UVES Spectroscopy with Astropy

# from : https://learn.astropy.org/rst-tutorials/UVES.html?highlight=filtertutorials


'''
 As the accretion rate decreases, the impact on the central star must change, which means that
the accretion rate is depended on the center of the star. Furthermore, the emission frequence is alse
depended on the center of the star. (The accretion disk of protostar emits the infrared radiation, and
                                     the nuetron star and blackhole emit the X-ray)
from : https://ui.adsabs.harvard.edu/abs/2013ApJ...771...70G/abstract




'''

import matplotlib.pyplot as plt
%matplotlib inline

# download tar files and extract the files' data
import tarfile
from astropy.utils.data import download_file
#url = 'http://data.astropy.org/tutorials/UVES/data_UVES.tar.gz'
#f = tarfile.open(download_file(url, cache = True), mode = 'r|*')
#working_dir_path = '.' # change the path
#f.extractall(path = working_dir_path)

# analyze data from NM Lup, a T Tauri star in the Taurus-Auriga star forming region located at a distance of about 140pc
# MN Lup has been observed simultaneously with XMM-Newton and the UVES spectrograph on the VLT

'''
- Spatially resolving the accretion shocks on the rapidly-rotating M0 T-Tauri star MN Lupi
 The accretion mechanism is the cause of its rapid surface rotation because ongoing disk accretion.
 
from : https://ui.adsabs.harvard.edu/abs/2005A%26A...440.1105S/abstract
'''

# reading the data
from glob import glob
import os
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits

globpath = os.path.join(working_dir_path, 'UVES/*.fits') # os.path.join : join two directories

print(globpath)
filelist = glob(globpath) # search through directories similar to the Unix shell

filelist.sort()

# read the first FITS file in the list and check what is in there
sp = fits.open(filelist[0])
sp.info()

# extract the WCS from header to get the wavelength coordinate.
# Even though the warnings are invoked in this process about a non-standard RADECSYS, the WCS will still work
header = sp[0].header

wcs = WCS(header)
index = np.arange(header['NAXIS1']) # make the array

wavelength = wcs.wcs_pix2world(index[:, np.newaxis], 0)
wavelength.shape
wavelength = wavelength.flatten() # adjust the wrong dimension by utilizing flatten()

# the flux is contained in the primary image
flux = sp[0].data


# make the function to reuse the code
# input -> filename, returns -> wvaelength, flux arrays, and the time of the observation
#from spectra_utils import func
def read_spec(filename):
    ''' Read a UVES spectrum from the ESO pipeline
    
    Parameters
    -----------
    filename : string. name of the fits file with the data
    
    Returns
    -----------
    wavelength : np.ndarray. wavelength(in Ang)
    flux : np.ndarray. flux (in erg/s/cm**2)
    date_obs : string. time of observation
    '''
    sp = fits.open(filename)
    header = sp[0].header
    
    wcs = WCS(header)
    index = np.arange(header['NAXIS1']) # make index array
    
    wavelength = wcs.wcs_pix2world(index[:, np.newaxis], 0)
    wvaelength = wvaelength.flatten()
    flux = sp[0].data
    
    date_obs = header['Date-OBS']
    return wavelength, flux, date_obs

help(read_spec)