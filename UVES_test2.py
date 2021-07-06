import matplotlib.pyplot as plt
#%matplotlib inline

import tarfile
from astropy.utils.data import download_file

from glob import glob
import os
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits

url = 'http://data.astropy.org/tutorials/UVES/data_UVES.tar.gz'
f = tarfile.open(download_file(url, cache = True), mode = 'r|*')
working_dir_path = 'C:/Users/Administrator/Desktop/donghun/AA'
#f.extractall(path = working_dir_path)

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
    
    flux = sp[0].data
    
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


#for i in filelist:
 #   print(test(i))


flux = np.zeros((len(filelist), len(wavelength)))
date = np.zeros((len(filelist)), dtype = 'U23')

for i, fname in enumerate(filelist):
    w, f, date_obs = read_spec(fname)
    flux[i,:] = f
    date[i] = date_obs
    #print('flux{0} : {1}'.format(i, flux))
    
import astropy.units as u
from astropy.constants.si import c, G, M_sun, R_sun

wavelength = wavelength * u.AA

heliocentric = -23. * u.km/u.s
v_rad = -4.77 * u.km/ u.s
R_MN_Lup = 0.9 * R_sun
M_MN_Lup = 0.6 * M_sun
vsini = 74.6 * u.km / u.s
period = 0.439 * u.day

inclination = 45. * u.degree
incl = inclination.to(u.radian)
    
v_accr = (2. * G * M_MN_Lup / R_MN_Lup) ** 0.5
#print(v_accr)
#print(v_accr.cgs)
from astropy.units import imperial
#print(v_accr.to(imperial.yd / u.hour))

v_rot = vsini / np.sin(incl)
#print(v_accr / v_rot)
#print((v_accr / v_rot).decompose())

#print(wavelength)
#wavelength = wavelength * (1. + heliocentric / c)
#print(wavelength)
#wavelength = wavelength * (1. * u.dimensionless_unscaled + heliocentric / c)
#print(wavelength)

#energy = wavelength.to(u.keV, equivalencies=u.spectral())
#frequency = wavelength.to(u.Hz, equivalencies=u.spectral())
#print(energy)
#print(frequency)

gravi = np.log10((G*M_MN_Lup / R_MN_Lup**2) / u.cm*u.second**2)
print(gravi)