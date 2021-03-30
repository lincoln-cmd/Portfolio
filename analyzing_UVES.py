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
working_dir_path = '.' # change the path
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
    wavelength = wavelength.flatten()
    flux = sp[0].data
    
    date_obs = header['Date-OBS']
    return wavelength, flux, date_obs

help(read_spec)


# The dataset of UVES spectra should have been taken using all the same setup
# EXPTIME : exposure time
# CRVAL1 : wavelength zero point
# HIERARCH ESO INS PATH : arm used(UVES has a red and a blue arm)
def read_setup(filename):
    '''
    Get setup for UVES spectrum from the ESO pipeline
    
    Parameters
    ----------
    filename : string
    name of the fits file with the data
    
    Returns
    ----------
    exposure_time : float
    wavelength_zero_point : float
    optical_arm : string
    '''
    
    sp = fits.open(filelist[0])
    header = sp[0].header
    
    return header['EXPTIME'], header['CRVAL1'], header['HIERARCH ESO INS PATH']

for f in filelist:
    print(read_setup(f))
 
'''
 The UVES pipeline that was used to reduce the data employs a fixed wavelength grid, so the wavelength is the same for all apectra,
which allows us to define an array that can hold the fluxes of all observations easily.
'''
# loop over the list of all filenames and fill this array with data
flux = np.zeros((len(filelist), len(wavelength)))
date = np.zeros((len(filelist)), dtype='U23') # data comes as string with 23 characters (dytpe = 'S23')

for i, fname in enumerate(filelist):
    w, f, date_obs = read_spec(fname)
    flux[i, :] = f
    date[i] = date_obs

import astropy.units as u
from astropy.constants.si import c, G, M_sun, R_sun

# all the constants have to have their own units, so the constant values have to be multiplied with units
wavelength = wavelength * u.AA

heliocentric = -23. * u.km/u.s
v_rad = -4.77 * u.km / u.s
R_MN_Lup = 0.9 * R_sun
M_MN_Lup = 0.6 * M_sun
vsini = 74.6 * u.km / u.s
period = 0.439 * u.day

inclination = 45. * u.degree

incl = inclination.to(u.radian) # convert the degree unit to radian unit

'''
 MN Lup is T Tauri star(TTS), which is possibly surrounded by an accertion disk, so we can expect those
accretion signatures to appear close to the free-fall velocity(v) that a mass(m) reaches, when it hits the stellar surface.

 -> extract infall speed with energy conservation equation
 0.5mv^2 = gmM / r
 (kenetic energy = gravitational energy), (gravitational energy equation = U(r) = - GmM / r)
'''

v_accr = (2. * G * M_MN_Lup / R_MN_Lup)**0.5
print(v_accr) # si system
print(v_accr.cgs) # cgs system
from astropy.units import imperial
print(v_accr.to(imperial.yd / u.hour)) # imperial system(english measuring system)

# the relation accretion velocity with rotational velocity
v_rot = vsini / np.sin(incl)
print((v_accr / v_rot).decompose()) # (v_accr / v_rot) is composed different units, m and km

'''
 λ_heliocentric = λ_bariocentric * (1 + v_helio / c)
-> correct the wavelength scale to the heliocentric velocity scale by utilizing the 'astropy.units'
'''
# wavelength = wavelength * (1. + heliocentric / c) # all the constants have each different unit value, km, m, and just a number
wavelength = wavelength * (1. * u.dimensionless_unscaled + heliocentric / c)
print(wavelength.to(u.keV, equivalencies = u.spectral())) # 1eV = 1.60218 x 10^-19
print(wavelength.to(u.Hz, equivalencies = u.spectral()))

'''
 Spectroscopically, MN Lup is classified as spectral type M0 V, so the gravitational acceleration on the surface log(g) should be comparable to the sun.
log(g) can be computered with mass and radius values.
-> log(g) is consistent
'''
print(np.log10((G * M_MN_Lup / R_MN_Lup**2) / u.cm * u.second**2)) # cgs system


'''
# making the function to convert the wavelength scal into the velocity scale
 Convert a wavelength scale into a velocity scale. Input the wavelengths array and the rest wavelength of a spectral line.
This function can be employed to show the red and blueshift of the spectrum relative to the Ca II H line.
'''
waveclosetoHa = np.array([6562., 6563, 6565.]) * u.AA


