import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
from astropy.modeling import Fittable1DModel, Parameter
from astroquery.sdss import SDSS

spectrum = SDSS.get_spectra(plate = 1349, fiberID = 216, mjd = 52797)[0]

flux = spectrum[1].data['flux']
lam = 10**(spectrum[1].data['loglam'])

units_flux = spectrum[0].header['bunit']
units_wavelength_full = spectrum[0].header['WAT1_001']
units_wavelength = units_wavelength_full[36:]

'''
 - Exercise
 1. choose correct model
 2. the problem of fitter
 3. adjust the range of data
 4. redesign plotting of data

'''



# employing the RickerWavelet1D as a model
'''
rick_model = models.RickerWavelet1D()
fitter = fitting.LevMarLSQFitter()
rick_fit = fitter(rick_model, lam, flux)

plt.figure(figsize = (8, 5))
plt.plot(lam, flux, color = 'k')
plt.plot(lam, rick_fit(lam), color = 'r')
plt.xlim(6300, 6700)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux ({})'.format(units_flux))
plt.show()
'''

# using the Polynomial1D
'''
poly_model = models.Polynomial1D(1)
fitter = fitting.LevMarLSQFitter()
poly_fit = fitter(poly_model, lam, flux)

plt.figure(figsize = (8, 5))
plt.plot(lam, flux, color = 'k')
plt.plot(lam, poly_fit(lam), color = 'r')
plt.xlim(6300, 6700)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux ({})'.format(units_flux))
plt.show()
'''

# utilizing the SLSQPLSQFitter as the fitter
'''
gaussian_model = models.Gaussian1D()
fitter = fitting.SLSQPLSQFitter()
gaussian_fit = fitter(gaussian_model, lam, flux)

plt.figure(figsize = (8, 5))
plt.plot(lam, flux, color = 'k')
plt.plot(lam, gaussian_fit(lam), color = 'r')
plt.xlim(6300, 6700)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('flux ({})'.format(units_flux))
plt.show()
'''

# the usage of the SimplexLSQFitter() as the fitter
'''
gaussian_model = models.Gaussian1D()
fitter = fitting.SimplexLSQFitter()
gaussian_fit = fitter(gaussian_model, lam, flux)

plt.figure(figsize = (8, 5))
plt.plot(lam, flux, color = 'k')
plt.plot(lam, gaussian_fit(lam), color = 'r')
plt.xlim(6300, 6700)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux ({})'.format(units_flux))
plt.show()
'''

# adjusting the date range
gaussian_model = models.Gaussian1D(1, 6563, 10)
fitter = fitting.LevMarLSQFitter()
gaussian_fit = fitter(gaussian_model, lam, flux)

plt.figure(figsize = (8, 5))
plt.plot(lam, flux, color = 'k')
plt.plot(lam, gaussian_fit(lam), color = 'r')
plt.xlim(5500, 7000)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux ({})'.format(units_flux))
plt.show()




