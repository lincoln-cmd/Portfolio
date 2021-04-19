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



'''
 - Modify the value of delta to change the minimum and maximum values for the mean of the gaussian.
 
'''
compound_model_bounded = models.Gaussian1D(1, 6563, 10) + models.Polynomial1D(degree = 1)
delta = 4e100
compound_model_bounded.mean_0.max = 6563 + delta
compound_model_bounded.mean_0.min = 6563 - delta

fitter = fitting.LevMarLSQFitter()
compound_fit_bounded = fitter(compound_model_bounded, lam, flux)

plt.figure(figsize = (8, 5))
plt.plot(lam, flux, color = 'k')
plt.plot(lam, compound_fit_bounded(lam), color = 'darkorange')
plt.xlim(6300, 6700)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux ({})'.format(units_flux))
plt.show()


'''
 - Modify the initial conditions of the fit and check the relation between the best fit parameters
and the initial condition for the previous example. It is possible by utilizing the Reduced Chi-Squared equation and function.
If the result value is closer to the 1, the fit is better.

'''
# test1
x1 = np.linspace(0, 10, 100)

a = 10
b = -4
c = 3
y1 = a * np.cosh(x1 + np.tan(b) + c**2)
y1 += np.random.normal(0., 2., x1.shape)
y1_err = np.ones(x1.shape) ** 0.3336

plt.errorbar(x1, y1, yerr = y1_err, fmt = '.')
plt.show()

@custom_model
def custom_test(x, a = 5., b = -2., c = 2):
    return a * np.cosh(x1 + np.tan(b) + c**2)

test_model = custom_test(10., -4., 3.)
fitter = fitting.LevMarLSQFitter()
test_fit = fitter(test_model, x1, y1, weights = 1.0 / y1_err**2)

plt.errorbar(x1, y1, yerr = y1_err, fmt = '.')
plt.plot(x1, test_fit(x1))
plt.show()

# test2

x2 = np.linspace(0, 10, 100)

a = 10
b = -4
c = 3
y2 = a**4 * np.sqrt(np.cos(a * np.tanh(x2) + c))
y2_err = np.ones(x2.shape) ** 0.8

plt.errorbar(x2, y2, yerr = y2_err, fmt = '.')
plt.show()

def custom_test2(x, a = 10., b = -4., c = 3.):
    return a**4 * np.sqrt(np.cos(a * np.tanh(x2) + c))

test_model_2 = custom_test2(10., -4., 3.)
fitter = fitting.LevMarLSQFitter()
test_fit_2 = fitter(test_model_2, x2, y2, weights = 1.0 / y2_err**2)

plt.errorbar(x2, y2, yerr = y2_err, fmt = '.')
plt.plot(x2, test_model_2(x2))
plt.show()

# test3

x3 = np.linspace(0, 10, 100)

a = 10
b = -4
c = 3
y3 = a * (1.0 / ((np.e**x3 + e**(-x3)) / 2))
y2_err = np.ones(x3.shape) ** 0.2

def sech(x, a = 10., b = -4., c = 3.):
    return a * (1.0 / ((np.e**x + e**(-x)) / 2))

