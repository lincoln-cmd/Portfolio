# Modeling2: create a user defined model using astropy.modeling

# from : https://learn.astropy.org/rst-tutorials/User-Defined-Model.html?highlight=filtertutorials


'''
 - define a models with compound model and with custom model

'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
from astropy.modeling import Fittable1DModel, Parameter
from astroquery.sdss import SDSS


'''
 - Fit an emission line in a stellar spectrum
 The low mass stars do not behave opposite with the higher mass stars. For instance, they stay magnetically active longer than
higher mass stars.
 H-alpha is specific deep-red visible spectral line in the Balmer series which have wavelength of 656.28[nm] in air.
This occurs when a hydrogen electron falls from its third to second lowest energy level.
 Download the information of spectrum from the SDSS's website or utilize the SDSSClass from astroquery.sdss



 - about 'H-alpha' : https://en.wikipedia.org/wiki/H-alpha
 - about 'Emission nebula' : https://en.wikipedia.org/wiki/Emission_nebula
 - about 'H II region(전리수소영역)' : https://en.wikipedia.org/wiki/H_II_region
 - SDSS database(searching the spectra) : https://dr12.sdss.org/basicSpectra
     information of the 'M dwarf' -> Plate : 1349, Fiber : 216, MJD(Modified Julian Date) : 52797
 - about 'SDSS Queries' : https://astroquery.readthedocs.io/en/latest/sdss/sdss.html#module-astroquery.sdss
 - SDSS document's tutorial : https://www.sdss.org/dr12/tutorials/quicklook/#python
'''
# get the information of the spectrum of the 'M dwarf'
spectrum = SDSS.get_spectra(plate = 1349, fiberID = 216, mjd = 52797)[0]

# read the information of the inside of the fits file spectrum
print(spectrum[1].columns)

flux = spectrum[1].data['flux']
lam = 10**(spectrum[1].data['loglam']) # modify the value to eleminate the log scale value

# to get the unit of the flux and wavelength, employ the 'fitsfile[0].header'.
units_flux = spectrum[0].header['bunit'] # 'bunit' or 'BUNIT' contain the physical units of the array values
print(units_flux)

units_wavelength_full = spectrum[0].header['WAT1_001']
print(units_wavelength_full)

# select only the characters of the unit : Angstroms
units_wavelength = units_wavelength_full[36:]
print(units_wavelength)

'''
# plot the spectrum graph
plt.plot(lam, flux, color = 'k')
plt.xlim(6300, 6700)
plt.axvline(x = 6563, linestyle = '--')
plt.xlabel('Wavelength ({})'.format(units_wavelength))
plt.ylabel('Flux ({})'.format(units_flux))
plt.show()
'''

'''
 - Fit an Emission Line with a Gaussian Model
 The blue dashed line marks the Hα emission line.
 Measure the height of this line by utilizingthe 'astropy.modeling' to fit a gaussian to the Hα line.
Firstly, initialize a gaussian model at the position of the Hα.
'''
gaussian_model = models.Gaussian1D(1, 6563, 10)
fitter = fitting.LevMarLSQFitter()
gaussian_fit = fitter(gaussian_model, lam, flux)

plt.figure(figsize = (8, 5))
plt.plot(lam, flux, color = 'k')
plt.plot(lam, gaussian_fit(lam), color = 'darkorange')
plt.xlim(6300, 6700)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux ({})'.format(units_flux))
plt.show()

print(gaussian_fit)


'''
 - Compound models
 Merely utilizing one model in above plot to measure the height of the Hα line is not enough to make this fit work, so
the combination of the models is required, which is the compound model in astropy. This allows to add, divide or multiply models that
already exist in 'astropy.modeling'.
 Combine the gaussian model with a polynomial one of degree 1 to account for the background spectrum close to the Hα line.
 
'''
# make the compound model
compound_model = models.Gaussian1D(1, 6563, 10) + models.Polynomial1D(degree = 1)
fitter = fitting.LevMarLSQFitter()
compound_fit = fitter(compound_model, lam, flux)

plt.figure(figsize = (8, 5))
plt.plot(lam, flux, color = 'k')
plt.plot(lam, compound_fit(lam), color = 'darkorange')
plt.xlim(6300, 6700)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux ({})'.format(units_flux))
plt.show()

print(compound_fit)

for x, y in zip(compound_fit.param_names, compound_fit.parameters): print(x, y)
# this shows the not only the fit parameters from the gaussian(mean, std, and amplitude), but also shows the two coefficients from the polynomial of degree 1(c0_1, c0_2)

print(compound_fit.amplitude_0)

'''
 - about fixed parameters : https://docs.astropy.org/en/stable/api/astropy.modeling.Parameter.html#astropy.modeling.Parameter.fixed
'''

compound_model_fixed = models.Gaussian1D(1, 6563, 10) + models.Polynomial1D(degree = 1)
compound_model_fixed.mean_0.fixed = True # utilizing the fixed parameter to fit the data
fitter = fitting.LevMarLSQFitter()
compound_fit_fixed = fitter(compound_model_fixed, lam, flux)


plt.figure(figsize = (8, 5))
plt.plot(lam, flux, color = 'k')
plt.plot(lam, compound_fit_fixed(lam), color = 'darkorange')
plt.xlim(6300, 6700)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux ({})'.format(units_flux))
plt.show()
# in the above plot, the height of the fit does not match the Hα line's height. This is because the too strict with the mean value cannot make good fit. However,
#the location of the height of the Hα line is matched.

'''
 - about minimum and maximum value : https://docs.astropy.org/en/stable/api/astropy.modeling.Parameter.html#astropy.modeling.Parameter.max
'''

print(compound_fit_fixed)

# let's loosen this condition a little

compound_model_bounded = models.Gaussian1D(1, 6563, 10) + models.Polynomial1D(degree = 1)
delta = 0.5
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

print(compound_fit_bounded)

'''
 - Custom model
 The Astropy provides the another tool to create the model when we face the models astropy.modeling does not offer.
 
 - about custom model : https://docs.astropy.org/en/stable/api/astropy.modeling.custom_model.html
'''
# basic custom model (with regard to the exponential model)
x1 = np.linspace(0, 10, 100)

a = 3
b = -2
c = 0
y1 = a * np.exp(b * x1 + c)
y1 += np.random.normal(0., 0.2, x1.shape)
y1_err = np.ones(x1.shape) * 0.2

plt.errorbar(x1, y1, yerr = y1_err, fmt = '.')
plt.show()

@custom_model
def exponential(x, a = 1., b = 1., c = 1.):
    '''
    f(x) = a * exp(b * x + c)
    '''
    return a * np.exp(b * x + c)

exp_model = exponential(1., -1., 1.)
fitter = fitting.LevMarLSQFitter()
exp_fit = fitter(exp_model, x1, y1, weights = 1.0 / y1_err**2)

plt.errorbar(x1, y1, yerr = y1_err, fmt = '.')
plt.plot(x1, exp_fit(x1))
plt.show()

print(exp_fit)

'''
 - about reduced chi-squared statistic(Wikipedia) : https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
 - about chi-squared distribution : https://ko.wikipedia.org/wiki/%EC%B9%B4%EC%9D%B4%EC%A0%9C%EA%B3%B1_%EB%B6%84%ED%8F%AC
'''
# check the parameters and the reduced chi-squared value, which provides the information about the goodness of the fit
def calc_reduced_chi_square(fit, x, y, yerr, N, n_free):
    '''
    fit(array) : values for the fit
    x, y, yerr (array) : data
    N : total number of points
    n_free : the number of parameters we are fitting
    '''
    return 1.0 / (N - n_free) * sum(((fit - y) / yerr)**2)

print('the goodness of the fit : {}'.format(calc_reduced_chi_square(exp_fit(x1), x1, y1, y1_err, len(x1), 3)))
# !!! The fits of non-linear parameters are extremely dependent on initial conditions
# the closer the value of the Reduced Chi-Squared value to 1, the higher the fit's accuracy





