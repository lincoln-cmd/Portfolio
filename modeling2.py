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

# plot the spectrum graph
plt.plot(lam, flux, color = 'k')
plt.xlim(6300, 6700)
plt.axvline(x = 6563, linestyle = '--')
plt.xlabel('Wavelength ({})'.format(units_wavelength))
plt.ylabel('Flux ({})'.format(units_flux))
plt.show()
