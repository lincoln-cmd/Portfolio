# Analyzing interstellar reddening and calculating synthetic photometry

# from : https://learn.astropy.org/rst-tutorials/color-excess.html?highlight=filtertutorials
# install the subsidary modules synphot, dus_extinction, astroquery
import matplotlib.pyplot as plt
%matplotlib inline

import numpy as np
import astropy.units as u
from astropy.table import Table
from dust_extinction.parameter_averages import CCM89, F99
from synphot import units, config
from synphot import SourceSpectrum, SpectralElement, Observation, ExtinctionModel1D
from synphot.models import BlackBodyNorm1D
from synphot.spectrum import BaseUnitlessSpectrum
from synphot.reddening import ExtinctionCurve
from astroquery.simbad import Simbad
from astroquery.mast import Observations
import astropy.visualization

# 1. Investigate Extinction Models
# create wavelengths array
wav = np.arange(0.1, 3.0, 0.001)*u.micron

for model in [CCM89, F99]:
    for R in (2.0, 3.0, 4.0):
        # initialize the extinction model
        ext = model(Rv = R)
        plt.plot(1/wav, ext(wav), label = model.name + 'R = ' + str(R))
        
plt.xlabel('$\lambda^{-1}$ ($\mu$m$^{-1}$')
plt.ylabel('A($\lambda$) / A(V)')
plt.legend(loc = 'best')
plt.title('Some Extinction Laws')
plt.show()


# 2. Deredden a Spectrum
# fetch the archival by utilizing the astroquery
obsTable = Observations.query_object('HD 147933', radius = '1 arcsec')
obsTable_spec = obsTable[obsTable['dataproduct_type'] == 'spectrum']
obsTable_spec.pprint()

obsids = ['3000022829']
dataProductsByID = Observations.get_product_list(obsids)
manifest = Observations.download_products(dataProductsByID)

# read the downloaded files into an astropy table
t_lwr = Table.read('./mastDownload/IUE/lwr05639/lwr05639mxlo_vo.fits')
print(t_lwr)

wav_UV = t_lwr['WAVE'][0,].quantity
UVflux = t_lwr['FLUX'][0,].quantity

custom_query = Simbad()
custom_query.add_votable_fields('fluxdata(U)', 'fluxdata(B)', 'fluxdata(V)')
phot_table = custom_query.query_object('HD 147933')
Umag = phot_table['FLUX_U']
Bmag = phot_table['FLUX_B']
Vmag = phot_table['FLUX_V']

wav_U = 0.3660 * u.micron
zeroflux_U_nu = 1.81E-23 * u.Watt / (u.m*u.m*u.Hz)
wav_B = 0.4400 * u.micron
zeroflux_B_nu = 4.26E-23 * u.Watt / (u.m*u.m*u.Hz)
wav_V = 0.5530 * u.micron
zeroflux_V_nu = 3.64E-23 * u.Watt / (u.m*u.m*u.Hz)

zeroflux_U = zeroflux_U_nu.to(u.erg/u.AA/u.cm/u.cm/u.s,
                              equivalencies = u.spectral_density(wav_U))
zeroflux_B = zeroflux_B_nu.to(u.erg/u.AA/u.cm/u.cm/u.s,
                              equivalencies = u.spectral_density(wav_B))
zeroflux_V = zeroflux_V_nu.to(u.erg/u.AA/u.cm/u.cm/u.s,
                              equivalencies = u.spectral_density(wav_V))

# convert from photometry to flux using the definition of magnuitude
Uflux = zeroflux_U * 10.**(-0.4*Umag)
Bflux = zeroflux_B * 10.**(-0.4*Bmag)
Vflux = zeroflux_V * 10.**(-0.4*Vmag)

# drwaing the plot
astropy.visualization.quantity_support()

plt.plot(wav_UV, UVflux, 'm', label = 'UV')
plt.plot(wav_V, Vflux, 'ko', label = 'U, B, V')
plt.plot(wav_B, Bflux, 'ko')
plt.plot(wav_U, Uflux, 'ko')
plt.legend(loc = 'best')
plt.ylim(0, 3E-10)
plt.title('rho Oph')
plt.show()

# initialize the extinction model
Rv = 5.0 # usually around 3, but about 5 for this star
Ebv = 0.5
ext = F99(Rv = Rv)

# to extinguish(redden) a spectrum, multiply by the 'ext.extinguish' function
# to unextinguish (deredden), divide by the same 'ext.extinguish' function
plt.semilogy(wav_UV, UVflux, 'm', label = 'UV')
plt.semilogy(wav_V, Vflux, 'ko', label = 'U, B, V')
plt.semilogy(wav_B, Bflux, 'ko')
plt.semilogy(wav_U, Uflux, 'ko')

plt.semilogy(wav_UV, UVflux / ext.extinguish(wav_UV, Ebv = Ebv), 'b',
             label = 'dereddened: EBV = 0.5, RV = 5')
plt.semilogy(wav_V, Vflux / ext.extinguish(wav_V, Ebv = Ebv), 'ro',
             label = 'dereddened: EBV = 0.5, RV = 5')
plt.semilogy(wav_B, Bflux / ext.extinguish(wav_B, Ebv = Ebv), 'ro')
plt.semilogy(wav_U, Uflux / ext.extinguish(wav_U, Ebv = Ebv), 'ro')

plt.legend(loc = 'best')
plt.title('rho Oph')
plt.show()


# 3. Calculate Color Excess with 'synphot'
