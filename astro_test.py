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

wav = np.arange(0.1, 3.0, 0.001)*u.micron
'''
for model in [CCM89, F99]:
    for R in (2.0, 3.0, 4.0):
        # initialize the extinction model
        ext = model(Rv = R)
        plt.plot(1/wav, ext(wav), label = model.name + 'R = ' + str(R))
        
plt.xlabel('$\lambda^{-1}$ ($\mu$m$^{-1}$)')
plt.ylabel('A($\lambda$) / A(V)')
plt.legend(loc = 'best')
plt.title('Some Extinction Laws')
plt.show()
'''
for model in [CCM89, F99]:
    print(model)