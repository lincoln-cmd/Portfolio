# Astrophysics

# modules : Matplotlib, IPython, NetworkX, PyGraphviz, pyfits, asciitable, ATpy, pidly ...
# from : https://learn.astropy.org/rst-tutorials/redshift-plot.html?highlight=filtertutorials#exercise

'''
import astropy
from astropy import cosmology
import matplotlib.pyplot as plt
%matplotlib inline

from IPython.display import Image
Image(filename = 'ang_dist.png', width = 500)
'''

# plot with both redshift and universe age

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

cosmo = FlatLambdaCDM(H0 = 70*u.km/u.s/u.Mpc, Om0 = 0.3)

import numpy as np
zvals = np.arange(0, 6, 0.1)
dist = cosmo.angular_diameter_distance(zvals)

fig = plt.figure(figsize = (6, 4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist)

dist.unit

ages = np.array([13, 10, 8, 6, 5, 4, 3, 2, 1.5, 1.2, 1])*u.Gyr

from astropy.cosmology import z_at_value
ageticks = [z_at_value(cosmo.age, age) for age in ages]

fig = plt.figure(figsize = (6, 4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist)
ax2 = ax.twiny()
ax2.set_xticks(ageticks)