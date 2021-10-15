# -*- coding: utf-8 -*-

# https://learn.astropy.org/rst-tutorials/3-Coordinates-Velocities.html?highlight=filtertutorials

'''
 - Python for Astronomers: https://prappleizer.github.io
 - Python for Astronomers(UC Berkeley): http://ugastro.berkeley.edu/pydecal/textbook.pdf
 
 
 
 
'''

import warnings

import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord, Distance, Galactic
import astropy.coordinates as coord
from astropy.io import fits
from astropy.table import QTable
from astropy.time import Time
from astropy.utils.data import download_file
from astropy.wcs import WCS

from astroquery.gaia import Gaia

'''
 - working with velocities in Astropy.coordinates : https://docs.astropy.org/en/latest/coordinates/velocities.html
 
 
'''

# passing velocity data into skycoord

SkyCoord(ra = 10*u.deg, dec = 10 * u.deg, pm_ra_cosdec = 1 * u.mas/ u.yr, pm_dec = 2 * u.mas / u.yr)

SkyCoord(ra = np.linspace(0, 10, 5) * u.deg, dec = np.linspace(5, 20, 5) * u.deg, pm_ra_cosdec = np.linspace(-5, 5, 5) * u.mas / u.yr, pm_dec = np.linspace(-5, 5, 5) * u.mas / u.yr)

velocity_coord = SkyCoord(ra = 10 * u.deg, dec = 20 * u.deg, pm_ra_cosdec = 1 * u.mas / u.yr, pm_dec = 2 * u.mas / u.yr, radial_velocity = 100 * u.km / u.s)
print(velocity_coord)
print(velocity_coord.pm_ra_cosdec)
print(velocity_coord.radial_velocity)

velocity_coord_gal = velocity_coord.transform_to(Galactic)
print(velocity_coord_gal)
# ra -> l, dec -> b, pm_ra_cosdec -> pm_l_cosb, pm_dec -> pm_b
print(velocity_coord_gal.pm_l_cosb)
print(velocity_coord_gal.pm_b)






