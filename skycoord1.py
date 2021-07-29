# -*- coding: utf-8 -*-

# https://learn.astropy.org/rst-tutorials/1-Coordinates-Intro.html?highlight=filtertutorials#exercises

'''
 The astropy.coordinates provides the services which trnasforms the differnet systems of coordinates and
represents the coordinates of objects.

'''

import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.io import fits
from astropy.table import QTable
from astropy.utils.data import download_file

from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = 10000

'''
 - about SkyCoord : https://docs.astropy.org/en/stable/coordinates/skycoord.html
 - about 'NGC 188' : https://en.wikipedia.org/wiki/NGC_188





'''

ngc188_center = SkyCoord(12.11 * u.deg, 85.26*u.deg)
print(ngc188_center)





