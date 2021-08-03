# -*- coding: utf-8 -*-
'''
 Exercise
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

ngc188_center = SkyCoord.from_name('NGC 188')
print(ngc188_center.to_string(style = 'hmsdms', sep = ':', precision = 1))



