# -*- coding: utf-8 -*-
'''
 Exercise
'''

'''
 - PDF about parallax : https://arxiv.org/pdf/1507.02105.pdf
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

print(ngc188_gaia_coords.separation(ngc188_center))

ngc188_center_3d = SkyCoord(12.11 * u.deg, 85.26 * u.deg, distance = 1.96 * u.kpc)

parallax_snr = ngc188_table['parallax'] / ngc188_table['parallax_error']
ngc188_table_3d = ngc188_table[parallax_snr > 10]
print(len(ngc188_table_3d))

