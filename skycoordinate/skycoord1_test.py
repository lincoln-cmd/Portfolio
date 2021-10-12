# -*- coding: utf-8 -*-
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

a0620_center = SkyCoord.from_name('A0620')
#print('test : ', a0620_center)

a0620_center = SkyCoord('08h05m43.20s', '45d40m58.0s', frame = 'icrs')
#print(a0620_center)

a0620_center = SkyCoord('08:05:43.20', '45:40:58.0', unit = (u.hour, u.deg), frame = 'icrs')
print(a0620_center)

print('ra : ', a0620_center.ra, '/ dec : ', a0620_center.dec)
print('ra(hour) : ', a0620_center.ra.to(u.hourangle), '/ ra(rad) : ', a0620_center.ra.to(u.radian), '/ ra(deg) : ', a0620_center.ra.to(u.degree))
print('ra(hour) : ', a0620_center.ra.to_string(unit = u.hourangle, sep = ':', pad = True), '/ ra(deg) : ', a0620_center.ra.to_string(unit = u.degree, sep = ':', pad = True))

'''
a0620_table = QTable.read('gaia_results.fits')
print(len(a0620_table))

a0620_gaia_coords = SkyCoord(a0620_table['ra'], a0620_table['dec'])
print(a0620_gaia_coords)
print(len(a0620_gaia_coords))
#print(type(a0620_gaia_coords))
'''