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
 - about 'ICRS' : https://arxiv.org/abs/astro-ph/0602086
 - astronomical coordinate : https://en.wikipedia.org/wiki/Astronomical_coordinate_systems
 
 - search astronomical object : http://cdsweb.u-strasbg.fr/cgi-bin/Sesame
 - astronomical object : https://en.wikipedia.org/wiki/Astronomical_object
 



'''

ngc188_center = SkyCoord(12.11 * u.deg, 85.26*u.deg)
print(ngc188_center)

ngc188_center2 = SkyCoord(12.11 * u.deg, 85.26 * u.deg, frame = 'icrs')
print(ngc188_center2)

print(SkyCoord('00h48m26.4s', '85d15m36s', frame = 'icrs'))

print(SkyCoord('00:48:26.4 85:15:36', unit=(u.hour, u.deg), frame = 'icrs'))

ngc188_center = SkyCoord.from_name('NGC 188')
print(ngc188_center)

#new_center = SkyCoord.from_name('A0620')
#print('test : ', new_center)

print('ra : ', ngc188_center.ra, 'dec : ', ngc188_center.dec)

print('type of ra : ', type(ngc188_center.ra), 'type of dec : ', type(ngc188_center.dec))

print(ngc188_center.ra.to(u.hourangle), ngc188_center.ra.to(u.radian), ngc188_center.ra.to(u.degree))

print(ngc188_center.ra.hour, ngc188_center.ra.radian, ngc188_center.ra.degree)

'''
 formatting coordinate strings : https://docs.astropy.org/en/latest/coordinates/formatting.html
'''
print(ngc188_center.ra.to_string(unit = u.hourangle, sep = ':', pad = True))
# print(ngc188_center.ra.to_string(unit = u.radian, sep = ':', pad = True))
print(ngc188_center.ra.to_string(unit = u.degree, sep = ':', pad = True))

job = Gaia.cone_search_async(ngc188_center, radius = 0.5 * u.deg)
#ngc188_table = job.get_results()

#ngc188_table = ngc188_table[ngc188_table['phot_g_mean_mag'] < 19 * u.mag]


#cols = ['source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'radial_velocity', 'phot_g_mean_mag', 'phot_bp_mean_mag',
 #       'phot_rp_mean_mag']
#ngc188_table[cols].write('gaia_results.fits', overwrite = True)

#print(len(ngc188_table))


ngc188_table = QTable.read('gaia_results.fits')
print(len(ngc188_table))

print('ra : ', ngc188_table['ra'])
print('dec : ', ngc188_table['dec'])

ngc188_gaia_coords = SkyCoord(ngc188_table['ra'], ngc188_table['dec'])
print(ngc188_gaia_coords)


'''
 Exercise
'''
ngc188_center = SkyCoord.from_name('NGC 188')




