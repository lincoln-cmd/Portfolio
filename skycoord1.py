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

print(ngc188_center.ra.to_string(unit = u.hourangle, sep = ':', pad = True))

