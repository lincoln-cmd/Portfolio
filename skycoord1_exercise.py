# -*- coding: utf-8 -*-
'''
 Exercise
'''

'''
 - PDF about parallax : https://arxiv.org/pdf/1507.02105.pdf

 - angular separation(distance) : https://en.wikipedia.org/wiki/Angular_distance

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


'''
 - estimating distances from parallaxes PDF : https://arxiv.org/pdf/1507.02105.pdf
'''
parallax_snr = ngc188_table['parallax'] / ngc188_table['parallax_error']
ngc188_table_3d = ngc188_table[parallax_snr > 10]
print(len(ngc188_table_3d))

print(Distance(parallax = 1 * u.mas))

gaia_dist = Distance(parallax = ngc188_table_3d['parallax'])

ngc188_coords_3d = SkyCoord(ra = ngc188_table_3d['ra'], dec = ngc188_table_3d['dec'], distance = gaia_dist)
print(ngc188_coords_3d)



fig, ax = plt.subplots(figsize = (6.5, 5.2), constrained_layout = True)
cs = ax.scatter(ngc188_coords_3d.ra.degree, ngc188_coords_3d.dec.degree, c = ngc188_coords_3d.distance.kpc,
                s = 5, vmin = 1.5, vmax = 2.5, cmap = 'twilight')
cb = fig.colorbar(cs)
cb.set_label(f'distance [{u.kpc:latex_inline}]')

ax.set_xlabel('RA [deg]')
ax.set_ylabel('Dec [deg]')

ax.set_title('Gaia DR2 sources near NGC 188', fontsize = 18)


sep3d = ngc188_coords_3d.separation_3d(ngc188_center_3d)
print(sep3d)



# exercis2
ngc188_3d_mask = sep3d < 50 * u.pc
print(ngc188_3d_mask.sum())








