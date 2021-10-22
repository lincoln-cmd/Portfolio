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

'''
# cannot transform with only sky position and proper motion to frame with a positional or velocity
test_coord = SkyCoord(ra = 10 * u.deg, dec = 20 * u.deg, pm_ra_cosdec = 1 * u.mas / u.yr, pm_dec = 2 * u.mas / u.yr)

print(test_coord.transform_to(coord.Galactocentric()))
'''


# query the Gaia catalog
# connected the Internet

# more info: https://www.cosmos.esa.int/web/gaia/dr2
# http://vizier.u-strasbg.fr/viz-bin/VizieR-2?-source=I%2F345%2Fgaia2&-c=349.72687648%2B5.40561039&-c.rs=720.0
# https://github.com/cds-astro/cds.cdsclient

gaia_tbl = Gaia.query_object(SkyCoord.from_name('HD 219829'), radius = 1 * u.arcmin)
'''
with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)
    
    gaia_tbl = QTable.read('HD_219829_query_results.ecsv')
'''


hd219829_row = gaia_tbl[gaia_tbl['phot_g_mean_mag'].argmin()]
print(hd219829_row['source_id', 'pmra', 'pmdec'])

#print(hd219829_row['parallax'])
#print(type(hd219829_row['parallax']))

#hd219829_coord = SkyCoord(ra = hd219829_row['ra'], dec = hd219829_row['dec'], distance = Distance(parallax = hd219829_row['parallax']), pm_ra_cosdec = hd219829_row['pmra'], pm_dec = hd219829_row['pmdec'], obstime = Time(hd219829_row['ref_epoch'], format = 'jyear'))
hd219829_coord = SkyCoord(
    ra=hd219829_row['ra'],
    dec=hd219829_row['dec'],
    distance=Distance(parallax=hd219829_row['parallax']),
    pm_ra_cosdec=hd219829_row['pmra'],
    pm_dec=hd219829_row['pmdec'],
    obstime=Time(hd219829_row['ref_epoch'], format='jyear'))
print(hd219829_coord)


dss_cutout_filename = download_file( f'http://archive.stsci.edu/cgi-bin/dss_search?'
                                    f'f=FITS&ra={hd219829_coord.ra.degree}&dec={hd219829_coord.dec.degree}'
                                    f'&width=4&height=4')
dss_cutout_filename = 'dss_hd219829.fits'

hdu = fits.open(dss_cutout_filename)[0]
wcs = WCS(hdu.header)

fig, ax = plt.subplots(1, 1, figsize = (8, 8), subplot_kw = dict(projection=wcs))
ax.imshow(hdu.data, origin='lower', cmap = 'Greys_r')
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.set_autoscale_on(False)

ax.scatter(hd219829_coord.ra.degree, hd219829_coord.dec.degree, s = 500, transform = ax.get_transform('world'), facecolor = 'none', linewidth = 2, color = 'tab:red')


print(hd219829_coord.obstime)

with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)
    
    hd219829_coord_1950 = hd219829_coord.apply_space_motion(new_obstime = Time('J1950'))
    
fig, ax = plt.subplots(1, 1, figsize = (8, 8), subplot_kw = dict(projection = wcs))

ax.imshow(hdu.data, origin = 'lower', cmap = 'Greys_r')
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.set_autoscale_on(False)

ax.scatter(hd219829_coord.ra.degree, hd219829_coord.dec.degree, s = 500, transform = ax.get_transform('world'), facecolor = 'none', linewidth = 2, color = 'tab:red')

ax.scatter(hd219829_coord_1950.ra.degree, hd219829_coord_1950.dec.degree, s = 500, transform = ax.get_transform('world'), facecolor = 'none', linewidth = 2, color = 'tab:blue')




