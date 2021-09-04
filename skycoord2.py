# -*- coding: utf-8 -*-

# https://learn.astropy.org/rst-tutorials/2-Coordinates-Transforms

import matplotlib as mpl
import matplotlib.pyplot as plt
#%matplotlib inline
import numpy as np

from astropy import units as u
from astropy.coordinates import (SkyCoord, Distance, Galactic, EarthLocation, AltAz)

import astropy.coordinates as coord
from astropy.io import fits
from astropy.table import QTable
from astropy.time import Time
from astropy.utils.data import download_file

'''
 example of previous chapter : SkyCoord1
'''
c = SkyCoord(ra = 15.9935 * u.deg, dec = -10.52351344 * u.deg)
print(c.ra.hourangle)
print(c.to_string('hmsdms'))
print(c.dec.to_string(sep = ':', precision = 5))
'''
'''

print(c.represent_as('cartesian')) # (x, y, z)
print(c.represent_as(coord.CartesianRepresentation))

print([x for x in dir(coord) if x.endswith('Representation') and not x.startswith('Base')])

# exercise of other representations
'''
 about list of all coordinate representations : https://docs.astropy.org/en/stable/coordinates/representations.html#astropy-coordinates-representations
'''

represent_list = ['CartesianRepresentation', 'CylindricalRepresentation', 'PhysicsSphericalRepresentation', 'RadialRepresentation', 'SphericalRepresentation', 'UnitSphericalRepresentation']

#for i in represent_list:
    #print('{0} : {1}'.format(i, c.represent_as(coord. + i)))


#print(c.represent_as(coord.CylindricalRepresentation))
#print(c.represent_as(coord.PhysicsSphericalRepresentation))
#print(c.represent_as(coord.RadialRepresentation))
#print(c.represent_as(coord.SphericalRepresentation))
#print(c.represent_as(coord.UnitSphericalRepresentation))


# without paramerter 'distance' in the SkyCoord, return the dimensionless values, but with that, return the 3D positional units.

c2 = SkyCoord(ra = 15.9932 * u.deg, dec = -10.52351344 * u.deg, distance = 127.4 * u.pc)

represent_list2 = ['cartesian', 'unitspherical', 'radial', 'spherical', 'physicsspherical', 'cylindrical']

for i in represent_list2:
    print('{0} : {1}'.format(i, c2.represent_as(i)))

#print(c2.represent_as('cartesian'))
#print(c2.represent_as('cylindrical'))
#print(c2.represent_as('physicsspherical'))
#print(c2.represent_as('radial'))
#print(c2.represent_as('spherical'))
#print(c2.represent_as('unitspherical'))

c3 = SkyCoord(ra = 15.9932 * u.deg, dec = -10.52351344 * u.deg, distance = 127.4 * u.pc)
print(c3.representation_type)
print(c3)

c3.representation_type = coord.CylindricalRepresentation
#print(c3.representation_type)
print(c3)
print(c3.rho, c3.phi * u.deg, c3.z)

'''
 - about ICRS PDF : https://arxiv.org/pdf/astro-ph/0602086.pdf
'''

# Transforming betwwen coordinate frames

tb1 = QTable.read('C:/Users/Administrator/Desktop/donghun/Cantat-Gaudin-open-clusters.ecsv')
# the file path can be different
# file from : https://github.com/astropy/astropy-tutorials/blob/main/tutorials/notebooks/astropy-coordinates/Cantat-Gaudin-open-clusters.ecsv

open_cluster_c = SkyCoord(ra = tb1['ra'], dec = tb1['dec'], distance = tb1['distance'], frame = 'icrs')
print(len(open_cluster_c))

print(open_cluster_c[:4])

'''
 visualize the positions of all of these clusters
'''
def coordinates_aitoff_plot(coords):
    fig, ax = plt.subplots(figsize = (10, 4), subplot_kw = dict(projection = 'aitoff'))
    
    sph = coords.spherical
    cs = ax.scatter(-sph.lon.wrap_at(180 * u.deg).radian, sph.lat.radian, c = sph.distance.value)
    
    def fmt_func(x, pos):
        val = coord.Angle(-x * u.radian).wrap_at(360 * u.deg).degree
        return f'${val:.0f}' + r'^{\circ}$'
    
    ticker = mpl.ticker.FuncFormatter(fmt_func)
    ax.xaxis.set_major_formatter(ticker)
    
    ax.grid()
    
    cb = fig.colorbar(cs)
    cb.set_label('distance [pc]')
    
    return fig, ax

fig, ax = coordinates_aitoff_plot(open_cluster_c)
ax.set_xlabel('RA [deg]')
ax.set_ylabel('Dec [deg]')

# convert to the other coordinate system
open_cluster_gal = open_cluster_c.transform_to(Galactic())

open_cluster_gal = open_cluster_c.galactic

print(open_cluster_gal[:4])

print(open_cluster_gal.l[:4])
print(open_cluster_gal.b[:4])
print(open_cluster_gal.distance[:4])

print(len(open_cluster_gal.l), len(open_cluster_gal.b), len(open_cluster_gal.distance))
print(len(open_cluster_gal))

# plot the new coordinate in galactic coordinates
fig, ax = coordinates_aitoff_plot(open_cluster_gal)
ax.set_xlabel('Galactic longitude, $l$ [deg]')
ax.set_ylabel('Galactic latitude, $b$ [deg]')




