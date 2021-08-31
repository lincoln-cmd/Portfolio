# -*- coding: utf-8 -*-

# https://learn.astropy.org/rst-tutorials/2-Coordinates-Transforms

import matplotlib as mpl
import matplotlib.pyplot as plt
%matplotlib inline
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






