# Computing Galactic Orbits of Stars with Gala

# from : https://learn.astropy.org/rst-tutorials/gaia-galactic-orbits.html?highlight=filtertutorials

'''
 utilize the data from the 'Gaia mission' to get sky positions, distances(parallaxes), proper motions, and radial
velocities for a set of stars that are close to the Sun. Then ransform these observes, heliocentric
kinematic measurements to Galactocentric Cartesian coordinates and employ the positions and velocities as initial conditions
to compute the orbits of these stars in the galaxy using the 'gala' Python package.


 - about Gaia mission : https://www.cosmos.esa.int/web/gaia
 - about gala package : http://gala.adrian.pw/en/latest/
 
 - Installing Dependencies
  : pip install astro-gala astroquery
'''

import astropy.coordinates as coordinates
from astropy.table import QTable
import astropy.units as u
from astroquery.gaia import Gaia

# Third-party imports
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
%matplotlib inline

# gala imports
import gala.coordinates as gc
import gala.dynamics as gd
import gala.potential as gp
from gala.units import galactic

# conda-forge / astro-gala and astropy / astro-gala packages only offer the services to linux-64 and osx-64 platforms
