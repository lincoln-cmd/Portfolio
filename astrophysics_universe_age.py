# Astrophysics

# modules : Matplotlib, IPython, NetworkX, PyGraphviz, pyfits, asciitable, ATpy, pidly ...
# from : https://learn.astropy.org/rst-tutorials/redshift-plot.html?highlight=filtertutorials#exercise

'''
import astropy
from astropy import cosmology
import matplotlib.pyplot as plt
%matplotlib inline

from IPython.display import Image
Image(filename = 'ang_dist.png', width = 500)
'''

# plot with both redshift and universe age

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

# hubble parameter(H0) at z = 0 : 70km/s/Mpc
# redefine the cosmo variable
cosmo = FlatLambdaCDM(H0 = 70*u.km/u.s/u.Mpc, Om0 = 0.3)
# FlatLambdaCDM(H0 : Hubble constant(float or Quantity),
# Om0 : Omega matter : non-relativistic matter's density in units of the critical density at z = 0,
# Tcmb0 : Temperature of the CMB z = 0(float, scalar Quantity, optional),
# Neff : Effective number of Neutrino species. Default = 3.04(float, optional),
# m_nu : Mass of each neutrino species(Quantitiy, optional),
# Ob0 : Omega baryons(float, None, optional),
# name : Name for this cosmological object(str, optional))

import matplotlib.pyplot as plt
import numpy as np

# step1 : calculate the angular diameter distance for a range of redshifts
zvals = np.arange(0, 6, 0.1)
dist = cosmo.angular_diameter_distance(zvals)

fig = plt.figure(figsize = (6, 4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist)

# step2

dist.unit # check the units of the angular diameter distance[Mpc]

ages = np.array([13, 10, 8, 6, 5, 4, 3, 2, 1.5, 1.2, 1])*u.Gyr

from astropy.cosmology import z_at_value
# z_at_value : link the redshift and age axes.
ageticks = [z_at_value(cosmo.age, age) for age in ages]

fig = plt.figure(figsize = (6, 4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist)
ax2 = ax.twiny()
# ax.twiny : create a new axes with an invisible y-axis and an independent
# x-axis positioned opposite to the original one
ax2.set_xticks(ageticks)
# the top of x-axis is created by ax2 by utilizing ax.twiny()
# the botton of x-axis is originated from original plot

# step3
# modify the top axis at the correct ages, not redshift
fig = plt.figure(figsize = (6, 4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist)
ax2 = ax.twiny()
ax2.set_xticks(ageticks)
ax2.set_xticklabels(['{:g}'.format(age) for age in ages.value])

# step4
# make sure the top and bottom axes have the same sedshift limits
fig = plt.figure(figsize = (6, 4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist)
ax2 = ax.twiny()
ax2.set_xticks(ageticks)
ax2.set_xticklabels(['{:g}'.format(age) for age in ages.value])
zmin, zmax = 0.0, 5.9
ax.set_xlim(zmin, zmax)
ax2.set_xlim(zmin, zmax)

# step5
# set the labels and some minor ticks
fig = plt.figure(figsize = (6, 4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist)
ax2 = ax.twiny()
ax2.set_xticks(ageticks)
ax2.set_xticklabels(['{:g}'.format(age) for age in ages.value])
zmin, zmax = 0, 5.9
ax.set_xlim(zmin, zmax)
ax2.set_xlim(zmin, zmax)
ax2.set_xlabel('Time since Big Bang (Gyr)')
ax.set_xlabel('redshift')
ax.set_ylabel('Angular diameter distance(Mpc)')
ax.set_ylim(0, 1890)
ax.minorticks_on()

# final
# add the angular diameter distance for a different cosmology, from the Planck 2013
from astropy.cosmology import Planck13
dist2 = Planck13.angular_diameter_distance(zvals)

fig = plt.figure(figsize = (6, 4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist2, label = 'Planck 2013')
ax.plot(zvals, dist, label = '$h = 0.7,\ \Omega_M = 0.3,\ \Omega_\Lambda = 0.7$')
# ΩM + ΩA + Ωk = 1
# Ωk : curvature density = 0
# ΩM : matter density = 0.3
# ΩA : dark energy density = 0.7
'''
https://astro.kasi.re.kr/learning/pageView/6383
https://ko.wikipedia.org/wiki/%ED%94%84%EB%A6%AC%EB%93%9C%EB%A7%8C_%EB%B0%A9%EC%A0%95%EC%8B%9D - 프리드만 방정식
'''
ax.legend(frameon = 0, loc = 'lower right')
ax2 = ax.twiny()
ax2.set_xticks(ageticks)
ax2.set_xticklabels(['{:g}'.format(age) for age in ages.value])
zmin, zmax = 0.0, 5.9
ax.set_xlim(zmin, zmax)
ax2.set_xlim(zmin, zmax)
ax2.set_xlabel('Time sincs Big Bang (Gyr)')
ax.set_xlabel('Redshift')
ax.set_ylabel('Angular diameter distance (Mpc)')
ax.minorticks_on()
ax.set_ylim(0, 1890)
fig.savefig('ang_dist.png', dpi = 200, bbox_inches = 'tight')