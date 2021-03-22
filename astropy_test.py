from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

cosmo = FlatLambdaCDM(H0 = 70*u.km/u.s/u.Mpc, Om0 = 0.3)


import matplotlib.pyplot as plt
import numpy as np

# step1
zvals = np.arange(0, 6, 0.1)
dist = cosmo.angular_diameter_distance(zvals)

fig = plt.figure(figsize = (6, 4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist)

# step2
dist.unit

ages = np.array([13, 10, 8, 6, 5, 4, 3, 2, 1.5, 1.2, 1])*u.Gyr

from astropy.cosmology import z_at_value
ageticks = [z_at_value(cosmo.age, age) for age in ages]

from astropy.cosmology import Planck13
dist2 = Planck13.angular_diameter_distance(zvals)

fig = plt.figure(figsize = (6, 4))
ax = fig.add_subplot(111)
ax.plot(zvals, dist2, label = 'Planck 2013')
ax.plot(zvals, dist, label = '$h = 0.7,\ \Omega_M = 0.3,\ \Omega_\Lambda = 0.7$')
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

print(cosmo.shape())