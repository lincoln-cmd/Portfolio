# Coords2 : Transforming between coordinate systems

# from : https://learn.astropy.org/rst-tutorials/Coordinates-Transform.html?highlight=filtertutorials

'''
 By utilizing the 'astropy.coordinates', 'frame' classes, define the astronomical coordinates,
then transform between the different built-in coordinate frames such as from ICRS(RA, Dec) to Galactic(l, b). Finally, compute
altitude and azimuth from a specific observing site.

    - the astroplan affiliated package : https://astroplan.readthedocs.io/en/latest/
'''

from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

# set up matplotlib and use a nicer set of plot parameters
from astropy.visualization import astropy_mpl_style
import matplotlib.pyplot as plt
plt.style.use(astropy_mpl_style)
%matplotlib inline


'''
 - Quickstart
 make the ICRS frame by passing the degrees explicitly or by passing in strings.
'''

hcg7_center = SkyCoord(9.81625 * u.deg, 0.88806 * u.deg, frame = 'icrs')

# hcg7_center = Skycoord('0h39m15.9s', '0d53m17.016s', frame = 'icrs')
# above objects are equivalent

print(hcg7_center.ra)
print(hcg7_center.dec)

'''
 - Introducing frame transformations
 All coordinates in Astropy are in particular 'frames'.
Transform center of HCG 7 from ICRS to Galactic coordinates

 - about 'astropy coordinates scheme' : http://astropy.readthedocs.org/en/latest/coordinates/index.html#overview-of-astropy-coordinates-concepts
'''

# transforming coordinates using attributes
print(hcg7_center.galactic)

# using the transform_to() method and other coordinate object
# It is possible to pass in an empty coordinate class to specify what coordinate system to transform into
from astropy.coordinates import Galactic
print(hcg7_center.transform_to(Galactic()))

# employing the transform_to() method and as string
# string is the name of a built-in coordinate system
print(hcg7_center.transform_to('galactic'))
# this method offers plenty of options to transform the coordinates frames and equinoxes

# alternate to FK5
hcg7_center_fk5 = hcg7_center.transform_to('fk5')
print(hcg7_center_fk5)

'''
! almost of the default values have different units systems.
i.e.) Galactic coordinate has (l, b)
'''

print(hcg7_center.galactic.l, hcg7_center.galactic.b)

'''
 - Trnasform frames to get to altitude-azimuth('AltAz')
 For observability, it is essential requesite to convert to a frame local to an on-earth observer.
horizontal altitude-azimuth(AltAz) is general choice.
'''
# firstly, specify both where and when we want to observe
from astropy.coordinates import EarthLocation
from astropy.time import Time

# Kitt Peak, Arizona
kitt_peak = EarthLocation(lat = '31d57.5m', lon = '-111d35.8m', height = 2096 * u.m)
kitt_peak = EarthLocation.of_site('Kitt Peak')
#print(EarthLocation.get_site_names()) # list of the observing points

observing_time = Time('2010-12-21 1:00')

# 'AltAz' frame has some other information about the atmosphere, which can be used to correct for atmospheric refraction.
from astropy.coordinates import AltAz

aa = AltAz(location = kitt_peak, obstime = observing_time)
print(aa)

# transform ICRS SkyCoord to AltAz to get the location in the sky over Kitt Peak at the requested time
print(hcg7_center.transform_to(aa))

print(hcg7_center.transform_to(aa).alt)


# creat the airmass plot
# this gives a Time object with an array of times
delta_hours = np.linspace(0, 6, 100) * u.hour
full_night_times = observing_time + delta_hours
full_night_aa_frames = AltAz(location = kitt_peak, obstime = full_night_times)
full_night_aa_coos = hcg7_center.transform_to(full_night_aa_frames)

'''
plt.plot(delta_hours, full_night_aa_coos.secz)
plt.xlabel('Hours from 6pm AZ time')
plt.ylabel('Airmass [Sec(z)]')
plt.ylim(0.9, 3)
plt.tight_layout()
'''
# utilize the 'get_sun' function for proper dark sky to observe
from astropy.coordinates import get_sun

full_night_sun_coos = get_sun(full_night_times).transform_to(full_night_aa_frames)
plt.plot(delta_hours, full_night_sun_coos.alt.deg)
plt.axhline(-18, color = 'k')
plt.xlabel('Hours from 6pm AZ time')
plt.ylabel('Sun altitude')
plt.tight_layout()

# object's altitude at the present time and date
now = Time.now()
hcg7_center = SkyCoord(9.81625 * u.deg, 0.88806 * u.deg, frame = 'icrs')
kitt_peak_aa = AltAz(location = kitt_peak, obstime = now)
print(hcg7_center.transform_to(kitt_peak_aa))