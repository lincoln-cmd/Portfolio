# Coords 1. : Getting STarted with astropy.coordinates

# from : https://learn.astropy.org/rst-tutorials/Coordinates-Intro.html?highlight=filtertutorials

'''
 Investigate the area of sky around the picturesque grop of galaxies names 'Hickson Compact Group 7', download and image,
and do something with its coordinagtes.
'''

from urllib.parse import urlencode
from urllib.request import urlretrieve

from astropy import units as u
from astropy.coordinates import SkyCoord
from IPython.display import Image

'''
 - Describing on-sky locations with 'coordinates'
 'SkyCoord' class package is used to represent celestial coordinates.
 Make the SkyCoord object based on object's name, 'Hickson Compact Group 7' or 'HCG 7'.
'SESAME',a service which queries Simbad, NED, and VizieR and returns the object's type and its J2000 position, offers the most astronomical object names.
This service can be used via the 'SkyCoord.from_name()'.

 - about SESAME : http://cdsweb.u-strasbg.fr/cgi-bin/Sesame
 - about SkyCoord.from_name : https://julien.danjou.info/blog/2013/guide-python-static-class-abstract-methods
'''

# initialize a SkyCood object named hcg7_center at the location of HCG 7
hcg7_center = SkyCoord.from_name('HCG 7')
# this object is made by the information which comes from the Internet, so above command requires an Internet connection.
# If not, following this : 
    # hcg7_center = SkyCoord(9.81625*u.deg, 0.88806*u.deg, frame = 'icrs')

print(type(hcg7_center))

print(dir(hcg7_center))

print(hcg7_center.ra, hcg7_center.dec)
print(hcg7_center.ra.hour, hcg7_center.dec)
# HCG7 is located at ra = 9.849 degree and dec = 0.878 degree
# The ra and dec attributes are specialized Quantitiy objects, which, actually, a subclass called 'Angle', which
#in turn is subclassed by 'Latitude' and 'Longitude'.

print(type(hcg7_center.ra), type(hcg7_center.dec))
print(hcg7_center.ra, hcg7_center.dec)
print(hcg7_center) # -> longitude and latitude
print(hcg7_center.ra.hour)

# SkyCoord accepts string-fromatted coordinates either as separate strings for RA/Dec or as single string
# If they are not part of the string itself, it needs to be given the units
# i.e)
print(SkyCoord('0h39m15.9s', '0d53m17.016s', frame = 'icrs'))

'''
 - Download an image
 Download the image as an infromation to utilize the object to access data from 'Sloan Digitial Sky Survey(SDSS)'
If there is not available the Internet connection, employ the local file 'HCG7_SDSS_cutout.jpg'

 - about Sloan Digitial Sky Survey(SDSS) : https://www.sdss.org/
'''
# define how big of a cutout is needed
im_size = 12 * u.arcmin # get a 12 arcmin square
im_pixels = 1024
cutoutbaseurl = 'http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx'
query_string = urlencode(dict(ra = hcg7_center.ra.deg, dec = hcg7_center.dec.deg,
                              width = im_pixels, height = im_pixels,
                              scale = im_size.to(u.arcsec).value / im_pixels))
url = cutoutbaseurl + '?' + query_string

# download the image in our disks
urlretrieve(url, 'HCG7_SDSS_cutout.jpg')

Image('HCG7_SDSS_cutout.jpg')