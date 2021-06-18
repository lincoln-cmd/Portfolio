# Analyzing UVES Spectroscopy with Astropy

# from : https://learn.astropy.org/rst-tutorials/UVES.html?highlight=filtertutorials


'''
 As the accretion rate decreases, the impact on the central star must change, which means that
the accretion rate is depended on the center of the star. Furthermore, the emission frequence is alse
depended on the center of the star. (The accretion disk of protostar emits the infrared radiation, and
                                in the form of heat     the nuetron star and blackhole emit the X-ray)
from : https://ui.adsabs.harvard.edu/abs/2013ApJ...771...70G/abstract

- about accretion disk(ko) : https://ko.wikipedia.org/wiki/%EA%B0%95%EC%B0%A9%EC%9B%90%EB%B0%98
                            https://naverkpsdictionary.miraheze.org/wiki/%EA%B0%95%EC%B0%A9%EC%9B%90%EB%B0%98
- about accretion(eng) : https://en.wikipedia.org/wiki/Accretion_disk
- about T-tauri(ko) : https://ko.wikipedia.org/wiki/%ED%99%A9%EC%86%8C%EC%9E%90%EB%A6%AC_T%ED%98%95_%ED%95%AD%EC%84%B1
- about T-tauri(eng) : https://en.wikipedia.org/wiki/T_Tauri_star
- about protostar(ko) : https://ko.wikipedia.org/wiki/%EC%9B%90%EC%8B%9C%EB%B3%84
- acrretion disk pdf file : file:///C:/Users/Administrator/Downloads/review-accretion.pdf
- acrretion disk pdf file2 : http://www.mso.anu.edu.au/~geoff/HEA/11_Accretion_Disks_II.pdf
- jeans instability(ko) : https://ko.wikipedia.org/wiki/%EC%A7%84%EC%8A%A4_%EB%B6%88%EC%95%88%EC%A0%95%EC%84%B1#%EC%A7%84%EC%8A%A4_%EC%A7%88%EB%9F%89
- sort of main-sequence star : https://steamcommunity.com/sharedfiles/filedetails/?id=1460382055
- proton-proton chain reaction(ko) : https://ko.wikipedia.org/wiki/%EC%96%91%EC%84%B1%EC%9E%90-%EC%96%91%EC%84%B1%EC%9E%90_%EC%97%B0%EC%87%84_%EB%B0%98%EC%9D%91
- proton-proton chain reaction(eng) : https://en.wikipedia.org/wiki/Proton%E2%80%93proton_chain





- words of the mathematics and physics : https://blog.daum.net/williamockham/83

'''

import matplotlib.pyplot as plt
#from IPython.display import set_matplotlib_formats
# %matplotlib inline

# download tar files and extract the files' data
import tarfile
from astropy.utils.data import download_file
url = 'http://data.astropy.org/tutorials/UVES/data_UVES.tar.gz'
f = tarfile.open(download_file(url, cache = True), mode = 'r|*')
working_dir_path = 'C:/Users/Administrator/Desktop/donghun' # change the path
f.extractall(path = working_dir_path)

# analyze data from NM Lup, a T Tauri star in the Taurus-Auriga star forming region located at a distance of about 140pc
# MN Lup has been observed simultaneously with XMM-Newton and the UVES spectrograph on the VLT

'''
- Spatially resolving the accretion shocks on the rapidly-rotating M0 T-Tauri star MN Lupi
 The accretion mechanism is the cause of its rapid surface rotation because ongoing disk accretion.
 
from : https://ui.adsabs.harvard.edu/abs/2005A%26A...440.1105S/abstract
'''

# reading the data
from glob import glob

import os
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits

globpath = os.path.join(working_dir_path, 'UVES/*.fits') # os.path.join : join two directories

print(globpath)
filelist = glob(globpath) # search through directories similar to the Unix shell

filelist.sort()

# read the first FITS file in the list and check what is in there
sp = fits.open(filelist[0])
sp.info()

# extract the WCS from header to get the wavelength coordinate.
# Even though the warnings are invoked in this process about a non-standard RADECSYS, the WCS will still work
header = sp[0].header

wcs = WCS(header)
index = np.arange(header['NAXIS1']) # make the array

wavelength = wcs.wcs_pix2world(index[:, np.newaxis], 0)
wavelength.shape
wavelength = wavelength.flatten() # adjust the wrong dimension by utilizing flatten()

# the flux is contained in the primary image
flux = sp[0].data


# make the function to reuse the code
# input -> filename, returns -> wvaelength, flux arrays, and the time of the observation
#from spectra_utils import func
def read_spec(filename):
    ''' Read a UVES spectrum from the ESO pipeline
    
    Parameters
    -----------
    filename : string. name of the fits file with the data
    
    Returns
    -----------
    wavelength : np.ndarray. wavelength(in Ang)
    flux : np.ndarray. flux (in erg/s/cm**2)
    date_obs : string. time of observation
    '''
    sp = fits.open(filename)
    header = sp[0].header
    
    wcs = WCS(header)
    index = np.arange(header['NAXIS1']) # make index array
    
    wavelength = wcs.wcs_pix2world(index[:, np.newaxis], 0)
    wavelength = wavelength.flatten()
    flux = sp[0].data
    
    date_obs = header['Date-OBS']
    return wavelength, flux, date_obs

help(read_spec)


# The dataset of UVES spectra should have been taken using all the same setup
# EXPTIME : exposure time
# CRVAL1 : wavelength zero point
# HIERARCH ESO INS PATH : arm used(UVES has a red and a blue arm)
def read_setup(filename):
    '''
    Get setup for UVES spectrum from the ESO pipeline
    
    Parameters
    ----------
    filename : string
    name of the fits file with the data
    
    Returns
    ----------
    exposure_time : float
    wavelength_zero_point : float
    optical_arm : string
    '''
    
    sp = fits.open(filelist[0])
    header = sp[0].header
    
    return header['EXPTIME'], header['CRVAL1'], header['HIERARCH ESO INS PATH']

for f in filelist:
    print(read_setup(f))
 
'''
 The UVES pipeline that was used to reduce the data employs a fixed wavelength grid, so the wavelength is the same for all apectra,
which allows us to define an array that can hold the fluxes of all observations easily.
'''
# loop over the list of all filenames and fill this array with data
flux = np.zeros((len(filelist), len(wavelength)))
date = np.zeros((len(filelist)), dtype='U23') # data comes as string with 23 characters (dytpe = 'S23')

for i, fname in enumerate(filelist):
    w, f, date_obs = read_spec(fname)
    flux[i, :] = f
    date[i] = date_obs

import astropy.units as u
from astropy.constants.si import c, G, M_sun, R_sun

# all the constants have to have their own units, so the constant values have to be multiplied with units
wavelength = wavelength * u.AA

heliocentric = -23. * u.km/u.s
v_rad = -4.77 * u.km / u.s
R_MN_Lup = 0.9 * R_sun
M_MN_Lup = 0.6 * M_sun
vsini = 74.6 * u.km / u.s
period = 0.439 * u.day

inclination = 45. * u.degree

incl = inclination.to(u.radian) # convert the degree unit to radian unit

'''
 MN Lup is T Tauri star(TTS), which is possibly surrounded by an accertion disk, so we can expect those
accretion signatures to appear close to the free-fall velocity(v) that a mass(m) reaches, when it hits the stellar surface.

 -> extract infall speed with energy conservation equation
 0.5mv^2 = gmM / r
 (kenetic energy = gravitational energy), (gravitational energy equation = U(r) = - GmM / r)
'''

v_accr = (2. * G * M_MN_Lup / R_MN_Lup)**0.5
print(v_accr) # si system
print(v_accr.cgs) # cgs system
from astropy.units import imperial
print(v_accr.to(imperial.yd / u.hour)) # imperial system(english measuring system)

# the relation accretion velocity with rotational velocity
v_rot = vsini / np.sin(incl)
print((v_accr / v_rot).decompose()) # (v_accr / v_rot) is composed different units, m and km

'''
 λ_heliocentric = λ_bariocentric * (1 + v_helio / c)
-> correct the wavelength scale to the heliocentric velocity scale by utilizing the 'astropy.units'
'''
wavelength = wavelength * (1. + heliocentric / c) # all the constants have each different unit value, km, m, and just a number
wavelength = wavelength * (1. * u.dimensionless_unscaled + heliocentric / c)
print(wavelength.to(u.keV, equivalencies = u.spectral())) # 1eV = 1.60218 x 10^-19
print(wavelength.to(u.Hz, equivalencies = u.spectral()))

'''
 Spectroscopically, MN Lup is classified as spectral type M0 V, so the gravitational acceleration on the surface log(g) should be comparable to the sun.
log(g) can be computered with mass and radius values.
-> log(g) is consistent
'''
print(np.log10((G * M_MN_Lup / R_MN_Lup**2) / u.cm * u.second**2)) # cgs system


'''
# making the function to convert the wavelength scal into the velocity scale
 Convert a wavelength scale into a velocity scale. Input the wavelengths array and the rest wavelength of a spectral line.
This function can be employed to show the red and blueshift of the spectrum relative to the Ca II H line.
'''
waveclosetoHa = np.array([6562., 6563., 6565.]) * u.AA

# make the function to utilize the 'Doppler equivalency' between wavelength and velocity
import astropy.units as u
def wave2doppler(w, w0):
    w0_equiv = u.doppler_optical(w0)
    w_equiv = w.to(u.km/u.s, equivalencies = w0_equiv)
    return w_equiv

print(wave2doppler(waveclosetoHa, 656.489 * u.nm).to(u.km/u.s))


'''
 Make the function which has the wavelength array and rest wavelength of spectral line as input values and Doppler shift as return value.
This function subtracts the radial velocity of MN Lup(4.77 km/s) and expresses the resulting velocity in units of vsini.
Furthermore, this function can be employed to show the red- and blueshift of the spectrum relative to the Ca 2 H line.

def my_w2vsini(wavelength_array, wavelength_line):
    w1 = u.doppler_optical(wavelength_line)
    w2 = wavelength_array.to(u.km/u.s, equivalencies = w1)
    array_of_shift_in_vsini = w2 - 4.77 * u.km/u.s
    return array_of_shift_in_vsini / vsini

print('my test', end = ' : ')
print(my_w2vsini(waveclosetoHa, 656.489 * u.nm))
'''

def w2vsini(w, w0):
    v = wave2doppler(w, w0) - 4.77 * u.km/u.s
    return v / vsini


'''
 - Converting times
 Utilize the 'astropy.time' to convert time and dates between different systems and formats.
If unless the format is unambiguous, the format needs to be specified(i.e. a number colud mean JD, MJD, or year).
Moreover, the time system needs to be given (i.e. UTC).
 It can be just read the keywords in the time system. This is because the 'ESO FITS headers' already contain the time of the observation in different systems.
'''
# initialized from different header keywords
from astropy.time import Time
t1 = Time(header['MJD-Obs'], format = 'mjd', scale = 'utc')
t2 = Time(header['Date-Obs'], scale = 'utc')

# times can be expressed in different formats
print(t1)
t1
print(t1.isot) # scale = 'utc', format = 'isot'
t1.isot
print(t2)
t2

print(t1.tt) # convert to a different time system. scale = 'tt', format = 'mjd'
t1.tt

# times can be initialized from arrays, and the time differences can be calculated
obs_times = Time(date, scale = 'utc')
delta_t = obs_times - Time(date[0], scale = 'utc') # the unit of delta_t is day

# to express the time difference between the individual spectra of MN Lup in rotational periods, to begin with, the 'delta_t' object has to be converted
#This is because this object has days as an unit. Howevert, the 'astropy.time.Time' and 'astropy.units.Quantity' do not offer the service with that units,
#so this has to be converted from one to the other explicityly
delta_p = delta_t.value * u.day /period

'''
 - Normalize the flux to the local continuum
 To estimate the equivalent width or make reasonable plots, the flux should to be normalize to the local continuum.
In this case, the emission line is bright and the continuum can be described reasonably by a second-order polynomial.
 Define two parts of regions left and right of the emission line where are fit the polynomial.
'''
def region_around_line(w, flux, cont):
    '''
    cut out and normalize flux around a line
    
    Parameters
    ----------
    w : one-dimension np.nadarray. (array of wavelengths)
    flux : np.ndarray of shape (N, len(w)). array of flux values for different spectra in the series
    cont : list of lists. wavelengths for continuum normalization [[low1, up1], [low2, up2]] that described two areas on both sides of the line
    '''
    # index is true in the region where we fit the polynomial
    indcont = ((w > cont[0][0]) & (w < cont[0][1])) | ((w > cont[1][0]) & (w < cont[1][1]))
    # index of the region we want to return
    indrange = (w > cont[0][0]) & (w < cont[1][1])
    # make a flux array of shape
    # (number of spectra, number of points in indrange)
    f = np.zeros((flux.shape[0], indrange.sum()))
    for i in range(flux.shape[0]):
        # fit polynomial of second order to the continuum region
        linecoeff = np.polyfit(w[indcont], flux[i, indcont], 2)
        # divide the flux by the polynomial and put the result in our new flux array
        f[i, :] = flux[i, indrange] / np.polyval(linecoeff, w[indrange].value)
    return w[indrange], f

wcaII, fcaII = region_around_line(wavelength, flux, 
                                  [[3925 * u.AA, 3930 * u.AA], 
                                   [3938 * u.AA, 3945 * u.AA]])

# calculate the equivalent width in Angstroms of the emission line for the first spectrum
ew = fcaII[0, :] -1. # ?? how can subtract the number from array??
ew = ew[:-1] * np.diff(wcaII.to(u.AA).value) # ?? metrix complication or just complication??
print(ew.sum()) # not 19.37760714447237... -> 20.21238214515653


# process all spectra at once by employing numpy array notation
delta_lam = np.diff(wcaII.to(u.AA).value)
ew = np.sum((fcaII - 1.)[:, :-1] * delta_lam[np.newaxis, :], axis = 1)

# generate a LaTeX table of the observation times, period, and equivalent width which can directyl paste into manuscript.
# for that, firstly, collect all the columns and make an astropy.table.Table object
from astropy.table import Column, Table
from astropy.io import ascii

datecol = Column(name = 'Obs Date', data = date)
pcol = Column(name = 'phase', data = delta_p, format = '{:.1f}')
ewcol = Column(name = 'EW', data = ew, format = '{:.1f}', unit = '\\AA')
tab = Table((datecol, pcol, ewcol))
# latexdicts['AA'] contains the style specifics for A&A (\hline etc..)
tab.write(os.path.join(working_dir_path, 'EWtab.tex'), latexdict = ascii.latexdicts['AA'])
# tab.write(os.path.join(working_dir_path, 'Ewtab2.py'), latexdict = ascii.latexdicts['AA'])

'''
 - Plot
 The x-axis shows the Doppler shift expressed in units of the rotational velocity. In this process, the features that are rotationally modulated
will stick out between -1 and +1
'''
x = w2vsini(wcaII, 393.366 * u.nm).decompose()

# set reasonable figsize for 1-column figures
# this plot shows only a single spectrum
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(x, fcaII[0,:], marker = '', drawstyle = 'steps-mid')
ax.set_xlim([-3, +3])
ax.set_xlabel('line shift [v sin(i)]')
ax.set_ylabel('flux')
ax.set_title('Ca II H line in MN Lup')
plt.draw()
fig.savefig('Ca_II_H_line_in_MN_Lup.png', dpi = 200, bbox_inches = 'tight')

# make the plot which shows all spectra into a single plot and introduce a sensibel offset between them. To do so, following the time evolution of the line is essential
yshift = np.arange((fcaII.shape[0])) * 0.5
yshift[:] += 1.5 # shift the second night up by a little more
yshift[13:] += 1

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for i in range(25):
    ax.plot(x, fcaII[i, :] + yshift[i], 'k')

# separately show the mean line profile in a different color
ax.plot(x, np.mean(fcaII, axis = 0))
ax.set_xlim([-2.5, +2.5])
ax.set_xlabel('line shift [$v \\sin i$]')
ax.set_ylabel('flux')
ax.set_title('Ca II H line in MN Lup')
fig.subplots_adjust(bottom = 0.15)
plt.draw()

fmean = np.mean(fcaII, axis = 0)
fdiff = fcaII - fmean[np.newaxis, :]

# the axis scales are not right, the gap between both nights is not visible and there is no proper labeling
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
im = ax.imshow(fdiff, aspect = 'auto', origin = 'lower')

# plot the spectra from both nights separately and pass the extent keyword to ax.imshow which takes care of the axis
ind1 = delta_p < 1 * u.dimensionless_unscaled
ind2 = delta_p > 1 * u.dimensionless_unscaled

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for ind in [ind1, ind2]:
    im = ax.imshow(fdiff[ind, :], extent = (np.min(x), np.max(x), np.min(delta_p[ind]), np.max(delta_p[ind])), aspect = 'auto', origin = 'lower')
    
ax.set_ylim([np.min(delta_p), np.max(delta_p)])
ax.set_xlim([-1.9, 1.9])
plt.draw()


'''
 more enhanced plot than above one
 
-> introduce an offset on the y-axis to reduce the amount of white space
-> above plot has wrong scale strictly because the extent keyword gives the edges of the image shown, byt x and delta_p contain the bin mid-points
-> utilize the gray scale instead of color to save publication charges
-> add labels to the axis
'''
# shift a little for plotting purposes
pplot = delta_p.copy().value
pplot[ind2] -= 1.5
# image goes from x1 to x2, but really x1 should be middle of first pixel
delta_t = np.median(np.diff(delta_p)) / 2. # numpy.median : calculate the median along the specified axis and return the median of the array elements
delta_x = np.median(np.diff(x)) / 2.
# normalization for plotting by hand to ensure it goes -1 to +1, not merely employ the imshow
fdiff = fdiff / np.max(np.abs(fdiff))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for ind in [ind1, ind2]:
    im = ax.imshow(fdiff[ind, :],
                   extent = (np.min(x) - delta_x, np.max(x) + delta_x,
                             np.min(pplot[ind]) - delta_t, np.max(pplot[ind]) + delta_t),
                   aspect = 'auto', origin = 'lower', cmap = plt.cm.Greys_r)
    
ax.set_ylim([np.min(pplot) - delta_t, np.max(pplot) + delta_t])
ax.set_xlim([-1.9, 1.9])
ax.set_xlabel('vel in $v\\sin i$')
ax.xaxis.set_major_locator(plt.MaxNLocator(4))

def pplot(y, pos):
    #'the two args are the value and tick position'
    #'Function to make tick labels look good.'
    if y < 0.5:
        yreal = y
    else:
        yreal = y + 1.5
    return yreal

formatter = plt.FuncFormatter(pplot)
ax.yaxis.set_major_formatter(formatter)
ax.set_ylabel('period')
fig.subplots_adjust(left = 0.15, bottom = 0.15, right = 0.99, top = 0.99)
plt.draw()