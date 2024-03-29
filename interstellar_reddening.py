# Analyzing interstellar reddening and calculating synthetic photometry

# from : https://learn.astropy.org/rst-tutorials/color-excess.html?highlight=filtertutorials
# install the subsidary modules synphot, dus_extinction, astroquery

import matplotlib.pyplot as plt
%matplotlib inline

'''
 - about Synthetic Photometry(synphot) document : https://synphot.readthedocs.io/en/latest/#installation-and-setup
 - about Bessell & Murphu : https://ui.adsabs.harvard.edu/abs/2012PASP..124..140B/abstract
 - 천문학 및 천체물리학 15장(pdf) : file:///C:/Users/sbs/Downloads/%EC%B2%9C%EB%AC%B8%ED%95%99%20%EB%B0%8F%20%EC%B2%9C%EC%B2%B4%EB%AC%BC%EB%A6%AC%ED%95%99%2015%EC%9E%A5_%EC%88%98%EC%A0%95.pdf
 - interstellar dust and extiction(pdf) : http://w.astro.berkeley.edu/~ay216/08/NOTES/Lecture05-08.pdf
 - model options 'dust_extinction'(Model Flavors) : https://dust-extinction.readthedocs.io/en/latest/dust_extinction/model_flavors.html
 - Dust attenuation, dust emission, and dust temperature in galaxies : https://watermark.silverchair.com/stz1324.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAqIwggKeBgkqhkiG9w0BBwagggKPMIICiwIBADCCAoQGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMVpvRcJeeGG4S-xQvAgEQgIICVQEpa-nXMwsw5JhgVske_zzmhe1tVAh9FpdUN6gNGxE8M2qcFcYNYrrzwqybceOklVy1gx-DXyrYdM8mFuIFl4kmuOGXgrApIiCNTFSaQ2QzQnd4_oLDqC1y8CN51AqM5XLslwdDnb8atdcoOMgwK2zn0fQvYTt4LfVanzWPh6EVepg2YZn4CRbnuW8KsFI8oHcrKK5hR154M7y6IFw__l1Y7rPQDYR_apZ8bF55OecJ4Em9R19odXzO1kZ6_g3xaJImJCif2MwX7ieh1BC3ULmQ097JeGV2hYdL9bZ_aNq3KrB2L75W1ryYpMDNF1GeRP1_q7bIrTOWrvvQ59ojWqjtkLZbOCarVfyFVFt1A4ZSnIFDIcws3-bXEYr3As6wCd8KzjRXrEl2Rl4IzbjbzW0SZB5Y4RtB5OzU2lzAwlvIfcnCLIKdevfxSPKiVYcyh8MdSPnVKEBFqDR5v6OSCUUCMbw9TTCfe9zei6tiwMxQLOkU2S1P5w4Vr8FauvnRpGvMb-1LGvZwr5xzAmSoNfyIPnUAjgSP-E3c2jT06obkZoBQIi-1LANJhgakiO-3kGQFCQ5-TGQ-P9Tr1IecankfHLUFPJpFX_JMQX-H3yDzD-qWSVftwEd4E4IsZrz31tfrs7GktYD43q4rTfE1X6WHQA3Ej3gask7G3T1ocBKEDQaELRKZwhey_aqQkilzXcIlS5a4hyhHVVCPaiOAFNElchCazwRYY9Z2NQHp4YHGyq4ggSd5NJpWBFxpPB5bL00OT0_-8eRHbjZ209YnzdPePGKMSA
 
 
'''
 

'''
성간 물질 : interstellar medium(ISM)
양자 역학 모듈 : QuTiP
상대성이론 모듈 : EinsteinPy

A_λ > 0 : 시선을 따라 흡수나 산란 때문에 올라가는 등급
A_λ = 1.086τ_λ -> 소광에 의한 겉보기 등급의 변화는 시선 방향으로 광학 두께와 거의 비슷하다.


τ_λ : 광학 두께
빛 세기의 상대적 변화 : I_λ / I_λ,0 = e^-τ_λ = T(transmittance)



'''

import numpy as np
import astropy.units as u
from astropy.table import Table
from dust_extinction.parameter_averages import CCM89, F99
from synphot import units, config
from synphot import SourceSpectrum, SpectralElement, Observation, ExtinctionModel1D
from synphot.models import BlackBodyNorm1D
from synphot.spectrum import BaseUnitlessSpectrum
from synphot.reddening import ExtinctionCurve
from astroquery.simbad import Simbad
from astroquery.mast import Observations
import astropy.visualization

# 1. Investigate Extinction Models
# create wavelengths array
wav = np.arange(0.1, 3.0, 0.001)*u.micron

for model in [CCM89, F99]:
    for R in (2.0, 3.0, 4.0):
        # initialize the extinction model
        ext = model(Rv = R)
        plt.plot(1/wav, ext(wav), label = model.name + 'R = ' + str(R))
        
plt.xlabel('$\lambda^{-1}$ ($\mu$m$^{-1}$)')
plt.ylabel('A($\lambda$) / A(V)')
plt.legend(loc = 'best')
plt.title('Some Extinction Laws')
plt.show()


# 2. Deredden a Spectrum
# fetch the archival by utilizing the astroquery
obsTable = Observations.query_object('HD 147933', radius = '1 arcsec')
obsTable_spec = obsTable[obsTable['dataproduct_type'] == 'spectrum']
obsTable_spec.pprint()

obsids = ['3000022829']
dataProductsByID = Observations.get_product_list(obsids)
manifest = Observations.download_products(dataProductsByID)
print('manifestf :', manifest)

# read the downloaded files into an astropy table
t_lwr = Table.read('./mastDownload/IUE/lwr05639/lwr05639mxlo_vo.fits')
#t_lwr = Table.read('./lwr05639mxlo_vo.fits')
print(t_lwr)

wav_UV = t_lwr['WAVE'][0,].quantity
UVflux = t_lwr['FLUX'][0,].quantity

custom_query = Simbad()
custom_query.add_votable_fields('fluxdata(U)', 'fluxdata(B)', 'fluxdata(V)')
phot_table = custom_query.query_object('HD 147933')
Umag = phot_table['FLUX_U']
Bmag = phot_table['FLUX_B']
Vmag = phot_table['FLUX_V']

'''
 - about astroquery.simbad : https://astroquery.readthedocs.io/en/latest/simbad/simbad.html
 - about convert from photometry to flux(Pogson equation) : http://ganymede.nmsu.edu/holtz/a535/ay535notes/node2.html
'''

wav_U = 0.3660 * u.micron
zeroflux_U_nu = 1.81E-23 * u.Watt / (u.m*u.m*u.Hz)
wav_B = 0.4400 * u.micron
zeroflux_B_nu = 4.26E-23 * u.Watt / (u.m*u.m*u.Hz)
wav_V = 0.5530 * u.micron
zeroflux_V_nu = 3.64E-23 * u.Watt / (u.m*u.m*u.Hz)

zeroflux_U = zeroflux_U_nu.to(u.erg/u.AA/u.cm/u.cm/u.s,
                              equivalencies = u.spectral_density(wav_U))
zeroflux_B = zeroflux_B_nu.to(u.erg/u.AA/u.cm/u.cm/u.s,
                              equivalencies = u.spectral_density(wav_B))
zeroflux_V = zeroflux_V_nu.to(u.erg/u.AA/u.cm/u.cm/u.s,
                              equivalencies = u.spectral_density(wav_V))

# convert from photometry to flux using the definition of magnuitude
Uflux = zeroflux_U * 10.**(-0.4*Umag)
Bflux = zeroflux_B * 10.**(-0.4*Bmag)
Vflux = zeroflux_V * 10.**(-0.4*Vmag)

# drwaing the plot
astropy.visualization.quantity_support()

plt.plot(wav_UV, UVflux, 'm', label = 'UV')
plt.plot(wav_V, Vflux, 'ko', label = 'U, B, V')
plt.plot(wav_B, Bflux, 'ko')
plt.plot(wav_U, Uflux, 'ko')
plt.legend(loc = 'best')
plt.ylim(0, 3E-10)
plt.title('rho Oph')
plt.show()

# initialize the extinction model
Rv = 5.0 # usually around 3, but about 5 for this star
Ebv = 0.5
ext = F99(Rv = Rv)

# to extinguish(redden) a spectrum, multiply by the 'ext.extinguish' function
# to unextinguish (deredden), divide by the same 'ext.extinguish' function
plt.semilogy(wav_UV, UVflux, 'm', label = 'UV')
plt.semilogy(wav_V, Vflux, 'ko', label = 'U, B, V')
plt.semilogy(wav_B, Bflux, 'ko')
plt.semilogy(wav_U, Uflux, 'ko')

plt.semilogy(wav_UV, UVflux / ext.extinguish(wav_UV, Ebv = Ebv), 'b',
             label = 'dereddened: EBV = 0.5, RV = 5')
plt.semilogy(wav_V, Vflux / ext.extinguish(wav_V, Ebv = Ebv), 'ro',
             label = 'dereddened: EBV = 0.5, RV = 5')
plt.semilogy(wav_B, Bflux / ext.extinguish(wav_B, Ebv = Ebv), 'ro')
plt.semilogy(wav_U, Uflux / ext.extinguish(wav_U, Ebv = Ebv), 'ro')

plt.legend(loc = 'best')
plt.title('rho Oph')
plt.show()


# 3. Calculate Color Excess with 'synphot'
# synthetic photometry
'''
download the data files

# Optional, for when the STScI ftp server is not answering:
config.conf.vega_file='http://ssb.stsci.edu/cdbs/calspec/alpha_lyr_stis_008.fits'
config.conf.johnson_u_file='http://ssb.stsci.edu/cdbs/comp/nonhst/johnson_u_004_syn.fits'
config.conf.johnson_b_file='http://ssb.stsci.edu/cdbs/comp/nonhst/johnson_b_004_syn.fits'
config.conf.johnson_v_file='http://ssb.stsci.edu/cdbs/comp/nonhst/johnson_v_004_syn.fits'
config.conf.johnson_r_file='http://ssb.stsci.edu/cdbs/comp/nonhst/johnson_r_003_syn.fits'
config.conf.johnson_i_file='http://ssb.stsci.edu/cdbs/comp/nonhst/johnson_i_003_syn.fits'
config.conf.bessel_j_file='http://ssb.stsci.edu/cdbs/comp/nonhst/bessell_j_003_syn.fits'
config.conf.bessel_h_file='http://ssb.stsci.edu/cdbs/comp/nonhst/bessell_h_004_syn.fits'
config.conf.bessel_k_file='http://ssb.stsci.edu/cdbs/comp/nonhst/bessell_k_003_syn.fits'
'''

# getting the filter transmission curves
u_band = SpectralElement.from_filter('johnson_u')
b_band = SpectralElement.from_filter('johnson_b')
v_band = SpectralElement.from_filter('johnson_v')
r_band = SpectralElement.from_filter('johnson_r')
i_band = SpectralElement.from_filter('johnson_i')
j_band = SpectralElement.from_filter('bessel_j')
h_band = SpectralElement.from_filter('bessel_h')
k_band = SpectralElement.from_filter('bessel_k')

# making a background flux
# making 10,000K blackbody
sp = SourceSpectrum(BlackBodyNorm1D, temperature = 10000)
# sp.plot(left = 1, right = 15000, flux_unit = 'flam', title = 'Blackbody')

# get the vega spectrum as the zero point flux
vega = SourceSpectrum.from_vega()
# vega.plot(left = 1, right = 15000, title = 'Vega')


# normalize the blackbody to some chosen magnitude, say V = 10
vmag = 10
v_band = SpectralElement.from_filter('johnson_v')
sp_norm = sp.normalize(vmag * units.VEGAMAG, v_band, vegaspec = vega)
# sp_norm.plot(left = 1, right = 15000, flux_unit = 'flam', title = 'Normed Blackbody')

# initialize the extinction model and choose an extinction of A_v = 2
# to get the 'dust_extinction' model working with 'synphot', create a wavelength array and spectral element with the extinction model

# initialize the extinction model and choose the extinction. A_v = 2
ext = CCM89(Rv = 3.1)
Av = 2.

# create a wavelength array
wav = np.arange(0.1, 3, 0.001) * u.micron

# make the extinction model in synphot using a lookup table
ex = ExtinctionCurve(ExtinctionModel1D, points = wav, lookup_table = ext.extinguish(wav, Av = Av))
sp_ext = sp_norm * ex
# sp_ext.plot(left = 1, right = 15000, flux_unit = 'flam', title = 'Normed Blackbody with Extinction')

# !! synthetic photometry refers to modeling an observation of a star by multiplying the theoretical model for the
# astronomical flux through a certain filter response function, then integrating !!

# observe the star through the filter and integrate to get photometric mag
sp_obs = Observation(sp_ext, v_band)
sp_obs_before = Observation(sp_norm, v_band)
# sp_obs.plot(left = 1, right = 1500, flux_unit = 'flam', title = 'Normed Blackbody with Extinction through V Filter')

# integration and computes magnitudes in the Vega system by 'synphot'
sp_stim_before = sp_obs_before.effstim(flux_unit = 'vegamag', vegaspec = vega)
sp_stim = sp_obs.effstim(flux_unit = 'vegamag', vegaspec = vega)
print('before dust, V =', np.round(sp_stim_before, 1))
print('after dust, V =', np.round(sp_stim, 1))

# calculate extinction and compare to chosen value
Av_calc = sp_stim - sp_stim_before
print('$A_V$ =', np.round(Av_calc, 1))


# reddening with the model extinction curve for comparison in the plot
bands = [u_band, b_band, v_band, r_band, i_band, j_band, h_band, k_band]

for band in bands:
    # calculate photometry with dust
    sp_obs = Observation(sp_ext, band, force = 'extrap')
    obs_effstim = sp_obs.effstim(flux_unit = 'vegamag', vegaspec = vega)
    #calculate photometry without dust
    sp_obs_i = Observation(sp_norm, band, force = 'extrap')
    obs_i_effstim = sp_obs_i.effstim(flux_unit = 'vegamag', vegaspec = vega)
    
    # extinction = mag with dust - mag without dust
    # color excess = extinction at lambda - extinction at V
    color_excess = obs_effstim - obs_i_effstim - Av_calc
    plt.plot(sp_obs_i.effective_wavelength(), color_excess, 'or')
    print(np.round(sp_obs_i.effective_wavelength(), 1), ',', np.round(color_excess, 2))
    
# plot the model extinction curve for comparison
plt.plot(wav, Av * ext(wav) - Av, '--k')
plt.ylim([-2, 2])
plt.xlabel('$\lambda$ (Angstrom)')
plt.ylabel('E($\lambda$-V)')
plt.title('Reddening of T = 10,000K Background Source with Av = 2')
plt.show()

# cannot convert flam unit to vegamag
# -> if you draw the plot employing flam unit, you cannot draw the plot utilizing the vegamag unit
# -> if you want to draw the plot with vegamag unit, you have to eliminate all of the plots which are composing the flam unit
