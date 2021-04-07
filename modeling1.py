# Modeling1 : Make a quick fit using astropy.modeling

# from : https://learn.astropy.org/rst-tutorials/Models-Quick-Fit.html?highlight=filtertutorials

'''
 - about 'astropy.modeling' : https://docs.astropy.org/en/stable/modeling/
'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astroquery.vizier import Vizier
import scipy.optimize
%matplotlib inline

'''
 - Fit a Linear model : Three steps to fit data using astropy.modeling
 This is a catalog of 'Type II Cepheids', which is a type of variable stars that pulsate with a period between
1 and 50 days.
 Measure the Cepheids Period-Luminosity relation using astropy.modeling. This relation states that
if a star has a longer period, the luminosity is higher.

 - about linear fit to real data(from Bhardwaj et al. 2017) : https://ui.adsabs.harvard.edu/abs/2017A%26A...605A.100B/abstract
'''

catalog = Vizier.get_catalogs('J/A+A/605/A100')
# above catalog has tons of infromation, but only utilize the periods and magnitudes with keywords 'Period' and '__Ksmag_'
period = np.array(catalog[0]['Period'])
log_period = np.log10(period)
k_mag = np.array(catalog[0]['__Ksmag_'])
k_mag_err = np.array(catalog[0]['e__Ksmag_']) # 'e__Ksmag_' : refers to the error bars in the magnitude measurements

# the magnitude measurements as a function of period
plt.errorbar(log_period, k_mag, k_mag_err, fmt = 'k.')
plt.xlabel(r'$\log_{10}$(Period [days])')
plt.ylabel('Ks')

# plot : there si a linear relationship between log period and magnitudes. To probe it, make a fit to the data


'''
 - Models in Astropy
 Models in Astropy are known parametrized functions. With this format, it is easier to define and to employ, which means that we do not need
to write the function expression whenever we utilize the model, if we know name. They can be linear or non-linear in the variables.


 - Fitters in Astropy
 Fitters in Astropy are the classes resposable for making the fit. They can be linear or non-linear in the parameters, not the variable like models.

 - about 'Models and Fitting' : https://docs.astropy.org/en/stable/modeling/#using-models
'''
# step1 : model
# choose the model to use for data
model = models.Linear1D()

# step2 : fitter
# choose the fitter to use to fit the model for the data
fitter = fitting.LinearLSQFitter() # LinearLSQFitter : a class performing a linear least square fitting.

# step3 : fit data
# ??? fitter is given the model and the data to perform the fit.
# In this process, the weights are included, which means that values with higher error will have smaller weight, but the data has smaller errors,
# in the fit. In contrast with it, the values with samller error will have higher weights, but the data has higher errors.
# This means that the values are inverse proportion with weights, and the data is direct proportion with weights.
# This is the 'Weighted Linear Least Squares'
# - about that : https://www.mathworks.com/help/curvefit/least-squares-fitting.html
# - : https://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares
best_fit = fitter(model, log_period, k_mag, weights = 1.0 / k_mag_err ** 2)
print(best_fit)

plt.errorbar(log_period, k_mag, k_mag_err, fmt = 'k.')
plt.plot(log_period, best_fit(log_period), color = 'g', linewidth = 3)
plt.xlabel(r'$\log_{10}$(Period [days])')
plt.ylabel('Ks')

'''
 It is possible to fit the data with three lines of code
 1. choose a model
 2. choose a fitter
 3. pass to the fitter the model and the data to perform fit
'''






