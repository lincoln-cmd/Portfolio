# exercise1 of modeling1

'''
 - exercis 1)
 - Fit a Polynomial model : Choose fitter wisely
 Utilize the model 'Polynomial1D(degree = 1)' to fit the same data and compare the results
 Create the fake data to make the fit by adding gaussian noise to the data with the function 'np.random.normal(0, 2)'
 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astroquery.vizier import Vizier
import scipy.optimize
%matplotlib inline

N = 100
x1 = np.linspace(0, 4, N) # make an array from 0 to 4 of N elements
y1 = x1**3 - 6*x1**2 + 12*x1 - 9
y1 += np.random.normal(0, 2, size = len(y1)) # add some gaussian noise
sigma = 1.5
y1_err = np.ones(N)*sigma

plt.errorbar(x1, y1, yerr = y1_err, fmt = 'k.')
plt.xlabel('$x_1$')
plt.ylabel('$y_1$')

# fit this data following three steps : model, fitter, and perform fit
model_poly = models.Polynomial1D(degree = 3)
fitter_poly = fitting.LinearLSQFitter()
best_fit_poly = fitter_poly(model_poly, x1, y1, weights = 1.0 / y1_err**2)
print(best_fit_poly)

# utilize the 'SimplexLSQFitter' as fitter
# LinearLSQFitter() : a class performing a linear least square fitting
# SimplexLSQFitter() : simplex algorithm and least squares statistic
fitter_poly_2 = fitting.SimplexLSQFitter()
#best_fit_poly_2 = fitter_poly_2(model_poly, x1, y1, weights = 1.0 / y1_err**2)
print(best_fit_poly_2)
# model : y = c0 + c1 * x + c2 * x^2 + c3 * x^3 ... -> linear in the parameters c_i
# 'SimplexLSQFitter' works better with models that are not linear in the parameters.

'''
 - reference sites
 check which model parameters are a better fit for calculationg : https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
 ordinary least squares : https://en.wikipedia.org/wiki/Ordinary_least_squares#Reduced_chi-squared
 weighted arithmetic mean#Correcting for over- or under-dispersion : https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Correcting_for_over-_or_under-dispersion
'''


def calc_reduced_chi_square(fit, x, y, yerr, N, n_free):
    '''
    fit(array) : values for the fit
    x, y, yerr(array) : data
    N : total number of points
    n_free : number of parameters we are fitting
    '''
    return 1.0 / (N - n_free) * sum(((fit - y) / yerr)**2)

reduced_chi_squared = calc_reduced_chi_square(best_fit_poly(x1), x1, y1, y1_err, N, 4)
print('Reduced Chi Squared with LinearLSQFitter : {0}'.format(reduced_chi_squared))
reduced_chi_squared = calc_reduced_chi_square(best_fit_poly_2(x1), x1, y1, y1_err, N, 4)
print('Reduced Chi Squared with SimplexLSQFitter : {}'.format(reduced_chi_squared))
# the result of first fit is closer to 1, which means that this fit is better.

# compared the two fits, LinearLSQFitter and SimplexLSQFitter, with plot
plt.errorbar(x1, y1, yerr = y1_err, fmt = 'k.')
plt.plot(x1, best_fit_poly(x1), color = 'r', linewidth = 3, label = 'LinearLSQFitter()')
plt.plot(x1, best_fit_poly_2(x1), color = 'g', linewidth = 3, label = 'SimplexLSQFitter()')
plt.xlabel(r'$=log_{10}$(Period [days])')
plt.ylabel('Ks')
plt.legend()

