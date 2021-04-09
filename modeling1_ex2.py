# Fit a Gaussian : Let's compare to scipy

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astroquery.vizier import Vizier
import scipy.optimize
#%matplotlib inline

'''
 - make the fake data in the shape of a Gaussian Distribution
 
 
 - about Gaussian Distribution : https://docs.scipy.org/doc/scipy-1.0.0/reference/generated/scipy.optimize.curve_fit.html
'''
mu, sigma, amplitude = 0.0, 10.0, 10.0
N = 100
N2 = 100
x2 = np.linspace(-30, 30, N2)
y2 = amplitude * np.exp(-(x2 - mu)**2 / (2*sigma**2))
y2 = np.array([y_point + np.random.normal(0, 1) for y_point in y2])
sigma = 1
y2_err = np.ones(N) * sigma

plt.errorbar(x2, y2, yerr = y2_err, fmt = 'k.')
plt.xlabel('$x_2$')
plt.ylabel('$y_2$')

# utilize the non-linear fitter, 'LevMarLSQFitter'. This is because the model 'Gaussian1D' is non-linear in the parameters
# LevMarLSQFitter() : Levenberg-Marquardt algorithm and least squares statistic
model_gauss = models.Gaussian1D()
fitter_gauss = fitting.LevMarLSQFitter()
best_fit_gauss = fitter_gauss(model_gauss, x2, y2, weights = 1 / y2_err**2)
print(best_fit_gauss) # LevMarLSQFitter provide the covariance matrix
# this also provides fit parameters error by doing 'fitter.fit_info['param_cov'] because the elements in the diagonal of this matrix are the square of the errors
# Covariance Matrix : https://mathworld.wolfram.com/CovarianceMatrix.html

# check the order of the parameters
print(model_gauss.param_names)
cov_diag = np.diag(fitter_gauss.fit_info['param_cov'])
print(cov_diag)

print('Amplitude : {} +\- {}'.format(best_fit_gauss.amplitude.value, np.sqrt(cov_diag[0])))
print('Mean : {} +\- {}'.format(best_fit_gauss.mean.value, np.sqrt(cov_diag[1])))
print('Standard Deviation : {} +\- {}'.format(best_fit_gauss.stddev.value, np.sqrt(cov_diag[2])))

# create the comparable thing by employing the 'scipy.optimize.curve_fit'
def f(x, a, b, c):
    return a * np.exp(-(x - b)**2 / (2.0 * c**2))
p_opt, p_cov = scipy.optimize.curve_fit(f, x2, y2, sigma = y1_err)
a, b, c = p_opt
best_fit_gauss_2 = f(x2, a, b, c)
print(p_opt)

print('Amplitude : {} +\- {}'.format(p_opt[0], np.sqrt(p_cov[0, 0])))
print('Mean : {} +\- {}'.format(p_opt[1], np.sqrt(p_cov[1, 1])))
print('Standard Deviation : {} +\- {}'.format(p_opt[2], np.sqrt(p_cov[2, 2])))

def calc_reduced_chi_square(fit, x, y, yerr, N, n_free):
    '''
    fit(array) : values for the fit
    x, y, yerr(array) : data
    N : Total number of points
    n_free : number of parameters we are fitting
    '''
    return 1.0 / (N - n_free) * sum(((fit - y) / yerr)**2)
reduced_chi_squared = calc_reduced_chi_square(best_fit_gauss(x2), x2, y2, y2_err, N2, 3)
print('Reduced Chi Squared using astropy.modeling : {}'.format(reduced_chi_squared))

reduced_chi_squared = calc_reduced_chi_square(best_fit_gauss_2, x2, y2, y2_err, N2, 3)
print('Reduced Chi Squared using scipy : {}'.format(reduced_chi_squared))
'''
 Above values have small difference in the Reduced Chi Squared. While the astropy.modeling needs to change the name of the fitter
and the model to perform a completely different fit, scipy requires to remember the expression of the function
'''

plt.errorbar(x2, y2, yerr = y2_err, fmt='k.')
plt.plot(x2, best_fit_gauss(x2), 'g-', linewidth = 6, label = 'astropy.modeling')
plt.plot(x2, best_fit_gauss_2, 'r-', linewidth = 2, label = 'scipy')
plt.xlabel('$x_2$')
plt.ylabel('$y_2$')
plt.legend()

