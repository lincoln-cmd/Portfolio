# modeling2_ex2

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
from astropy.modeling import Fittable1DModel, Parameter
from astroquery.sdss import SDSS

'''
 - Custom models for unusual function
 Custom models are useful when we want to fit an unusual function to our data
'''
x3 = np.linspace(-2, 0, 100)
y3 = x3**2 * np.exp(-0.5 * (x3)**3 / 2**2)
# a = 2, b = -0.5, c = 3, d = 2
y3 += np.random.normal(0., 0.5, x3.shape)
y3_err = np.ones(x3.shape) * 0.5

plt.errorbar(x3, y3, yerr = y3_err, fmt = '.')
plt.show()

class expnew(Fittable1DModel):
    a = Parameter(default = 1.)
    b = Parameter(default = 1.)
    c = Parameter(default = 1.)
    d = Parameter(default = 1.)
    
    @staticmethod
    def evaluate(x, a, b, c, d):
        return x**a * np.exp(b * (x)**3 / d**2)
    
    @staticmethod
    def fit_deriv(x, a, b, c, d):
        d_a = np.exp(b * (x)**3 / d**2)
        d_b = x**a * np.exp(b * (x)**3)
        d_c = x**a * np.exp((x)**3 / d**2)
        d_d = np.ones(x.shape)
        return [d_a, d_b, d_c, d_d]

exp_model = expnew(a = 2., b = -0.5, c = 3., d = 2.)
fitter = fitting.LevMarLSQFitter()
exp_fit = fitter(exp_model, x3, y3, weights = 1.0 / y3_err**2)

fitter_2 = fitting.SLSQPLSQFitter()
exp_fit_2 = fitter_2(exp_model_2, x3, y3, weights = 1.0 / y3_err**2)

fitter_3 = fitting.SimplexLSQFitter()
exp_fit_3 = fitter_3(exp_model, x3, y3, weights = 1.0/ y3_err*2)

gaussian_model = models.Gaussian1D()
fitter_ga = fitting.LevMarLSQFitter()
gaussian_fit = fitter_ga(gaussian_model, x3, y3, weights = 1.0 / y3_err**2)

sine_model = models.Sine1D()
fitter_sine = fitting.LevMarLSQFitter()
sine_fit = fitter_sine(sine_model, x3, y3, weights = 1.0 / y3_err**2)

#poly_model = models.Polynomial1D(degree = 1)
#fitter_poly = fitting.LevMarLSQFitter()
#poly_fit = fitter_poly(poly_model, x3, y3, weights = 1.0 / y3_err**2)

print('exp_fit : {}'.format(exp_fit))
print('exp_fit_2 : {}'.format(exp_fit_2))
print('exp_fit_3 : {}'.format(exp_fit_3))
print('gaussian : {}'.format(gaussian_fit))
print('sine : {}'.format(sine_fit))
#print('poly : {}'.format(poly_fit))

def calc_reduced_chi_square(fit, x, y, yerr, N, n_free):
    return 1.0 / (N - n_free) * sum(((fit - y) / yerr)**2)

print('exp_fit : {}'.format(calc_reduced_chi_square(exp_fit(x3), x3, y3, y3_err, len(x3), 3)))
print('exp_fit_2 : {}'.format(calc_reduced_chi_square(exp_fit_2(x3), x3, y3, y3_err, len(x3), 3)))
print('exp_fit_3 : {}'.format(calc_reduced_chi_square(exp_fit_3(x3), x3, y3, y3_err, len(x3), 3)))
print('gaussian : {}'.format(calc_reduced_chi_square(gaussian_fit(x3), x3, y3, y3_err, len(x3), 3)))
print('sine : {}'.format(calc_reduced_chi_square(sine_fit(x3), x3, y3, y3_err, len(x3), 3)))
#print('poly : {}'.format(calc_reduced_chi_square(poly_fit(x3), x3, y3, y3_err, len(x3), 3)))



plt.errorbar(x3, y3, yerr = y3_err, fmt = '.')
plt.plot(x3, exp_fit(x3), color = 'r',
         label = 'Custom Model')
#plt.plot(x3, exp_fit_2(x3), color = 'b', linewidth = 4,
 #        label = 'SLSQPLSQFitter')
#plt.plot(x3, exp_fit_3(x3), color = 'y', linewidth = 2,
 #        label = 'SimplexLSQFitter')
plt.plot(x3, gaussian_fit(x3), color = 'g', 
         label = 'Gaussian')
plt.plot(x3, sine_fit(x3), color = 'y',
         label = 'Sine')
#plt.plot(x3, poly_fit(x3), color = 'b', 
 #        label = 'Poly')
plt.title('Graph construction range with in -2 ~ 3')
plt.legend()
plt.show()