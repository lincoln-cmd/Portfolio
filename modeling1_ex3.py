# modeling1_ex3

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astroquery.vizier import Vizier
import scipy.optimize

'''
 1. model
 2. fitter
 3. perform fit to data
'''



'''
 - about expectation, variance, standard deviation, and covariance Matrix : http://matrix.skku.ac.kr/math4ai-intro/W11/
     https://blogs.ubc.ca/math105/continuous-random-variables/expected-value-variance-standard-deviation/
 - about 
'''
# model : Sine1D
# fitter : 

N3 = 100
x3 = np.linspace(0, 3, N3)
y3 = 5.0 * np.sin(2 * np.pi * x3)
y3 = np.array([y_point + np.random.normal(0, 1) for y_point in y3])
sigma = 1.5
y3_err = np.ones(N3) * sigma


plt.errorbar(x3, y3, yerr = y3_err, fmt = 'k.')
plt.xlabel('$x_3$')
plt.ylabel('$y_3$')


model_sine = models.Sine1D()
fitter_sine = fitting.LevMarLSQFitter()
best_fit_sine = fitter_sine(model_sine, x3, y3, weights = 1.0 / y3_err**2)
print(best_fit_sine)

print(model_sine.param_names)
cov_diag = np.diag(fitter_sine.fit_info['param_cov'])
print(cov_diag)

print('Amplitude : {} +\- {}'.format(best_fit_sine.amplitude.value, np.sqrt(cov_diag[0])))
print('Frequency : {} +\- {}'.format(best_fit_sine.frequency.value, np.sqrt(cov_diag[1])))
print('Phase : {} +\- {}'.format(best_fit_sine.phase.value, np.sqrt(cov_diag[2])))

def calc_reduced_chi_square(fit, x, y, yerr, N, n_free):
    return 1.0 / (N - n_free) * sum(((fit - y)/ yerr)**2)

reduced_chi_squared = calc_reduced_chi_square(best_fit_sine(x3), x3, y3, y3_err, N3, 4)
print('Reduced Chi Squared with LevMarLSQFitter : {}'.format(reduced_chi_squared))

#fitter_sine_2 = fitting.SimplexLSQFitter()
#best_fit_sine_2 = fitter_sine_2(model_sine, x3, y3, weights = 1.0 / y3_err**2)

#reduced_chi_squared_2 = calc_reduced_chi_square(best_fit_sine_2(x3), x3, y3, y3_err, N3, 4)
#print('Reduced Chi Squared with SimplexLSQFitter : {}'.format(reduced_chi_squared_2))

fitter_sine_3 = fitting.SLSQPLSQFitter()
best_fit_sine_3 = fitter_sine_3(model_sine, x3, y3, weights = 1.0 / y3_err**2)

reduced_chi_squared_3 = calc_reduced_chi_square(best_fit_sine_3(x3), x3, y3, y3_err, N3, 4)
print('Reduced Chi Squared with SLSQPLSQFitter : {}'.format(reduced_chi_squared_3))

print(best_fit_sine)
print(best_fit_sine_2)
print(best_fit_sine_3)
'''
 If the consquence of the Reduced Chi Square is closer to one, the fit is better.
 

'''

plt.errorbar(x3, y3, yerr = y3_err, fmt = 'k.')
plt.plot(x3, best_fit_sine(x3), 'r', linewidth = 6, label = 'LevMarLSQFitter')
#plt.plot(x3, best_fit_sine_2(x3), color = 'g', linewidth = 5, label = 'SimplexLSQFitter')
plt.plot(x3, best_fit_sine_3(x3), color = 'b', linewidth = 3, label = 'SLSQPLSQFitter')
plt.xlabel(r'$\log_{10}$(Period [days])')
plt.ylabel('Ks')
plt.legend()


