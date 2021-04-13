import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astroquery.vizier import Vizier
import scipy.optimize
#%matplotlib inline

N = 80
x = np.linspace(0, 3, N)
y = 5.0 * np.arcsinh(2 * np.pi * x)
y = np.array([y_point + np.random.normal(0, 1) for y_point in y])
sigma = 1.5
y_err = np.ones(N) * sigma

plt.errorbar(x, y, yerr = y_err, fmt = 'k.')
plt.xlabel('$x$')
plt.ylabel('$y$')


# model = Sine1D
# fitter = LevMarLSQFitter
model_arcsh = models.Sine1D()
fitter_arcsh = fitting.LevMarLSQFitter()
best_fit_arcsh = fitter_arcsh(model_arcsh, x, y, weights = 1.0 / y_err**2)

print(best_fit_arcsh)
print(model_arcsh.param_names)

cov_diag = np.diag(fitter_arcsh.fit_info['param_cov'])
print(cov_diag)

print('Amplitude : {} +/- {}'.format(best_fit_arcsh.amplitude.value, np.sqrt(cov_diag[0])))
print('Frequency : {} +/- {}'.format(best_fit_arcsh.frequency.value, np.sqrt(cov_diag[1])))
print('Phase : {} +/- {}'.format(best_fit_arcsh.phase.value, np.sqrt(cov_diag[2])))

'''
# fitter = LinearLSQFitter
fitter_arcsh_1 = fitting.LinearLSQFitter()
best_fit_arcsh_1 = fitter_arcsh_1(model_arcsh, x, y, weights = 1.0 / y_err**2)
'''


# model = Gaussian1D
# fitter = LevMarLSQFitter
model_arcsh_2 = models.Gaussian1D()
fitter_arcsh_2 = fitting.LevMarLSQFitter()
best_fit_arcsh_2 = fitter_arcsh_2(model_arcsh_2, x, y, weights = 1.0 / y_err**2)

print(best_fit_arcsh_2)
cov_diag_2 = np.diag(fitter_arcsh_2.fit_info['param_cov'])
print(cov_diag_2)

print('Amplitude : {} +/- {}'.format(best_fit_arcsh_2.amplitude.value, np.sqrt(cov_diag_2[0])))
print('Mean : {} +/- {}'.format(best_fit_arcsh_2.mean.value, np.sqrt(cov_diag_2[1])))
print('Stddev : {} +/- {}'.format(best_fit_arcsh_2.stddev.value, np.sqrt(cov_diag_2[2])))
'''
# fitter = SLSQPLSQFitter()
model_arcsh_2_1 = models.Gaussian1D()
fitter_arcsh_2_1 = fitting.SLSQPLSQFitter()
best_fit_arcsh_2_1 = fitter_arcsh_2_1(model_arcsh_2_1, x, y, weights = 1.0 / y_err**2)

#cov_diag_2_1 = np.diag(fitter_arcsh_2_1.fit_info['param_cov'])
print(best_fit_arcsh_2_1)


# model = KingProjectedAnalytic1D
# fitter = LevMarLSQFitter
model_arcsh_3 = models.KingProjectedAnalytic1D()
fitter_arcsh_3 = fitting.LevMarLSQFitter()
best_fit_arcsh_3 = fitter_arcsh_3(model_arcsh_3, x, y, weights = 1.0 /  y_err**2)

cov_diag_3 = np.diag(fitter_arcsh_3.fit_info['param_cov'])
print(best_fit_arcsh_3)

'''
# model = Polynomial1D
# fitter = SLSQPLSQFitter
model_arcsh_4 = models.Polynomial1D(degree = 10)
fitter_arcsh_4 = fitting.SLSQPLSQFitter()
best_fit_arcsh_4 = fitter_arcsh_4(model_arcsh_4, x, y, weights = 1.0 / y_err**2)

#cov_diag_4 = np.diag(fitter_arcsh_4.fit_info['param_cov'])
print(best_fit_arcsh_4)
'''
# fitter = LevMarLSQFitter
model_arcsh_4_2 = models.Polynomial1D(degree = 0)
fitter_arcsh_4_2 = fitting.LevMarLSQFitter()
best_fit_arcsh_4_2 = fitter_arcsh_4_2(model_arcsh_4_2, x, y, weights = 1.0 / y_err**2)

print(best_fit_arcsh_4_2)
'''



# plot
plt.errorbar(x, y, yerr = y_err, fmt = 'k.')
#plt.plot(x, best_fit_arcsh(x), color = 'r', linewidth = 3,
#         label = 'LevMarLSQFitter(Sine1D)')
#plt.plot(x, best_fit_arcsh(x), 'g', linewidth = 3,
 #        label = 'LinearLSQFitter(Sine1D)')
plt.plot(x, best_fit_arcsh_2(x), color = 'b', linewidth = 3,
         label = 'LevMarLSQFitter(Gaussian1D)')
#plt.plot(x, best_fit_arcsh_3(x), 'g', linewidth = 3, 
#         label = 'LevMarLSQFitter(KingProjectedAnalytic1D)')
#plt.plot(x, best_fit_arcsh_2_1(x), 'y', linewidth = 3,
 #        label = 'SLSQPLSQFitter(Gaussian1D)')
plt.plot(x, best_fit_arcsh_4(x), 'r', linewidth = 3,
         label = 'SLSQPLSQFitter(Polynomial1D)')
plt.xlabel(r'$\log_{10}$(Period [days])')
plt.ylabel('Ks')
plt.legend()