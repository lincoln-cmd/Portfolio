# exercise1 of modeling1

'''
 - exercis 1)
 - Fit a Polynomial model : Choose fitter wisely
 Utilize the model 'Polynomial1D(degree = 1)' to fit the same data and compare the results
 Create the fake data to make the fit by adding gaussian noise to the data with the function 'np.random.normal(0, 2)'
 
'''
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