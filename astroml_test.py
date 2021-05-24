'''
 <AstroML>
 4. Unsupervised Learning : Density Estimation
 4.1. Bayesian Blocks : Histograms the Rigth Way
 
 - from : https://www.astroml.org/user_guide/density_estimation.html#extreme-deconvolution
 - Extrem Deconvolution of Stellar Data : https://www.astroml.org/book_figures/chapter6/fig_stellar_XD.html#book-fig-chapter6-fig-stellar-xd
'''

import numpy as np
from astroML.plotting import hist
x = np.random.normal(size = 1000)
hist(x, bins = 'blocks')

