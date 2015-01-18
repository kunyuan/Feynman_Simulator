from numpy import linspace,exp
import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt
from scipy.interpolate import LSQUnivariateSpline
from scipy.signal import wiener, filtfilt, butter, gaussian, freqz
from scipy.ndimage import filters

x = linspace(-3, 3, 32)
sigma=1.0/50
y = exp(-x**2) + randn(32)*sigma

#t = [-1,-0.5, 0,0.5, 1]
t = [-1, 0, 1]
s = LSQUnivariateSpline(x, y, t)
print s.get_knots()
print s.get_residual()
print np.sum((y-s(x))**2)/32/sigma**2
xs = linspace(-3, 3, 1000)
ys = s(xs)
plt.plot(x, y, '.-')
plt.plot(xs, ys)
plt.show()
