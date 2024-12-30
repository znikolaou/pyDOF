from include import *
from discreteFilter import *
from filters import *
import numpy as np
import matplotlib.pyplot as plt

k=Wavenumber()

#actual filter
filt=Gaussian(k, 4.0)
tFunc=filt.get()

#discrete filter
n=9
dfilt=DiscreteFilter(k, 9)
discreteTFunc=dfilt.get([0.1, 0.8, 0.1])

plt.plot(k.get(), tFunc)
plt.plot(k.get(), discreteTFunc)
plt.show()

