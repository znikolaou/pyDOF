from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import *

#define waveno object
wavenumber=Wavenumber()

#define Gaussian filter having Delta/h 4
flt=Gaussian(wavenumber, 16.0e0)

#define a 9-point stencil discrete filter
dflt=DiscreteFilter(wavenumber)

#run optimiser
opt=Optimiser(flt, dflt)
opt.run()
print("Optimised discrete filter coeffs=", dflt.coeffs)
plt.plot(wavenumber.get(), flt.getTransferFunction(), 'ko', label="$\hat{G}$")
plt.plot(wavenumber.get(), dflt.getTransferFunction(), label="$\hat{G}_d$")
plt.grid("on")
plt.xlabel("k")
plt.legend()
plt.show()


