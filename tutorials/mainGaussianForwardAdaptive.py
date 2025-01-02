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
flt=Gaussian(wavenumber, 4.0e0)

#define 9-point stencil discrete filter
dflt=DiscreteFilter(wavenumber)

#run optimiser
opt=Optimiser(flt, dflt, isZeroDampingAtPi=False, isPositive=True, isAdaptive=True)
opt.run()
coeffs=dflt.coeffs
print("Coeffs=", coeffs)

out=Output("output/", "coeffs", flt, dflt)
out.write()

plt.plot(wavenumber.get(), flt.getTransferFunction(), 'ko')
plt.plot(wavenumber.get(), flt.getTransferFunction())
plt.show()

