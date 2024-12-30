from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *

#define waveno object
wavenumber=Wavenumber(100)

#define Implicit filter having a=0.0
flt=Implicit(wavenumber, -0.49)

#define discrete filter
dflt=DiscreteFilter(wavenumber)

#run optimiser
opt=Optimiser(flt, dflt, isZeroDampingAtPi=True, isPositive=False, isAdaptive=True)
opt.run()
coeffs=opt.get()
print(coeffs)
plt.plot(wavenumber.get(), flt.getTransferFunction(), 'ko')
plt.plot(wavenumber.get(), opt.getSolution())
plt.show()

