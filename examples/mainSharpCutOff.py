from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *

#define waveno object
wavenumber=Wavenumber()

#define sharp cut-off filter at pi/2.0
flt=SharpCutOff(wavenumber, np.pi/2.0)

#define discrete filter
dflt=DiscreteFilter(wavenumber,nStencil=31)

#run optimiser
opt=Optimiser(flt, dflt, isZeroDampingAtPi=False, isPositive=True, isAdaptive=False, lamda=2.0)
opt.run()
print("Optimised discrete filter coeffs=", dflt.coeffs)
plt.plot(wavenumber.get(), flt.getTransferFunction(), 'ko', label="$\hat{G}$")
plt.plot(wavenumber.get(), dflt.getTransferFunction(), label="$\hat{G}_d$")
plt.grid("on")
plt.xlabel("k")
plt.legend()
plt.show()

#write to file
#output=Output("./output/", coeffs, flt, dflt)
#output.writeToFile("filterCoeffs")
