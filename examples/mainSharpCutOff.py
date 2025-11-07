from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import *

wavenumber=Wavenumber()

flt=SharpCutOff(wavenumber, np.pi/2.0)
dflt=DiscreteFilter(wavenumber, nStencil=33)

opt=Optimiser(flt, dflt, isPositive=True, isMonotonic=True, isAdaptive=False)
opt.run()

out=Output(flt, dflt)
out.write()

plt.plot(wavenumber.get(), flt.getTransferFunction(), 'ko', label="$G$")
plt.plot(wavenumber.get(), dflt.getTransferFunction(), label="$G_d$")
plt.grid("on")
plt.xlabel("$kh$")
plt.legend()
plt.show()

