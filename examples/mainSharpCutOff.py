from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import *

wavenumber=Wavenumber()

flt=SharpCutOff(wavenumber, np.pi/2.0)
dflt=DiscreteFilter(wavenumber,nStencil=31)

opt=Optimiser(flt, dflt, isMonotonic=False, isAdaptive=False)
opt.run()

out=Output("./output/", "coeffs", flt, dflt)
out.write()

plt.plot(wavenumber.get(), flt.getTransferFunction(), 'ko', label="$G$")
plt.plot(wavenumber.get(), dflt.getTransferFunction(), label="$G_d$")
plt.grid("on")
plt.xlabel("$kh$")
plt.legend()
plt.show()

