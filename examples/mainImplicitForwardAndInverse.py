from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import * 

#define waveno object
wavenumber=Wavenumber()

#define Implicit filter
flt=Implicit(wavenumber, 0.45) #strict case

#define discrete filter
dflt=DiscreteFilter(wavenumber)

#run optimiser for forward filter
opt=Optimiser(flt, dflt)
opt.run()
ns=dflt.getNs()
exit()
#define target for inverse discrete filter
tdflt=ReconstructedTransferFunction(dflt, flt, nVanCittert=50)

#define inverse discrete filter
vdflt=DiscreteFilter(wavenumber, isForward=False)

#run optimiser for inverse filter
iopt=Optimiser(tdflt, vdflt)
iopt.run()

VGTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()

#Output
out=Output("./output/", "coeffs", flt, dflt, tdflt, vdflt)
out.write()

#Plot
k=wavenumber.get()
plt.plot(k, flt.getTransferFunction(), 'ko', label="Filter")
plt.plot(k, dflt.getTransferFunction(), label="Discrete filter")
plt.plot(k, VGTfunc, label="Product of inverse and discrete filters")
plt.plot(k, tdflt.getTransferFunction(), 'ro', label="Target rec transfer function")
plt.legend()
plt.grid('on')
plt.show()


