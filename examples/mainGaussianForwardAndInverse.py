from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import * 

#define waveno object
wavenumber=Wavenumber()

#define Gaussian filter having Delta/h 8
flt=Gaussian(wavenumber, 8.0)
nsNominal=flt.getNominalStencil()

#define discrete filter
dflt=DiscreteFilter(wavenumber, nStencil=nsNominal)

#run optimiser for forward filter
opt=Optimiser(flt, dflt)
opt.run()

#define target for inverse discrete filter
tdflt=ReconstructedTransferFunction(dflt, nVanCittert=100)

#define inverse discrete filter
vdflt=DiscreteFilter(wavenumber, isForward=False)

#run optimiser for inverse filter
iopt=Optimiser(tdflt, vdflt)
iopt.run()

VGTfunc=iopt.getSolution()

#Write to file
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


