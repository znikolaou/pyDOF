from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import * 
from fortranGenerator import *

#define waveno object
wavenumber=Wavenumber()

#define Gaussian filter having Delta/h 4
flt=Implicit(wavenumber, -0.45)

#define discrete filter
dflt=DiscreteFilter(wavenumber)

#run optimiser for forward filter
opt=Optimiser(flt, dflt)
opt.run()
nsForward=dflt.getNs()

#define target for inverse discrete filter
tdflt=ReconstructedTransferFunction(dflt, flt, nVanCittert=100) #70

#define inverse discrete filter
vdflt=DiscreteFilter(wavenumber, nStencil=nsForward, isForward=False)

#run optimiser for inverse filter
iopt=Optimiser(tdflt, vdflt)
iopt.run()

VGTfunc=iopt.getSolution()

#Fortran code
fortranOutDir='../tests/fortran/'
fCode=FortranGenerator(fortranOutDir, 'f90', [dflt, vdflt])
fCode.generateF90Module()

#Plot
k=wavenumber.get()
plt.plot(k, flt.getTransferFunction(), 'ko', label="Filter")
plt.plot(k, dflt.getTransferFunction(), label="Discrete filter")
plt.plot(k, VGTfunc, label="Product of inverse and discrete filters")
plt.plot(k, tdflt.getTransferFunction(), 'ro', label="Target rec transfer function")
plt.legend()
plt.show()


