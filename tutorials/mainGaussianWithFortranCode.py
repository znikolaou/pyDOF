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
flt=Gaussian(wavenumber, 4.0)
nsNominal=flt.getNominalStencil()

#define discrete filter
dflt=DiscreteFilter(wavenumber, nStencil=11)

#run optimiser for forward filter
opt=Optimiser(flt, dflt, isZeroDampingAtPi=True, isAdaptive=False)
opt.run()
coeffs=dflt.coeffs
print(coeffs)

#define target for inverse discrete filter
tdflt=ReconstructedTransferFunction(dflt, nVanCittert=70)

#define inverse discrete filter
vdflt=DiscreteFilter(wavenumber, nStencil=2*nsNominal+1, isForward=False)

#run optimiser for inverse filter
iopt=Optimiser(tdflt, vdflt, isZeroDampingAtPi=True, isAdaptive=True)
iopt.run()

VGTfunc=iopt.getSolution()
print(vdflt.coeffs)

#Output
#out=Output("./fortran/", "coeffs", flt, dflt, tdflt, vdflt)
#out.write()

#Fortran code
fortranOutDir='../tests/fortran/'
fCode=FortranGenerator(fortranOutDir, [dflt, vdflt, vdflt, dflt])
fCode.generateHeader()


#Plot
k=wavenumber.get()
plt.plot(k, flt.getTransferFunction(), 'ko', label="Filter")
plt.plot(k, dflt.getTransferFunction(), label="Discrete filter")
plt.plot(k, VGTfunc, label="Product of inverse and discrete filters")
plt.plot(k, tdflt.getTransferFunction(), 'ro', label="Target rec transfer function")
plt.legend()
plt.show()


