from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *

#define waveno object
wavenumber=Wavenumber()

#define Gaussian filter having Delta/h 4
fltG=Gaussian(wavenumber, 4.0e0)
#define Helmholtz, Implicit, and sharp cut-off filters
fltH=Helmholtz(wavenumber, 1.0/fltG.getCutOff())
fltI=Implicit(wavenumber, -0.5*np.cos(fltG.getCutOff()))
fltS=SharpCutOff(wavenumber, fltG.getCutOff())

k=wavenumber.get()
plt.plot(k, fltG.getTransferFunction(), label="Gaussian")
plt.plot(k, fltH.getTransferFunction(), label="Helmholtz")
plt.plot(k, fltI.getTransferFunction(), label="3-point Implicit")
plt.plot(k, fltS.getTransferFunction(), label="Sharp cut-off")
plt.xlabel("kh")
plt.ylabel("$\hat{G}$")
plt.grid("on")
plt.legend()
plt.show()

