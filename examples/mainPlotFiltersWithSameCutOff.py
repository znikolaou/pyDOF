from include import *
from utils import *
from discreteFilter import *
from filters import *

#define waveno object
wavenumber=Wavenumber()

#define Gaussian filter having Delta/h 4
fltG=Gaussian(wavenumber, 4.0e0)
kc=fltG.getCutOff()
fltH=Helmholtz(wavenumber, 1.0/kc)
fltI=Implicit(wavenumber, -0.5*np.cos(kc))
fltS=SharpCutOff(wavenumber, kc)
fltB=Butterworth(wavenumber, kc, 6)
lovh=Box.getParamForBoxFilterCutOff(kc)
fltBox=Box(wavenumber, lovh)

k=wavenumber.get()
plt.plot(k, fltG.getTransferFunction(), label="Gaussian")
plt.plot(k, fltH.getTransferFunction(), label="Helmholtz")
plt.plot(k, fltI.getTransferFunction(), label="3-point Implicit")
plt.plot(k, fltS.getTransferFunction(), label="Sharp cut-off")
plt.plot(k, fltB.getTransferFunction(), label="Butterworth")
plt.plot(k, fltBox.getTransferFunction(), label="Box")
plt.xlabel("kh")
plt.ylabel("$\hat{G}$")
plt.grid("on")
plt.legend()
plt.show()

