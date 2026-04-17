from include import *
from utils import *
from discreteFilter import *
from filters import *

nButterworth=3
wavenumber=Wavenumber()

fltG=Gaussian(wavenumber, 4.0e0)
kc=fltG.getCutOff()
print('Filter cut-off wavernumber (kc)=', kc)
fltI=Implicit(wavenumber, -0.5*np.cos(kc))
fltS=Sharp(wavenumber, kbounds=[kc])
fltB=Butterworth(wavenumber, kc, nButterworth, fType='lowpass')

k=wavenumber.get()
plt.plot(k, fltG.getTransferFunction(), label="Gaussian")
plt.plot(k, fltI.getTransferFunction(), label="3-point Implicit")
plt.plot(k, fltS.getTransferFunction(), label="Sharp cut-off")
plt.plot(k, fltB.getTransferFunction(), label="Butterworth (n="+str(nButterworth)+")")
plt.xlabel("kh")
plt.ylabel("$\hat{G}$")
plt.grid("on")
plt.legend(frameon=False)
plt.show()

