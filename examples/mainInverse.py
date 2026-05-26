from include import *
from utils import *
from discreteFilter import *
from filters import *

nButterworth=10
nV=100
kc=np.pi/2.0
wavenumber=Wavenumber()
fltB=Butterworth(wavenumber, kc, nButterworth, fType='lowpass')
k=wavenumber.get()

fltI=Inverse(fltB, 100)
plt.plot(k, fltB.getTransferFunction(), label="Butterworth-F (n="+str(nButterworth)+")")
plt.plot(k, fltI.getTransferFunction(), label="Butterworth-I (nV="+str(nV)+")")
plt.xlabel("kh")
plt.grid("on")
plt.legend(frameon=False)
plt.show()

