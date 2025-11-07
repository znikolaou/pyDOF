from include import *
from filters import *
import numpy as np

k=Wavenumber()
kc=np.pi/2.0
scutoff=SharpCutOff(k, kc)

x=np.linspace(-8.0*np.pi/kc, 8.0*np.pi/kc, 33)
print(np.trapz(scutoff.get(x),x))
plt.plot(x*kc/np.pi,scutoff.get(x))
plt.grid('on')
plt.show()
