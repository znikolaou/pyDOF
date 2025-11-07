from include import *
from filters import *
import numpy as np

k=Wavenumber()
kc=np.pi/2.0
flt=ErrorFunction(k, kc, 0.1)

plt.plot(k.get(), flt.getTransferFunction())
plt.grid('on')
plt.show()
